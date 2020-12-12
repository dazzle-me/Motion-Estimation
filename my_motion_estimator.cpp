#include "my_motion_estimator.h"

namespace py = pybind11;

template<typename T>
std::pair<T, T> operator+(const std::pair<T, T>& a, const std::pair<T, T>& b) {
    return std::make_pair(a.first + b.first, a.second + b.second);
}

MotionEstimator::MotionEstimator(
    int width, 
    int height,
    int quality,
    bool use_halfpixel
) : _width(width),
    _height(height),
    _quality(quality),
    _use_halfpixel(use_halfpixel),
    SEARCH_MODE(MODE::DiamondSearch),
    _3DRS_offset_index(0),
    _brute_force_stride(1),
    _brute_force_width(16),
    _brute_force_height(16),
    _cross_search_side(8),
    _cross_search_error_threshold(100),
    _cross_search_split_threshold(1000),
    _orthonormal_search_step_size(9), 
    _three_step_search_side(8),
    _static_threshold(450),
    is_first(true),
    border_size(16),
    new_height(2 * border_size + height),
    new_width(2 * border_size + width),
    candidate_threshold(450),
    _stop_threshold(450) {
        if (_use_halfpixel) {
            this -> previous_up = new unsigned char[_height * _width];
            this -> previous_left = new unsigned char[_height * _width];
            this -> previous_up_left = new unsigned char[_height * _width];
        }
        previous_extended = new unsigned char[new_height * new_width];
        if (quality == 0) {
            _error_threshold = std::numeric_limits<int>::max();
        } else if (quality == 20) {
            _error_threshold = 75'000;
        } else if (quality == 40) {
            _error_threshold = 50'000;
        } else if (quality == 60) {
            _error_threshold = 25'000;
        } else if (quality == 80) {
            _error_threshold = 15'000;
        } else if (quality == 100) {
            _error_threshold = 10'000;
        }
        this -> large_diamond = {{
                            {-2, 0},
                    {-1, -1},        {-1, 1},
            {0, -2},         {0,  0}        ,{0, 2},
                    {1,  -1},        {1,  1},
                            {2,  0}
        }};
        this -> small_diamond = {{ {0, -1}, {-1, 0}, {0, 1}, {-1, 0} }};

        this -> large_diamond_shifted = {{
                            {-4, 0},
                    {-2, -2},        {-2, 2},
            {0, -4},        {0,  0}        ,{0, 4},
                    {2,  -2},        {2,  2},
                            {4,  0}
        }};
        this -> small_diamond_shifted = {{ {0, -1}, {-1, 0}, {0, 1}, {-1, 0}, 
                                           {0, -2}, {-2, 0}, {0, 2}, {-2, 0} }};
    }

MotionEstimator::~MotionEstimator() {
    delete[] previous_extended;
    if (_use_halfpixel) {
        delete[] previous_up;
        delete[] previous_left;
        delete[] previous_up_left;
    }
}
int MotionEstimator::ComputeAbsDifference(
    const Matrix& domain, 
    int domain_h,
    int domain_w,
    const Matrix& rank,
    int rank_h,
    int rank_w,
    int block_size, 
    int error
)  {
    // Rank blocks are always under control, so just check if domain 
    // block lays inside the picture.
    if (domain_h < 0 || domain_h + block_size - 1 >= domain.getHeight() || 
        domain_w < 0 || domain_w + block_size - 1 >= domain.getWidth()) 
    {
           return std::numeric_limits<int>::max();
    }
    // if (block_size == 16) {
    //     int domain_abs_sum = this -> previous_frame_precomputed[(domain_h + domain_w) >> 8];
    //     int rank_abs_sum = this -> current_frame_precomputed[(rank_w + rank_h) >> 8];
        
    //     if (std::abs(domain_abs_sum - rank_abs_sum) >= error) {
    //         return std::numeric_limits<int>::max();
    //     }
    // }

    int sum = 0;
    for (int h = 0; h < block_size; h++) {
        for (int w = 0; w < block_size; w++) {
            int16_t value = domain.get(h + domain_h, w + domain_w) - rank.get(h + rank_h, w + rank_w);
            sum += value * value;
            // sum += std::abs(domain.get(h + domain_h, w + domain_w) - rank.get(h + rank_h, w + rank_w)); 
            if (sum >= error) {
                return std::numeric_limits<int>::max();
            }
        }
    }
    return sum;
}

inline MotionVector MotionEstimator::FindBlock_BruteForce(
    const Matrix& previous_frame,
    const Matrix& current_frame,
    int h,
    int w,
    int shifted_h,
    int shifted_w
) {
    int error = ComputeAbsDifference(previous_frame, h, w, current_frame, h, w);
    int found_h = 0, found_w = 0;
    for (int dh = -_brute_force_height; dh <= _brute_force_height; dh += this -> _brute_force_stride) {
        for (int dw = -_brute_force_width; dw <= _brute_force_width; dw += this -> _brute_force_stride) {
            int current_error = ComputeAbsDifference(previous_frame, dh + shifted_h, dw + shifted_w, current_frame, h, w, 16, error);
            if (current_error < error) {
                error = current_error;
                found_h = dh;
                found_w = dw;
            }
        }
    }
    return MotionVector(h + found_h, w + found_w, error);
}

inline MotionVector MotionEstimator::FindBlock_CrossSearch(
    const Matrix& previous_frame,
    const Matrix& current_frame,
    int dh,
    int dw,
    size_t side,
    int shifted_h,
    int shifted_w,
    int error,
    int block_size
) {
    if (side <= 1) {
        if (error >= this -> _cross_search_split_threshold && block_size == 16) {
            block_size >>= 1;
            std::vector<MotionVector> subvectors;
            std::array<std::pair<int, int>, 4> shifts = {
                {{0, 0},         {0, block_size},
                 {block_size, block_size},{block_size, 0}}
            };
            for (size_t i = 0; i < 4; i++) {
                subvectors.push_back(FindBlock_CrossSearch(
                    previous_frame, 
                    current_frame,
                    dh + shifts[i].first,
                    dw + shifts[i].second,
                    this -> _cross_search_side,
                    dh + shifts[i].first,
                    dw + shifts[i].second,
                    std::numeric_limits<int>::max(),
                    block_size
                ));
            }
            return MotionVector(subvectors);
        } else { 
            return MotionVector(shifted_h, shifted_w, error);
        }
    }
    size_t found_h = 0, found_w = 0;
    size_t halfside = side >> 1;
    std::array<std::pair<int, int>, 5> candidates{
        {{0, 0}, 
        {-halfside, halfside}, {halfside, halfside},
        {-halfside, -halfside}, {halfside, -halfside}}
    };
    for (const auto&[offset_h, offset_w] : candidates) {
        int current_error = ComputeAbsDifference(previous_frame, offset_h + shifted_h, offset_w  + shifted_w, current_frame, dh, dw, block_size, error);
        if (current_error < error) {
            error = current_error;
            found_h = offset_h;
            found_w = offset_w;
        }
        if (error < this -> _cross_search_error_threshold) {
            return MotionVector(shifted_h + found_h, shifted_w + found_w, error);
        }
    } 
    // Reference point stays the same, but offset updates and side halfs every iteration.
    return FindBlock_CrossSearch(previous_frame, current_frame, dh, dw, halfside, shifted_h + found_h, shifted_w + found_w, error, block_size);
}
inline MotionVector MotionEstimator::FindBlock_OrthonormalSearch(
    const Matrix& previous_frame,
    const Matrix& current_frame,
    int dh,
    int dw,
    size_t step_size,
    int shifted_h,
    int shifted_w,
    int error,
    bool is_horizontal
) {
    if (step_size == 0) {
        return MotionVector(shifted_h, shifted_w, error);
    }
    size_t found_h = 0, found_w = 0;
    std::vector<std::pair<int, int>> candidates{{0, 0}};
    if (is_horizontal) {
        candidates.push_back({0, -step_size});
        candidates.push_back({0, step_size});
    } else {
        candidates.push_back({-step_size, 0});
        candidates.push_back({step_size, 0});
    }
    for (const auto&[offset_h, offset_w] : candidates) {
        int current_error = ComputeAbsDifference(previous_frame, offset_h + shifted_h, offset_w  + shifted_w, current_frame, dh, dw, 16, error);
        if (current_error < error) {
            error = current_error;
            found_h = offset_h;
            found_w = offset_w;
        }
    }
    return FindBlock_OrthonormalSearch(previous_frame, current_frame, dh, dw, step_size >> 1, shifted_h + found_h, shifted_w + found_w, error, is_horizontal ^ true); 
}

inline MotionVector MotionEstimator::FindBlock_ThreeStepSearch(
    const Matrix& previous_frame,
    const Matrix& current_frame,
    int dh,
    int dw,
    size_t side,
    int shifted_h,
    int shifted_w,
    int error
) {
    if (side <= 1) {
        return MotionVector(shifted_h, shifted_w, error);
    }  
    int found_h = 0, found_w = 0;
    size_t halfside = side >> 1;
    std::array<std::pair<int, int>, 9> candidates{
        {{0, 0}, {-halfside, -halfside}, {-halfside, 0}, {-halfside, halfside},
                 {0, -halfside}                        , {0, halfside},
                 {halfside, -halfside},  {halfside, 0} , {halfside, halfside}}
    };
    for (const auto&[offset_h, offset_w] : candidates) {
        int current_error = ComputeAbsDifference(previous_frame, offset_h + shifted_h, offset_w  + shifted_w, current_frame, dh, dw, 16, error);
        if (current_error < error) {
            error = current_error;
            found_h = offset_h;
            found_w = offset_w;
        }
        // Static background
        if (error < this -> _cross_search_error_threshold) {
            return MotionVector(shifted_h + found_h, shifted_w + found_w, error);
        }
    } 
    // Reference point stays the same, but offset updates and side halfs every iteration.
    return FindBlock_ThreeStepSearch(previous_frame, current_frame, dh, dw, halfside, shifted_h + found_h, shifted_w + found_w, error);
}

inline MotionVector MotionEstimator::FindBlock_3DRS(
    const Matrix& previous_frame,
    const Matrix& current_frame,
    int dh,
    int dw,
    int shifted_h,
    int shifted_w,
    int error
) {
    int found_h = 0, found_w = 0;
    for (const auto& const_candidate : _3DRS_current_frame_offset) {
        auto candidate = const_candidate + _3DRS_random_fluctuations[_3DRS_offset_index++];
        _3DRS_offset_index %= _3DRS_random_fluct_size;
        int current_error = compute_abs_difference(previous_frame, candidate.first, candidate.second, current_frame, dh, dw);
        if (current_error < error) {
            error = current_error;
            found_h = candidate.first;
            found_w = candidate.second;
        }
    }
    return MotionVector(0, 0, error);
}

inline MotionVector MotionEstimator::CheckIfStatic(
    const Matrix& previous_frame,
    const Matrix& current_frame,
    int dh,
    int dw,
    int shifted_h,
    int shifted_w,
    int block_size
) {
    int error = this -> ComputeAbsDifference(
        previous_frame,
        shifted_h,
        shifted_w,
        current_frame,
        dh,
        dw,
        block_size
    );
    if (error < this -> _static_threshold) {
        return MotionVector(shifted_h, shifted_w, error);
    }
    return MotionVector(-1, -1, error);
}

inline MotionVector MotionEstimator::GetCandidates(
    const Matrix& previous_frame,
    const Matrix& current_frame,
    int dh,
    int dw
) {
    if (is_first) {
        return MotionVector(0, 0, std::numeric_limits<int>::max());
    }
    // Block size shift
    // ~ value /= 16;
    int error = std::numeric_limits<int>::max();
    int found_h = dh, found_w = dw;
    
    dh >>= 4;
    dw >>= 4;
    // Position of blocks, where candidates are located
    std::array<std::pair<int, int>, 4> previous_frame_c = {{
        {dh + 1, dw - 1}, {dh + 1, dw + 1}, {dh + 2, dw - 2}, {dh + 2, dw + 2}
    }};
    std::array<std::pair<int, int>, 2> current_frame_c = {{
        {dh - 1, dw - 1}, {dh - 1, dw + 1}
    }};
    dh <<= 4;
    dw <<= 4;

    // Previous frame candidates
    // For now implementation has only block-16 size (should I fix this now?)
    for (const auto&[candidate_h, candidate_w] : previous_frame_c) {
        if (candidate_h < 0 || candidate_h >= ((this -> _height) >> 4) ||
            candidate_w < 0 || candidate_w >= ((this -> _width) >> 4)) {
                continue;
        }
        MotionVector candidate = this -> previous_storage.at(candidate_h * (this -> _width / this -> _block_size) + candidate_w);
        
        // Doesn't support splitted version for now
        if (candidate._splitted) continue;

        // Subtract reference point
        candidate._h -= candidate_h * this -> _block_size;
        candidate._w -= candidate_w * this -> _block_size;
        int current_error = ComputeAbsDifference(previous_frame, dh + candidate._h, dw + candidate._w, current_frame, dh, dw, this -> _block_size, error);
        if (current_error < error) {
            error = current_error;
            found_h = candidate._h;
            found_w = candidate._w;
        }
    }
    // Current frame candidates
    for (const auto&[candidate_h, candidate_w] : current_frame_c) {
        if (candidate_h < 0 || candidate_h >= ((this -> _height) >> 4) ||
            candidate_w < 0 || candidate_w >= ((this -> _width) >> 4)) {
                continue;
        }
        
        MotionVector candidate = this -> current_storage.at(candidate_h * (this -> _width >> 4) + candidate_w);
        if (candidate._splitted) continue;
        // Subtract reference point
        candidate._h -= candidate_h * this -> _block_size;
        candidate._w -= candidate_w * this -> _block_size;

        int current_error = ComputeAbsDifference(previous_frame, dh + candidate._h, dw + candidate._w, current_frame, dh, dw, this -> _block_size, error);
        if (current_error < error) {
            error = current_error;
            found_h = candidate._h;
            found_w = candidate._w;
        }
    }
    return MotionVector(dh + found_h, dw + found_w, error);
}

inline MotionVector MotionEstimator::FindBlock_DiamondSearch(
    const Matrix& previous_frame,
    const Matrix& current_frame,
    int dh,
    int dw,
    int shifted_h,
    int shifted_w,
    int error, 
    int block_size,
    int shift_dir
) { 
    MotionVector not_moving = CheckIfStatic(previous_frame, current_frame, dh, dw, shifted_h, shifted_w, block_size);
    not_moving.shift_dir = shift_dir;
    if (not_moving._error < this -> _static_threshold) {
        return not_moving;
    }

    int found_h = 0, found_w = 0;
    for (const auto&[offset_h, offset_w] : this -> large_diamond) {
        // if (_iteration_count >= 35) {
        //     return MotionVector(shifted_h + found_h, shifted_w + found_w, error, shift_dir);
        // }
        int current_error = ComputeAbsDifference(previous_frame, offset_h + shifted_h, offset_w  + shifted_w, current_frame, dh, dw, block_size, error);
        if (current_error < error) {
            error = current_error;
            found_h = offset_h;
            found_w = offset_w;
        }
        if (error <= this -> _stop_threshold) {
            return MotionVector(shifted_h + found_h, shifted_w + found_w, error, shift_dir);
        }
        // _diamond_search_error_map[key] = current_error;
    }
    if (found_h == 0 && found_w == 0) {
        for (const auto&[offset_h, offset_w] : this -> small_diamond) {
            int current_error = ComputeAbsDifference(previous_frame, offset_h + shifted_h, offset_w  + shifted_w, current_frame, dh, dw, block_size, error);
            if (current_error < error) {
                error = current_error;
                found_h = offset_h;
                found_w = offset_w;
            }
        }
        if (error >= _error_threshold && block_size >= 16) {
            std::vector<MotionVector> subvectors;
            std::array<std::pair<int, int>, 4> shifts = {
                {{0, 0},         {0, block_size >> 1},
                 {block_size >> 1, block_size >> 1},{block_size >> 1, 0}}
            };

            int new_error = 0;
            for (int i = 0; i < 4; i++) {
                subvectors.push_back(FindBlock_DiamondSearch(
                    previous_frame, 
                    current_frame,
                    dh + shifts[i].first,
                    dw + shifts[i].second,
                    dh + shifts[i].first,
                    dw + shifts[i].second,
                    std::numeric_limits<int>::max(),
                    block_size >> 1,
                    shift_dir
                ));
                new_error += subvectors[i]._error;
            }
            if (new_error < error) {
                return MotionVector(subvectors, new_error);
            } 
        }
        return MotionVector(found_h + shifted_h, found_w + shifted_w, error, shift_dir);   
    }
    return FindBlock_DiamondSearch(previous_frame, current_frame, dh, dw, found_h + shifted_h, found_w + shifted_w, error, block_size, shift_dir);
}

inline MotionVector MotionEstimator::FindBlock_HexagonSearch(
    const Matrix& previous_frame,
    const Matrix& current_frame,
    int dh,
    int dw,
    int shifted_h,
    int shifted_w,
    int error, 
    int block_size
) {
    std::array<std::pair<int, int>, 9> large_hexagon = {{
                {-2, -1},       {-2, 1},
        
        {0, -2},         {0,  0},       {0, 2},
                
                {2,  -1},       {2, 1}
    }};
    int found_h = 0, found_w = 0;
    for (const auto&[offset_h, offset_w] : large_hexagon) {
        if (_iteration_count >= 50) {
            return MotionVector(shifted_h + found_h, shifted_w + found_w, error);
        }

        int current_error = ComputeAbsDifference(previous_frame, offset_h + shifted_h, offset_w  + shifted_w, current_frame, dh, dw, block_size, error);
        if (current_error < error) {
            error = current_error;
            found_h = offset_h;
            found_w = offset_w;
        }
        if (error <= this -> _static_threshold) {
            return MotionVector(shifted_h + found_h, shifted_w + found_w, error);
        }
        // _diamond_search_error_map[key] = current_error;
    }
    if (found_h == 0 && found_w == 0) {
        std::array<std::pair<int, int>, 4> small_hexagon = {{
            {0, -1}, {-1, 0}, {0, 1}, {-1, 0}
        }};
        for (const auto&[offset_h, offset_w] : small_hexagon) {
            // Make use of previously computed stuff
            int current_error = ComputeAbsDifference(previous_frame, offset_h + shifted_h, offset_w  + shifted_w, current_frame, dh, dw, block_size, error);
            if (current_error < error) {
                error = current_error;
                found_h = offset_h;
                found_w = offset_w;
            }
        }
        if (error >= _error_threshold && block_size == 16) {
            block_size >>= 1;
            std::vector<MotionVector> subvectors;
            std::array<std::pair<int, int>, 4> shifts = {
                {{0, 0},         {0, block_size},
                 {block_size, block_size},{block_size, 0}}
            };
            for (size_t i = 0; i < 4; i++) {
                subvectors.push_back(FindBlock_HexagonSearch(
                    previous_frame, 
                    current_frame,
                    dh + shifts[i].first,
                    dw + shifts[i].second,
                    dh + shifts[i].first,
                    dw + shifts[i].second,
                    std::numeric_limits<int>::max(),
                    block_size
                ));
            }
            return MotionVector(subvectors);
        }
        return MotionVector(found_h + shifted_h, found_w + shifted_w, error);   
    }
    return FindBlock_HexagonSearch(previous_frame, current_frame, dh, dw, found_h + shifted_h, found_w + shifted_w, error, block_size);
}

int MotionEstimator::ComputeSum(
    const Matrix& frame,
    int h, 
    int w
) {
    int sum = 0;
    for (int dh = 0; dh < this -> _block_size; dh++) {
        for (int dw = 0; dw < this -> _block_size; dw++) {
            sum += std::abs(frame.get(dh + h, dw + w));
        }
    }
    return sum;
}

void MotionEstimator::GenerateSubpixelArrays(
    unsigned char* input,
    unsigned char* output_up,
    unsigned char* output_left,
    unsigned char* output_up_left,
    int height,
    int width
) {
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            int cur_pixel_pos = y * width + x;
            int left_pixel_pos = y * width + x - 1;
            int left_up_pixel_pos = (y - 1) * width + x - 1;
            int up_pixel_pos = (y - 1) * width + x;
            
            if (y > 0) {
                output_up[cur_pixel_pos] = (int(input[cur_pixel_pos]) + input[up_pixel_pos]) / 2;
            } else {
                output_up[cur_pixel_pos] = input[cur_pixel_pos];
            }
            if (x > 0) {
                output_left[cur_pixel_pos] = (int(input[cur_pixel_pos]) + input[left_pixel_pos]) / 2;
            } else {
                output_left[cur_pixel_pos] = input[cur_pixel_pos];
            }

            if (x > 0 && y > 0) {   
                output_up_left[cur_pixel_pos] = (
                        int(input[cur_pixel_pos]) + 
                        input[left_pixel_pos] + 
                        input[left_up_pixel_pos] + 
                        input[up_pixel_pos]
                ) / 4;
            } else if (y == 0) {
                output_up_left[cur_pixel_pos] = output_left[cur_pixel_pos];
            } else {
                output_up_left[cur_pixel_pos] = output_up[cur_pixel_pos];
            }
        }
    }
}

void MotionEstimator::ExtendBorders(
    unsigned char* input,
    unsigned char* output
) {
    // Copy frame to center of new
    auto p_output = output + new_width * border_size + border_size;
    auto p_input = input;
    for (size_t y = 0; y < _height; ++y, p_output += 2 * border_size) {
        for (size_t x = 0; x < _width; ++x, ++p_output, ++p_input) {
            *p_output = *p_input;
        }
    }

    // Left and right borders.
    p_output = output + new_width * border_size;
    for (size_t y = 0; y < _height; ++y) {
        memset(p_output, p_output[border_size], border_size);
        p_output += border_size + _width;
        memset(p_output, p_output[-1], border_size);
        p_output += border_size;
    }

    // Top and bottom borders.
    p_output = output;
    auto p_output_reference_row = p_output + new_width * border_size;

    for (size_t y = 0; y < border_size; ++y) {
        memcpy(p_output, p_output_reference_row, new_width);
        p_output += new_width;
    }
    p_output = output + new_width * (_height + border_size);
    p_output_reference_row = p_output_reference_row - new_width;

    for (size_t y = 0; y < border_size; ++y) {
        memcpy(p_output, p_output_reference_row, new_width);
        p_output += new_width;
    }
}

void MotionEstimator::Estimate(
    py::array_t<unsigned char> _previous_frame,
    py::array_t<unsigned char> _current_frame
) {   
    previous_storage = std::move(current_storage);
    this -> current_storage.clear();
    
    // For every block in current_frame we have to find corresponding (the closest)
    // block in the previous_frame
    unsigned char* previous_frame_ptr = static_cast<unsigned char*>(_previous_frame.request().ptr);
    unsigned char* current_frame_ptr = static_cast<unsigned char*>(_current_frame.request().ptr);

    // If we want to extend borders

    // ExtendBorders(previous_frame_ptr, previous_extended);
    // previous_frame_ptr = previous_extended; 
    // int first_row_offset = new_width * border_size + border_size; 
    // frames = {Matrix(previous_extended, this -> new_height, this -> new_width, first_row_offset, this -> border_size)};
 
    frames = {Matrix(previous_frame_ptr, this -> _height, this -> _width)};
    if (_use_halfpixel) {
        GenerateSubpixelArrays(
            previous_frame_ptr,
            this -> previous_up,
            this -> previous_left,
            this -> previous_up_left,
            this -> _height,
            this -> _width
        );
        this -> frames.push_back(Matrix(this -> previous_up, this -> _height, this -> _width));
        this -> frames.push_back(Matrix(this -> previous_left, this -> _height, this -> _width));
        this -> frames.push_back(Matrix(this -> previous_up_left, this -> _height, this -> _width));
    }

    Matrix current_frame = Matrix(current_frame_ptr, this -> _height, this -> _width);

    for (int h = 0; h < this -> _height; h += this -> _block_size) {
        for (int w = 0; w < this -> _width; w += this -> _block_size) {
            MotionVector found_motion_vector = MotionVector(0, 0, std::numeric_limits<int>::max(), 0);
            for (int shift_dir = 0; shift_dir < frames.size(); shift_dir++, this -> _iteration_count = 0) {
                MotionVector candidate = GetCandidates(frames[shift_dir], current_frame, h, w);
                if (candidate._error < this -> candidate_threshold) {
                    candidate.shift_dir = shift_dir;
                    found_motion_vector = candidate;
                    break;
                }
                MotionVector motion_vector = this -> FindBlock_DiamondSearch(frames[shift_dir], current_frame, h, w, h, w, std::numeric_limits<int>::max(), this -> _block_size, shift_dir);
                // if (this -> SEARCH_MODE == MODE::BruteForce) {
                //     motion_vector = this -> FindBlock_BruteForce(frames[shift_dir], current_frame, h, w, h, w);
                // } else if (this -> SEARCH_MODE == MODE::CrossSearch) {
                //     motion_vector = this -> FindBlock_CrossSearch(frames[shift_dir], current_frame, h, w, this -> _cross_search_side, h, w, std::numeric_limits<int>::max(), this -> _block_size);
                // } else if (this -> SEARCH_MODE == MODE::OrthonormalSearch) {
                //     motion_vector = this -> FindBlock_OrthonormalSearch(frames[shift_dir], current_frame, h, w, this -> _orthonormal_search_step_size, h, w, std::numeric_limits<int>::max(), true);
                // } else if (this -> SEARCH_MODE == MODE::_3DRS) {
                //     motion_vector = this -> FindBlock_3DRS(frames[shift_dir], current_frame, h, w, h, w, std::numeric_limits<int>::max());
                // } else if (this -> SEARCH_MODE == MODE::ThreeStepSearch) {
                //     motion_vector = this -> FindBlock_ThreeStepSearch(frames[shift_dir], current_frame, h, w, this -> _three_step_search_side, h, w, std::numeric_limits<int>::max());
                // } else if (this -> SEARCH_MODE == MODE::DiamondSearch) {
                //     motion_vector = this -> FindBlock_DiamondSearch(frames[shift_dir], current_frame, h, w, h, w, std::numeric_limits<int>::max(), this -> _block_size, shift_dir);
                // } else if(this -> SEARCH_MODE == MODE::HexagonSearch) {
                //     motion_vector = this -> FindBlock_HexagonSearch(frames[shift_dir], current_frame, h, w, h, w, std::numeric_limits<int>::max(), this -> _block_size);
                // }
                if (motion_vector._error < found_motion_vector._error) {
                    found_motion_vector = motion_vector;
                }
            }
            this -> current_storage.emplace_back(found_motion_vector);
        }
    }
    return; 
}

py::array_t<unsigned char> MotionEstimator::Remap(
    py::array_t<unsigned char> _previous_frame
) {
    py::array_t<unsigned char> result(this -> _height * this -> _width);
    unsigned char* result_ptr = static_cast<unsigned char*>(result.request().ptr);
    
    int index = 0;
    for (int h = 0; h < this -> _height; h += this -> _block_size) {
        for (int w = 0; w < this -> _width; w += this -> _block_size, index++) {
            if (current_storage[index]._splitted) {
                int block_size = 8;
                std::array<std::pair<int, int>, 4> shifts = {
                    {{0, 0},         {0, block_size},
                    {block_size, block_size},{block_size, 0}}
                };
                for (int i = 0; i < 4; i++) {
                    int shift_dir = current_storage[index]._subvectors[i].shift_dir;
                    AssignBlock(result_ptr, h + shifts[i].first, w + shifts[i].second, current_storage[index]._subvectors[i], this -> frames[shift_dir], block_size);
                }
            } else {
                int shift_dir = current_storage[index].shift_dir;
                AssignBlock(result_ptr, h, w, current_storage[index], this -> frames[shift_dir], this -> _block_size);
            }
            // AssignBlock(result_ptr, h, w, current_storage[index], this -> _block_size);
        }
    }
    result.resize({this -> _height, this -> _width});
    return result;
}

void MotionEstimator::AssignBlock(
    unsigned char* result_ptr,
    int dh,
    int dw,
    MotionVector& motion_vector,
    const Matrix& previous_frame,
    int block_size
) {
    // if (motion_vector._splitted) {
    //     std::array<std::pair<int, int>, 4> shifts = {
    //         {{0, 0},         {0, block_size},
    //         {block_size, block_size},{block_size, 0}}
    //     };
    //     for (size_t i = 0; i < 4; i++) {
    //         AssignBlock(result_ptr, dh + shifts[i].first, dw + shifts[i].second, motion_vector._subvectors[i], block_size >> 1);
    //     }
    //     return;
    // } // else it is not splitted
    for (int h = 0; h < block_size; h++) {
        for (int w = 0; w < block_size; w++) {
            result_ptr[(dh + h) * this -> _width + w + dw] = previous_frame.get(h + motion_vector._h, w + motion_vector._w);
        }
    }
}

inline int MotionEstimator::GetKey(int h, int w) const {
    return this -> _width * h + w;
}

void MotionEstimator::set_SearchMethod(py::array_t<int> value) {
    this -> SEARCH_MODE = *(int*)value.request().ptr;
}

void MotionEstimator::set_CrossSearch_Side(py::array_t<int> value) {
    this -> _cross_search_side = *(int*)value.request().ptr;
}
void MotionEstimator::set_CrossSearch_ErrorThreshold(py::array_t<int> value) {
    this -> _cross_search_error_threshold = *(int*)value.request().ptr;
}