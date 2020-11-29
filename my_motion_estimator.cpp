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
    SEARCH_MODE(MODE::BruteForce),
    _3DRS_offset_index(0),
    _brute_force_stride(1),
    _brute_force_width(16),
    _brute_force_height(16),
    _cross_search_side(8),
    _cross_search_error_threshold(100),
    _cross_search_split_threshold(1000),
    _orthonormal_search_step_size(9), 
    _three_step_search_side(8),
    _static_threshold(750),
    _error_threshold(10000),
    is_first(true),
    border_size(16),
    new_height(2 * border_size + height),
    new_width(2 * border_size + width) {}


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
    if (domain_h < 0 || domain_h + block_size >= domain.getHeight() + 1 || 
        domain_w < 0 || domain_w + block_size >= domain.getWidth() + 1) 
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
            int value = domain.get(h + domain_h, w + domain_w) - rank.get(h + rank_h, w + rank_w);
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
    int w
) {
    int error = ComputeAbsDifference(previous_frame, h, w, current_frame, h, w);
    int found_h = 0, found_w = 0;
    for (int dh = -_brute_force_height; dh <= _brute_force_height; dh += this -> _brute_force_stride) {
        for (int dw = -_brute_force_width; dw <= _brute_force_width; dw += this -> _brute_force_stride) {
            int current_error = ComputeAbsDifference(previous_frame, dh + h, dw + w, current_frame, h, w, 16, error);
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
        is_first = false;
        return MotionVector(0, 0, std::numeric_limits<int>::max());
    }
    // Block size shift
    // ~ value /= 16;
    dh >>= 4;
    dw >>= 4; 

    int error = std::numeric_limits<int>::max();
    int found_h = 0, found_w = 0;
    // Position of blocks, where candidates located
    std::array<std::pair<int, int>, 4> previous_frame_c = {{
        {dh + 1, dw - 1}, {dh + 1, dw + 1}, {dh + 2, dw - 2}, {dh + 2, dw + 2}
    }};
    std::array<std::pair<int, int>, 2> current_frame_c = {{
        {dh - 1, dw - 1}, {dh - 1, dw + 1}
    }};
    
    // Previous frame candidates
    // For now implementation has only block-16 size (should I fix this now?)
    for (const auto&[offset_h, offset_w] : previous_frame_c) {
        MotionVector candidate = this -> previous_storage.at((dh + offset_h) * this -> _width + offset_w);
        // Doesn't support splitted version for now
        if (candidate._splitted) continue;

        int current_error = ComputeAbsDifference(previous_frame, candidate._h, candidate._w, current_frame, dh, dw, 16, error);
        if (current_error < error) {
            error = current_error;
            found_h = candidate._h;
            found_w = candidate._w;
        }
    }
    // Current frame candidates
    for (const auto&[offset_h, offset_w] : current_frame_c) {
        MotionVector candidate = this -> current_storage.at((dh + offset_h) * this -> _width + offset_w);
        if (candidate._splitted) continue;

        int current_error = ComputeAbsDifference(current_frame, candidate._h, candidate._w, current_frame, dh, dw, 16, error);
        if (current_error < error) {
            error = current_error;
            found_h = candidate._h;
            found_w = candidate._w;
        }
    }
    return MotionVector(found_h, found_w, error);
}

inline MotionVector MotionEstimator::FindBlock_DiamondSearch(
    const Matrix& previous_frame,
    const Matrix& current_frame,
    int dh,
    int dw,
    int shifted_h,
    int shifted_w,
    int error, 
    int block_size
) {
    MotionVector not_moving = this -> CheckIfStatic(previous_frame, current_frame, dh, dw, shifted_h, shifted_w, block_size);
    if (not_moving._h != -1) return not_moving;

    std::array<std::pair<int, int>, 9> large_diamond = {{
                         {-2, 0},
                {-1, -1},        {-1, 1},
        {0, -2},         {0,  0}        ,{0, 2},
                {1,  -1},        {1,  1},
                         {2,  0}
    }};
    int found_h = 0, found_w = 0;
    for (const auto&[offset_h, offset_w] : large_diamond) {
        // Then is has been previosly calculated, skip this because
        // 1) If it was the minimum, it has been added to map and stored in (shifted_h, shifted_w) and error
        // 2) Otherwise, it's useless piece of block, and hence we don't need it
        // std::pair<int, int> key = std::make_pair(GetKey(offset_h + shifted_h, offset_w  + shifted_w), GetKey(dh, dw));
        // if (_diamond_search_error_map.count(key)) {
        //     continue;
        // }
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
        std::array<std::pair<int, int>, 4> small_diamond = {{
            {0, -1}, {-1, 0}, {0, 1}, {-1, 0}
        }};
        for (const auto&[offset_h, offset_w] : small_diamond) {
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
                subvectors.push_back(FindBlock_DiamondSearch(
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
    return FindBlock_DiamondSearch(previous_frame, current_frame, dh, dw, found_h + shifted_h, found_w + shifted_w, error, block_size);
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

py::list MotionEstimator::Estimate(
    py::array_t<unsigned char> _previous_frame,
    py::array_t<unsigned char> _current_frame
) {   
    previous_storage = current_storage;
    this -> current_storage.clear();
    
    // Since we were storing all of the previous Motion Vectors, we have to clear 'em.
    // this -> _diamond_search_error_map.clear();

    // For every block in current_frame we have to find corresponding (the closest)
    // block in the previous_frame
    unsigned char* previous_frame_ptr = static_cast<unsigned char*>(_previous_frame.request().ptr);
    unsigned char* current_frame_ptr = static_cast<unsigned char*>(_current_frame.request().ptr);
    
    previous_frame_with_borders_ptr = new unsigned char[new_height * new_width];
    ExtendBorders(previous_frame_ptr, previous_frame_with_borders_ptr);

    Matrix previous_frame = Matrix(previous_frame_with_borders_ptr, this -> new_height, this -> new_width);
    Matrix current_frame = Matrix(current_frame_ptr, this -> _height, this -> _width);

    // for (int h = 0; h < this -> _height; h += this -> _block_size) {
    //     for (int w = 0; w < this -> _width; w += this -> _block_size) {
    //         previous_frame_precomputed.push_back(this -> ComputeSum(previous_frame, h, w));
    //         current_frame_precomputed.push_back(this -> ComputeSum(current_frame, h, w));
    //     }
    // }

    MotionVector motion_vector = MotionVector(0, 0);
    for (int h = 0; h < this -> _height; h += this -> _block_size) {
        for (int w = 0; w < this -> _width; w += this -> _block_size, this -> _iteration_count = 0) {
            // MotionVector candidate = GetCandidates(
            //     previous_frame, current_frame, h, w
            // );
            if (this -> SEARCH_MODE == MODE::BruteForce) {
                motion_vector = this -> FindBlock_BruteForce(previous_frame, current_frame, h, w);
            } else if (this -> SEARCH_MODE == MODE::CrossSearch) {
                motion_vector = this -> FindBlock_CrossSearch(previous_frame, current_frame, h, w, this -> _cross_search_side, h, w, std::numeric_limits<int>::max(), this -> _block_size);
            } else if (this -> SEARCH_MODE == MODE::OrthonormalSearch) {
                motion_vector = this -> FindBlock_OrthonormalSearch(previous_frame, current_frame, h, w, this -> _orthonormal_search_step_size, h, w, std::numeric_limits<int>::max(), true);
            } else if (this -> SEARCH_MODE == MODE::_3DRS) {
                motion_vector = this -> FindBlock_3DRS(previous_frame, current_frame, h, w, h, w, std::numeric_limits<int>::max());
            } else if (this -> SEARCH_MODE == MODE::ThreeStepSearch) {
                motion_vector = this -> FindBlock_ThreeStepSearch(previous_frame, current_frame, h, w, this -> _three_step_search_side, h, w, std::numeric_limits<int>::max());
            } else if (this -> SEARCH_MODE == MODE::DiamondSearch) {
                motion_vector = this -> FindBlock_DiamondSearch(previous_frame, current_frame, h, w, h + border_size, w + border_size, std::numeric_limits<int>::max(), this -> _block_size);
            }
            this -> current_storage.emplace_back(motion_vector);
        }
    }
    std::vector<std::pair<std::pair<int, int>, int>> castable_storage;
    for (const auto& motion_vector : this -> current_storage) {
        castable_storage.push_back(std::make_pair(std::make_pair(motion_vector._h, motion_vector._w), motion_vector._error));
    }
    return py::cast(castable_storage); 
}

py::array_t<unsigned char> MotionEstimator::Remap(
    py::array_t<unsigned char> _previous_frame
) {
    unsigned char* previous_frame_ptr = static_cast<unsigned char*>(_previous_frame.request().ptr);
    Matrix previous_frame = Matrix(previous_frame_ptr, this -> _height, this -> _width);
    size_t index = 0;

    py::array_t<unsigned char> result(this -> _height * this -> _width);
    unsigned char* result_ptr = static_cast<unsigned char*>(result.request().ptr);
    
    for (size_t h = 0; h < this -> _height; h += this -> _block_size) {
        for (size_t w = 0; w < this -> _width; w += this -> _block_size, index++) {
            if (current_storage[index]._splitted) {
                int block_size = 8;
                std::array<std::pair<int, int>, 4> shifts = {
                    {{0, 0},         {0, block_size},
                    {block_size, block_size},{block_size, 0}}
                };
                for (size_t i = 0; i < 4; i++) {
                    AssignBlock(result_ptr, h + shifts[i].first, w + shifts[i].second, current_storage[index]._subvectors[i], previous_frame, block_size);
                }
            } else {
                AssignBlock(result_ptr, h, w, current_storage[index], previous_frame, this -> _block_size);
            }
        }
    }
    result.resize({this -> _height, this -> _width});
    return result;
}

void MotionEstimator::AssignBlock(
    unsigned char* result_ptr,
    size_t dh,
    size_t dw,
    MotionVector& motion_vector,
    Matrix& previous_frame,
    int block_size
) {
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