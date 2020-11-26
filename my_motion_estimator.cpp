#include "my_motion_estimator.h"

namespace py = pybind11;

template<typename T>
std::pair<T, T> operator+(const std::pair<T, T>& a, const std::pair<T, T>& b) {
    return std::make_pair(a.first + b.first, a.second + b.second);
}

MotionEstimator::MotionEstimator(
    size_t width, 
    size_t height,
    size_t quality,
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
    _orthonormal_search_step_size(9)  {}


inline MotionVector MotionEstimator::FindBlock_BruteForce(
    const Matrix& previous_frame,
    const Matrix& current_frame,
    int h,
    int w
) {
    int error = compute_abs_difference(previous_frame, h, w, current_frame, h, w);
    int found_h = 0, found_w = 0;
    for (int dh = -_brute_force_height; dh <= _brute_force_height; dh += this -> _brute_force_stride) {
        for (int dw = -_brute_force_width; dw <= _brute_force_width; dw += this -> _brute_force_stride) {
            int current_error = compute_abs_difference(previous_frame, dh + h, dw + w, current_frame, h, w);
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
    int error
) {
    if (side <= 1) {
        return MotionVector(shifted_h, shifted_w, error);
    }
    size_t found_h = 0, found_w = 0;
    size_t halfside = side >> 1;
    std::array<std::pair<int, int>, 5> candidates{
        {{0, 0}, 
        {-halfside, halfside}, {halfside, halfside},
        {-halfside, -halfside}, {halfside, -halfside}}
    };
    for (const auto&[offset_h, offset_w] : candidates) {
        int current_error = compute_abs_difference(previous_frame, offset_h + shifted_h, offset_w  + shifted_w, current_frame, dh, dw);
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
    return FindBlock_CrossSearch(previous_frame, current_frame, dh, dw, halfside, shifted_h + found_h, shifted_w + found_w, error);
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
    std::vector<std::pair<int, int>> candidates;
    if (is_horizontal) {
        candidates.push_back({0, -step_size});
        candidates.push_back({0, step_size});
    } else {
        candidates.push_back({-step_size, 0});
        candidates.push_back({step_size, 0});
    }
    for (const auto&[offset_h, offset_w] : candidates) {
        int current_error = compute_abs_difference(previous_frame, offset_h + shifted_h, offset_w  + shifted_w, current_frame, dh, dw);
        if (current_error < error) {
            error = current_error;
            found_h = offset_h;
            found_w = offset_w;
        }
    }
    return FindBlock_OrthonormalSearch(previous_frame, current_frame, dh, dw, step_size >> 1, shifted_h + found_h, shifted_w + found_w, error, is_horizontal ^ true); 
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

py::list MotionEstimator::Estimate(
    py::array_t<unsigned char> _previous_frame,
    py::array_t<unsigned char> _current_frame
) {   
    // Since we were storing all of the previous Motion Vectors, we have to clear 'em.
    this -> storage.clear();

    // For every block in current_frame we have to find corresponding (the closest)
    // block in the previous_frame
    unsigned char* previous_frame_ptr = static_cast<unsigned char*>(_previous_frame.request().ptr);
    unsigned char* current_frame_ptr = static_cast<unsigned char*>(_current_frame.request().ptr);
    
    Matrix previous_frame = Matrix(previous_frame_ptr, this -> _height, this -> _width);
    Matrix current_frame = Matrix(current_frame_ptr, this -> _height, this -> _width);
    MotionVector motion_vector = MotionVector(0, 0);
    for (size_t h = 0; h < this -> _height; h += this -> _block_size) {
        for (size_t w = 0; w < this -> _width; w += this -> _block_size) {
            if (this -> SEARCH_MODE == MODE::BruteForce) {
                motion_vector = this -> FindBlock_BruteForce(previous_frame, current_frame, h, w);
            } else if (this -> SEARCH_MODE == MODE::CrossSearch) {
                motion_vector = this -> FindBlock_CrossSearch(previous_frame, current_frame, h, w, this -> _cross_search_side, h, w, std::numeric_limits<int>::max());
            } else if (this -> SEARCH_MODE == MODE::OrthonormalSearch) {
                motion_vector = this -> FindBlock_OrthonormalSearch(previous_frame, current_frame, h, w, this -> _orthonormal_search_step_size, h, w, std::numeric_limits<int>::max(), true);
            } else if (this -> SEARCH_MODE == MODE::_3DRS) {
                motion_vector = this -> FindBlock_3DRS(previous_frame, current_frame, h, w, h, w, std::numeric_limits<int>::max());
            }
            this -> storage.emplace_back(motion_vector);
        }
    }
    return py::list();
    std::vector<std::pair<std::pair<int, int>, int>> castable_storage;
    for (const auto& motion_vector : this -> storage) {
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
            AssignBlock(result_ptr, h, w, storage[index], previous_frame);
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
    Matrix& previous_frame
) {
    for (int h = 0; h < this -> _block_size; h++) {
        for (int w = 0; w < this -> _block_size; w++) {
            result_ptr[(dh + h) * this -> _width + w + dw] = previous_frame.get(h + motion_vector._h, w + motion_vector._w);
        }
    }
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