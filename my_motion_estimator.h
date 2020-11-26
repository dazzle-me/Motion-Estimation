#pragma once

#include <limits>
#include <vector>
#include <array>

#include "matrix.h"
#include "my_metric.h"
#include "MotionVector.h"

#include <pybind11/functional.h>
#include <pybind11/stl.h> 
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

class MotionEstimator {
public:
    MotionEstimator(
        size_t width, 
        size_t height,
        size_t quality,
        bool use_halfpixel
    );
    py::list Estimate(
        py::array_t<unsigned char> _previous_frame,
        py::array_t<unsigned char> _current_frame
    );

    MotionVector FindBlock_BruteForce(
        const Matrix& previous_frame, 
        const Matrix& current_frame,
        int dh,
        int dw
    );
    MotionVector FindBlock_CrossSearch(
        const Matrix& previous_frame, 
        const Matrix& current_frame,
        int dh,
        int dw,
        size_t halfside, 
        int shifted_h,
        int shifted_w, 
        int error
    );
    MotionVector FindBlock_OrthonormalSearch(
        const Matrix& previous_frame, 
        const Matrix& current_frame,
        int dh,
        int dw,
        size_t step_size, 
        int shifted_h,
        int shifted_w, 
        int error, 
        bool is_horizontal
    );
    MotionVector FindBlock_3DRS(
        const Matrix& previous_frame,
        const Matrix& current_frame,
        int dh,
        int dw,
        int shifted_h,
        int shifted_w,
        int error
    );
    py::array_t<unsigned char> Remap(py::array_t<unsigned char> _previous_frame);
    void AssignBlock(unsigned char* result_ptr, size_t dh, size_t dw, MotionVector& motion_vector, Matrix& previous_frame);

    // python setters
    void set_SearchMethod(py::array_t<int> value);
    void set_CrossSearch_Side(py::array_t<int> value);
    void set_CrossSearch_ErrorThreshold(py::array_t<int> value);
    
private:
    enum MODE {
        BruteForce = 0,
        CrossSearch,
        OrthonormalSearch,
        _3DRS
    };
    std::vector<MotionVector> storage; 
    // @params
    const size_t _width;
    const size_t _height;
    const size_t _quality;
    const bool _use_halfpixel;

    static constexpr int _block_size = 16;
    size_t SEARCH_MODE;
    
    int _brute_force_stride; // 1
    int _brute_force_height; // 16
    int _brute_force_width; // 16

    int _cross_search_error_threshold; // 100
    uint32_t _cross_search_side; // 8

    uint32_t _orthonormal_search_step_size;

    static constexpr std::array<std::pair<int, int>, 2> _3DRS_current_frame_offset{
        {{-_block_size, -_block_size}, {-_block_size, _block_size}}
    };
    static constexpr std::array<std::pair<int, int>, 4> _3DRS_previous_frame_offset{
        {{2 * _block_size, 2 * _block_size}, {2 * _block_size, 2 * _block_size},
         {_block_size, -_block_size},         {_block_size, _block_size}}
    };
    static constexpr std::array<std::pair<int, int>, 9> _3DRS_random_fluctuations{
        {{0, 0}, 
         {0, _block_size}, {_block_size, 0}, 
         {0, -_block_size}, {-_block_size, 0},
         {0, 3 * _block_size}, {0, - 3 * _block_size},
         {2 * _block_size, -2 * _block_size}}
    };
    static const size_t _3DRS_random_fluct_size = 9;
    size_t _3DRS_offset_index;
};