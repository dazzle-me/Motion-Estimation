#pragma once

#include <limits>

#include "matrix.h"
#include "my_metric.h"
#include "MotionVector.h"

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
    void Estimate(
        py::array_t<unsigned char> _current_frame,
        py::array_t<unsigned char> _previous_frame
    );
    MotionVector FindBlock_BadBruteForce(const Matrix& previous_frame, 
                                         const Matrix& current_frame, size_t dh, size_t dw);
    MotionVector FindBlock_CrossSearch(const Matrix& previous_frame, 
                                      const Matrix& current_frame, size_t dh, size_t dw, size_t halfside, 
                                                                   size_t shifted_h, size_t shifted_w);
    py::array_t<unsigned char> Remap(py::array_t<unsigned char> _previous_frame);
    void AssignBlock(unsigned char* result_ptr, size_t dh, size_t dw, MotionVector& motion_vector, Matrix& previous_frame);
private:
    std::vector<MotionVector> storage;
    
    const size_t _block_size;
    const size_t _width;
    const size_t _height;
    const size_t _quality;
    const bool _use_halfpixel;

    // @params
    static const uint32_t _cross_search_halfside = 8;
};