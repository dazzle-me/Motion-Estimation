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
    MotionVector FindBlock(const Matrix& previous_frame, const Matrix& current_frame, size_t dh, size_t dw);
    py::array_t<unsigned char> Remap(py::array_t<unsigned char> _previous_frame);
    void AssignBlock(unsigned char* result_ptr, size_t dh, size_t dw, MotionVector& motion_vector, Matrix& previous_frame);
    // ISO CPP tells us that if we define function inside the class, eventually
    // compiler makes it inline
private:
    std::vector<MotionVector> storage;
    
    const size_t _block_size;

    const size_t _width;
    const size_t _height;
    const size_t _quality;
    const bool _use_halfpixel;
};