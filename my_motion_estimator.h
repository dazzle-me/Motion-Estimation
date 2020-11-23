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
        py::array_t<unsigned char> current_frame,
        py::array_t<unsigned char> previous_frame
    );
    MotionVector FindBlock(const Matrix& rank_block, unsigned char* previous_frame_ptr);
    py::array_t<unsigned char> Remap(py::array_t<unsigned char> previous_frame);
    void AssignBlock(unsigned char* result_ptr, const Matrix& domain_block, size_t h, size_t w);
private:
    std::vector<MotionVector> storage;
    
    const size_t _block_size;

    const size_t _width;
    const size_t _height;
    const size_t _quality;
    const bool _use_halfpixel;
};