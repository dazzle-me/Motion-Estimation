#include <pybind11/pybind11.h>

#include "my_motion_estimator.h"

namespace py = pybind11;

PYBIND11_MODULE(me_estimator, m) {
    py::class_<MotionEstimator>(m, "MotionEstimator")
        .def(py::init<size_t, size_t, size_t, bool>())
        .def("Estimate", &MotionEstimator::Estimate)
        .def("Remap", &MotionEstimator::Remap)
        .def("AssignBlock", &MotionEstimator::AssignBlock);
    py::class_<Matrix>(m, "Matrix")
        .def(py::init<unsigned char*, size_t, size_t>())
        .def("getHeight", &Matrix::getHeight)
        .def("getWidth", &Matrix::getWidth)
        .def("get", &Matrix::get);
};