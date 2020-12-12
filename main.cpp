#include <pybind11/pybind11.h>

#include "my_motion_estimator.h"

namespace py = pybind11;

PYBIND11_MODULE(me_estimator, m) {
    py::class_<MotionEstimator>(m, "MotionEstimator")
        .def(py::init<size_t, size_t, size_t, bool>())
        .def("Estimate", &MotionEstimator::Estimate)
        .def("Remap", &MotionEstimator::Remap)
        .def("AssignBlock", &MotionEstimator::AssignBlock)
        .def("set_SearchMethod", &MotionEstimator::set_SearchMethod)
        .def("set_CrossSearch_ErrorThreshold", &MotionEstimator::set_CrossSearch_ErrorThreshold)
        .def("set_CrossSearch_Side", &MotionEstimator::set_CrossSearch_Side);
    py::class_<Matrix>(m, "Matrix")
        .def(py::init<unsigned char*, size_t, size_t>())
        .def("getHeight", &Matrix::getHeight)
        .def("getWidth", &Matrix::getWidth)
        .def("get", &Matrix::get);
    py::class_<MotionVector>(m, "MotionVector")
        .def(py::init<int, int>())
        .def(py::init<int, int, int>())
        .def(py::init<const std::vector<MotionVector>&, int>())
        .def(py::init<int, int, int, int>())
        .def("getHeight", &MotionVector::getHeight)
        .def("getWidth", &MotionVector::getWidth)
        .def("getError", &MotionVector::getError)
        .def("is_splitted", &MotionVector::is_splitted)
        .def("getSubvectors", &MotionVector::getSubvectors);
};