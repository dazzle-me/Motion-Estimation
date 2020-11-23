import pybind11
from distutils.core import setup, Extension

ext_modules = [
    Extension(
        'me_estimator',
        ['my_motion_estimator.cpp', 'matrix.cpp', 'my_metric.cpp', 'main.cpp'],
        include_dirs=[pybind11.get_include()],
        language='c++',
        extra_compile_args=['-std=c++2a', '-Wall']
    ),
]

setup(
    name='library',
    version='0.0.1',
    author='Martynov Eduard, 210',
    author_email='git.mart.eduard@gmail.com',
    description='ME estimator template for MSU VideoCourse',
    ext_modules=ext_modules,
    requires=['pybind11']
)
