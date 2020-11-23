#pragma once

#include <stdlib.h>
#include <stdexcept>
#include <fstream>
#include <iostream>

#include "MotionVector.h"

class Matrix {
public:
    Matrix(
        unsigned char* vector, 
        size_t big_height,
        size_t big_width,
        size_t offset_height = 0,
        size_t offset_width = 0
    );
    Matrix(
        unsigned char* vector,
        size_t big_height, 
        size_t big_width,
        const MotionVector& motion_vector
    );
    size_t getHeight() const;
    size_t getWidth() const;
    int get(size_t h, size_t w) const;
    
    // friend std::ostream& operator<<(
    //     std::ostream& in,
    //     const Matrix& matrix
    // );
private:
    size_t _big_height;
    size_t _big_width;

    size_t _offset;
    unsigned char* _vector;
};