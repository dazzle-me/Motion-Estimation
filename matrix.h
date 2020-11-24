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
        size_t height,
        size_t width
    );
    uint32_t getHeight() const;
    uint32_t getWidth() const;

    int get(size_t h, size_t w) const;
private:
    uint32_t _height;
    uint32_t _width;
    uint32_t _total;
    unsigned char* _vector;
};