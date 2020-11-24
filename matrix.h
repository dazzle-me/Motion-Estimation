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
    // ISO CPP tells us that if we define function inside the class, eventually
    // compiler makes it inline
    uint32_t getHeight() const  {
        return this -> _height;
    };
    uint32_t getWidth() const {
        return this -> _width;;
    };

    int get(size_t h, size_t w) const {
        return static_cast<int>(_vector[h * getWidth() + w]);
    };
private:
    uint32_t _height;
    uint32_t _width;
    uint32_t _total;
    unsigned char* _vector;
};