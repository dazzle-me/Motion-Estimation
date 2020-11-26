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
        int height,
        int width
    );
    // ISO CPP tells us that if we define function inside the class, eventually
    // compiler makes it inline
    int getHeight() const  {
        return this -> _height;
    };
    int getWidth() const {
        return this -> _width;;
    };

    int get(size_t h, size_t w) const {
        return static_cast<int>(_vector[h * getWidth() + w]);
    };
private:
    int _height;
    int _width;
    int _total;
    unsigned char* _vector;
};