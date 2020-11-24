#include "matrix.h"

Matrix::Matrix(unsigned char* vector,
               size_t height,
               size_t width) :
               _vector(vector),
               _height(height),
               _width(width) {
                   this -> _total = height * width;
               }

int Matrix::get(size_t h, size_t w) const {
    return static_cast<int>(_vector[h * getWidth() + w]);
}

uint32_t Matrix::getHeight() const {
    return this -> _height;
}

uint32_t Matrix::getWidth() const {
    return this -> _width;;
}