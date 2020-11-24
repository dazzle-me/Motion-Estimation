#include "matrix.h"

Matrix::Matrix(unsigned char* vector,
               size_t height,
               size_t width) :
               _vector(vector),
               _height(height),
               _width(width) {
                   this -> _total = height * width;
               }