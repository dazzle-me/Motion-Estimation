#include "matrix.h"

Matrix::Matrix(unsigned char* vector,
               int height,
               int width) :
               _vector(vector),
               _height(height),
               _width(width) {
                   this -> _total = height * width;
               }