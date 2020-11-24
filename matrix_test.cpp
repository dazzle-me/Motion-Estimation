#include <iostream>

#include "matrix.h"

int main() {
    unsigned char* a = new unsigned char[50];
    for (unsigned char i = 0; i < 50; i++) {
        a[i] = i;
    }
    Matrix matrix = Matrix(a, 5, 10);
    std::cout << matrix;
}