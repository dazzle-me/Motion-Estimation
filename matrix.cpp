#include "matrix.h"

Matrix::Matrix(const unsigned char* vector,
               size_t big_height, 
               size_t big_width,
               size_t offset_height,
               size_t offset_width) :
               _vector(vector),
               _big_height(big_height),
               _big_width(big_width),
               _offset(big_width * offset_height + offset_width) {}

Matrix::Matrix(const unsigned char* vector,
               size_t big_height,
               size_t big_width,
               const MotionVector& motion_vector) :
               _vector(vector),
               _big_height(big_height),
               _big_width(big_width),
               _offset(motion_vector._h * big_width + motion_vector._w) {}

int Matrix::get(size_t h, size_t w) const {
    int real_offset = this -> _offset + h * getWidth() + w;
    if (real_offset >= this -> _big_height * this -> _big_width) {
        throw std::invalid_argument("Real offset : " + std::to_string(real_offset) + " is too big, total : " + std::to_string(this -> _big_height * this -> _big_width) + 
                                    "\nheight : " + std::to_string(h) + 
                                    "\nwidth : " + std::to_string(w));
    }
    return static_cast<int  >(_vector[this -> _offset + h * getWidth() + w]);
}

size_t Matrix::getHeight() const {
    return this -> _big_height;
}

size_t Matrix::getWidth() const {
    return this -> _big_width;;
}

// std::ostream& operator<<(std::ostream& out, const Matrix& matrix) {
//     for (size_t h = 0; h < matrix.getHeight(); h++) {
//         out << "{";
//         for (size_t w = 0; w < matrix.getWidth(); w++) {
//             out << matrix.get(h, w);
//             if (w + 1 != matrix.getWidth()) {
//                 out << ", ";
//             }
//         }
//         out << "}\n";
//     }
//     return out;
// }