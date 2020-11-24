#include "my_motion_estimator.h"

MotionEstimator::MotionEstimator(
    size_t width, 
    size_t height,
    size_t quality,
    bool use_halfpixel
) : _width(width),
    _height(height),
    _quality(quality),
    _use_halfpixel(use_halfpixel),
    _block_size(16) {}


MotionVector MotionEstimator::FindBlock(const Matrix& previous_frame,
                                        const Matrix& current_frame,
                                        size_t dh,
                                        size_t dw)  
{
    // Brute-force search
    int error = std::numeric_limits<int>::max();
    size_t found_h = 0, found_w = 0;
    for (int h = 0; h < this -> _height; h += this -> _block_size) {
        for (int w = 0; w < this -> _width; w += this -> _block_size) {
            int current_error = compute_abs_difference(previous_frame, h, w, current_frame, dh, dw);
            if (current_error < error) {
                error = current_error;
                found_h = h;
                found_w = w;
            }
        }
    }
    return MotionVector(found_h, found_w);
}

void MotionEstimator::Estimate(py::array_t<unsigned char> _current_frame,
                               py::array_t<unsigned char> _previous_frame) 
{   
    // Since we were storing all of the previous Motion Vectors, we have to clear 'em.
    this -> storage.clear();

    // For every block in current_frame we have to find corresponding (the closest)
    // block in the previous_frame
    unsigned char* previous_frame_ptr = static_cast<unsigned char*>(_previous_frame.request().ptr);
    unsigned char* current_frame_ptr = static_cast<unsigned char*>(_current_frame.request().ptr);
    
    Matrix previous_frame = Matrix(previous_frame_ptr, this -> _height, this -> _width);
    Matrix current_frame = Matrix(current_frame_ptr, this -> _height, this -> _width);
    for (size_t h = 0; h < this -> _height; h += this -> _block_size) {
        for (size_t w = 0; w < this -> _width; w += this -> _block_size) {
            MotionVector move_vector = this -> FindBlock(previous_frame, current_frame, h, w);
            this -> storage.emplace_back(move_vector);
        }
    } 
}

py::array_t<unsigned char> MotionEstimator::Remap(py::array_t<unsigned char> _previous_frame) {
    unsigned char* previous_frame_ptr = static_cast<unsigned char*>(_previous_frame.request().ptr);
    Matrix previous_frame = Matrix(previous_frame_ptr, this -> _height, this -> _width);
    size_t index = 0;

    py::array_t<unsigned char> result(this -> _height * this -> _width);
    unsigned char* result_ptr = static_cast<unsigned char*>(result.request().ptr);
    
    for (size_t h = 0; h < this -> _height; h += this -> _block_size) {
        for (size_t w = 0; w < this -> _width; w += this -> _block_size, index++) {
            AssignBlock(result_ptr, h, w, storage[index], previous_frame);
        }
    }
    result.resize({this -> _height, this -> _width});
    return result;
}

void MotionEstimator::AssignBlock(unsigned char* result_ptr, size_t dh, size_t dw,
                                         MotionVector& motion_vector, Matrix& previous_frame) 
{
    for (size_t h = 0; h < this -> _block_size; h++) {
        for (size_t w = 0; w < this -> _block_size; w++) {
            result_ptr[(dh + h) * this -> _width + w + dw] = previous_frame.get(h + motion_vector._h, w + motion_vector._w);
        }
    }
}