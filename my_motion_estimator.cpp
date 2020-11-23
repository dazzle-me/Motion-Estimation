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

MotionVector MotionEstimator::FindBlock(const Matrix& rank_block,
                                         unsigned char* previous_frame_ptr)  
{
    // Brute-force search
    uint32_t error = std::numeric_limits<uint32_t>::max();

    size_t found_h = 0, found_w = 0;
    for (size_t h = 0; h < this -> _height; h += this -> _block_size) {
        for (size_t w = 0; w < this -> _width; w += this -> _block_size) {
            Matrix domain_block = Matrix(previous_frame_ptr, this -> _height, this -> _width, h, w);
            uint32_t current_error = compute_abs_difference(rank_block, domain_block);
            if (current_error < error) {
                error = current_error;
                found_h = h;
                found_w = w;
            }
        }
    }
    return MotionVector(found_h, found_w);
}

void MotionEstimator::Estimate(py::array_t<unsigned char> current_frame,
                               py::array_t<unsigned char> previous_frame) 
{   
    // Since we were storing all of the previous Motion Vectors, we have to clear 'em.
    this -> storage.clear();

    // For every block in current_frame we have to find corresponding (the closest)
    // block in the previous_frame
    unsigned char* previous_frame_ptr = static_cast<unsigned char*>(previous_frame.request().ptr);
    unsigned char* current_frame_ptr = static_cast<unsigned char*>(current_frame.request().ptr);
    
    for (size_t h = 0; h < this -> _height; h += this -> _block_size) {
        for (size_t w = 0; w < this -> _width; w += this -> _block_size) {
            Matrix rank_block = Matrix(current_frame_ptr, this -> _height, this -> _width, h, w);

            MotionVector move_vector = this -> FindBlock(rank_block, previous_frame_ptr);
            this -> storage.emplace_back(move_vector);
        }
    } 
}

py::array_t<unsigned char> MotionEstimator::Remap(py::array_t<unsigned char> previous_frame) {
    unsigned char* previous_frame_ptr = static_cast<unsigned char*>(previous_frame.request().ptr);
    size_t index = 0;

    py::array_t<unsigned char> result(this -> _height * this -> _width);
    unsigned char* result_ptr = static_cast<unsigned char*>(result.request().ptr);

    for (size_t h = 0; h < this -> _height; h += this -> _block_size) {
        for (size_t w = 0; w < this -> _width; w += this -> _block_size, index++) {
            Matrix domain_block = Matrix(previous_frame_ptr, this -> _height, this -> _width, storage.at(index));
            AssignBlock(result_ptr, domain_block, h, w);
        }
    }
    result.resize({this -> _height, this -> _width});
    return result;
}

void MotionEstimator::AssignBlock(unsigned char* result_ptr, const Matrix& domain_block, 
                                  size_t offset_h, size_t offset_w) 
{
    size_t offset = this -> _width * offset_h + offset_w; 
    for (size_t h = 0; h < 16; h++) {
        for (size_t w = 0; w < 16; w++) {
            result_ptr[offset + h * this -> _width + w] = domain_block.get(h, w);
        }
    }
}