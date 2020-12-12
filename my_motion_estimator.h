#pragma once

#include <limits>
#include <vector>
#include <array>
#include <unordered_map>
#include <stdexcept>

#include "matrix.h"
#include "my_metric.h"
#include "MotionVector.h"

#include <pybind11/functional.h>
#include <pybind11/stl.h> 
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

class MotionEstimator {
public:
    MotionEstimator(
        int width, 
        int height,
        int quality,
        bool use_halfpixel
    );
    ~MotionEstimator();

    void Estimate(
        py::array_t<unsigned char> _previous_frame,
        py::array_t<unsigned char> _current_frame
    );

    MotionVector FindBlock_BruteForce(
        const Matrix& previous_frame, 
        const Matrix& current_frame,
        int dh,
        int dw,
        int shifted_h,
        int shifted_w
    );
    MotionVector FindBlock_CrossSearch(
        const Matrix& previous_frame, 
        const Matrix& current_frame,
        int dh,
        int dw,
        size_t halfside, 
        int shifted_h,
        int shifted_w, 
        int error,
        int block_size
    );
    MotionVector FindBlock_OrthonormalSearch(
        const Matrix& previous_frame, 
        const Matrix& current_frame,
        int dh,
        int dw,
        size_t step_size, 
        int shifted_h,
        int shifted_w, 
        int error, 
        bool is_horizontal
    );
    MotionVector FindBlock_ThreeStepSearch(
        const Matrix& previous_frame,
        const Matrix& current_frame,
        int dh,
        int dw,
        size_t side,
        int shifted_h,
        int shifted_w,
        int error
    );
    MotionVector FindBlock_3DRS(
        const Matrix& previous_frame,
        const Matrix& current_frame,
        int dh,
        int dw,
        int shifted_h,
        int shifted_w,
        int error
    );
    MotionVector FindBlock_DiamondSearch(
        const Matrix& previous_frame,
        const Matrix& current_frame,
        int dh,
        int dw,
        int shifted_h,
        int shifted_w,
        int error,
        int block_size,
        int shift_dir
    );
    MotionVector FindBlock_HexagonSearch(
        const Matrix& previous_frame,
        const Matrix& current_frame,
        int dh,
        int dw,
        int shifted_h,
        int shifted_w,
        int error, 
        int block_size
    );
    py::array_t<unsigned char> Remap(
        py::array_t<unsigned char> _previous_frame
    );
    void AssignBlock(
        unsigned char* result_ptr, 
        int dh, 
        int dw, 
        MotionVector& motion_vector, 
        const Matrix& previous_frame,
        int block_size
    );
    int ComputeSum(const Matrix& frame, int h, int w);
    int ComputeAbsDifference(
        const Matrix& domain, 
        int domain_h,
        int domain_w,
        const Matrix& rank,
        int rank_h,
        int rank_w,
        int block_size = 16, 
        int error = std::numeric_limits<int>::max()
    );
    MotionVector CheckIfStatic(
        const Matrix& previous_frame,
        const Matrix& current_frame,
        int dh,
        int dw,
        int shifted_h,
        int shifted_w,
        int block_size
    );

    int clip(int pos, int total);
    
    MotionVector GetCandidates(
        const Matrix& preivous_frame,
        const Matrix& current_frame,
        int dh,
        int dw
    );
    void ExtendBorders(
        unsigned char* frame,
        unsigned char* new_frame
    );
    void GenerateSubpixelArrays(
        unsigned char* input,
        unsigned char* output_up,
        unsigned char* output_left,
        unsigned char* output_up_left,
        int height,
        int width
    );
    int GetKey(int h, int w) const;
    // python setters
    void set_SearchMethod(py::array_t<int> value);
    void set_CrossSearch_Side(py::array_t<int> value);
    void set_CrossSearch_ErrorThreshold(py::array_t<int> value);
private:
    enum MODE {
        BruteForce = 0,
        CrossSearch,
        OrthonormalSearch,
        _3DRS,
        ThreeStepSearch,
        DiamondSearch,
        HexagonSearch
    };
    std::vector<MotionVector> previous_storage;
    std::vector<MotionVector> current_storage;

    std::vector<int> previous_frame_precomputed, current_frame_precomputed;
    // Global params
    const int _width;
    const int _height;
    const int _quality;
    const bool _use_halfpixel;

    static constexpr int _block_size = 16;

    size_t SEARCH_MODE;

    int _static_threshold;
    int _error_threshold;
    int _stop_threshold;
    int candidate_threshold;
    
    // _use_halfpixel == true
    std::vector<Matrix> frames;
    unsigned char* previous_up;
    unsigned char* previous_up_left;
    unsigned char* previous_left;

    // If we extend borders, we need this
    unsigned char* previous_extended;
    int border_size;
    int new_width;
    int new_height;
    
    // Candidates search
    bool is_first;

    // Brute-force params
    int _brute_force_stride;
    int _brute_force_height;
    int _brute_force_width;

    // Cross-search params
    int _cross_search_error_threshold;
    int _cross_search_split_threshold;
    uint32_t _cross_search_side;

    // Orthonormal search params
    uint32_t _orthonormal_search_step_size;

    // TSS params
    int _three_step_search_side;

    // Diamond-search params
    // This map maps pair of  coordinates ([h, w], [h', w']) -> error.
    // How can we make sure, that the key is unique for some value (h_1, w_1)?
    // Well, it can be just flatten coordinate of given element (h_1 * this -> _width + w_1)
    // This value has to be unique because we have [0, _height *_width) elements in the matrix
    int _iteration_count;
    std::map<std::pair<int, int>, int> _diamond_search_error_map;
    std::map<std::pair<int, int>, int> _diamon_search_error_map_old;
    std::array<std::pair<int, int>, 4> small_diamond;
    std::array<std::pair<int, int>, 9> large_diamond;

    std::array<std::pair<int, int>, 8> small_diamond_shifted;
    std::array<std::pair<int, int>, 9> large_diamond_shifted;
    

    // 3DRS params
    static constexpr std::array<std::pair<int, int>, 2> _3DRS_current_frame_offset{
        {{-_block_size, -_block_size}, {-_block_size, _block_size}}
    };
    static constexpr std::array<std::pair<int, int>, 4> _3DRS_previous_frame_offset{
        {{2 * _block_size, 2 * _block_size}, {2 * _block_size, 2 * _block_size},
         {_block_size, -_block_size},         {_block_size, _block_size}}
    };
    static constexpr std::array<std::pair<int, int>, 9> _3DRS_random_fluctuations{
        {{0, 0}, 
         {0, _block_size}, {_block_size, 0}, 
         {0, -_block_size}, {-_block_size, 0},
         {0, 3 * _block_size}, {0, - 3 * _block_size},
         {2 * _block_size, -2 * _block_size}}
    };
    static const size_t _3DRS_random_fluct_size = 9;
    size_t _3DRS_offset_index;
};

// motion_vector = GetCandidates(frames[shift_dir], current_frame, h, w);
// motion_vector.shift_dir = shift_dir;
// if (motion_vector._error < this -> candidate_threshold) {
//     found_motion_vector = motion_vector;
//     if (motion_vector.getHeight() != h && motion_vector.getWidth() != w) {
//         std::cout << "w : " << w << " h : " << h << " motion_vector.getHeight() " << motion_vector.getHeight() << " motion_vector.getWidth() " << motion_vector.getWidth() << std::endl;
//         std::cout << "Found candidate, on " + std::to_string(h) + " " + std::to_string(w) << std::endl;
//         break;
//     }
// }
// if (this -> SEARCH_MODE == MODE::BruteForce) {
//     motion_vector = this -> FindBlock_BruteForce(frames[shift_dir], current_frame, h, w, h, w);
// } else if (this -> SEARCH_MODE == MODE::CrossSearch) {
//     motion_vector = this -> FindBlock_CrossSearch(frames[shift_dir], current_frame, h, w, this -> _cross_search_side, h, w, std::numeric_limits<int>::max(), this -> _block_size);
// } else if (this -> SEARCH_MODE == MODE::OrthonormalSearch) {
//     motion_vector = this -> FindBlock_OrthonormalSearch(frames[shift_dir], current_frame, h, w, this -> _orthonormal_search_step_size, h, w, std::numeric_limits<int>::max(), true);
// } else if (this -> SEARCH_MODE == MODE::_3DRS) {
//     motion_vector = this -> FindBlock_3DRS(frames[shift_dir], current_frame, h, w, h, w, std::numeric_limits<int>::max());
// } else if (this -> SEARCH_MODE == MODE::ThreeStepSearch) {
//     motion_vector = this -> FindBlock_ThreeStepSearch(frames[shift_dir], current_frame, h, w, this -> _three_step_search_side, h, w, std::numeric_limits<int>::max());
// } else if (this -> SEARCH_MODE == MODE::DiamondSearch) {
//     motion_vector = this -> FindBlock_DiamondSearch(frames[shift_dir], current_frame, h, w, h, w, std::numeric_limits<int>::max(), this -> _block_size);
// } else if(this -> SEARCH_MODE == MODE::HexagonSearch) {
//     motion_vector = this -> FindBlock_HexagonSearch(frames[shift_dir], current_frame, h, w, h, w, std::numeric_limits<int>::max(), this -> _block_size);
// }