#pragma once

#include <limits>
#include <vector>
#include <array>
#include <unordered_map>

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
    py::list Estimate(
        py::array_t<unsigned char> _previous_frame,
        py::array_t<unsigned char> _current_frame
    );

    MotionVector FindBlock_BruteForce(
        const Matrix& previous_frame, 
        const Matrix& current_frame,
        int dh,
        int dw
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
        int block_size
    );
    py::array_t<unsigned char> Remap(py::array_t<unsigned char> _previous_frame);
    void AssignBlock(unsigned char* result_ptr, size_t dh, size_t dw, MotionVector& motion_vector, Matrix& previous_frame, int block_size);
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
        DiamondSearch
    };
    std::vector<MotionVector> previous_storage;
    std::vector<MotionVector> current_storage;

    std::vector<int> previous_frame_precomputed, current_frame_precomputed;
    // @params
    const int _width;
    const int _height;
    const int _quality;
    const bool _use_halfpixel;

    int border_size;
    int new_width;
    int new_height;
    unsigned char* previous_frame_with_borders_ptr;

    int _static_threshold;

    static constexpr int _block_size = 16;
    size_t SEARCH_MODE;
    

    int _brute_force_stride;
    int _brute_force_height;
    int _brute_force_width;

    int _cross_search_error_threshold;
    int _cross_search_split_threshold;
    uint32_t _cross_search_side;

    uint32_t _orthonormal_search_step_size;

    int _three_step_search_side;

    int _error_threshold;
    int _iteration_count;
    bool is_first;
    // This map maps pair of  coordinates ([h, w], [h', w']) -> error.
    // How can we make sure, that the key is unique for some value (h_1, w_1)?
    // Well, it can be just flatten coordinate of given element (h_1 * this -> _width + w_1)
    // This value has to be unique because we have [0, _height *_width) elements in the matrix
    std::map<std::pair<int, int>, int> _diamond_search_error_map;
    
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