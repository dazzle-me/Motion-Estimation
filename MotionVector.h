#pragma once

#include <vector>

struct MotionVector {
    MotionVector(int h, int w) {
        _h = h;
        _w = w;
        _error = -1;
        _splitted = false;
    }
    MotionVector(int h, int w, int error) {
        _h = h;
        _w = w;
        _error = error;
        _splitted = false;
    }
    MotionVector(const std::vector<MotionVector>& subvectors) {
        _h = 0;
        _w = 0;
        _splitted = true;
        _subvectors = subvectors;
    }
    std::vector<MotionVector> _subvectors;
    bool _splitted;
    int _error;
    int _h;
    int _w;
};