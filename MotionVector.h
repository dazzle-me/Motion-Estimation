#pragma once

#include <vector>

class MotionVector {
public:
    MotionVector() = default;

    MotionVector(int h, int w) {
        _h = h;
        _w = w;
        _error = -1;
        _splitted = false;
        shift_dir = 0;
    }
    MotionVector(int h, int w, int error) {
        _h = h;
        _w = w;
        _error = error;
        _splitted = false;
        shift_dir = 0;
    }
    MotionVector(const std::vector<MotionVector>& subvectors, int error = 0) {
        _h = 0;
        _w = 0;
        _splitted = true;
        _subvectors = subvectors;
        shift_dir = 0;
        _error = error;
    }
    MotionVector(int h, int w, int error, int dir) {
        _h = h;
        _w = w;
        _error = error;
        _splitted = false;
        shift_dir = dir;
    }
    int getHeight() {
        return this -> _h;
    }
    int getWidth() {
        return this -> _w;
    }
    int is_splitted() {
        return this -> _splitted;
    }
    int getError() {
        return this -> _error;
    }
    std::vector<MotionVector> getSubvectors() {
        return this -> _subvectors;
    }
    int shift_dir;
    std::vector<MotionVector> _subvectors;
    bool _splitted;
    int _error;
    int _h;
    int _w;
};