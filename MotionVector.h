#pragma once

struct MotionVector {
    MotionVector(int h, int w) {
        _h = h;
        _w = w;
        _error = -1;
    }
    MotionVector(int h, int w, int error) {
        _h = h;
        _w = w;
        _error = error;
    }
    int _error;
    int _h;
    int _w;
};