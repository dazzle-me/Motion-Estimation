#pragma once

#include <stdlib.h>

struct MotionVector {
    MotionVector(size_t h, size_t w) {
        _h = h;
        _w = w;
    }
    int _h;
    int _w;
};