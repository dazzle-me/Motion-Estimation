#include "my_metric.h"

uint32_t compute_abs_difference(const Matrix& a, 
                                const Matrix& b) 
{
    uint32_t sum = 0;
    for (size_t h = 0; h < 16; h++) {
        for (size_t w = 0; w < 16; w++) {
            sum += std::abs(a.get(h, w) - b.get(h, w)); 
        }
    }
    return sum;
}