#include "my_metric.h"

int compute_abs_difference(const Matrix& a, 
                            int dh_a,
                            int dw_a,
                            const Matrix& b,
                            int dh_b,
                            int dw_b) 
{
    int sum = 0;
    for (int h = 0; h < 16; h++) {
        for (int w = 0; w < 16; w++) {
            sum += std::abs(a.get(h + dh_a, w + dw_a) - b.get(h + dh_b, w + dw_b)); 
        }
    }
    return sum;
}