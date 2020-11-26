#include "my_metric.h"

int compute_abs_difference(
    const Matrix& domain, 
    int domain_h,
    int domain_w,
    const Matrix& rank,
    int rank_h,
    int rank_w,
    int block_size, 
    int error
)  {
    // Rank blocks are always under control, so just check if domain 
    // block lays inside the picture.
    if (domain_h < 0 || domain_h + block_size >= domain.getHeight() + 1 || 
        domain_w < 0 || domain_w + block_size >= domain.getWidth() + 1) 
    {
           return std::numeric_limits<int>::max();
    }
    int sum = 0;
    for (int h = 0; h < block_size; h++) {
        for (int w = 0; w < block_size; w++) {
            sum += std::abs(domain.get(h + domain_h, w + domain_w) - rank.get(h + rank_h, w + rank_w)); 
            if (sum >= error) {
                return std::numeric_limits<int>::max();
            }
        }
    }
    return sum;
}