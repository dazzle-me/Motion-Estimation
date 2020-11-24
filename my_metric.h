#pragma once

#include "matrix.h"

// Computes the absolute difference between given slices
int compute_abs_difference(const Matrix& a,
                            int dh_a,
                            int dw_a,
                            const Matrix& b,
                            int dh_b,
                            int dw_b);
