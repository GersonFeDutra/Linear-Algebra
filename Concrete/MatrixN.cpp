#include "Matrix.hpp"

namespace algebra {


MatrixN& MatrixN::transpose()
{
    for (std::size_t i = 0; i < get_m(); i++)
        for (std::size_t j = i + 1; j < get_n(); j++)
            std::swap(get(j, i), get(i, j));

    return *this;
}


struct laplace_res {
    double det = 0;
    bool valid = true;
};


static laplace_res laplace_det(const Matrix& M) noexcept
{
    // base cases
    if (M.get_n() == 1) 
        return {M.get(0, 0)};

    if (M.get_n() == 2)
        return {M.get(0, 0) * M.get(1, 1) - M.get(1, 0) * M.get(0, 1)};

    // recursive case, matrix_size > 2
    // expansion indiscriminately along the first row
    //
    //   minor  matrix with i-th row and j-th column
    //          removed for the purpose of calculating
    //          the minor.
    //   r, c   row and column index variables
    //   col    current column in expansion
    //   d      determinant accumulator
    //
    double det = 0;

    // TODO -> Add memory errors throws
    MatrixMxN minor {M.get_n() - 1, M.get_n() - 1};

    // expand along first row
    for(std::size_t col = 0; col < M.get_n(); col++) {
        // copy into minor matrix, left side then right side
        for(std::size_t r = 1; r < M.get_n(); r++) {
            for(std::size_t c = 0; c < col; c++)
                minor.get(r - 1, c) = minor.get(r, c);
            for(std::size_t c = col+1; c < M.get_n(); c++)
                minor.get(r - 1, c - 1) = minor.get(r, c);
        }

        // calculate minor
        laplace_res tmp = laplace_det(minor);
        if(!tmp.valid)
            return {0, tmp.valid};

        det += (col & 1 ? -1 : 1) * M.get(0, col) * tmp.det;
    }

    return {det, true};
}


double MatrixN::determinant() const
{
    // TODO -> Use a more efficient algorithm
    return laplace_det(*this).det;
}


} // namespace algebra
