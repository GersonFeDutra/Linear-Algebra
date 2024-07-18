#include "Matrix.hpp"

namespace algebra {


Matrix2::Matrix2(std::initializer_list<std::initializer_list<double>> v)
{
    double *p = data;
    for (auto l : v)
        for (auto e : l)
            *p++ = e;
}


Matrix2& Matrix2::operator=(const Matrix2&copy)
{
    double *q = const_cast<double *>(copy.data);
    for (double *p = data; p < data + 4; p++)
        *p = *(q++);

    return *this;
}


void Matrix2::invert()
{
    double scalar = determinant();

    if (scalar == 0)
        throw MatrixN::matrix_determinant_error();

    // swap with scalar
    double tmp = _x1 * scalar;
    _x1 = _y2 * scalar;
    _y2 = tmp;

    _x2 *= -scalar;
    _y1 *= -scalar;
}


}
