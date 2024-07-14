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
    std::swap(_x1, _y2);
    _x2 *= -1;
    _y1 *= -1;
}


}
