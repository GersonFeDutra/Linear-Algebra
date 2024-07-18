#include "Vector.hpp"

namespace algebra {


Vector& Vector::operator=(const Vector&copy)
{
    if (copy.size() != size())
        throw std::range_error("Vectors sizes differ.");

    for (std::size_t i = 0; i < size(); i++)
        (*this)[i] = copy[i];

    return *this;
}


Vector& Vector::operator+=(const Vector& b)
{
    Vector &a = *this;

    if (a.size() != size())
        throw std::range_error("Vectors sizes differ.");

    for (std::size_t i = 0; i < size(); i++)
        a[i] += b[i];

    return a;
}


Vector& Vector::operator-=(const Vector& b)
{
    Vector &a = *this;

    if (a.size() != size())
        throw std::range_error("Vectors sizes differ.");

    for (std::size_t i = 0; i < size(); i++)
        a[i] -= b[i];

    return a;
}


Vector& Vector::operator*=(double scalar)
{
    Vector &a = *this;

    if (a.size() != size())
        throw std::range_error("Vectors sizes differ.");

    for (auto& el : *this)
        el *= scalar;

    return a;
}


Vector& Vector::operator/=(double scalar)
{
    Vector &a = *this;

    if (a.size() != size())
        throw std::range_error("Vectors sizes differ.");

    for (auto& el : *this)
        el /= scalar;

    return a;
}


Vector& Vector::operator-()
{
    for (auto& el : *this)
        el *= -1;

    return *this;
}


std::ostream& operator<<(std::ostream& os, const Vector& v)
{
    os << '(';
    const double *it;
    if (v.size() > 0)
        os << *(it = v.begin());
    for (++it; it < v.end(); it++)
        os << ", " << *it;
    os << ')';

    return os;
}


} // namespace algebra
