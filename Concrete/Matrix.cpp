#include "Matrix.hpp"

namespace algebra {


Matrix& Matrix::operator=(const Matrix&copy)
{
    if (get_m() != copy.get_m() || get_n() != copy.get_n())
        throw std::range_error("Matrix size differ");

    for (std::size_t i = 0; i < get_m(); i++)
        for (std::size_t j = 0; j < get_n(); j++)
            get(i, j) = copy.get(i, j);

    return *this;
}


Matrix& Matrix::operator+=(const Matrix &b)
{
    Matrix& a = *this;
    if (a.get_m() != b.get_m() || a.get_n() != b.get_n())
        throw std::range_error("Matrix sizes differ");

    for (std::size_t i = 0; i < get_m(); i++)
        for (std::size_t j = 0; j < get_n(); j++)
            a.get(i, j) += b.get(i, j);

    return *this;
}


Matrix& Matrix::operator-=(const Matrix &b)
{
    Matrix& a = *this;
    if (a.get_m() != b.get_m() || a.get_n() != b.get_n())
        throw std::range_error("Matrix sizes differ");

    for (std::size_t i = 0; i < get_m(); i++)
        for (std::size_t j = 0; j < get_n(); j++)
            a.get(i, j) -= b.get(i, j);

    return *this;
}


Matrix& Matrix::operator*=(double scalar)
{
    for (std::size_t i = 0; i < get_m(); i++)
        for (std::size_t j = 0; j < get_n(); j++)
            this->get(i, j) *= scalar;

    return *this;
}


Matrix& Matrix::operator/=(double scalar)
{
	for (std::size_t i = 0; i < get_m(); i++)
		for (std::size_t j = 0; j < get_n(); j++)
			this->get(i, j) /= scalar;

    return *this;
}


std::ostream& operator<<(std::ostream& os, const Matrix& M)
{
    os << "[\n";
    for (std::size_t i = 0; i < M.get_m(); i++) {
        MatrixVector row = M[i];
        os << M[i] << '\n';
    }
    os << ']';

    return os;
}


} // namespace algebra
