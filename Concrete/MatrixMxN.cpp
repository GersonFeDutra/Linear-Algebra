#include "Matrix.hpp"
#include <algorithm>
#include <iostream>

namespace algebra {


MatrixMxN::MatrixMxN(std::size_t m, std::size_t n)
    : m{m}, n{n}, data{new double[m * n]}
{
    for (std::size_t i = 0; i < m; i++)
        for (std::size_t j = 0; j < n; j++)
            *(data + i * n + j) = 0;
}


MatrixMxN::MatrixMxN(std::size_t m, std::size_t n,
        double (* init)(std::size_t, std::size_t))
    : m{m}, n{n}, data{new double[m * n]}
{
    double *p = data;
    for (std::size_t i = 0; i < m; i++)
        for (std::size_t j = 0; j < n; j++)
            *(p++) = init(i, j);
}


MatrixMxN::MatrixMxN(std::initializer_list<std::initializer_list<double>> il)
    : m{il.size()}, n{il.begin()->size()}
{
    data = new double[m * n];
    double *p = data;

    for (auto l : il) {
        assert(l.size() == n);
        for (auto e : l)
            *(p++) = e;
    }
}


MatrixMxN::MatrixMxN(std::initializer_list<std::reference_wrapper<Vector>> il)
    : m{il.begin()->get().size()}, n{il.size()}
{
    data = new double[m * n];
    std::size_t j = 0;

    for (auto vr : il) {
        Vector &v = vr.get();
        assert(v.size() == m);

        for (std::size_t i = 0; i < v.size(); i++)
            *(data + i * n + j) = v[i];

        j++;
    }
}


MatrixMxN::MatrixMxN(const MatrixMxN&copy)
    : m{copy.m}, n{copy.n}, data{new double[copy.m * copy.n]}
{
    for (std::size_t i = 0; i < m; i++)
        for (std::size_t j = 0; j < n; j++)
            get(i, j) = copy.get(i, j);
}


MatrixMxN::MatrixMxN(MatrixMxN&& move)
    : m{move.m}, n{move.n}, data{move.data}
{
    move.n = move.m = 0;
    move.data = nullptr;
}


MatrixMxN& MatrixMxN::operator=(const MatrixMxN&copy)
{
    if (copy.m * copy.n != m * n) {
        delete[] data;
        m = copy.m;
        n = copy.n;
        data = new double[m * n];
    }

    for (std::size_t i = 0; i < copy.m; i++)
        for (std::size_t j = 0; j < copy.n; j++)
            get(i, j) = copy.get(i, j);

    return *this;
}


MatrixMxN& MatrixMxN::operator=(MatrixMxN&&move)
{
    m = move.m;
    n = move.n;
    data = move.data;
    move.m = move.n = 0;
    move.data = nullptr;

    return *this;
}


MatrixMxN& MatrixMxN::transpose()
{
    // TODO -> Parallelism
    if (m == n) {
        for (std::size_t i = 0; i < m; i++)
            for (std::size_t j = i + 1; j < n; j++)
                std::swap(get(j, i), get(i, j));
    }
    else {
        double *new_data = new double[m * n];

        for(std::size_t k = 0; k < m * n; k++) {
            std::size_t i = k / m;
            std::size_t j = k % m;

            new_data[k] = data[n * j + i];
        }

        delete[] data;
        data = new_data;
    }
    std::swap(m, n);

    return *this;
}


MatrixMxN transposed(const Matrix& M)
{
    MatrixMxN N{M.get_n(), M.get_m()};

    for(std::size_t k = 0; k < M.get_m() * M.get_n(); k++) {
        std::size_t i = k / M.get_m();
        std::size_t j = k % M.get_m();
        N.get(i, j) = M.get(j, i);
    }

    return N;
}


MatrixMxN operator+(const Matrix& M, const Matrix& N)
{
    if (M.get_m() != N.get_m() || M.get_n() != N.get_n())
        throw std::range_error("Matrices dimensions differ");

    std::size_t m = M.get_m();
    std::size_t n = M.get_n();
    MatrixMxN result{m, n};

    for (std::size_t i = 0; i < m; i++)
        for (std::size_t j = 0; j < n; j++)
            result.get(i, j) = M.get(i, j) + N.get(i, j);

    return result;
}


MatrixMxN operator-(const Matrix& M, const Matrix& N)
{
    if (M.get_m() != N.get_m() || M.get_n() != N.get_n())
        throw std::range_error("Matrices dimensions differ");

    std::size_t m = M.get_m();
    std::size_t n = M.get_n();
    MatrixMxN result{m, n};

    for (std::size_t i = 0; i < m; i++)
        for (std::size_t j = 0; j < n; j++)
            result.get(i, j) = M.get(i, j) - N.get(i, j);

    return result;
}


MatrixMxN operator*(const Matrix& M, double scalar)
{
    std::size_t m = M.get_m();
    std::size_t n = M.get_n();
    MatrixMxN result{m, n};

    for (std::size_t i = 0; i < m; i++)
        for (std::size_t j = 0; j < n; j++)
            result.get(i, j) = M.get(i, j) * scalar;

    return result;
}


MatrixMxN operator/(const Matrix& M, double scalar)
{
    std::size_t m = M.get_m();
    std::size_t n = M.get_n();
    MatrixMxN result{m, n};

    for (std::size_t i = 0; i < m; i++)
        for (std::size_t j = 0; j < n; j++)
            result.get(i, j) = M.get(i, j) / scalar;

    return result;
}


MatrixMxN multiply(const Matrix& A, const Matrix& B)
{
    if (A.get_n() != B.get_m())
        throw std::range_error(
            "Matrix dimensions incompatible for multiplication.");

    std::size_t n = A.get_m();
    std::size_t m = A.get_n(); // B.get_m()
    std::size_t p = B.get_n();
    MatrixMxN C{n, p};

    // TODO -> Add parallelism
    // TODO -> Use a more efficient algorithm
    for (std::size_t i = 0; i < n; i++)
        for (std::size_t j = 0; j < p; j++)
            for (std::size_t k = 0; k < m; k++)
                C.get(i, j) += A.get(i, k) * B.get(k, j);

    return C;
}


} // namespace algebra
