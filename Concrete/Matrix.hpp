#pragma once
#include "Vector.hpp"
#include <algorithm>
#include <cstdlib>
#include <cassert>
#include <stdexcept>
#include <functional>
#include <iostream>

namespace algebra {


class MatrixVector : public Vector {
public:
    // Send nullptr to construct a null-vector slice (for constants only)
    MatrixVector(double *p, std::size_t n) : data(p), _size(n) {}
    ~MatrixVector() override = default;

    MatrixVector(const MatrixVector& copy)
        : data{copy.data}, _size{copy._size} {}
    MatrixVector(MatrixVector&& move) = default;
    MatrixVector& operator=(const MatrixVector& copy) = delete;
    MatrixVector& operator=(MatrixVector&& move) = delete;

    double &operator[](std::size_t i) override {
        if (i >= _size)
            throw std::range_error("Index out of bounds");

        return data[i];
    }
    const double& operator[](std::size_t i) const override {
        if (i >= _size)
            throw std::range_error("Index out of bounds");

        return data[i];
    }

    std::size_t size() const override { return _size; }

    friend bool operator==(const MatrixVector& a, const MatrixVector& b) { return a.data == b.data; }
    friend bool operator!=(const MatrixVector& a, const MatrixVector& b) { return !(a.data == b.data); }
    friend bool operator<=(const MatrixVector& a, const MatrixVector& b) { return a.data <= b.data; }
    friend bool operator>=(const MatrixVector& a, const MatrixVector& b) { return a.data >= b.data; }
    friend bool operator<(const MatrixVector& a, const MatrixVector& b) { return a.data < b.data; }
    friend bool operator>(const MatrixVector& a, const MatrixVector& b) { return a.data > b.data; }
protected:
    mutable double* data;
    std::size_t _size;
};


class MatrixIterator : public MatrixVector {
public:
    MatrixIterator(double *p, std::size_t n) : MatrixVector{p, n} {}

    // Iterator operators
    const MatrixIterator& operator++() const
    {
        data += _size;
        return *this;
    }
    const MatrixIterator operator++(int) const
    {
        MatrixIterator iterator = *this;
        data += _size;
        return iterator;
    }
    const MatrixIterator& operator--() const
    {
        data -= _size;
        return *this;
    }
    const MatrixIterator operator--(int) const
    {
        MatrixIterator iterator = *this;
        data -= _size;
        return iterator;
    }
    const MatrixIterator operator+(std::size_t offset) const
    {
        return {data + (offset * _size), _size};
    }
    const MatrixIterator operator-(std::size_t offset) const
    {
        return {data - (offset * _size), _size};
    }
    double& operator*() { return *data; }
    const double& operator*() const { return *data; }
};

/*
 * M_(m* x n*) = [m_ij];
 * m*: row, n*: col;
 * m: row-1, n: row-1
 * | 0'0 0'1 0'2 ... 0'j ... 0'n | ⇾ row vector
 * | 1'0 1'1 1'2 ... 1'j ... 1'n |
 * | 2'0 2'1 2'2 ... 2'j ... 2'n |
 * |  :   :   :   *   :   *   :  |
 * | i'0 i'1 i'2 ... i'j ... i'n |
 * |  :   :   :   *   :       :  |
 * | m'0 m'1 m'2 ... m'j ... m'n |
 *    ↳ col vector
 * Inicialização numérica por vetor-linhas (Row-major order)
 * {               [
 *   {, , },         * * *
 *   {, , }          * * *
 * }               ]
 * Inicialização vetorial por vetor-colunas
 * { u, v, w }     [ [] [] [] ]
 */
class Matrix {
public:
    using Iterator = MatrixIterator;
public:
    Matrix() {}
    virtual ~Matrix() {}
    Matrix(const Matrix&copy) = delete;

    Matrix& operator=(const Matrix& copy);
    Matrix& operator+=(const Matrix& b);
    Matrix& operator-=(const Matrix& b);
    Matrix& operator*=(double scalar);
    Matrix& operator/=(double scalar);

    virtual MatrixVector operator[](std::size_t i) = 0;
    virtual const MatrixVector operator[](std::size_t i) const = 0;

    virtual Iterator begin() = 0;
    virtual Iterator end() = 0;
    virtual const Iterator begin() const = 0;
    virtual const Iterator end() const = 0;
    std::size_t size() const { return get_m(); }; // Used by iterators

    virtual std::size_t get_m() const = 0;
    virtual std::size_t get_n() const = 0;

    virtual double& get(std::size_t i, std::size_t j) = 0;
    virtual double get(std::size_t i, std::size_t j) const = 0;

    friend std::ostream& operator<<(std::ostream& os, const Matrix& M);
};


class MatrixMxN : public Matrix {
private:
    std::size_t m = 0; // rows
    std::size_t n = 0; // cols
    double *data = nullptr;
public:
    MatrixMxN() {}
    MatrixMxN(std::size_t m, std::size_t n);
    MatrixMxN(std::size_t m, std::size_t n,
        double (* init)(std::size_t, std::size_t));
    MatrixMxN(std::initializer_list<std::initializer_list<double>>);
    MatrixMxN(std::initializer_list<std::reference_wrapper<Vector>>);

    ~MatrixMxN() override {
        delete[] data;
        data = nullptr;
        m = n = 0;
    }

    MatrixMxN(const MatrixMxN& copy);
    MatrixMxN(MatrixMxN&& move);
    MatrixMxN& operator=(const MatrixMxN& copy);
    MatrixMxN& operator=(MatrixMxN&& move);

    MatrixVector operator[](std::size_t i) override { return {data + i * m, n}; }
    const MatrixVector operator[](std::size_t i) const override
    {
        return {data + i * n, m};
    }

    Iterator begin() override { return {data, n}; };
    Iterator end() override { return {data + m * n, n}; };
    const Iterator begin() const override { return {data, n}; };
    const Iterator end() const override { return {data + m * n, n}; };

    MatrixMxN& transpose();

    std::size_t get_m() const override { return m; }
    std::size_t get_n() const override { return n; }

    double& get(std::size_t i, std::size_t j) override
    {
        return *(data + (i * n + j));
    }
    double get(std::size_t i, std::size_t j) const override
    {
        return *(data + (i * n + j));
    }
};

MatrixMxN operator+(const Matrix&, const Matrix&);
MatrixMxN operator-(const Matrix&, const Matrix&);
MatrixMxN operator*(const Matrix&, double);
MatrixMxN operator/(const Matrix&, double);
MatrixMxN multiply(const Matrix&, const Matrix&);
MatrixMxN transposed(const Matrix&);


class MatrixN : public Matrix {
public:
    virtual ~MatrixN() {}

    virtual double determinant() const;
    virtual void invert() = 0;

    virtual MatrixN& transpose();

	class matrix_determinant_error : public std::exception {};
};


class Matrix2 : public MatrixN {
public:
    Matrix2() {}
    Matrix2(std::initializer_list<std::initializer_list<double>> v);
    Matrix2(Vector2& u, Vector2& v) : data{u.x(), v.x(), u.y(), v.y()} {}

    ~Matrix2() override = default;

    Matrix2(const Matrix2& copy)
        : data{copy._x1, copy._y1,
               copy._x2, copy._y2} {}
    Matrix2& operator=(const Matrix2& copy);

    Matrix2(Matrix2&& move) = default;
    Matrix2 &operator=(Matrix2&& move) = default;

    MatrixVector operator[](std::size_t i) override
    {
        return {data + (2 * i), 2};
    }
    const MatrixVector operator[](std::size_t i) const override
    {
        return {const_cast<double *>(data + 2 * i), (std::size_t)2};
    }

    Iterator begin() override { return {data, 2}; };
    Iterator end() override { return {data + 4, 2}; };
    const Iterator begin() const override { return {const_cast<double *>(data), 2}; };
    const Iterator end() const override { return {const_cast<double *>(data) + 4, 2}; };

    constexpr std::size_t get_m() const override { return 2; }
    constexpr std::size_t get_n() const override { return 2; }

    double& get(std::size_t i, std::size_t j) override
    {
        return *(data + (2 * i) + j);
    }
    double get(std::size_t i, std::size_t j) const override
    {
        return *(data + (2 * i) + j);
    }

    double determinant() const override { return _x1 * _y2 - _x2 * _y1; }
    void invert() override;
    inline Matrix2 adjoint() { return {{_y2, -_y1}, {-_x2, _x1}}; }
    inline Matrix2 inverse()
    {
        double det = determinant();
        if (det == 0)
            throw MatrixN::matrix_determinant_error();
        Matrix2 adj = adjoint();
        adj *= 1 / det;
        return adj;
    }
    inline double dot() { return _x1 * _y1 + _x2 * _y2; }

    friend std::ostream& operator<<(std::ostream& os, const Matrix2& v)
    {
        return os << "[\n(" << v._x1 << ", " << v._y1 << ")\n" \
                << v._x2 << ", " << v._y2 << ")\n]";
    }

    inline double& x1() { return _x1; }
    inline double& x2() { return _x2; }
    inline double& y1() { return _y1; }
    inline double& y2() { return _y2; }
    inline double x1() const { return _x1; }
    inline double x2() const { return _x2; }
    inline double y1() const { return _y1; }
    inline double y2() const { return _y2; }
private:
    union {
        double data[4]{0, 0, 0, 0};
        struct {
            double _x1, _y1, _x2, _y2;
        };
    };
};

inline Matrix2 operator+(const Matrix2& a, const Matrix2& b)
{
    return {
        {a.x1() + b.x1(), a.y1() + b.y1()},
        {a.x2() + b.x2(), a.y2() + b.y2()},
    };
}
inline Matrix2 operator-(const Matrix2& a, const Matrix2& b)
{
    return {
        {a.x1() - b.x1(), a.y1() - b.y1()},
        {a.x2() - b.x2(), a.y2() - b.y2()},
    };
}
inline Matrix2 operator*(const Matrix2& M, double scalar)
{
    return {
        {M.x1() * scalar, M.y1() * scalar},
        {M.x2() * scalar, M.y2() * scalar},
    };
}
inline Matrix2 operator/(const Matrix2& M, double scalar)
{
    return {
        {M.x1() / scalar, M.y1() / scalar},
        {M.x2() / scalar, M.y2() / scalar},
    };
}

inline Matrix2 operator*(double scalar, const Matrix2& M) { return M * scalar; }
inline Matrix2 operator/(double scalar, const Matrix2& M) { return M / scalar; }

inline Matrix2 multiply(const Matrix2& A, const Matrix2& B)
{
    return {
        { A.x1() * B.x1() + A.y1() * B.x2(),
          A.x2() * B.x1() + A.y2() * B.x2() },
        { A.x1() * B.y1() + A.y1() * B.y2(),
          A.x2() * B.y1() + A.y2() * B.y2(), }
    };
}


} // namespace algebra
