#pragma once
#include "Vector.hpp"
#include <stdlib.h>
#include <algorithm>
#include <stdexcept>
namespace algebra {

class Matrix {
public:
    virtual ~Matrix() = 0;
    virtual MatrixVector operator[](size_t i) = 0;
    virtual Matrix &transpose() = 0;
    virtual size_t get_M() const = 0;
    virtual size_t get_N() const = 0;

    Matrix& operator+=(const Matrix &b)
    {
        Matrix& a = *this;
        if (a.get_M() != b.get_M() || a.get_N() != b.get_N())
            throw std::range_error("Matrix sizes differ");

        for (size_t i = 0; i < get_M(); i++)
            for (size_t j = 0; j < get_N(); j++)
                a.get(i, j) += b.get(i, j);
    }

    Matrix& operator-=(const Matrix &b)
    {
        Matrix& a = *this;
        if (a.get_M() != b.get_M() || a.get_N() != b.get_N())
            throw std::range_error("Matrix sizes differ");

        for (size_t i = 0; i < get_M(); i++)
            for (size_t j = 0; j < get_N(); j++)
                a.get(i, j) += b.get(i, j);
    }

    Matrix& operator*=(double scalar) {
        for (size_t i = 0; i < get_M(); i++)
            for (size_t j = 0; j < get_N(); j++)
                this->get(i, j) *= scalar;
    }

    Matrix& operator/=(double scalar) {
        for (size_t i = 0; i < get_M(); i++)
            for (size_t j = 0; j < get_N(); j++)
                this->get(i, j) /= scalar;
    }

    virtual double &get(size_t i, size_t j) = 0;
    virtual double &get(size_t i, size_t j) const = 0;
};


class MatrixN : public Matrix {
public:
    virtual double determinant() = 0;

    virtual void invert() = 0;

    MatrixN &transpose() {
        size_t n = get_N();
        for (size_t j = 0; j < n; j++)
            for (size_t i = j + 1; i < n; i++)
                std::swap(this->get(j, i), this->get(i, j));
    }
};

class Matrix2 : public MatrixN {
private:
    double *data;
public:
    Matrix2() : data(new double[4]) {}
    Matrix2(std::initializer_list<std::initializer_list<double>> v) : data(new double[4]) {
        double *p;
        for (auto l : v)
            for (auto e : l)
                *p++ = e;
    }

    ~Matrix2() override { delete[] data; }

    MatrixVector operator[](size_t i) override
    {
        return MatrixVector(data + (2 * i), 2);
    }

    size_t get_M() const override { return 2; }
    size_t get_N() const override { return 2; }

    double &get(size_t i, size_t j) override { return *(data + (2 * i) + j); }
    double &get(size_t i, size_t j) const override { return *(data + (2 * i) + j); }

    double determinant() override
    {
        return this->get(0, 0) * this->get(1, 1) - this->get(0, 1) * this->get(1, 0);
    }

    virtual void invert() override
    {
        std::swap(this->get(0, 0), this->get(1, 1));
        this->get(0, 1) *= -1;
        this->get(1, 0) *= -1;
    }

    Matrix2 inverse()
    {
        return {{this->get(1, 1), -this->get(0, 1)}, {-this->get(1, 0), this->get(0, 0)}};
    }

    double dot(const Matrix &b)
    {
        return this->get(0, 0) * this->get(1, 0) + this->get(0, 1) * this->get(1, 0);
    }
};

Matrix2 &operator+(const Matrix2& a, const Matrix2& b)
{
    Matrix2 c{a};
    c += b;
    return c;
}

Matrix2 &operator-(const Matrix2& a, const Matrix2& b)
{
    Matrix2 c{a};
    c += b;
    return c;
}

Matrix2 &operator*(const Matrix2& a, double b)
{
    Matrix2 c{a};
    c *= b;
    return c;
}

Matrix2 &operator/(const Matrix2& a, double b)
{
    Matrix2 c{a};
    c /= b;
    return c;
}

}
