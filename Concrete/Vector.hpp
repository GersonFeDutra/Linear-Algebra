#pragma once
#include <initializer_list>
#include <stdlib.h>
#include <stdexcept>
namespace algebra {

class Vector {
public:
    Vector() {}
    Vector(const Vector&copy) = delete;
    Vector operator=(const Vector&copy) = delete;

    virtual double &operator[](size_t) = 0;

    virtual size_t get_Size() = 0;
};

class Vector2 : public Vector {
public:
    Vector2() {
        data[1] = data[0] = 0.0;
    }

    Vector2(double x, double y) : data{{x}, {y}} {}

    Vector2(const Vector2&copy) = default;
    Vector2(const Vector2&&move) = default;
    Vector2 &operator=(const Vector2&copy) = default;
    Vector2 &operator=(Vector2&&move) = default;

    constexpr size_t get_Size() { return 2; }

    double &operator[](size_t i) override { return data[i]; }
private:
    double data[2];
};

class MatrixVector : public Vector {
public:
    MatrixVector(double *p, size_t n) : data(p), size(n) {}

    MatrixVector(const MatrixVector&copy) = default;
    MatrixVector(const MatrixVector&&move) = default;
    MatrixVector &operator=(const MatrixVector&copy) = default;
    MatrixVector &operator=(MatrixVector&&move) = default;

    double &operator[](size_t i) override {
        if (i >= size)
            throw std::range_error("Index out of bounds");

        return data[i];
    }

    size_t get_Size() override { return size; }

private:
    double *const data;
    size_t size;
};

}
