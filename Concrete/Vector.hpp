#pragma once
#include <initializer_list>
#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <vector>
#include <iostream>

namespace algebra {


class Vector {
public:
    Vector() {}
    Vector(const Vector&copy) = delete;
    virtual ~Vector() {}

    Vector& operator=(const Vector&copy);
    Vector& operator+=(const Vector&);
    Vector& operator-=(const Vector&);
    Vector& operator*=(double);
    Vector& operator/=(double);
	Vector& operator+() { return *this; };
	Vector& operator-();

    virtual double &operator[](std::size_t) = 0;
    virtual const double& operator[](std::size_t) const = 0;

    virtual std::size_t size() const = 0;

    double* begin() { return {&(*this)[0]}; }
    double* end() { return {&(*this)[0] + size()}; }

    const double* begin() const { return &((*this)[0]); }
    const double* end() const { return &((*this)[0]) + size(); }

    friend std::ostream& operator<<(std::ostream& os, const Vector& v);
};


class VectorN : public Vector
{
private:
    std::size_t n = 0;
    double *data = nullptr;
public:
    VectorN(std::size_t n) : n{n}, data{new double[n]}
    {
        for (std::size_t i = 0; i < n; i++)
            data[i] = 0;
    }
    VectorN(std::initializer_list<double> il) : n{il.size()}, data{new double[il.size()]}
    {
        std::size_t i = 0;
        for (double el : il)
            data[i++] = el;
    }
    VectorN(std::size_t n, double (* f)(std::size_t i)) : n{n}, data{new double[n]}
    {
        for (std::size_t i = 0; i < n; i++)
            data[i] = f(i);
    }

    ~VectorN() override
    {
        delete[] data;
        data = nullptr;
    }

    double &operator[](std::size_t i) override { return *(data + i); }
    const double& operator[](std::size_t i) const override { return *(data + i); }

    std::size_t size() const override { return n; }
};


class Vector2 : public Vector {
public:
    Vector2() {
        data[1] = data[0] = 0.0;
    }
    ~Vector2() override = default;

    Vector2(double x, double y) : data{{x}, {y}} {}
    Vector2(const double v[2]) : _x{v[0]}, _y{v[1]} {}

    Vector2(const Vector2&copy) : data{{copy._x}, {copy._y}} {}
    Vector2(Vector2&&move) = default;
    Vector2 &operator=(const Vector2&copy)
    {
        _x = copy._x;
        _y = copy._y;
        return *this;
    }
    Vector2 &operator=(Vector2&&move) = default;

    constexpr std::size_t size() const { return 2; }

    double& operator[](std::size_t i) override { return data[i]; }
    const double& operator[](std::size_t i) const override { return data[i]; }

    inline double dot(const Vector2 &v) const { return _x * v._x + _y * v._y; }
    inline double norm() const { return sqrt(_x * _x + _y * _y); }

    Vector2& normalized();

    inline double angle() const { return atan2(_y, _x); }
    inline double angle_to(const Vector2 &v) const
    {
        return acos(dot(v) / (norm() * v.norm()));
    }

    friend std::ostream& operator<<(std::ostream& os, const Vector2& v)
    {
        return os << '(' << v._x << ", " << v._y << ')';
    }

    inline double& x() { return _x; }
    inline double& y() { return _y; }
    inline double x() const { return _x; }
    inline double y() const { return _y; }
private:
    union {
        double data[2];
        struct {
            double _x;
            double _y;
        };
    };
};


inline Vector2 operator+(const Vector2& u, const Vector2& v)
{
    return {u.x() + v.x(), u.y() + v.y()};
}
inline Vector2 operator-(const Vector2& u, const Vector2& v)
{
    return {u.x() - v.x(), u.y() - v.y()};
}
inline Vector2 operator*(const Vector2& u, double scalar)
{
    return {u.x() * scalar, u.y() * scalar};
}
inline Vector2 operator/(const Vector2& u, double scalar)
{
    return {u.x() / scalar, u.y() / scalar};
}

inline Vector2 operator*(double scalar, const Vector2& u) { return u * scalar; }
inline Vector2 operator/(double scalar, const Vector2& u) { return u / scalar; }


} // namespace algebra
