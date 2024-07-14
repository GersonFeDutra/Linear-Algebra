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

    virtual double &operator[](std::size_t) = 0;
    virtual double operator[](std::size_t) const = 0;

    virtual std::size_t size() const = 0;

    friend std::ostream& operator<<(std::ostream& os, const Vector& v);
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

    double &operator[](std::size_t i) override { return data[i]; }
    double operator[](std::size_t i) const override { return data[i]; }

    inline double dot(const Vector2 &v) const
    {
        return _x * v._x + _y * v._y;
    }

    inline double norm() const {
        return sqrt(_x * _x + _y * _y);
    }
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
