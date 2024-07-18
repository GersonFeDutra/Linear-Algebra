#pragma once
#include <cstdlib>
#include <iostream>
#include <functional>
#include <cassert>

namespace algebra {


template<typename T, std::size_t N> class Vector
{
private:
    T data[N];

public: // Special methods
    constexpr std::size_t size() const { return N; }

    T* begin() { return const_cast<T *>(data); }
    T* end() { return const_cast<T *>(data) + size(); }
    const T* begin() const { return const_cast<T *>(data); }
    const T* end() const { return const_cast<T *>(data) + size(); }

public: // Constructors
    Vector()
    {
        for (T& el : *this)
            el = 0;
    }

    // TODO -> Experiment with type castings
    // TODO -> try expansion expressions
    Vector(const T arr[N])
    {
        for (std::size_t i = 0; i < N; i++)
            data[i] = arr[i];
    }

    Vector(const std::initializer_list<T>& il)
    {
        assert(il.size() == N);
        T *p = data;
        for (auto& el : il)
            *(p++) = el;
    }

    Vector(const std::function<T(std::size_t)>& f)
    {
        T *p = data;
        for (std::size_t i = 0; i < N; i++)
            *(p++) = f(i);
    }

public: // Operators
    Vector& operator+=(const Vector& other)
    {
        for (std::size_t i = 0; i < N; i++)
            (*this)[i] += other[i];

        return *this;
    }
    Vector& operator-=(const Vector& other)
    {
        for (std::size_t i = 0; i < N; i++)
            (*this)[i] -= other[i];

        return *this;
    }
    Vector& operator*=(T scalar)
    {
        for (T& el : *this)
            el *= scalar;

        return *this;
    }
    Vector& operator/=(T scalar)
    {
        for (T& el : *this)
            el /= scalar;

        return *this;
    }
    Vector& operator+()
    {
        return *this;
    }
    Vector operator-()
    {
        Vector copy{*this};
        for (T& el : copy)
            el *= -1;

        return copy;
    }

    T& operator[](std::size_t i) { return data[i]; }
    const T& operator[](std::size_t i) const { return data[i]; }

    friend std::ostream& operator<<(std::ostream& os, const Vector& v)
    {
        const T *el = v.begin();
        os << '(';
        if (v.size() > 0)
            os << *el;
        for (++el; el < v.end(); el++)
            os << ", " << *el;
        os << ')';
        return os;
    }
};


} // namespace algebra
