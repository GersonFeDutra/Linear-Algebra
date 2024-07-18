#pragma once
#include <array>
#include <cassert>
#include <iostream>
#include <memory>
#include <functional>
#include "vector.hpp"

namespace algebra {


template<typename T, std::size_t N> class MatrixVector
{
// Vector members
private:
    T data[N];

public:
    constexpr std::size_t size() const { return N; }

    T* begin() { return const_cast<T *>(data); }
    T* end() { return const_cast<T *>(data) + size(); }
    const T* begin() const { return const_cast<T *>(data); }
    const T* end() const { return const_cast<T *>(data) + size(); }

public:
    MatrixVector()
    {
        for (T& el : *this)
            el = 0;
    }

    MatrixVector& operator+=(const MatrixVector& other)
    {
        for (std::size_t i = 0; i < N; i++)
            (*this)[i] += other[i];
    }
    MatrixVector& operator-=(const MatrixVector& other)
    {
        for (std::size_t i = 0; i < N; i++)
            (*this)[i] -= other[i];
    }
    MatrixVector& operator*=(T scalar)
    {
        for (T& el : *this)
            el *= scalar;
    }
    MatrixVector& operator/=(T scalar)
    {
        for (T& el : *this)
            el /= scalar;
    }

    T& operator[](std::size_t i) { return data[i]; }
    const T& operator[](std::size_t i) const { return data[i]; }

    friend std::ostream& operator<<(std::ostream& os, const MatrixVector& v)
    {
        const T *el = v.begin();
        os << '(';
        for (++el; el < v.end(); el++)
            os << ", " << el;
        os << ')';

		return os;
    }
// End Vector members
// Matrix Vector members
public:
    friend bool operator==(const MatrixVector& a, const MatrixVector& b) { return a.data == b.data; }
    friend bool operator!=(const MatrixVector& a, const MatrixVector& b) { return !(a.data == b.data); }
    friend bool operator<=(const MatrixVector& a, const MatrixVector& b) { return a.data <= b.data; }
    friend bool operator>=(const MatrixVector& a, const MatrixVector& b) { return a.data >= b.data; }
    friend bool operator<(const MatrixVector& a, const MatrixVector& b) { return a.data < b.data; }
    friend bool operator>(const MatrixVector& a, const MatrixVector& b) { return a.data > b.data; }
// End Matrix Vector members
};


template<typename T, std::size_t N> class MatrixIterator : public MatrixVector<T, N>
{
// Vector members
private:
    T *data;

public:
    MatrixIterator()
    {
        for (T& el : *this)
            el = 0;
    }

public:
	constexpr std::size_t size() { return N; }

public: // Operators
    MatrixIterator& operator+=(const MatrixIterator& other)
    {
        for (std::size_t i = 0; i < N; i++)
            (*this)[i] += other[i];
    }
    MatrixIterator& operator-=(const MatrixIterator& other)
    {
        for (std::size_t i = 0; i < N; i++)
            (*this)[i] -= other[i];
    }
    MatrixIterator& operator*=(T scalar)
    {
        for (T& el : *this)
            el *= scalar;
    }
    MatrixIterator& operator/=(T scalar)
    {
        for (T& el : *this)
            el /= scalar;
    }

    T& operator[](std::size_t i)
	{
		static_assert(i < N, "index out of bounds");
		return data[i];
	}
    const T& operator[](std::size_t i) const
	{
		static_assert(i < N, "index out of bounds");
		return data[i];
	}

    constexpr std::size_t size() const { return N; }

    T* begin() { return const_cast<T *>(data); }
    T* end() { return const_cast<T *>(data) + size(); }
    const T* begin() const { return const_cast<T *>(data); }
    const T* end() const { return const_cast<T *>(data) + size(); }

    friend std::ostream& operator<<(std::ostream& os, const MatrixIterator& v)
    {
        const T *el = v.begin();
        os << '(';
        for (++el; el < v.end(); el++)
            os << ", " << el;
        os << ')';

		return os;
    }
// End Vector members
// Matrix Vector members
public:
    friend bool operator==(const MatrixVector<T, N>& a, const MatrixVector<T, N>& b) { return a.data == b.data; }
    friend bool operator!=(const MatrixVector<T, N>& a, const MatrixVector<T, N>& b) { return !(a.data == b.data); }
    friend bool operator<=(const MatrixVector<T, N>& a, const MatrixVector<T, N>& b) { return a.data <= b.data; }
    friend bool operator>=(const MatrixVector<T, N>& a, const MatrixVector<T, N>& b) { return a.data >= b.data; }
    friend bool operator<(const MatrixVector<T, N>& a, const MatrixVector<T, N>& b) { return a.data < b.data; }
    friend bool operator>(const MatrixVector<T, N>& a, const MatrixVector<T, N>& b) { return a.data > b.data; }
// End Matrix Vector members
// Matrix Iterator members
public:
    // Iterator operators
    const MatrixIterator& operator++() const
    {
        data += N;
        return *this;
    }
    const MatrixIterator operator++(int) const
    {
        MatrixIterator iterator = *this;
        data += N;
        return iterator;
    }
    const MatrixIterator& operator--() const
    {
        data -= N;
        return *this;
    }
    const MatrixIterator operator--(int) const
    {
        MatrixIterator iterator = *this;
        data -= N;
        return iterator;
    }
    const MatrixIterator operator+(std::size_t offset) const
    {
        return {data + (offset * N), N};
    }
    const MatrixIterator operator-(std::size_t offset) const
    {
        return {data - (offset * N), N};
    }
    T& operator*() { return *data; }
    const T& operator*() const { return *data; }
	T *operator->() { return data; }
	const T *operator->() const { return data; }
// End Matrix Iterator members
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
 *
 * Atenção: Os operadores binários +,-,* criam e copiam novos dados
 */
template <class T, std::size_t M, std::size_t N> class Matrix
{
protected: // Attributes
	// TODO -> Store it as an Vector array container itself
	// Use union to store 2 data formats (single-vector) and (raw array)
	// Use casting to access individual vectors
	T matrix[M * N];
public:
	constexpr std::size_t get_m() const { return M; }
	constexpr std::size_t get_n() const { return N; }

public: // Constructors
	Matrix()
	{
		for (T *p = matrix; p < matrix + M * N; p++)
			*p = 0;
	}
	Matrix(std::initializer_list<std::initializer_list<T>> list)
	{
		if (M == 0)
			return;
		assert(list.size() <= M);
		assert(list.begin()->size() <= N);
		assert(list.size() >= M);
		assert(list.begin()->size() >= N);

		T *p = const_cast<T *>(matrix);

		for (auto &row : list) {
			if (row.size() > N)
				throw std::out_of_range{"Incompatible initializer list length."};

			for (auto &el : row)
				*(p++) = el;
		}
	}
	Matrix(std::initializer_list<Vector<T, N>> il)
	{
		T *m = matrix;
		assert(il.size() == M);
		for (const auto& v : il) {
			assert(v.size() == N);

			std::size_t i = 0;
			for (auto el : v)
				*(m + i++ * N) = el; // init by vector-column
			m++;
		}
	}

	Matrix(const std::function<T(std::size_t, std::size_t)>& f) { fill(f); }

public: // Constexpr's
	constexpr bool is_square() const { return M == N; }

public: // Helpers
	Matrix& fill(const std::function<T(std::size_t, std::size_t)>& f)
	{
		for (std::size_t i = 0; i < M; i++)
			for (std::size_t j = 0; j < N; j++)
				get(i, j) = f(i, j);

		return *this;
	}

public: // Accessors
	T& get(std::size_t i, std::size_t j) { return *(matrix + i * N + j); }
	T get(std::size_t i, std::size_t j) const { return *(matrix + i * N + j); }

	T *begin() { return matrix; }
	T *end() { return matrix + M * N; }
	const T *begin() const { return matrix; }
	const T *end() const { return matrix + M * N; }

public: // Features
	Matrix<T, N, M> *transposed() const { return transposed(*this); }

public: // Operators
	T *operator[](const std::size_t i) { return (T[N])(matrix + i * N); }

	Matrix &operator+=(const Matrix<T, M, N> &mat)
	{
		for (std::size_t i = 0; i < M; i++)
			for (std::size_t j = 0; j < N; j++)
				get(i, j) += mat.get(i, j);
		return *this;
	}

	Matrix &operator-=(const Matrix<T, M, N> &mat)
	{
		for (std::size_t i = 0; i < M; i++)
			for (std::size_t j = 0; j < N; j++)
				get(i, j) -= mat.get(i, j);
		return *this;
	}

	Matrix &operator*=(const T scalar)
	{
		std::size_t v = M * N;
		for (T *p = matrix; p < matrix + v; p++)
			*p *= scalar;
		return *this;
	}

	Matrix& operator/=(const T scalar)
	{
		for (T *p = matrix; p < matrix + M * N; p++)
			*p *= scalar;
		return *this;
	}

	Matrix& operator+() { return *this; }
	Matrix operator-()
	{
		Matrix copy{*this};
		for (T *p = copy.matrix; p < p + M * N; p++)
			*p *= -1;
		return copy;
	}

	friend std::ostream &operator<<(std::ostream &out, const Matrix &mat)
	{
		out << '[' << std::endl;
		for (std::size_t i = 0; i < M; i++) {
			out << '[';
			if (N > 0)
				out << mat.get(i, 0);
			for (std::size_t j = 1; j < N; j++)
				out << ", " << mat.get(i, j);
			out << ']' << std::endl;
		}
		out << ']';

		return out;
	}
};

// TODO -> Helper function to create a matrix and deduce M and N.
// template <typename T>
// auto make_matrix(std::initializer_list<std::initializer_list<T>> init)
// {
//     constexpr std::size_t M = init.size();
//     constexpr std::size_t N = init.begin()->size();
//     return Matrix<T, M, N>(init);
// }

// TODO -> Use concepts to avoid type repetitions on specifications
template<typename T, std::size_t M, std::size_t N>
Matrix<T, M, N> operator+(const Matrix<T, M, N>& a, const Matrix<T, M, N>& b)
{
	return {
		[&](std::size_t i, std::size_t j) {
			return a.get(i, j) + b.get(i, j);
		}
	};
}
template<typename T, std::size_t M, std::size_t N>
Matrix<T, M, N> operator-(const Matrix<T, M, N>& a, const Matrix<T, M, N>& b)
{
	return {
		[&](std::size_t i, std::size_t j) {
			return a.get(i, j) - b.get(i, j);
		}
	};
}
template<typename T, std::size_t M, std::size_t N>
Matrix<T, M, N> operator*(const Matrix<T, M, N>& m, T scalar)
{
	return {
		[&](std::size_t i, std::size_t j){
			return m.get(i, j) * scalar;
		}
	};
}
template<typename T, std::size_t M, std::size_t N>
Matrix<T, M, N> operator*(T scalar, const Matrix<T, M, N>& m)
{
	return m * scalar;
}
template<typename T, std::size_t M, std::size_t N>
Matrix<T, M, N> operator/(const Matrix<T, M, N>& m, T scalar)
{
	return {
		[&](std::size_t i, std::size_t j){
			return m.get(i, j) / scalar;
		}
	};
}
template<typename T, std::size_t M, std::size_t N>
Matrix<T, M, N> operator/(T scalar, const Matrix<T, M, N>& m)
{
	return {
		[&](std::size_t i, std::size_t j){
			return scalar / m.get(i, j);
		}
	};
}
template<typename T, std::size_t M, std::size_t N>
Matrix<bool, M, N> operator==(const Matrix<T, M, N>& a, const Matrix<T, M, N>& b)
{
	return {
		[&](std::size_t i, std::size_t j){
			return a.get(i, j) == b.get(i, j);
		}
	};
}
template<typename T, std::size_t M, std::size_t N>
Matrix<bool, M, N> operator!=(const Matrix<T, M, N>& a, const Matrix<T, M, N>& b)
{
	return {
		[&](std::size_t i, std::size_t j){
			return a.get(i, j) != b.get(i, j);
		}
	};
}

template <typename T, std::size_t M, std::size_t N, std::size_t P>
Matrix<T, M, P> multiply(const Matrix<T, M, N>& m, const Matrix<T, N, P>& n)
{
	// DEBUG -> Check rightness
	// TODO -> usar algoritmo otimizado
	return {
		[&m, &n](std::size_t i, std::size_t j) {
			T e = 0;
			for (std::size_t k = 0; k < N; k++)
				e += m.get(i, k) * n.get(k, j);
			return e;
		}
	};
}

template<typename T, std::size_t M, std::size_t N>
Matrix<T, N, M> transposed(const Matrix<T, M, N>& m)
{
	return {
		[&m](std::size_t i, std::size_t j){
			return m.get(M - i - 1, N - j - 1);
		}
	};
}


namespace __helper_methods {

struct __laplace_res {
	double det = 0;
	bool valid = true;
};

#define __get(M, i, j, n) (*(M + i * n + j))

template<typename T>
static __laplace_res __laplace_det(const T *M, std::size_t n) noexcept
{
    // base cases
    if (n == 1) 
        return {*M};

    if (n == 2)
        return {*M * __get(M, 1, 1, n) - __get(M, 1, 0, n) * __get(M, 0, 1, n)};

    // recursive case, matrix_size > 2
    // expansion indiscriminately along the first row
    //
    //   minor  matrix with i-th row and j-th column
    //          removed for the purpose of calculating
    //          the minor.
    //   r, c   row and column index variables
    //   col    current column in expansion
    //   d      determinant accumulator

    T det = 0;
    auto minor = new std::unique_ptr<T[(n - 1) * (n - 1)]>;

    // expand along first row
    for(std::size_t col = 0; col < n; col++) {
        // copy into minor matrix, left side then right side
        for(std::size_t r = 1; r < n; r++) {
            for(std::size_t c = 0; c < col; c++)
				__get(minor, r - 1, c, n - 1) = __get(minor, r, c, n - 1);
            for(std::size_t c = col + 1; c < n; c++)
				__get(minor, r - 1, c - 1, n - 1) = __get(minor, r, c, n - 1);
        }

        // calculate minor
        __laplace_res tmp = __laplace_det(minor);
        if(!tmp.valid)
            return {0, tmp.valid};

        det += (col & 1 ? -1 : 1) * M.get(0, col) * tmp.det;
    }

    return {det, true};
}

} // namespace __helper_methods


template <class T, std::size_t N> class SquareMatrix : public Matrix<T, N, N>
{
public:
	constexpr std::size_t get_m() const { return N; }
	constexpr std::size_t get_n() const { return N; }
	constexpr std::size_t size() const { return N * N; }
	constexpr bool is_square() const { return true; }

public:
	SquareMatrix() : Matrix<T, N, N>() {}
	SquareMatrix(std::initializer_list<std::initializer_list<T>> il) : Matrix<T, N, N>(il) {}
	SquareMatrix(std::initializer_list<Vector<T, N>> il) : Matrix<T, N, N>(il) {}
	SquareMatrix(const std::function<T(std::size_t, std::size_t)>& f) : Matrix<T, N, N>(f) {}

public: // Features
	SquareMatrix& transpose()
	{
		for (std::size_t i = 0; i < N; i++)
			for (std::size_t j = i + 1; j < N; j++)
				std::swap(Matrix<T, N, N>::get(i, j), Matrix<T, N, N>::get(j, i));
		return *this;
	}
	SquareMatrix& invert()
	{
		*this = get_inverse();
		return *this;
	}

	// TODO -> Implement
	T get_determinant() const
	{
		// TODO -> Use a more efficient algorithm
		// FIXME -> use dynamic array here
		//return __helper_methods::__laplace_det(Matrix<T, N, N>::matrix, N).det;
		return 0;
	}

	SquareMatrix<T, N> get_adjoint()
	{
		return {
			[this](std::size_t i, std::size_t j) {
				// TODO -> Implement
				// FIXME
				//return (((i + j) & 0b1) ? -1 : 1) * this->get(j, i);
				return 0;
			}
		};
	}

	class matrix_determinant_error : public std::exception {};

	// TODO -> Usar algoritmo otimizado
	SquareMatrix get_inverse()
	{
		T det = get_determinant();
		if (det == 0)
			throw matrix_determinant_error();

		SquareMatrix adj = get_adjoint();
		adj *= 1 / det;

		return adj;
	}
};

template <typename T, std::size_t M, std::size_t N>
constexpr Matrix<T, M, N> operator*(T scalar, Matrix<T, M, N> &m)
{
	return m * scalar;
}

template <typename T, std::size_t N>
SquareMatrix<T, N> multiply(const SquareMatrix<T, N>& m, const SquareMatrix<T, N>& n)
{
	// DEBUG -> Check rightness
	// TODO -> usar algoritmo otimizado
	return {
		[&m, &n](std::size_t i, std::size_t j) {
			T e = 0;
			for (std::size_t k = 0; k < N; k++)
				e += m.get(i, k) * n.get(k, j);
			return e;
		}
	};
}

template <typename T, std::size_t N, std::size_t P>
SquareMatrix<T, N> multiply(const Matrix<T, N, P>& m, const Matrix<T, P, N>& n)
{
	// DEBUG -> Check rightness
	// TODO -> usar algoritmo otimizado
	return {
		[&m, &n](std::size_t i, std::size_t j) {
			T e = 0;
			for (std::size_t k = 0; k < N; k++)
				e += m.get(i, k) * n.get(k, j);
			return e;
		}
	};
}

template<typename T, std::size_t N>
SquareMatrix<T, N> transposed(const SquareMatrix<T, N>& m)
{
	return {
		[&m](std::size_t i, std::size_t j){
			return m.get(N - i - 1, N - j - 1);
		}
	};
}

template<typename T, std::size_t N>
auto inverted(const SquareMatrix<T, N>& m)
{
	return SquareMatrix<T, N>{
		// TODO -> Implement
	};
}


// TODO -> Use concepts to avoid type repetitions on specifications
template<typename T, std::size_t N>
SquareMatrix<T, N> operator+(const SquareMatrix<T, N>& a, const SquareMatrix<T, N>& b)
{
	return {
		[&](std::size_t i, std::size_t j) {
			return a.get(i, j) + b.get(i, j);
		}
	};
}
template<typename T, std::size_t N>
SquareMatrix<T, N> operator-(const SquareMatrix<T, N>& a, const SquareMatrix<T, N>& b)
{
	return {
		[&](std::size_t i, std::size_t j) {
			return a.get(i, j) - b.get(i, j);
		}
	};
}
template<typename T, std::size_t N>
SquareMatrix<T, N> operator*(const SquareMatrix<T, N>& m, T scalar)
{
	return {
		[&](std::size_t i, std::size_t j){
			return m.get(i, j) * scalar;
		}
	};
}
template<typename T, std::size_t N>
SquareMatrix<T, N> operator*(T scalar, const SquareMatrix<T, N>& m)
{
	return m * scalar;
}
template<typename T, std::size_t N>
SquareMatrix<T, N> operator/(const SquareMatrix<T, N>& m, T scalar)
{
	return {
		[&](std::size_t i, std::size_t j){
			return m.get(i, j) / scalar;
		}
	};
}
template<typename T, std::size_t N>
SquareMatrix<T, N> operator/(T scalar, const SquareMatrix<T, N>& m)
{
	return {
		[&](std::size_t i, std::size_t j){
			return scalar / m.get(i, j);
		}
	};
}
template<typename T, std::size_t N>
SquareMatrix<bool, N> operator==(const SquareMatrix<T, N>& a, const SquareMatrix<T, N>& b)
{
	return {
		[&](std::size_t i, std::size_t j){
			return a.get(i, j) == b.get(i, j);
		}
	};
}
template<typename T, std::size_t N>
SquareMatrix<bool, N> operator!=(const SquareMatrix<T, N>& a, const SquareMatrix<T, N>& b)
{
	return {
		[&](std::size_t i, std::size_t j){
			return a.get(i, j) != b.get(i, j);
		}
	};
}


template <typename T> class Matrix2 : public SquareMatrix<T, 2>
{
public:
	Matrix2() : SquareMatrix<T, 2>::matrix{ {0, 0}, {0, 0} } {}
	Matrix2(std::initializer_list<std::initializer_list<T>> list)
		: SquareMatrix<T, 2>(list) {}
	Matrix2(std::initializer_list<Vector<T, 2>> il)
		: SquareMatrix<T, 2>(il) {}
	Matrix2(const std::function<T(std::size_t, std::size_t)>& f)
		: SquareMatrix<T, 2>(f) {}
public:
	Matrix2& transpose()
	{
		std::swap(Matrix<T, 2, 2>::get(0, 1), Matrix<T, 2, 2>::get(1, 0));
	}
	Matrix2& invert()
	{
		T scalar = get_determinant();

		if (scalar == 0)
			throw SquareMatrix<T, 2>::matrix_determinant_error();

		scalar = 1 / scalar;
		// swap multiply
		T tmp = Matrix<T, 2, 2>::get(0, 0) * scalar;
		Matrix<T, 2, 2>::get(0, 0) = Matrix<T, 2, 2>::get(1, 1) * scalar;
		Matrix<T, 2, 2>::get(1, 1) = tmp;

		Matrix<T, 2, 2>::get(0, 1) *= -scalar;
		Matrix<T, 2, 2>::get(1, 0) *= -scalar;
	}
	inline T get_determinant() const
	{
		return Matrix<T, 2, 2>::get(0, 0) * Matrix<T, 2, 2>::get(1, 1) \
			 - Matrix<T, 2, 2>::get(1, 0) * Matrix<T, 2, 2>::get(0, 1);
	}
	inline Matrix2 get_adjoint()
	{
		return {
			{ Matrix<T, 2, 2>::get(1, 1), -Matrix<T, 2, 2>::get(0, 1) },
			{ -Matrix<T, 2, 2>::get(1, 0), Matrix<T, 2, 2>::get(0, 0) },
		};
	}
	Matrix2 get_inverse()
	{
		T det = get_determinant();
		if (det == 0)
			throw SquareMatrix<T, 2>::matrix_determinant_error();
		Matrix2 adj = get_adjoint();
		adj *= 1 / det;

		return adj;
	}

    inline double dot()
	{
		return Matrix<T, 2, 2>::get(0, 0) * Matrix<T, 2, 2>::get(0, 1) \
		 	 + Matrix<T, 2, 2>::get(1, 0) * Matrix<T, 2, 2>::get(1, 1);
	}
};

template <class T> using Matrix3 = SquareMatrix<T, 3>;
template <class T> using Matrix4 = SquareMatrix<T, 4>;

namespace Matrices {
	template <typename T, std::size_t N>
	SquareMatrix<T, N> I{
		[](std::size_t i, std::size_t j){
			return i == j;
		}
	};
}




} // namespace algebra
