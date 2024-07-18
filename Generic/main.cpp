#include "matrix.hpp"
#include "vector.hpp"
#include <algorithm>
#include <complex>
#include <iostream>
#include <memory>

using namespace algebra;

namespace std {
    template<typename T>
    std::ostream& operator<<(std::ostream& os, const std::complex<T>& c) {
        os << "(" << c.imag() << "i + " << c.real() << ")";
        return os;
    }
}

template <typename T, std::size_t N>
inline void print_sized(const Vector<T, N> &v, const char *name)
{
    std::cout << name << '[' << v.size() << "]: " << v << '\n';
}


template <typename T, std::size_t M, std::size_t N>
inline void print_sized(const Matrix<T, M, N> &m, const char *name)
{
    std::cout << name << '[' << m.get_m() << 'x' << m.get_n() << "]: " << m << '\n';
}

template <typename T, std::size_t N>
inline void print_sized(const SquareMatrix<T, N> &m, const char *name)
{
    std::cout << name << '[' << m.get_n() << "]: " << m << '\n';
}

inline void print(void)
{
	std::cout << '\n';
}
inline void print(const char *s)
{
	std::cout << s << '\n';
}

template <typename T>
inline void print(const T& obj, const char *name)
{
	std::cout << name << ": " << obj << '\n';
}

template<typename T, std::size_t N>
inline void print(const SquareMatrix<T, N>& mat, const char *name)
{
	std::cout << name << ": " << mat
			// FIX -> determinant function
			//<< ", det(" << name << ") = " << mat.get_determinant() << '\n';
			<< '\n';
}

template<typename T>
inline void print(const Matrix2<T>& mat, const char *name)
{
	std::cout << name << ": " << mat << ", det(" << name << ") = " \
			<< mat.get_determinant() << '\n';
}

template <typename T, int M, int N>
inline void print_square(Matrix<T, M, N> mat)
{
	std::cout << mat << " is" << (mat.is_square() ? " " : " not ") << "square\n";
}


template <typename T, std::size_t N>
void use(Vector<T, N> &u, Vector<T, N> &v)
{
	print();

	const char *s;
	try
	{
		s = "u = u + v";
		u += v;
		print(u, s);

		s = "u = u - v";
		u -= v;
		print(u, s);

		s = "v = -2.5v";
		v *= -2.5;
		print(v, s);

		s = "v = v/-2.5";
		v /= -2.5;
		print(v, s);
	}
	catch(const std::exception& e)
	{
		std::cerr << s << ": " << e.what() << '\n';
	}
}

template <typename T, std::size_t M, std::size_t N>
void use(Matrix<T, M, N>& a, Matrix<T, M, N> b)
{
	print();
	print_sized(a, "A");
	print_sized(b, "B");

	const char *s;
	try
	{
		s = "A = A + B";
		a += b;
		print(a, s);

		s = "A = A - B";
		a -= b;
		print(a, s);

		s = "B = -2.5B";
		b *= -2.5;
		print(b, s);

		s = "B = B/-2.5";
		b /= -2.5;
		print(b, s);
	}
	catch (std::range_error err) {
		std::cerr << s << ": " << err.what() << '\n';
	}

	Matrix<T, N, M> a_t = transposed(a);
	print_sized(a_t, "A^t");

	Matrix<T, M, M> m = multiply(b, a_t);
	print_sized(m, "BA^t");

	Matrix<T, N, N> n = multiply(a_t, b);
	print_sized(n, "A^tB");
}


int main()
{
	using namespace std::complex_literals;

	// Initializations
	Vector<int, 2>coord2d;
	coord2d[0] = -4;
	coord2d[1] = 13;
	
	// TODO -> Conversion in-between value types / make the methods coerce to correct type
	auto coords2d_ops{-coord2d};

	// TODO -> support compatible types arithmetics
	double pair[2] = {M_PI, M_E};
	Vector<double, 2> u{pair};
	Vector<double, 2> v{u};
	v[0] = -1;

	Vector<float, 3> w {
		[](std::size_t i){
			return i + 1;
		}
	};
	Vector<std::complex<double>, 3>v3 { 2.0i + 1., 1.0i + 2., 3 };

	// TODO -> Try adding deduction guides.
	SquareMatrix<int, 2> S { coord2d, coords2d_ops };
	SquareMatrix<int, 3> transform {
		{4, 5, 0},
		{-4, -2, 0},
		{8, 7, 13}
	};
	SquareMatrix<float, 3> M {
		[](std::size_t i, std::size_t j){
			return i * 3 + j + 1;
		}
	};
	SquareMatrix<double, 2> O;

	Matrix<int, 2, 3> A;
	A.get(0, 0) = -1;
	A.get(0, 1) = 1;
	A.get(1, 1) = 3;
	A.get(1, 2) = -5;
	Matrix<int, 2, 3> O2x3;

    Matrix<double, 3, 1> Q{ {1}, {2}, {3} };
    Matrix<double, 3, 1> R{ {-3}, {-2}, {} };
    auto& I = Matrices::I<double, 4>;

	// TODO -> Make method for type deduction
	//auto N = make_matrix({{1, 2}, {3, 4}});

	// Prints
	print("Vectors");
	print(coord2d, "2D coords");
	print(pair, "pair");
	print(u, "u");
	print(v, "v");
	print_sized(w, "w");
	print_sized(v3, "complex 3d vector");

	print("\nMatrices");
	print(O, "O2");
	print(S, "S");
	print(transform, "Transform");
	print_sized(M, "M");

	//print(N, "N");

	print(O2x3, "O2x3");
	print_sized(A, "A");
	print_sized(I, "I");

	// Uses
	print("\nUsing Vectors");
	use(u, v);

	print("\nUsing Matrices");
	use(Q, R);

	// Matrix<int, 4, 5> mat;
	print_sized(mat, "MAT");
	mat.fill([](int i, int j) { return 2 * (i + 1) + 3 * (j + 1); });
	print(mat, "MAT");
	
	Matrix<int, 4, 4> mat4;
	mat4.fill([](int i, int j) { return -(i + 1) + 4 * (j + 1); });
	std::cout << mat4 << std::endl;

	Matrix2<double> a = {{2, 3}, {-5, 9}}, b = {{4, -5}, {3, 1}},
	             c = {{6, 2}, {7, -8}};
	print(a, "A");
	print(b, "B");
	print(c, "C");
	print(a + b, "A + B");
	print(a + c, "A + B");
	print(c - a, "A + B");
	print(a * -2.0, "-2A");
	print(-6.0 * b, "-6B");
	print(c / 3.0, "C/3");
	print(1.0 / c, "1/C");
	print(multiply(a, b), "AB");
}
