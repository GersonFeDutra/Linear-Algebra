#include <iostream>
#include "Vector.hpp"
#include "Matrix.hpp"
#include <cmath>

#define rad2deg(RAD) (RAD * (180.0 / M_PI))

using namespace algebra;


inline void print_sized(const Vector &v, const char *name)
{
    std::cout << name << '[' << v.size() << "]: " << v << '\n';
}

inline void print_sized(const Matrix &M, const char *name)
{
    std::cout << name << '[' << M.get_m() << 'x' << M.get_n() << ']' \
            << ": " << M << '\n';
}

inline void print(const Vector &v, const char *name)
{
    std::cout << name << ": " << v << '\n';
}

inline void print(const Matrix &M, const char *name)
{
    std::cout << name << ": [\n" << M << '\n';
}

inline void print(const Vector2 v, const char *name)
{
    std::cout << name << ": " << v << '\n';
}

inline void print(const Matrix2& M, const char *name)
{
    std::cout << name << ": " << M << \
            ", det(" << name << ")=" << M.determinant() << "\n";
}


void use(Vector &u, Vector &v, Vector &w)
{
    print_sized(u, "u");
    print_sized(v, "v");
    print_sized(w, "w");

    const char *s;
    try {
        s = "u = u + v";
        u += v;
        print(u, s);

        s = "u = u - v";
        u -= v;
        print(u, s);

        s = "w = -2.5w";
        w *= -2.5;
        print(w, s);

        s = "w = w/-2.5";
        w /= -2.5;
        print(w, s);
    }
    catch (std::range_error err) {
        std::cerr << s << ": " << err.what() << '\n';
    }
}


void use(Matrix &M, Matrix &N, Matrix &P)
{
    print_sized(M, "M");
    print_sized(N, "N");
    print_sized(P, "P");

    const char *s;
    try {
        s = "M = M + N";
        M += N;
        print(M, s);

        s = "M = M - N";
        M -= N;
        print(M, s);

        s = "P = -2.5P";
        P *= -2.5;
        print(P, s);

        s = "P = P/-2.5";
        P /= -2.5;
        print(P, s);
    }
    catch (std::range_error err) {
        std::cerr << s << ": " << err.what() << '\n';
    }
}


void use(Vector2 u, Vector2 v, Vector2 w)
{
    print(u + v, "u + v");
    print(u - v, "u - v");
    print(2 * w, "2 * w");
    print((u + v) / 1.5 - w, "(u + v) / 1.5 - w");

    std::cout << "<u, v>: " << u.dot(v) << '\n';
    std::cout << "angle(u, v): " << u.angle_to(v) << " = " << rad2deg(w.angle_to(v)) << "°\n";
    std::cout << "angle(w): " << w.angle() << " = " << rad2deg(w.angle()) << "°\n";
    std::cout << "||w||: " << w.norm() << '\n';
    u = w.normalized();
    print(u, "w -> u (unitary)");
}


void use(Matrix2 M, Matrix2 N, Matrix2 P)
{
    print(M + N, "M + N");
    print(M - N, "M - N");
    print(2 * P, "2P");
    print((M + N) / 1.5 - P, "(M + N) / 1.5 - P");
    print(M.inverse(), "M^-1");
    Matrix2 P_inv = P;
    P_inv.invert();
    print(P_inv, "P^-1");
    print(multiply(P, P_inv), "PP^-1 = I");
}


constexpr const std::size_t BIG_M = 10, BIG_N = 5;
double f(std::size_t i, std::size_t j)
{
    return i * BIG_N + j;
}


int main()
{
    double arr[] = {1, 2};
    Vector2 u{arr};
    Vector2 v{u};
    v[0] = -1;

    auto w = Vector2{3, 4};
    Vector2 o;

    Matrix2 M{ {1, 0}, {2, 1} };
    Matrix2 N{u, w};
    Matrix2 P{ {cos(M_PI), sin(M_PI)}, {-sin(M_PI), cos(M_PI)} };
    Matrix2 O;

    MatrixMxN O2x3(2, 3);
    MatrixMxN Q{ {1}, {2}, {3} };
    MatrixMxN R{u, v, w};
    MatrixMxN I{ { 1, 0, 0, 0 }, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1} };

    MatrixMxN Big{BIG_M, BIG_N, f};

    std::cout << "Vectors\n";
    print_sized(o, "o");
    use((Vector&)u, (Vector&)v, (Vector&)w);
    use(u, v, w);

    std::cout << "\nMatrices\n";
    print_sized(O, "O");

    use((Matrix&)M, (Matrix&)N, (Matrix&)P);
    use(M, N, P);

    print_sized(O2x3, "O'");
    print_sized(Q, "Q");

    MatrixMxN Q_t = transposed(Q);
    print_sized(Q_t, "Q^t");

    print_sized(R, "R");

    try {
        MatrixMxN res = multiply(R, Q);
        print_sized(res, "RQ");
    }
    catch (std::range_error err) {
        std::cerr << "RQ: " << err.what() << '\n';
    }

    print_sized(R.transpose(), "R^t");

    print_sized(I, "I");
    print(I.transpose(), "I^t == I");

    print_sized(Big, "Big");
    print_sized(Big.transpose(), "Big^t");
}
