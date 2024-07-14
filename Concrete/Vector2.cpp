#include "Vector.hpp"

namespace algebra {


Vector2& Vector2::normalized()
{
    double _norm = norm();
    _x /= _norm;
    _y /= _norm;

    return *this;
}


} // namespace algebra
