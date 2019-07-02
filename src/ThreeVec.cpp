/*
  Three-vector class including mathematical operations and IO
*/

#include "ThreeVec.h"

ThreeVec::ThreeVec() : ThreeVec(0.0, 0.0, 0.0)
{
}

// Cartesian constructor
ThreeVec::ThreeVec(double x, double y, double z)
{
    coord_[0] = x;
    coord_[1] = y;
    coord_[2] = z;
}

// Access function for x coordinate
double ThreeVec::get_x() { return coord_[0]; }

// Access function for y coordinate
double ThreeVec::get_y() { return coord_[1]; }

// Access function for z coordinate
double ThreeVec::get_z() { return coord_[2]; }

// Access function for ith coordinate
double ThreeVec::get(int i) { return coord_[i]; }

// Modifier method for x coordinate
void ThreeVec::set_x(double value) { coord_[0] = value; }

// Modifier method for y coordinate
void ThreeVec::set_y(double value) { coord_[1] = value; }

// Modifier method for z coordinate
void ThreeVec::set_z(double value) { coord_[2] = value; }

// Modifier method for ith coordinate -> INSERT
void ThreeVec::set(int i, double value) { coord_[i] = value; }

void ThreeVec::set_all(double x, double y, double z)
{
    coord_[0] = x;
    coord_[1] = y;
    coord_[2] = z;
}

// Alternative modifier method for ith coordinate -> ADD
void ThreeVec::inc(int i, double value) { coord_[i] += value; }

// Square the threevector
double ThreeVec::square()
{
    double answer = 0.0;
    for (int i = 0; i < 3; ++i)
    {
        answer += coord_[i] * coord_[i];
    }
    return answer;
}

// Magnitude of the threevector
double ThreeVec::mag() { return sqrt(square()); }

/*
Overload the operators +,-,* and ^ to represent vector operations
*/
// Addition
ThreeVec ThreeVec::operator+(ThreeVec vec)
{
    ThreeVec ans(vec.get_x() + coord_[0],
                 vec.get_y() + coord_[1],
                 vec.get_z() + coord_[2]);
    return ans;
}

// Subtraction
ThreeVec ThreeVec::operator-(ThreeVec vec)
{
    ThreeVec ans(coord_[0] - vec.get_x(),
                 coord_[1] - vec.get_y(),
                 coord_[2] - vec.get_z());
    return ans;
}

// Increment
ThreeVec ThreeVec::operator+=(ThreeVec vec)
{
    coord_[0] += vec.get_x();
    coord_[1] += vec.get_y();
    coord_[2] += vec.get_z();
    return *(this); //not necessary but returns modified ThreeVec
}

// Decrement
ThreeVec ThreeVec::operator-=(ThreeVec vec)
{
    coord_[0] -= vec.get_x();
    coord_[1] -= vec.get_y();
    coord_[2] -= vec.get_z();
    return *(this); //not necessary but returns modified ThreeVec
}

// Scalar multiplication
ThreeVec ThreeVec::operator*(double value)
{
    ThreeVec ans(coord_[0] * value,
                 coord_[1] * value,
                 coord_[2] * value);
    return ans;
}

// Scalar division
ThreeVec ThreeVec::operator/(double value)
{
    ThreeVec ans(coord_[0] / value,
                 coord_[1] / value,
                 coord_[2] / value);
    return ans;
}

// Vector multiplication - dot product
double ThreeVec::operator*(ThreeVec vec)
{
    double ans = 0.0;
    for (int i = 0; i < 3; ++i)
    {
        ans += coord_[i] * vec.get(i);
    }
    return ans;
}

// Vector multiplication - cross product
ThreeVec ThreeVec::operator^(ThreeVec vec)
{
    ThreeVec ans(coord_[1] * vec.get_z() - coord_[2] * vec.get_y(),
                 coord_[2] * vec.get_x() - coord_[0] * vec.get_z(),
                 coord_[0] * vec.get_y() - coord_[1] * vec.get_x());
    return ans;
}

// Friend function to e.g. print to screen via cout
// std::ostream &operator<<(std::ostream &outstream, const ThreeVec &v)
// {
//     for (int i = 0; i < 3; ++i)
//     {
//         outstream << v.coord_[i] << '\t';
//     }
//     return outstream;
// }
