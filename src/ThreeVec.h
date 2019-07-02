/*
  Three-vector class including mathematical operations and IO
*/

#ifndef THREEVEC_H
#define THREEVEC_H

#include <cmath>
#include <iostream>

class ThreeVec
{
    private:
        double coord_[3]; // Private data members e.g. x,y,z

    public:
        // Default constructor
        ThreeVec();

        // Cartesian constructor
        ThreeVec(double x, double y, double z);

        // Access function for x coordinate
        double get_x();

        // Access function for y coordinate
        double get_y();

        // Access function for z coordinate
        double get_z();

        // Access function for ith coordinate
        double get(int i);

        // Modifier method for x coordinate
        void set_x(double value);

        // Modifier method for y coordinate
        void set_y(double value);

        // Modifier method for z coordinate
        void set_z(double value);

        // Modifier method for ith coordinate -> INSERT
        void set(int i, double value);

        void set_all(double x, double y, double z);

        // Alternative modifier method for ith coordinate -> ADD
        void inc(int i, double value);

        // Square the threevector
        double square();

        // Magnitude of the threevector
        double mag();

        /*
        Overload the operators +,-,* and ^ to represent vector operations
        */
        // Addition
        ThreeVec operator+(ThreeVec vec);

        // Subtraction
        ThreeVec operator-(ThreeVec vec);

        // Increment
        ThreeVec operator+=(ThreeVec vec);

        // Decrement
        ThreeVec operator-=(ThreeVec vec);

        // Scalar multiplication
        ThreeVec operator*(double value);

        // Scalar division
        ThreeVec operator/(double value);

        // Vector multiplication - dot product
        double operator*(ThreeVec vec);

        // Vector multiplication - cross product
        ThreeVec operator^(ThreeVec vec);

        // friend std::ostream &operator<<(std::ostream &outstream, const ThreeVec &v);
};

#endif
