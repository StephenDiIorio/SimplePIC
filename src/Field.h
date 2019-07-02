#ifndef FIELD_H
#define FIELD_H

#include <iostream>
#include <vector>
#include <functional>

#include "FFT.h"

typedef unsigned int uint;

class Field
{
    private:
        std::vector<double> K2;
        std::vector<double> kappa;

        void init_field(std::function<void(Field &, uint)> init_fcn);

    public:
        uint size;
        std::vector<double> f1;
        std::vector<double> f2;
        std::vector<double> f3;

        double total_U;

        /**********************************************************
        CONSTRUCTORS/DESTRUCTORS
        ***********************************************************/
        Field(); //TODO: see if this can be removed
        Field(uint nx, double dx, std::function<void(Field &, uint)> init_fcn);
        ~Field();
        //-----------------------------------------

        int solve_field(std::vector<double> re, std::vector<double> im);
        void print_field();
};

#endif
