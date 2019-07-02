#ifndef SPECIES_H
#define SPECIES_H

#include <iostream>
#include <vector>
#include <functional>

#include "Particle.h"
#include "Field.h"
#include "ThreeVec.h"

class Species
{
    private:
        std::vector<Particle> parts;

        void init_species(std::function<void(Species &, uint)> init_fcn);
        void apply_bc(ThreeVec &pos, const double L_sys, const double dx);

      public:
        // uint ppc;
        uint npar;

        double density;
        std::vector<double> density_arr;

        double Qpar;

        // Diagnostics
        double total_KE;

        /**********************************************************
        CONSTRUCTORS/DESTRUCTORS
        ***********************************************************/
        Species(uint ppc, uint nx, double Qpar, double density, std::function<void(Species &, uint)> init_fcn);

        ~Species();
        //-----------------------------------------

        void add_particle(double x_pos, double y_pos, double z_pos, double x_mom, double y_mom, double z_mom, double Wpar);
        void add_particle(ThreeVec pos, ThreeVec mom, double Wpar);
        void add_particle(Particle p);

        int deposit_charge(const double dx, const double L_sys, const uint Nx);
        int map_field_to_part(const Field &f, const double dx, const double L_sys, const uint Nx);
        int push_particles(const unsigned long n_iter, const double L_sys, const double dt, const double dx);
        int boris_push(const unsigned long n_iter, const double L_sys, const double dt, const double dx);

        std::vector<double> get_x_phasespace();
        std::vector<double> get_px_phasespace();
        std::vector<double> get_py_phasespace();
        void print_density();
};

#endif
