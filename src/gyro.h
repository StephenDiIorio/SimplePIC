#include "Simulation.h"
#include "Species.h"
#include "Field.h"

const uint ndump = 1;
const uint Nx = 64;           // number of grid points
const double L_sys = 2.0 * M_PI; // system length
const double tmax = 20.0 * M_PI;
const double dt = 0.025;


void gyro(Species &spec, uint Npar)
{
    const double Wpar = spec.Qpar / (double(Npar) / double(Nx));
    double x_pos, x_mom;

    for (uint i = 0; i < Npar; ++i)
    {
        x_pos = L_sys / 2.0;

        x_mom = 0.1;

        spec.add_particle(x_pos, 0.0, 0.0, x_mom, 0.0, 0.0, Wpar);
    }
}

void init_e_field(Field &f, uint size)
{
    f.f1 = std::vector<double>(size, 0.0);
    f.f2 = std::vector<double>(size, 0.0);
    f.f3 = std::vector<double>(size, 0.0);
}

void init_b_field(Field &f, uint size)
{
    f.f1 = std::vector<double>(size, 0.0);
    f.f2 = std::vector<double>(size, 0.0);
    f.f3 = std::vector<double>(size, 0.0);
}

void Simulation::init_simulation()
{
    const double Qpar = -1.0;
    const uint ppc = 1;
    const uint npar = 1;//ppc * Nx;
    const double density = 1.0;

    this->add_species(npar, Qpar, density, gyro);

    // Initialize fields
    this->add_e_field(init_e_field);
    this->add_b_field(init_b_field);
}
