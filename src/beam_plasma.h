#include "Simulation.h"
#include "Species.h"
#include "Field.h"

const uint ndump = 25;
const uint Nx = 64;           // number of grid points
const double L_sys = 8.0 * M_PI * sqrt(3.0);//2.0 * M_PI; // system length
const double tmax = 30.0 * M_PI;
const double dt = 0.025;

const uint nspec = 2;

// Gaussian random distribution function
double gaussian()
{
    double gaussian_width = 2.0; // +/-2 sigma range
    double x, y;

    do
    {
        x = (2.0 * rand() / RAND_MAX - 1.0) * gaussian_width;
        y = exp(-x * x * 0.5);
    } while (1.0 * rand() / RAND_MAX > y);

    return x;
}

void beam(Species &spec, uint Npar)
{
    const double Wpar = spec.Qpar / (double(Npar) / double(Nx));
    double x_pos, x_mom;
    double dx = L_sys / double(Npar);

    for (uint i = 0; i < Npar; ++i)
    {
        x_pos = ((L_sys / double(Npar)) * double(i));
        if (x_pos < -dx / 2.0)
        {
            x_pos += L_sys;
        }
        else if (x_pos >= (L_sys - (dx / 2.0)))
        {
            x_pos -= L_sys;
        }

        x_mom = 3.0;

        spec.parts.push_back(Particle(x_pos, 0.0, 0.0, x_mom, 0.0, 0.0, Wpar));
    }
    //TODO: call species BC here?
}

void plasma(Species &spec, uint Npar)
{
    const double Wpar = spec.Qpar / (double(Npar) / double(Nx));
    double x_pos, x_mom;
    double dx = L_sys / double(Npar);

    for (uint i = 0; i < Npar; ++i)
    {
        x_pos = ((L_sys / double(Npar)) * double(i));
        if (x_pos < -dx / 2.0)
        {
            x_pos += L_sys;
        }
        else if (x_pos >= (L_sys - (dx / 2.0)))
        {
            x_pos -= L_sys;
        }

        x_mom = 0.1 * gaussian();

        spec.parts.push_back(Particle(x_pos, 0.0, 0.0, x_mom, 0.0, 0.0, Wpar));
    }
    //TODO: call species BC here?
}

void single_plasma(Species &spec, uint Npar)
{
    const double Wpar = spec.Qpar / (double(Npar) / double(Nx));
    double x_pos, x_mom;
    double dx = L_sys / double(Npar);

    for (uint i = 0; i < Npar; ++i)
    {
        x_pos = ((L_sys / double(Npar)) * double(i));
        if (x_pos < -dx / 2.0)
        {
            x_pos += L_sys;
        }
        else if (x_pos >= (L_sys - (dx / 2.0)))
        {
            x_pos -= L_sys;
        }

        x_mom = 0.0;
        if (i % 10 == 0)
        {
            x_mom = 3.0;
        }
        else
        {
            x_mom = 0.1 * gaussian();
        }

        spec.parts.push_back(Particle(x_pos, 0.0, 0.0, x_mom, 0.0, 0.0, Wpar));
    }
    //TODO: call species BC here?
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
    // Attributes for beam
    const double Qpar_b = -0.1;   // charge of particle
    const uint ppc_b = 256;
    const uint npar_b = ppc_b * Nx;
    const double density_b = 1.0;
    this->spec.push_back(Species(npar_b, this->Nx, Qpar_b, density_b, beam));

    // Attributes for plasma
    const double Qpar_p = -1.0;
    const uint ppc_p = 256;
    const uint npar_p = ppc_p * Nx;
    const double density_p = 1.0;
    this->spec.push_back(Species(npar_p, this->Nx, Qpar_p, density_p, plasma));

    // Initialize fields
    this->e_field = Field(this->Nx, init_e_field);
    this->b_field = Field(this->Nx, init_b_field);
}
