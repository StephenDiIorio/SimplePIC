#include "Simulation.h"
#include "Species.h"
#include "Field.h"

const uint ndump = 1;
const uint Nx = 64;           // number of grid points
const double L_sys = 2.0 * M_PI; // system length
const double tmax = 10.5 * M_PI;
const double dt = 0.025;

void uh_2a(Species &spec, uint Npar)
{
    const double Wpar = spec.Qpar / (double(Npar) / double(Nx));
    const double B0 = sqrt(3.0);

    double x_pos, x_pert, x_tot, x_mom;
    double pert = 0.01 * (1 + (B0 * B0));
    double k = 2.0 * M_PI / L_sys;
    double dx = L_sys / double(Npar);

    double wp = 1.0;

    for (uint i = 0; i < Npar; ++i)
    {
        x_pos = ((L_sys / double(Npar)) * double(i));
        x_pert = 0.0;
        x_tot = x_pos + x_pert + 0.01;
        if (x_tot < -dx / 2.0)
        {
            x_tot += L_sys;
        }
        else if (x_tot >= (L_sys - (dx / 2.0)))
        {
            x_tot -= L_sys;
        }

        x_mom = ((wp / k) * (pert / spec.Qpar) * cos(k * x_tot));

        spec.add_particle(x_tot, 0.0, 0.0, x_mom, 0.0, 0.0, Wpar);
    }
}

void uh_2b(Species &spec, uint Npar)
{
    const double Wpar = spec.Qpar / (double(Npar) / double(Nx));

    double x_pos, x_pert, x_tot, x_mom;
    double pert = 0.01;
    double k = 2.0 * M_PI / L_sys;
    double dx = L_sys / double(Npar);

    for (uint i = 0; i < Npar; ++i)
    {
        x_pos = ((L_sys / double(Npar)) * double(i));
        x_pert = (pert / (spec.Qpar * k)) * sin(k * x_pos);
        x_tot = x_pos + x_pert + 0.01;
        if (x_tot < -dx / 2.0)
        {
            x_tot += L_sys;
        }
        else if (x_tot >= (L_sys - (dx / 2.0)))
        {
            x_tot -= L_sys;
        }

        x_mom = 0.0;

        spec.add_particle(x_tot, 0.0, 0.0, x_mom, 0.0, 0.0, Wpar);
    }
}

void uh_2c(Species &spec, uint Npar)
{
    const double Wpar = spec.Qpar / (double(Npar) / double(Nx));
    const double B0 = sqrt(3.0);

    double x_pos, x_pert, x_tot, y_mom;
    double pert = 0.01;
    double k = 2.0 * M_PI / L_sys;
    double dx = L_sys / double(Npar);

    double wp = 1.0;

    for (uint i = 0; i < Npar; ++i)
    {
        x_pos = ((L_sys / double(Npar)) * double(i));
        x_pert = (pert / (spec.Qpar * k)) * cos(k * x_pos);
        x_tot = x_pos + x_pert + 0.01;
        if (x_tot < -dx / 2.0)
        {
            x_tot += L_sys;
        }
        else if (x_tot >= (L_sys - (dx / 2.0)))
        {
            x_tot -= L_sys;
        }

        y_mom = (pert / (spec.Qpar * k * B0)) * cos(k * x_tot);

        spec.add_particle(x_tot, 0.0, 0.0, 0.0, y_mom, 0.0, Wpar);
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
    const uint ppc = 64;
    const uint npar = ppc * Nx;
    const double density = 1.0;

    this->add_species(npar, Qpar, density, uh_2b);

    // Initialize fields
    this->add_e_field(init_e_field);
    this->add_b_field(init_b_field);
}
