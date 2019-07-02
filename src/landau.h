#include "Simulation.h"
#include "Species.h"
#include "Field.h"

const uint ndump = 25;
const uint Nx = 2048;           // number of grid points
const double L_sys = (1.0/30.0) * M_PI;// / 10.0;//8.0 * M_PI * sqrt(3.0); //M_PI / 5.0; //2.0 * M_PI; // system length
const double tmax = 4.0 * M_PI;
const double dt = 0.01;//0.025;

const uint nspec = 1;

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

void landau(Species &spec, uint Npar)
{
    const double Wpar = spec.Qpar / (double(Npar) / double(Nx)); //double(spec.ppc); // weighting (should be / NPPC)
    double x_pos, x_pert, x_tot, x_mom, vth;
    double pert = 0.05;
    double k = 60; //1.0 / (2.0 * sqrt(3.0));
    double dx = L_sys / double(Npar);

    double wp = 1.0;
    vth = 0.0125;//0.0125;

    for (uint i = 0; i < Npar; ++i)
    {
        x_pos = ((L_sys / double(Npar)) * double(i));
        x_pert = (pert / (spec.Qpar * k)) * sin(k * x_pos);
        x_tot = x_pos + x_pert + 0.01; //TODO: do I need to account for Kaiser-Wilhelm effect here? pg 91 B-L. Dependent on ppc and Nx
        if (x_tot < -dx / 2.0)
        {
            x_tot += L_sys;
        }
        else if (x_tot >= (L_sys - (dx / 2.0)))
        {
            x_tot -= L_sys;
        }

        x_mom = (gaussian() * vth) + ((wp / k) * (pert / spec.Qpar) * cos(k * x_tot));

        spec.add_particle(x_tot, 0.0, 0.0, x_mom, 0.0, 0.0, Wpar);
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
    const uint ppc = 256;
    const uint npar = ppc * Nx;
    const double density = 1.0;
    this->add_species(npar, Qpar, density, landau);

    // Initialize fields
    this->add_e_field(init_e_field);
    this->add_b_field(init_b_field);
}
