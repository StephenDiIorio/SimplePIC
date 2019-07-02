#include "Simulation.h"
#include "Species.h"
#include "Field.h"

const uint ndump = 1;
const uint Nx = 256;           // number of grid points
const double L_sys = 8.0 * M_PI * sqrt(3.0);//2.0 * M_PI; // system length
const double tmax = 10.0 * M_PI;
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

void two_stream(Species &spec, uint Npar)
{
    const double Wpar = spec.Qpar / (double(Npar) / double(Nx)); //double(spec.ppc); // weighting (should be / NPPC)
    double x_pos, x_pert, x_tot, x_mom;
    double pert = 0.001;
    double k = 1.0 / (2.0 * sqrt(3.0));//2.0 * 2.0 * M_PI / L_sys;
    double dx = L_sys / double(Npar);

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

        if (i % 2 == 0)
        {
            x_mom = 3.0;
        }
        else
        {
            x_mom = -3.0;
        }

        spec.add_particle(x_tot, 0.0, 0.0, x_mom, 0.0, 0.0, Wpar);
    }
    //TODO: call species BC here?
}

void two_stream_plus(Species &spec, uint Npar)
{
    const double Wpar = spec.Qpar / (double(Npar) / double(Nx)); //double(spec.ppc); // weighting (should be / NPPC)
    double x_pos, x_pert, x_tot, x_mom;
    double pert = 0.001;
    double k = 1.0 / (2.0 * sqrt(3.0)); //2.0 * 2.0 * M_PI / L_sys;
    double dx = L_sys / double(Npar);

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

        x_mom = 3.0;

        spec.add_particle(x_tot, 0.0, 0.0, x_mom, 0.0, 0.0, Wpar);
    }
    //TODO: call species BC here?
}

void two_stream_minus(Species &spec, uint Npar)
{
    const double Wpar = spec.Qpar / (double(Npar) / double(Nx)); //double(spec.ppc); // weighting (should be / NPPC)
    double x_pos, x_pert, x_tot, x_mom;
    double pert = 0.001;
    double k = 1.0 / (2.0 * sqrt(3.0)); //2.0 * 2.0 * M_PI / L_sys;
    double dx = L_sys / double(Npar);

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

        x_mom = -3.0;

        spec.add_particle(x_tot, 0.0, 0.0, x_mom, 0.0, 0.0, Wpar);
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
    // Attributes for species
    // const double Qpar = 2.0;   // charge of particle
    // const uint ppc = 16;
    // const uint npar = ppc * Nx;
    // const double density = 2.0;
    // this->spec.push_back(Species(npar, this->Nx, Qpar, density, two_stream));

    // Attributes for species 1
    const double Qpar1 = 1.0; // charge of particle
    const uint ppc1 = 64;
    const uint npar1 = ppc1 * Nx;
    const double density1 = 1.0;
    this->add_species(npar1, Qpar1, density1, two_stream_plus);

    // Attributes for species 2
    const double Qpar2 = 1.0; // charge of particle
    const uint ppc2 = 64;
    const uint npar2 = ppc2 * Nx;
    const double density2 = 1.0;
    this->add_species(npar2, Qpar2, density2, two_stream_minus);

    // Initialize fields
    this->add_e_field(init_e_field);
    this->add_b_field(init_b_field);
}
