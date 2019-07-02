#include "Species.h"

/**********************************************************
CONSTRUCTORS/DESTRUCTORS
***********************************************************/
Species::Species(uint npar, uint nx, double Qpar, double density, std::function<void(Species &, uint)> init_fcn)
{
    this->npar = npar;
    this->parts.reserve(npar);

    this->density = density;
    this->density_arr = std::vector<double>(nx, 0.0);

    this->Qpar = Qpar;

    this->total_KE = 0.0;

    init_species(init_fcn);
}

Species::~Species()
{
}
//-----------------------------------------

void Species::add_particle(double x_pos, double y_pos, double z_pos, double x_mom, double y_mom, double z_mom, double Wpar)
{
    this->add_particle(Particle(x_pos, y_pos, z_pos, x_mom, y_mom, z_mom, Wpar));
}

void Species::add_particle(ThreeVec pos, ThreeVec mom, double Wpar)
{
    this->add_particle(Particle(pos, mom, Wpar));
}

void Species::add_particle(Particle p)
{
    this->parts.push_back(p);
}

//-----------------------------------------

int Species::deposit_charge(const double dx, const double L_sys, const uint Nx)
{
    const double idx = 1.0 / dx;
    unsigned int j;

    // Initialize
    for (auto &d : this->density_arr)
    {
        d = 0.0;
    }

    double weight, par_weight;

    // for (auto &p : this->parts)
    for (uint i = 0; i < this->npar; ++i)
    {
        par_weight = this->parts.at(i).get_weight();
        weight = this->parts.at(i).get_pos().get_x();

        // This is because I have chosen to start my boundary at -dx/2
        if (weight < 0.0)
        {
            weight += L_sys;
        }

        weight *= idx;
        j = (int)(weight); // cast float -> int

        weight -= (double)j; // now weight is difference between xi and xj

        // Increment rho
        this->density_arr.at(j) += (1.0 - weight) * par_weight;
        ++j;
        if (j >= Nx)
        {
            this->density_arr.at(0) += weight * par_weight;
        }
        else
        {
            this->density_arr.at(j) += weight * par_weight;
        }
    }
    return 0;
}

int Species::map_field_to_part(const Field &f, const double dx, const double L_sys, const uint Nx)
{
    const double idx = 1.0 / dx;
    unsigned long j;

    double weight;
    double loc_f_x1;
    double loc_f_x2 = 0.0;
    double loc_f_x3 = 0.0;
    for (auto &p : this->parts)
    {
        weight = p.get_pos().get_x();
        // This is because I have chosen to start my boundary at -dx/2
        if (weight < 0.0) //NOTE: Assumes L_sys min is 0.0
        {
            weight += L_sys;
        }

        weight *= idx;
        j = (int)(weight); // cast float -> int

        weight -= (double)j; // now weight is difference between xi and xj

        if ((j + 1) >= Nx)
        {
            loc_f_x1 = (1.0 - weight) * this->Qpar * f.f1.at(j) + weight * this->Qpar * f.f1.at(0);
        }
        else
        {
            loc_f_x1 = (1.0 - weight) * this->Qpar * f.f1.at(j) + weight * this->Qpar * f.f1.at(j + 1);
        }

        p.set_local_e_field(loc_f_x1, loc_f_x2, loc_f_x3);
    }
    return 0;
}

int Species::push_particles(const unsigned long n_iter, const double L_sys, const double dt, const double dx)
{
    // For total kinetic energy diagnostic
    double KE = 0.0;
    this->total_KE = 0.0;

    for (auto &p : this->parts)
    {
        ThreeVec pos = p.get_pos();
        ThreeVec mom = p.get_mom();

        // For total kinetic energy diagnostic
        KE = mom.mag();

        if (n_iter > 0)
        {
            mom -= (p.get_local_e_field() * dt);
        }
        else
        {
            mom -= (p.get_local_e_field() * (dt / 2.0));
        }
        p.set_mom(mom);

        // For total kinetic energy diagnostic
        KE *= mom.mag();
        this->total_KE += KE;

        pos += (mom * dt);
        this->apply_bc(pos, L_sys, dx);
        p.set_pos(pos);
    }

    // For total kinetic energy diagnostic
    this->total_KE *= (L_sys / (2.0 * double(this->npar)));

    return 0;
}

int Species::boris_push(const unsigned long n_iter, const double L_sys, const double dt, const double dx)
{
    // For total kinetic energy diagnostic
    double KE = 0.0;
    this->total_KE = 0.0;

    double B0 = sqrt(3.0);//1.0;
    ThreeVec t = ThreeVec(0.0, 0.0, B0 * dt / 2);
    ThreeVec s = t * (2 / (1 + t.square()));

    ThreeVec vstar, vperp;

    for (auto &p : this->parts)
    {
        ThreeVec pos = p.get_pos();
        ThreeVec mom = p.get_mom();

        // For total kinetic energy diagnostic
        KE = mom.mag();

        vperp = mom - (p.get_local_e_field() * (dt / 2.0));
        vstar = vperp + (vperp^t);
        vperp = vperp + (vstar^s);
        mom = vperp - (p.get_local_e_field() * (dt / 2.0));

        p.set_mom(mom);

        // For total kinetic energy diagnostic
        KE *= mom.mag();
        this->total_KE += KE;

        pos += (mom * dt);
        this->apply_bc(pos, L_sys, dx);
        p.set_pos(pos);
    }

    // For total kinetic energy diagnostic
    this->total_KE *= (L_sys / (2.0 * double(this->npar)));

    return 0;
}

std::vector<double> Species::get_x_phasespace()
{
    std::vector<double> to_ret = std::vector<double>(this->npar);

    for (uint i = 0; i < this->npar; ++i)
    {
        to_ret[i]  = this->parts[i].get_pos().get_x();
    }

    return to_ret;
}

std::vector<double> Species::get_px_phasespace()
{
    std::vector<double> to_ret = std::vector<double>(this->npar);

    for (uint i = 0; i < this->npar; ++i)
    {
        to_ret[i] = this->parts[i].get_mom().get_x();
    }

    return to_ret;
}

std::vector<double> Species::get_py_phasespace()
{
    std::vector<double> to_ret = std::vector<double>(this->npar);

    for (uint i = 0; i < this->npar; ++i)
    {
        to_ret[i] = this->parts[i].get_mom().get_y();
    }

    return to_ret;
}

void Species::print_density()
{
    for (auto &d : this->density_arr)
    {
        std::cout << d << '\t';
    }
    std::cout << std::endl;
}

/**********************************************************
PRIVATE FUNCTIONS
***********************************************************/
void Species::init_species(std::function<void(Species &, uint)> init_fcn)
{
    init_fcn(*this, this->npar);
}

void Species::apply_bc(ThreeVec &pos, const double L_sys, const double dx)
{
    double x1 = pos.get_x();
    if (pos.get_x() < -dx / 2.0)
    {
        pos.set_x(x1 + L_sys);
    }
    else if (pos.get_x() >= (L_sys - (dx / 2.0)))
    {
        pos.set_x(x1 - L_sys);
    }
    //TODO: add check to make sure that if the particle is farther out of bounds than L_sys that it actually goes back inside. put in while loop to keep subtracting while it's out of bounds.
}
