#include "Simulation.h"

/**********************************************************
CONSTRUCTORS/DESTRUCTORS
***********************************************************/
Simulation::Simulation(uint ndump, uint nx, double L, double dt, double tmax)
{
    this->err = 0;

    this->n_iter = 0;
    this->ndump = ndump;

    this->Nx = nx;
    this->L_sys = L;
    this->dx = L / double(nx);

    this->grid.reserve(nx);
    for (uint i = 0; i < nx; ++i)
    {
        this->grid.push_back(this->dx * double(i));
    }

    this->dt = dt;//0.99 * this->dx;
    this->tmax = tmax;

    this->spec.reserve(nspec);

    init_simulation();

    this->nspec = this->spec.size();

    // Initialize densities and fields after instantiation
    deposit_charge();
    solve_field();
}

Simulation::~Simulation()
{
}
//-----------------------------------------

void Simulation::add_species(uint npar, double Qpar, double density, std::function<void(Species &, uint)> init_fcn)
{
    this->spec.push_back(Species(npar, this->Nx, Qpar, density, init_fcn));
}

void Simulation::add_e_field(std::function<void(Field &, uint)> init_fcn)
{
    this->e_field = Field(this->Nx, this->dx, init_fcn);
}

void Simulation::add_b_field(std::function<void(Field &, uint)> init_fcn)
{
    this->b_field = Field(this->Nx, this->dx, init_fcn);
}

//-----------------------------------------

bool Simulation::dump_data()
{
    if (this->ndump > 0)
    {
        return !(this->n_iter % this->ndump);
    }
    else
    {
        return false;
    }

}

void Simulation::iterate()
{
    // Already deposited charge and fields on creation, so can immediately push
    map_field_to_species();
    push_species();
    deposit_charge();
    solve_field();

    ++this->n_iter;
}

std::vector<double> Simulation::get_total_density()
{
    std::vector<double> total_dens = std::vector<double>(this->Nx, 0.0);

    for (const auto &s : this->spec)
    {
        for (uint i = 0; i < this->Nx; ++i)
        {
            total_dens.at(i) += s.density_arr.at(i);
        }
    }

    return total_dens;
}

void Simulation::print_spec_density(uint i)
{
    spec.at(i).print_density();
}

/**********************************************************
PRIVATE FUNCTIONS
***********************************************************/
/**********************************************************
void Simulation::init_simulation() is defined for each run.
***********************************************************/

void Simulation::deposit_charge()
{
    for (auto &s : this->spec)
    {
        s.deposit_charge(this->dx, this->L_sys, this->Nx);
    }
}

void Simulation::solve_field()
{
    std::vector<double> total_dens_re = get_total_density();
    std::vector<double> total_dens_im = std::vector<double>(this->Nx, 0.0);

    this->e_field.solve_field(total_dens_re, total_dens_im); //TODO: having this function return a value and change state seems bad, maybe pass err as a parameter to also be changed?
}

void Simulation::map_field_to_species()
{
    for (auto &s : this->spec)
    {
        s.map_field_to_part(this->e_field, this->dx, this->L_sys, this->Nx);
    }
}

void Simulation::push_species()
{
    for (auto &s : this->spec)
    {
        // s.push_particles(n_iter, this->L_sys, this->dt, this->dx);
        s.boris_push(n_iter, this->L_sys, this->dt, this->dx);
    }
}
