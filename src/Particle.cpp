#include "Particle.h"

/**********************************************************
CONSTRUCTORS/DESTRUCTORS
***********************************************************/
Particle::Particle() : Particle(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0)
{
}

Particle::Particle(double weight) : Particle(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, weight)
{
}

Particle::Particle(double pos_x, double pos_y, double pos_z, double mom_x, double mom_y, double mom_z) : Particle(pos_x, pos_y, pos_z, mom_x, mom_y, mom_z, 1.0)
{
}

Particle::Particle(double pos_x, double pos_y, double pos_z, double mom_x, double mom_y, double mom_z, double weight)
{
    this->pos = ThreeVec(pos_x, pos_y, pos_z);
    this->mom = ThreeVec(mom_x, mom_y, mom_z);
    this->weight = weight;
    this->local_e_field = ThreeVec();
    this->local_b_field = ThreeVec();
}

Particle::Particle(ThreeVec pos, ThreeVec mom) : Particle(pos.get_x(), pos.get_y(), pos.get_z(), mom.get_x(), mom.get_y(), mom.get_z(), 1.0)
{
}

Particle::Particle(ThreeVec pos, ThreeVec mom, double weight) : Particle(pos.get_x(), pos.get_y(), pos.get_z(), mom.get_x(), mom.get_y(), mom.get_z(), weight)
{
}

Particle::~Particle()
{
}
//-----------------------------------------

/**********************************************************
GETTER FUNCTIONS
***********************************************************/
ThreeVec Particle::get_pos()
{
    return this->pos;
}
ThreeVec Particle::get_pos() const
{
    return this->pos;
}

ThreeVec Particle::get_mom()
{
    return this->mom;
}
ThreeVec Particle::get_mom() const
{
    return this->mom;
}

double Particle::get_weight()
{
    return this->weight;
}
double Particle::get_weight() const
{
    return this->weight;
}

ThreeVec Particle::get_local_e_field()
{
    return this->local_e_field;
}
ThreeVec Particle::get_local_e_field() const
{
    return this->local_e_field;
}

ThreeVec Particle::get_local_b_field()
{
    return this->local_b_field;
}
ThreeVec Particle::get_local_b_field() const
{
    return this->local_b_field;
}
//-----------------------------------------

/**********************************************************
SETTER FUNCTIONS
***********************************************************/
void Particle::set_pos(ThreeVec pos)
{
    this->set_pos(pos.get_x(), pos.get_y(), pos.get_z());
}
void Particle::set_pos(const double pos_x, const double pos_y, const double pos_z)
{
    this->pos.set_all(pos_x, pos_y, pos_z);
}

void Particle::set_mom(ThreeVec mom)
{
    this->set_mom(mom.get_x(), mom.get_y(), mom.get_z());
}
void Particle::set_mom(const double mom_x, const double mom_y, const double mom_z)
{
    this->mom.set_all(mom_x, mom_y, mom_z);
}

void Particle::set_weight(const double weight)
{
    this->weight = weight;
}

void Particle::set_local_e_field(ThreeVec field)
{
    this->set_local_e_field(field.get_x(), field.get_y(), field.get_z());
}
void Particle::set_local_e_field(const double x1, const double x2, const double x3)
{
    this->local_e_field.set_all(x1, x2, x3);
}

void Particle::set_local_b_field(ThreeVec field)
{
    this->set_local_b_field(field.get_x(), field.get_y(), field.get_z());
}
void Particle::set_local_b_field(const double x1, const double x2, const double x3)
{
    this->local_b_field.set_all(x1, x2, x3);
}
//-----------------------------------------

// std::ostream &
// operator<<(std::ostream &outstream, const Particle &p)
// {
//     outstream << p.getpos() << p.getmom() << '\t';
//     return outstream;
// }
