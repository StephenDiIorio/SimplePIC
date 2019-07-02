#ifndef PARTICLE_H
#define PARTICLE_H

#include <iostream>

#include "ThreeVec.h"

class Particle
{
    private:
        ThreeVec pos;
        ThreeVec mom;
        ThreeVec local_e_field;
        ThreeVec local_b_field;

        double weight;

    public:
        /**********************************************************
        CONSTRUCTORS/DESTRUCTORS
        ***********************************************************/
        Particle();

        Particle(double weight);

        Particle(double pos_x, double pos_y, double pos_z, double mom_x, double mom_y, double mom_z);

        Particle(double pos_x, double pos_y, double pos_z, double mom_x, double mom_y, double mom_z, double weight);

        Particle(ThreeVec pos, ThreeVec mom);

        Particle(ThreeVec pos, ThreeVec mom, double weight);

        ~Particle();
        //-----------------------------------------

        /**********************************************************
        GETTER FUNCTIONS
        ***********************************************************/
        ThreeVec get_pos();
        ThreeVec get_pos() const;

        ThreeVec get_mom();
        ThreeVec get_mom() const;

        double get_weight();
        double get_weight() const;

        ThreeVec get_local_e_field();
        ThreeVec get_local_e_field() const;

        ThreeVec get_local_b_field();
        ThreeVec get_local_b_field() const;
        //-----------------------------------------

        /**********************************************************
        SETTER FUNCTIONS
        ***********************************************************/
        void set_pos(ThreeVec pos);
        void set_pos(const double pos_x, const double pos_y, const double pos_z);

        void set_mom(ThreeVec mom);
        void set_mom(const double mom_x, const double mom_y, const double mom_z);

        void set_weight(const double weight);

        void set_local_e_field(ThreeVec field);
        void set_local_e_field(const double x1, const double x2, const double x3);

        void set_local_b_field(ThreeVec field);
        void set_local_b_field(const double x1, const double x2, const double x3);
        //-----------------------------------------

        // friend std::ostream &operator<<(std::ostream &outstream, const Particle &p);
};

#endif
