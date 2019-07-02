#include "Field.h"

/**********************************************************
CONSTRUCTORS/DESTRUCTORS
***********************************************************/
Field::Field() //TODO: See if this can be removed
{
}

Field::Field(uint nx, double dx, std::function<void(Field &, uint)> init_fcn)
{
    this->size = nx;

    this->total_U = 0.0;

    std::vector<double> k = get_k_vec(this->size, dx);
    this->K2 = get_K2_vec(k, this->size, dx);
    this->kappa = get_kappa_vec(k, this->size, dx);

    init_field(init_fcn);
}

Field::~Field()
{
}
//-----------------------------------------

int Field::solve_field(std::vector<double> re, std::vector<double> im)
{
    // For total electrostatic energy diagnostic
    this->total_U = 0.0;

    int err = 0;
    const uint fft = 1;   // to perform fft
    const uint ifft = -1; // to perform ifft

    err = FFT(re, im, this->size, fft);
    if (err)
    {
        return err;
    }

    re.at(0) = im.at(0) = 0.0; // set offset density pert to 0

    for (uint i = 1; i < this->size; ++i) // avoid divide by 0, start at 1
    {
        // For total electrostatic energy diagnostic
        // Can ignore at i=0 since the values are 0
        this->total_U += ((re.at(i) * re.at(i)) + (im.at(i) * im.at(i))) / this->K2.at(i);

        //calculate phi
        re.at(i) /= this->K2.at(i);
        im.at(i) /= this->K2.at(i);

        //calculate E
        re.at(i) *= -this->kappa.at(i);
        im.at(i) *= this->kappa.at(i); //multiply i by i so get extra factor of -1
    }
    // By the end of this, the real part has become the imaginary part,
    // and vice versa, so we switch the order in the inverse fft
    err = FFT(im, re, this->size, ifft);
    if (err)
    {
        return err;
    }
    this->f1 = im;

    // For total electrostatic energy diagnostic
    this->total_U *= 0.5;

    return err;
}

void Field::print_field()
{
    for (auto &f : this->f1)
    {
        std::cout << f << '\t';
    }
    std::cout << std::endl;
}

/**********************************************************
PRIVATE FUNCTIONS
***********************************************************/
void Field::init_field(std::function<void(Field &, uint)> init_fcn)
{
    init_fcn(*this, this->size);
}
