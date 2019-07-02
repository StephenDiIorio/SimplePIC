/*
****************************************************************************
*                                                                          *
*   FFT Header file, based on code from www.codeproject.com, originally    *
*   based on code from Numerical Recipes with refinement in speed.         *
*   Added functionality to load in complex and real arrays          	   *
*   separately. This requires a number of values which is a power of 2.    *
*   This version is also normalized so that the 1/N is taken into account. *
*   AGRT 2010                                                              *
*                                                                          *
****************************************************************************

    data_re -> float array that represent the real array of complex samples
    data_im -> float array that represent the imag array of complex samples
    NVALS -> length of real or imaginary arrays (N^2 order number)
    isign -> 1 to calculate FFT and -1 to calculate Reverse FFT

    The function returns an integer, 0 if FFT ran, 1 otherwise. It will be
    one if NVALS is not a power of 2
*/

#include "FFT.h"

#define SWAP(a, b) \
    tempr = (a);   \
    (a) = (b);     \
    (b) = tempr
//tempr is a variable from our FFT function

double sinc(const double x)
{
    if (x != 0.0)
    {
        return sin(x) / x;
    }
    else
    {
        return 1.0;
    }
}

int FFT(std::vector<double> &data_re, std::vector<double> &data_im, const unsigned long NVALS, const int isign)
{
    // TEST THAT NVALS IS A POWER OF 2
    const int err = !(NVALS && !(NVALS & (NVALS - 1)));
    if (!err)
    {
        //variables for trigonometric recurrences
        unsigned long mmax, m, j, istep, i;
        double wtemp, wr, wpr, wpi, wi, theta, tempr, tempi;
        /*
	    the complex array is real+complex so the array
	    as a size n = 2* number of samples
	    real part is data[index] and the complex part is data[index+1]
        */
        const unsigned long n = NVALS << 1; // bitwise multiply by 2

        /*
	    Pack components into double array of size n - the ordering
	    is data[0]=real[0], data[1]=imag[0] .......
        */
        std::vector<double> data(n);
        for (i = 0, j = 0; j < n; ++i, j += 2)
        {
            data.at(j) = data_re.at(i);
            data.at(j + 1) = data_im.at(i);
        }
        /*
	    binary inversion (note that
	    the indexes start from 1 which means that the
	    real part of the complex is on the odd-indexes
	    and the complex part is on the even-indexes
        */

        j = 0;
        for (i = 0; i < n / 2; i += 2)
        {
            if (j > i)
            {
                //swap the real part
                SWAP(data.at(j), data.at(i));

                //swap the complex part
                SWAP(data.at(j + 1), data.at(i + 1));

                // checks if the changes occurs in the first half
                // and use the mirrored effect on the second half
                if ((j / 2) < (n / 4))
                {
                    //swap the real part
                    SWAP(data.at((n - (i + 2))), data.at((n - (j + 2))));

                    //swap the complex part
                    SWAP(data.at((n - (i + 2)) + 1), data.at((n - (j + 2)) + 1));
                }
            }

            m = NVALS;
            while (m >= 2 && j >= m)
            {
                j -= m;
                m >>= 1;
            }
            j += m;
        }

        //Danielson-Lanzcos routine
        mmax = 2;
        while (n > mmax)
        {
            istep = mmax << 1;
            theta = isign * (2.0 * M_PI / mmax);
            wtemp = sin(0.5 * theta);
            wpr = -2.0 * wtemp * wtemp;
            wpi = sin(theta);
            wr = 1.0;
            wi = 0.0;
            //internal loops

            for (m = 1; m < mmax; m += 2)
            {
                for (i = m; i <= n; i += istep)
                {
                    j = i + mmax;
                    tempr = wr * data.at(j - 1) - wi * data.at(j);
                    tempi = wr * data.at(j) + wi * data.at(j - 1);
                    data.at(j - 1) = data.at(i - 1) - tempr;
                    data.at(j) = data.at(i) - tempi;
                    data.at(i - 1) += tempr;
                    data.at(i) += tempi;
                }
                wr = (wtemp = wr) * wpr - wi * wpi + wr;
                wi = wi * wpr + wtemp * wpi + wi;
            }
            mmax = istep;
        }
        /*
	    Return elements to real and complex components
        */
        if (isign > 0)
        {
            for (i = 0, j = 0; j < n; ++i, j += 2)
            {
                data_re.at(i) = data.at(j);
                data_im.at(i) = data.at(j + 1);
            }
        }
        else
        {
            double invNVALs = 1.0 / NVALS;
            for (i = 0, j = 0; j < n; ++i, j += 2)
            {
                data_re.at(i) = data.at(j) * invNVALs;
                data_im.at(i) = data.at(j + 1) * invNVALs;
            }
        }
    }
    return err;
}

std::vector<double> get_k_vec(const uint size, const double dx)
{
    std::vector<double> k = std::vector<double>(size);
    const double kmax = M_PI / dx;

    for (uint i = 0; i < size / 2; ++i)
    {
        k.at(i) = kmax * i / (size / 2);
    }
    for (uint i = size / 2; i < size; ++i)
    {
        k.at(i) = kmax * i / (size / 2) - 2 * kmax;
    }

    return k;
}

std::vector<double> get_K2_vec(const std::vector<double> k, const uint size, const double dx)
{
    std::vector<double> K2 = std::vector<double>(size);
    double val;

    for (uint i = 0; i < size; ++i)
    {
        val = k.at(i) * sinc(k.at(i) * dx / 2.0);
        K2.at(i) = val * val;
    }

    return K2;
}

std::vector<double> get_kappa_vec(const std::vector<double> k, const uint size, const double dx)
{
    std::vector<double> kappa = std::vector<double>(size);

    for (uint i = 0; i < size; ++i)
    {
        kappa.at(i) = k.at(i) * sinc(k.at(i) * dx);
    }

    return kappa;
}
