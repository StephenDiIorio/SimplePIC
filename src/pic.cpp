#include <cmath>
#include <iostream>
#include <fstream>
using namespace std;

#include "Simulation.h"
#include "two_stream.h"

int main()
{
    ofstream x_dom, t_dom, dens_out, dens_out_b, dens_out_p, field_out, KE_out, U_out, totE_out, part_x, part_px, part_py;

    x_dom.open("x_domain.txt");
    t_dom.open("t_domain.txt");
    dens_out.open("dens.txt");
    field_out.open("field.txt");

    double ke = 0.0, u = 0.0, tote = 0.0;
    KE_out.open("KE.txt");
    U_out.open("U.txt");
    totE_out.open("totE.txt");

    std::vector<double> phase;
    part_x.open("part_x.txt");
    part_px.open("part_px.txt");
    part_py.open("part_py.txt");

    Simulation sim(ndump, Nx, L_sys, dt, tmax);

    double t;
    for (t = 0.0; t < sim.tmax; t += sim.dt)
    {
        if (sim.dump_data())
        {
            for (auto &d : sim.spec.at(0).density_arr)
            {
                dens_out << d << '\t';
            }
            dens_out << '\n';

            for (auto &f : sim.e_field.f1)
            {
                field_out << f << '\t';
            }
            field_out << '\n';

            ke = sim.spec.at(0).total_KE;
            u = sim.e_field.total_U;
            tote = ke + u;
            KE_out << ke << '\t';
            U_out << u << '\t';
            totE_out << tote << '\t';

            for (uint i = 0; i < sim.nspec; ++i)
            {
                phase = sim.spec.at(i).get_x_phasespace();
                for (uint j = 0; j < sim.spec.at(i).npar; ++j)
                {
                    part_x << phase[j] << '\t';
                }
                phase = sim.spec.at(i).get_px_phasespace();
                for (uint j = 0; j < sim.spec.at(i).npar; ++j)
                {
                    part_px << phase[j] << '\t';
                }
                phase = sim.spec.at(i).get_py_phasespace();
                for (uint j = 0; j < sim.spec.at(i).npar; ++j)
                {
                    part_py << phase[j] << '\t';
                }
            }
            part_x << '\n';
            part_px << '\n';
            part_py << '\n';

            t_dom << t << '\t';
        }

        sim.iterate();
    }
    t_dom << '\n';
    KE_out << '\n';
    U_out << '\n';
    totE_out << '\n';


    for (uint i = 0; i < sim.Nx; ++i)
    {
        x_dom << sim.grid.at(i) << '\t';
    }
    x_dom << '\n';

    x_dom.close();
    t_dom.close();
    dens_out.close();
    field_out.close();
    KE_out.close();
    U_out.close();
    totE_out.close();
    part_x.close();
    part_px.close();
    part_py.close();

    return 0;
}
