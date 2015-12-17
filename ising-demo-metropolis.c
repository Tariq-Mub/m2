#include <stdio.h>
#include <math.h>

#include "ising_datastructure.h"
#include "ising.h"
#include "metropolis.h"
#include "progressbar.h"
#include "ising-demo-metropolis.h"

int main (int argc, char **argv)
{
    struct ising model;

    double t = 2.;                   // temperature (kT)
    double h = 0.;                   // magnetic field
    double j = 1.;                   // ferromagnetic coupling
    double tmin = 2., tmax = 3.;//, t;
    int lx = 32;
    int ly = lx;

    int mcsteps = 1024 * 128;        // 2^{17}
    int nt = 20;
    double dt = (tmax - tmin)/(nt - 1);
    t = tmin;

    if( ising_init (&model, lx, ly, j, h, t) != 0)
    {
        fprintf(stderr, "%s - cannot allocate spin array\n", argv[0]);
        argc = 1;
        exit(1);
    }

    print_header (&model, stderr, tmin, tmax);

    int therm_steps = mcsteps / 4;
    int report = therm_steps / 32;

    for(int i = 0; i < nt; i++){
    for (int s = 1; s <= therm_steps; s++)
    {
        (void) one_metropolis_step_per_spin (&model);
        if (s % report == 0)
        {
            text_progress_bar (s, therm_steps, 32, "Thermalization", stderr);
        }
    }
    fprintf (stderr, " Done\n");

    report = mcsteps / 32;
    double mav = 0., m2av = 0., eav = 0., e2av = 0.;
    double cv, xi;

    for (int s = 1; s <= mcsteps; s++)
    {
        (void) one_metropolis_step_per_spin (&model);
        double m = magnetization_per_spin (&model);
        double e = energy_per_spin (&model);

        mav += m;
        m2av += m * m;
        eav += e;
        e2av += e * e;

        if (s % report == 0)
        {
            text_progress_bar (s, mcsteps, 32, "Production", stderr);
        }

    }
    fprintf (stderr, " Done\n\n");

    mav /= mcsteps;
    eav /= mcsteps;
    m2av /= mcsteps;
    e2av /= mcsteps;

    // specific heat dE/dt
    cv = (e2av - eav*eav)/(t*t); 

    // magnetic susceptibility
    xi = (m2av - mav*mav)/t;

    printf ("%.3f   %.3f     % f  %f     % f  %f\n", t, h, mav, 
        xi, eav, cv);
    t += dt;
    ising_reinit(&model, h, t);
    }

    ising_free (model); 

    return 0;
}

void print_header (struct ising *model, FILE *fp, double tmin, double tmax)
{
    fprintf (fp, "\n");
    fprintf (fp, " Two-dimensional Ising model - Metropolis simulation\n");
    fprintf (fp, " ---------------------------------------------------\n");
    fprintf (fp, " Tmin/J = %.3f  Tmax/J = %.3f  H/J = %.3f  %dx%d\n\n", tmin/model->J, tmax/model->J,
        model->H/model->J, model->Lx, model->Ly);
}

