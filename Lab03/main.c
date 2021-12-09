#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <windows.h>

const double x_0 = 0.01;
const double v_0 = 0.;
const double dt_0 = 1.;
const double S = 0.75;
const double p = 2.;
const double t_max = 40.;
const double alfa = 5.;
const double delta = pow(10, -10);
const double TOL_1 = pow(10, -2);
const double TOL_2 = pow(10, -5);

struct data
{
    double x;
    double v;
};

struct data metoda_trapezy(double x_n, double v_n, double dt)
{
    struct data d_n_1;
    double x_n_1 = x_n;
    double v_n_1 = v_n;
    double dx, dv, F, G;

    do
    {
        double a_11 = 1.;
        double a_12 = -(dt / 2.);
        double a_21 = -(dt / 2.) * (-2 * alfa * x_n_1 * v_n_1 - 1);
        double a_22 = 1 - (dt / 2.) * alfa * (1 - pow(x_n_1, 2));
        F = x_n_1 - x_n - dt / 2. * (v_n + v_n_1);
        G = v_n_1 - v_n - dt / 2. * ((alfa * (1 - pow(x_n, 2)) * v_n - x_n) + (alfa * (1 - pow(x_n_1, 2)) * v_n_1 - x_n_1));

        dx = (-1. * F * a_22 - (-1. * G) * a_12) / (a_11 * a_22 - a_12 * a_21);
        dv = (a_11 * -1. * G - a_21 * -1. * F) / (a_11 * a_22 - a_12 * a_21);

        x_n_1 += dx;
        v_n_1 += dv;

    } while ((fabs(dv) > delta) || (fabs(dx) > delta));

    d_n_1.x = x_n_1;
    d_n_1.v = v_n_1;

    return d_n_1;
}

struct data RK2(double x_n, double v_n, double dt)
{
    struct data data_n_1;
    double k1x = v_n;
    double k1v = alfa * (1 - pow(x_n, 2)) * v_n - x_n;
    double k2x = v_n + dt * k1v;
    double k2v = alfa * (1- pow((x_n + dt * k1x), 2)) * (v_n + dt * k1v) - (x_n + dt * k1x);
    data_n_1.x = x_n + dt/2. * (k1x + k2x);
    data_n_1.v = v_n + dt/2. * (k1v + k2v);

    return data_n_1;
}

void alg_kroku_czasowego(FILE* fp, struct data (* schemat_numeryczny)(double x_n, double v_n, double dt), double TOL)
{
    double t = 0.;
    double dt = dt_0;
    double x_n = x_0;
    double v_n = v_0;

    fprintf(fp, "%f \t %f \t %f \t %f \n", t, dt, x_n, v_n);

    do
    {
        struct data d2_n_1 = schemat_numeryczny(x_n, v_n, dt);
        d2_n_1 = schemat_numeryczny(d2_n_1.x, d2_n_1.v, dt);
        struct data d1_n_1 = schemat_numeryczny(x_n, v_n, dt * 2);

        double E_x = (d2_n_1.x- d1_n_1.x)/ (pow(2., p) - 1.);
        double E_v = (d2_n_1.v - d1_n_1.v) / (pow(2., p) - 1.);

        if(fmax(fabs(E_x), fabs(E_v)) < TOL)
        {
            t += 2. * dt;
            x_n = d2_n_1.x;
            v_n = d2_n_1.v;
            fprintf(fp, "%f \t %f \t %f \t %f \n", t, dt, x_n, v_n);
        }

        dt *= pow(S * TOL/(fmax(fabs(E_x), fabs(E_v))), (1./(p + 1.)));

    }while(t < t_max);
}

int main()
{
    FILE *rk2_TOL1 = fopen("RK2_TOL1.txt", "w");
    alg_kroku_czasowego(rk2_TOL1, RK2, TOL_1);
    fclose(rk2_TOL1);

    FILE *rk2_TOL2 = fopen("RK2_TOL2.txt", "w");
    alg_kroku_czasowego(rk2_TOL2, RK2, TOL_2);
    fclose(rk2_TOL2);

    FILE* trapezy_TOL1 = fopen("trapezy_TOL1.txt", "w");
    alg_kroku_czasowego(trapezy_TOL1, metoda_trapezy, TOL_1);
    fclose(trapezy_TOL1);

    FILE* trapezy_TOL2 = fopen("trapezy_TOL2.txt", "w");
    alg_kroku_czasowego(trapezy_TOL2, metoda_trapezy, TOL_2);
    fclose(trapezy_TOL2);

    return 0;
}
