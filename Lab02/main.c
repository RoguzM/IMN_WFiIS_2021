#include <stdio.h>
#include <math.h>
#include <stdlib.h>


const double beta = 0.001;
const int N = 500;
const double gamm = 0.1;
const double t_max = 100.;
const double u_0 = 1.;
const double TOL = 10e-6;
const double a_11 = 0.25;
const double a_22 = 0.25;
const double a_12 = 0.25 - sqrt(3)/6.;
const double a_21 = 0.25 + sqrt(3)/6.;
const double c_1 = 0.5 - sqrt(3)/6.;
const double c_2 = 0.5 + sqrt(3)/6.;
const double b_1 = 0.5;
const double b_2 = 0.5;

double fun(double u)
{
    return (beta * N - gamm) * u - beta * u * u;
}

void Picard(FILE* fp)
{
    double u_n = u_0;
    double u_n_1 = u_n;
    double u_mi = 0;
    int it;
    for(double dt = 0; dt < t_max; dt += 0.1)
    {
        fprintf(fp, "%.1lf \t %lf \t %lf \n", dt, u_n, N - u_n);
        it = 0;
        u_n_1 = u_n;

        while((fabs(u_n - u_mi) > TOL) && it <= 20)
        {
            u_mi = u_n_1;
            u_n_1 = u_n + (0.1/2.0) * (fun(u_n) + fun(u_mi));
            it++;
        }
        u_n = u_n_1;
    }
}

void Newton(FILE* fp)
{
    double u_n = u_0;
    double u_n_1 = u_n;
    double u_mi = 0;
    int it;

    for(double dt = 0; dt < t_max; dt += 0.1)
    {
        it = 0;
        fprintf(fp, "%.1lf \t %lf \t %lf \n", dt, u_n, N - u_n);
        u_n_1 = u_n;
        while(fabs(u_n - u_mi) > TOL && it <= 20)
        {
            u_mi = u_n_1;
            u_n_1 = u_mi - (u_mi - u_n - 0.1/2.0*(fun(u_n) + fun(u_mi)))/ (1. - 0.1/2.0 * ((beta * N - gamm) - 2.0 * beta * u_mi));
            it++;
        }
        u_n = u_n_1;
    }
}

double F1(double u_n, double U_1, double U_2)
{
    return (U_1 - u_n - 0.1 * (a_11 * ((beta * N - gamm) * U_1 - beta * pow(U_1,2)) + a_12 * ((beta * N - gamm) * U_2 - beta * pow(U_2, 2))));
}

double F2(double u_n, double U_1, double U_2)
{
    return ( U_2 - u_n - 0.1 * (a_21 * ((beta * N - gamm) * U_1 - beta * pow(U_1,2)) + a_22 * ((beta * N - gamm) * U_2 - beta * pow(U_2, 2))));
}

void RK2(FILE* fp)
{
    double U_1, U_2, dU_1, dU_2, m_11, m_12, m_21, m_22;
    double u_n = u_0;
    double u_n_1 = u_n;
    int it;

    for(double dt = 0.; dt < t_max; dt+= 0.1)
    {
        it = 0;
        U_1 = u_n;
        U_2 = u_n;
        fprintf(fp, "%.1lf \t %lf \t %lf \n", dt, u_n, N - u_n);

        while(fabs(dU_1) > TOL && fabs(dU_2) > TOL && it <= 20)
        {
            m_11 = 1. - 0.1 * a_11 * ((beta * N - gamm) - 2. * beta * U_1);
            m_12 = -0.1 * a_12 * ((beta * N - gamm) - 2. * beta * U_2);
            m_22 = 1. - 0.1 * a_22 * ((beta * N - gamm) - 2. * beta * U_2);
            m_21 = -0.1 * a_21 * ((beta * N - gamm) - 2. * beta * U_1);
            dU_1 = ((F2(u_n, U_1, U_2) * m_12) - F1(u_n, U_1, U_2)* m_22)/ (m_11 * m_22 - m_12 * m_21);
            dU_2 = ((F1(u_n, U_1, U_2) * m_21) - F1(u_n, U_1, U_2)* m_11)/ (m_11 * m_22 - m_21 * m_12);
            U_1 += dU_1;
            U_2 += dU_2;
            it++;
        }
        u_n += 0.1 * (b_1 * fun(U_1) + b_2 * fun(U_2));
    }

}

int main()
{
    FILE* zad1 = fopen("zadanie1.txt", "w");
    Picard(zad1);
    fclose(zad1);

    FILE* zad2 = fopen("zadanie2.txt", "w");
    Newton(zad2);
    fclose(zad2);

    FILE* zad3 = fopen("zadanie3.txt", "w");
    RK2(zad3);
    fclose(zad3);

    return 0;
}
