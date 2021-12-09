#include <stdio.h>
#include <math.h>

const double y_0 = 1.;
const double lambda = -1.;
const double t0 = 0.;
const double t1 = 5.;
const double pi = 3.141592653589;
double t;
double y_dok;

double analytical(double t)
{
    return exp(lambda * t);
}

void euler(double dt, FILE *fp1, FILE *fp2)
{
    double y = y_0;
    for(t = t0; t < t1; t += dt)
    {
        y_dok = analytical(t);
        fprintf(fp1, "%lf \t %.10lf\n", t, y);
        fprintf(fp2, "%lf \t %.10lf\n", t, y - y_dok);
        y += dt * lambda * y;
    }
}

void RK2(double dt, FILE *fp1, FILE *fp2)
{
    double y = y_0;
    double k1, k2;
    for(t = t0; t < t1; t += dt)
    {
        y_dok = analytical(t);
        fprintf(fp1, "%lf \t %.10lf\n", t, y);
        fprintf(fp2, "%lf \t %.10lf\n", t, y - y_dok);
        k1 = lambda * y;
        k2 = lambda * (y + dt * k1);
        y += dt/2*(k1+k2);
    }
}

void RK4(double dt, FILE *fp1, FILE *fp2)
{
    double y = y_0;
    double k1, k2, k3, k4;

    for(t = t0; t < t1; t+=dt)
    {
        y_dok = analytical(t);
        fprintf(fp1, "%lf \t %.15lf\n", t, y);
        fprintf(fp2, "%lf \t %.15lf\n", t, y - y_dok);
        k1 = lambda * y;
        k2 = lambda * (y + dt/2 * k1);
        k3 = lambda * (y + dt/2 * k2);
        k4 = lambda * (y + dt * k3);
        y += dt/6 * (k1 + 2 * k2 + 2 * k3 + k4);
    }
}

double voltage(double factor, double omega, double t)
{
    return 10 * sin(factor * omega * t);
}

void RRZ2(double factor, FILE* fp1, FILE* fp2)
{
    double dt = 10e-4;
    double R = 100;
    double L = 0.1;
    double C = 0.001;
    double I = 0;
    double Q = 0;
    double omega = 1./sqrt(L*C);
    double t_0 = 2 * pi / omega;

    double k1Q, k1I, k2Q, k2I, k3Q, k3I, k4Q, k4I;

    for(t = 0.; t < 4 * t_0; t+=dt)
    {
        k1Q = I;
        k1I = voltage(factor, omega, t)/L - Q/L/C - R*I/L;
        k2Q = I + dt / 2.0 * k1I;
        k2I = voltage(factor, omega, t + dt/2.0)/L - (Q + dt/2.0 * k1Q)/L/C - R * (I + dt/2.0* k1I)/L;
        k3Q = I + dt/2.0 * k2I;
        k3I = voltage(factor, omega, t + dt/2.0)/L - (Q + dt/2.0 * k2Q)/L/C - R * (I + dt/2.0* k2I)/ L;
        k4Q = I + dt * k3I;
        k4I = voltage(factor, omega, t + dt)/L - (Q + dt * k3Q)/L/C - R * (I + dt * k3I)/L;
        fprintf(fp1, "%lf \t %lf\n", t, Q);
        fprintf(fp2, "%lf \t %lf\n", t, I);

        Q += dt/6.0*(k1Q + 2.0 * k2Q + 2.0* k3Q + k4Q);
        I += dt/6.0*(k1I + 2.0 *k2I + 2.0 *k3I + k4I);
    }
}

int main()
{
    FILE *analitic = fopen("analitic.txt", "w");

    for(t = t0; t < t1; t += 0.01)
    {
        fprintf(analitic, "%f \t %f\n", t, analytical(t));
    }
    fclose(analitic);

    FILE *euler001 = fopen("euler001.txt", "w");
    FILE *euler001_blad = fopen("euler001_blad.txt", "w");
    euler(0.01, euler001, euler001_blad);
    fclose(euler001);
    fclose(euler001_blad);

    FILE *euler01 = fopen("euler01.txt", "w");
    FILE *euler01_blad = fopen("euler01_blad.txt", "w");
    euler(0.1, euler01, euler01_blad);
    fclose(euler01);
    fclose(euler01_blad);

    FILE *euler1 = fopen("euler1.txt", "w");
    FILE *euler1_blad = fopen("euler1_blad.txt", "w");
    euler(1.0, euler1, euler1_blad);
    fclose(euler1);
    fclose(euler1_blad);

    FILE *rk2_001 = fopen("rk2_001.txt", "w");
    FILE *rk2_001_blad = fopen("rk2_001_blad.txt", "w");
    RK2(0.01, rk2_001, rk2_001_blad);
    fclose(rk2_001);
    fclose(rk2_001_blad);

    FILE *rk2_01 = fopen("rk2_01.txt", "w");
    FILE *rk2_01_blad = fopen("rk2_01_blad.txt", "w");
    RK2(0.1, rk2_01, rk2_01_blad);
    fclose(rk2_01);
    fclose(rk2_01_blad);

    FILE *rk2_1 = fopen("rk2_1.txt", "w");
    FILE *rk2_1_blad = fopen("rk2_1_blad.txt", "w");
    RK2(1, rk2_1, rk2_1_blad);
    fclose(rk2_1);
    fclose(rk2_1_blad);

    FILE *rk4_001 = fopen("rk4_001.txt", "w");
    FILE *rk4_001_blad = fopen("rk4_001_blad.txt", "w");
    RK4(0.01, rk4_001, rk4_001_blad);
    fclose(rk4_001);
    fclose(rk4_001_blad);

    FILE *rk4_01 = fopen("rk4_01.txt", "w");
    FILE *rk4_01_blad = fopen("rk4_01_blad.txt", "w");
    RK4(0.1, rk4_01, rk4_01_blad);
    fclose(rk4_01);
    fclose(rk4_01_blad);


    FILE *rk4_1 = fopen("rk4_1.txt", "w");
    FILE *rk4_1_blad = fopen("rk4_1_blad.txt", "w");
    RK4(1, rk4_1, rk4_1_blad);
    fclose(rk4_1);
    fclose(rk4_1_blad);

    FILE *rrz2_q05 = fopen("rrz2_q05.txt", "w");
    FILE *rrz2_i05 = fopen("rrz2_i05.txt", "w");
    RRZ2(0.5, rrz2_q05, rrz2_i05);
    fclose(rrz2_q05);
    fclose(rrz2_i05);

    FILE *rrz2_q08 = fopen("rrz2_q08.txt", "w");
    FILE *rrz2_i08 = fopen("rrz2_i08.txt", "w");
    RRZ2(0.8, rrz2_q08, rrz2_i08);
    fclose(rrz2_q08);
    fclose(rrz2_i08);

    FILE *rrz2_q1 = fopen("rrz2_q1.txt", "w");
    FILE *rrz2_i1 = fopen("rrz2_i1.txt", "w");
    RRZ2(1.0, rrz2_q1, rrz2_i1);
    fclose(rrz2_q1);
    fclose(rrz2_i1);

    FILE *rrz2_q12 = fopen("rrz2_q12.txt", "w");
    FILE *rrz2_i12 = fopen("rrz2_i12.txt", "w");
    RRZ2(1.2, rrz2_q12, rrz2_i12);
    fclose(rrz2_q12);
    fclose(rrz2_i12);

    return 0;
}
