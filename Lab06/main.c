#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "mgmres.h"
#include "mgmres.c"
#include <stdbool.h>

double delta = 0.1;
int e_1 = 1;
int e_2 = 1;
int n_x = 4;
int n_y = 4;
int V_1 = 10;
int V_2 = -10;
int V_3 = 10;
int V_4 = -10;

int reindex_j(int l, int n_x)
{
    return floor(l/(n_x+1.));
}

int reindex_i(int l, int n_x)
{
    return l - reindex_j(l, n_x) * (n_x + 1);
}

double ro(double x, double y, double x_max, double y_max, double sigma)
{
    return exp(-1. * pow(x - 0.25 * x_max, 2.)/pow(sigma, 2.) - pow(y - 0.5 * y_max, 2.)/pow(sigma, 2.)) - exp(-1. * pow(x - 0.75 * x_max, 2.)/pow(sigma, 2.) - pow(y - 0.5 * y_max, 2.)/pow(sigma, 2.));
}

int El(int n_x, int l, int epsilon_1, int epsilon_2)
{
    if(reindex_i(l, n_x) <= n_x/2)
    {
        return epsilon_1;
    }
    else
    {
        return epsilon_2;
    }
}

void poisson(double n_x, double n_y, double epsilon_1, double epsilon_2, double V_1, double V_2, double V_3, double V_4, bool is_ro, FILE* fp)
{
    int N = (n_x + 1)*(n_y + 1);

    double a[5 * N];
    int ja[5 * N];
    int ia[N + 1];

    double x_max = delta *n_x;
    double y_max = delta * n_y;
    double sigma = x_max/10.;

    for(int i = 0; i < 5 * N; i++)
    {
        a[i] = 0;
        ja[i] = 0;
    }

    double b[N];
    double V[N];

    for(int i=0; i< N + 1; i++)
    {
        ia[i] = -1;
    }

    int k = -1;

    for(int l = 0; l < N; l++)
    {
        int brzeg = 0;
        double vb = 0.;

        if(reindex_i(l, n_x) == 0)
        {
            brzeg = 1;
            vb = V_1;
        }

        if(reindex_j(l, n_x) == n_y)
        {
            brzeg = 1;
            vb = V_2;
        }

        if(reindex_i(l, n_x) == n_x)
        {
            brzeg = 1;
            vb = V_3;
        }

        if(reindex_j(l, n_x) == 0)
        {
            brzeg = 1;
            vb = V_4;
        }

        if(!is_ro)
        {
            b[l] = -ro(delta * reindex_i(l, n_x), delta * reindex_j(l, n_x), x_max, y_max, sigma);
        }
        else
        {
            b[l] = 0.;
        }

        if(brzeg == 1)
        {
            b[l] = vb;
        }

        ia[l] = -1;

        if(l - n_x - 1 >= 0 && brzeg == 0)
        {
            k++;

            if(ia[l] < 0)
            {
                ia[l] = k;
            }

            a[k] = El(n_x, l, epsilon_1, epsilon_2)/ (delta * delta);
            ja[k] = l - n_x - 1;
        }

        if(l - 1 >= 0 && brzeg == 0)
        {
            k++;

            if(ia[l] < 0)
            {
                ia[l] = k;
            }

            a[k] = El(n_x, l, epsilon_1, epsilon_2)/ (delta* delta);
            ja[k] = l - 1;
        }

        k++;

        if(ia[l] < 0)
        {
            ia[l] = k;
        }

        if(brzeg == 0)
        {
            a[k] = - (2 * El(n_x, l, epsilon_1, epsilon_2) + El(n_x, l + 1, epsilon_1, epsilon_2) + El(n_x, l + n_x + 1, epsilon_1, epsilon_2)) / (delta* delta);
        }
        else
        {
            a[k] = 1.;
        }

        ja[k] = l;

        if(l < N && brzeg == 0)
        {
            k++;

            a[k] = El(n_x, l + 1, epsilon_1, epsilon_2)/ (delta * delta);
            ja[k] = l + 1;
        }

        if(l < N - n_x - 1 && brzeg == 0)
        {
            k++;

            a[k] = El(n_x, l + n_x + 1,  epsilon_1, epsilon_2)/ (delta * delta);
            ja[k] = l + n_x + 1;
        }

        //fprintf(fp, "%d \t %d \t %d \t %f \n", l, reindex_i(l, n_x), reindex_j(l, n_x), b[l]);
    }

    int nz_num = k + 1;
    ia[N] = nz_num;

    int itr_max = 500;
    int mr = 500;
    double tol_abs = 1.0e-8;
    double tol_rel = 1.0e-8;


    /*for(int iter = 0; iter < nz_num; iter++)
    {
        fprintf(fp, "%d \t %0.f \n", iter, a[iter]);
    }*/

    pmgmres_ilu_cr(N, nz_num, ia, ja, a, V, b, itr_max, mr, tol_abs, tol_rel);

    double x= 0;

    for(int iter = 0; iter < N; iter++)
    {
        if(delta * reindex_i(iter, n_x) < x)
        {
            fprintf(fp, "\n");
        }
        fprintf(fp, "%f \t %f \t %lf \n", delta* reindex_i(iter, n_x), reindex_j(iter, n_x) * delta, V[iter]);
        x = delta *reindex_i(iter, n_x);
    }
}


int main()
{
    FILE* vect = fopen("vector.txt", "w");
    poisson(4, 4, 1, 1, 10, -10, 10, -10, true, vect);
    fclose(vect);

    FILE* matr = fopen("matrix.txt", "w");
    poisson(4, 4, 1, 1, 10, -10, 10, -10, true, matr);
    fclose(matr);

    FILE* no_ro_50 = fopen("no_ro_50.txt", "w");
    poisson(50, 50, 1, 1, 10, -10, 10, -10, true, no_ro_50);
    fclose(no_ro_50);

    FILE* no_ro_100 = fopen("no_ro_100.txt", "w");
    poisson(100, 100, 1, 1, 10, -10, 10, -10, true, no_ro_100);
    fclose(no_ro_100);

    FILE* no_ro_200 = fopen("no_ro_200.txt", "w");
    poisson(200, 200, 1, 1, 10, -10, 10, -10, true, no_ro_200);
    fclose(no_ro_200);

    FILE* ro_1_1 = fopen("ro_1_1.txt", "w");
    poisson(100, 100, 1, 1, 0, 0, 0, 0, false, ro_1_1);
    fclose(ro_1_1);

    FILE* ro_1_2 = fopen("ro_1_2.txt", "w");
    poisson(100, 100, 1, 2, 0, 0, 0, 0, false, ro_1_2);
    fclose(ro_1_2);

    FILE* ro_1_10 = fopen("ro_1_10.txt", "w");
    poisson(100, 100, 1, 10, 0, 0, 0, 0, false, ro_1_10);
    fclose(ro_1_10);

    return 0;
}
