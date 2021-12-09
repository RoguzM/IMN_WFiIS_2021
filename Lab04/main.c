#include <stdio.h>
#include <math.h>
#include <stdlib.h>

const double epsilon = 1.;
const double dt = 0.1;
const int n_x = 150;
const int n_y = 100;
const double V_1 = 10.;
const double V_2 = 0.;
const double x_max = 15.;
const double y_max = 10.;
const double sig_x = 1.5;
const double sig_y = 1.;
const double TOL = pow(10, -8);

double ro_1(double x, double y)
{
    return exp((-pow((x - 0.35 * x_max),2))/(pow(sig_x, 2)) - (pow(y - 0.5 * y_max, 2))/(pow(sig_y, 2)));
}

double ro_2(double x, double y)
{
    return -exp((-pow((x - 0.65 * x_max),2))/(pow(sig_x, 2)) - (pow(y - 0.5 * y_max, 2))/(pow(sig_y, 2)));
}

void global_relax(double omega, FILE* fp_S, FILE* fp_Vn, FILE* fp_E)
{
    double V_n[n_x + 1][n_y + 1];
    double V_s[n_x + 1][n_y + 1];
    double error[n_x + 1][n_y + 1];
    double ro[n_x + 1][n_y + 1];

    int iter = 0;

    for(int i = 0; i < n_x + 1; i++)
    {
        for(int j = 0; j < n_y + 1; j++)
        {
            V_n[i][j] = 0.;
            V_s[i][j] = 0.;
            error[i][j] = 0.;
            ro[i][j] = ro_1(i * dt, j * dt) + ro_2(i * dt, j * dt);
        }
    }

    for(int i = 0; i < n_x + 1; i++)
    {
        V_n[i][0] = V_1;
        V_s[i][0] = V_1;
    }

    double S_it, S_it_1;
    S_it = 0.;
    S_it_1 = 0.;

    do
    {
        for(int i = 1; i < n_x; i++)
        {
            for(int j = 1; j < n_y; j++)
            {
                V_n[i][j] = 0.25 * (V_s[i + 1][j] + V_s[i - 1][j] + V_s[i][j + 1] + V_s[i][j - 1] + pow(dt, 2)/epsilon * ro[i][j]);
            }
        }

        for(int j = 1;  j < n_y + 1; j++)
        {
            V_n[0][j] = V_n[1][j];
            V_n[n_x][j] = V_n[n_x - 1][j];
        }

        for(int i = 0; i < n_x + 1; i++)
        {
            for(int j = 0; j < n_y + 1; j++)
            {
                V_s[i][j] = (1. - omega) * V_s[i][j] + omega * V_n[i][j];
            }
        }

        S_it_1 = S_it;
        S_it = 0.;


        for(int i = 0; i < n_x; i++)
        {
            for(int j = 0; j < n_y; j++)
            {
                S_it +=  pow(dt, 2) * (0.5 * pow(((V_n[i+1][j] - V_n[i][j])/dt), 2) + 0.5 * pow((V_n[i][j+1] - V_n[i][j])/dt, 2) - ro[i][j] * V_n[i][j]);
            }
        }

        fprintf(fp_S, "%d \t %f \n", iter, S_it);
        iter++;
    }while(fabs((S_it - S_it_1) / S_it_1) > TOL);

    for(int i = 1; i < n_x; i++)
    {
        for(int j = 1; j < n_y; j++)
        {
            fprintf(fp_Vn, "%f \t %f \t %f \n", dt * i, dt * j, V_n[i][j]);
        }
    }

    for(int i = 1; i < n_x; i++)
    {
        for(int j = 1; j < n_y; j++)
        {
            error[i][j] = (V_n[i+1][j] - 2 * V_n[i][j] + V_n[i - 1][j])/ (dt* dt) + (V_n[i][j + 1] - 2 * V_n[i][j] + V_n[i][j - 1])/ (dt * dt) + ro[i][j]/epsilon;
        }
    }

        for(int i = 1; i < n_x; i++)
    {
        for(int j = 1; j < n_y; j++)
        {
            fprintf(fp_E, "%f \t %f \t %f \n", dt * i, dt * j, error[i][j]);
        }
    }

}

void local_relax(double omega, FILE* fp_S)
{
    double V_n[n_x + 1][n_y + 1];
    double ro[n_x + 1][n_y +1];

    int iter = 0;

    for(int i = 0; i < n_x + 1; i++)
    {
        for(int j = 0; j < n_y + 1; j++)
        {
            V_n[i][j] = 0.;
            ro[i][j] = ro_1(i * dt, j * dt) + ro_2(i * dt, j * dt);
        }
    }

    for(int i = 0; i < n_x; i++)
    {
        V_n[i][0] = V_1;
    }

    double S_it, S_it_1;
    S_it = 0.;
    S_it_1 = 0.;

    do
    {
        for(int i = 1; i < n_x; i++)
        {
            for(int j = 1; j < n_y; j++)
            {
                V_n[i][j] = 0.25 * omega * (V_n[i + 1][j] + V_n[i - 1][j] + V_n[i][j + 1] + V_n[i][j - 1] + pow(dt, 2)/epsilon * ro[i][j]) + (1. - omega) * V_n[i][j];
            }
        }

        for(int j = 1;  j < n_y + 1; j++)
        {
            V_n[0][j] = V_n[1][j];
            V_n[n_x][j] = V_n[n_x - 1][j];
        }

        S_it_1 = S_it;
        S_it = 0.;

        for(int i = 0; i < n_x; i++)
        {
            for(int j = 0; j < n_y; j++)
            {
                S_it +=  pow(dt, 2) * (0.5 * pow(((V_n[i+1][j] - V_n[i][j])/dt), 2) + 0.5 * pow((V_n[i][j+1] - V_n[i][j])/dt, 2) - ro[i][j] * V_n[i][j]);
            }
        }
        fprintf(fp_S, "%d \t %f \n", iter, S_it);
        iter++;
    }while(fabs((S_it - S_it_1)/S_it_1) > TOL);

}

int main()
{

    FILE* global_relax_S_06 = fopen("global_relax_S_06.txt", "w");
    FILE* global_relax_V_06 = fopen("global_relax_V_06.txt", "w");
    FILE* global_relax_E_06 = fopen("global_relax_E_06.txt", "w");
    global_relax(0.6, global_relax_S_06, global_relax_V_06, global_relax_E_06);
    fclose(global_relax_S_06);
    fclose(global_relax_V_06);
    fclose(global_relax_E_06);

    FILE* global_relax_S_1 = fopen("global_relax_S_1.txt", "w");
    FILE* global_relax_V_1 = fopen("global_relax_V_1.txt", "w");
    FILE* global_relax_E_1 = fopen("global_relax_E_1.txt", "w");
    global_relax(1., global_relax_S_1, global_relax_V_1, global_relax_E_1);
    fclose(global_relax_S_1);
    fclose(global_relax_V_1);
    fclose(global_relax_E_1);

    FILE* local_relax_1 = fopen("local_relax_S_1.txt", "w");
    local_relax(1., local_relax_1);
    fclose(local_relax_1);

    FILE* local_relax_14 = fopen("local_relax_S_14.txt", "w");
    local_relax(1.4, local_relax_14);
    fclose(local_relax_14);

    FILE* local_relax_18 = fopen("local_relax_S_18.txt", "w");
    local_relax(1.8, local_relax_18);
    fclose(local_relax_18);

    FILE* local_relax_19 = fopen("local_relax_S_19.txt", "w");
    local_relax(1.9, local_relax_19);
    fclose(local_relax_19);


    return 0;
}
