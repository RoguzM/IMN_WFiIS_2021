#include <stdio.h>
#include <math.h>
#include <stdlib.h>

const int n_x = 400;
const int n_y = 90;
const int i_1 = 200;
const int i_2 = 210;
const int j_1 = 50;
const double delta = 0.01;
const double sigma = 0.1;
const double x_A = 0.45;
const double y_A = 0.45;
const int IT_MAX = 10000;

void V_field(double v_x[n_x + 1][n_y + 1], double v_y[n_x + 1][n_y + 1], double psi[n_x + 1][n_y + 1], FILE* fp_1, FILE* fp_2)
{
    for(int i = 1; i <= n_x - 1; i++)
    {
        for(int j = 1; j <= n_y - 1; j++)
        {
            v_x[i][j] = (psi[i][j + 1] - psi[i][j - 1]) / (2. * delta);
            v_y[i][j] = -(psi[i + 1][j] - psi[i - 1][j]) / (2. * delta);
        }
    }

    for(int i = i_1; i <= i_2; i++)
    {
        for(int j = 0; j <= j_1; j++)
        {
            v_x[i][j] = 0.;
            v_y[i][j] = 0.;
        }
    }

    for(int i = 1; i <= n_x - 1; i++)
    {
        v_x[i][0] = 0.;
        v_y[i][n_y] = 0.;
    }

    for(int j = 0; j <= n_y; j++)
    {
        v_x[0][j] = v_x[1][j];
        v_x[n_x][j] = v_x[n_x - 1][j];
    }

    for(int i = 0; i <= n_x; i++)
    {
        for(int j = 0; j <= n_y; j++)
        {
            fprintf(fp_1, "%d \t %d \t %g \n", i, j, v_x[i][j]);
            fprintf(fp_2, "%d \t %d \t %g \n", i, j, v_y[i][j]);
        }
        fprintf(fp_1, "\n");
        fprintf(fp_2, "\n");
    }
}

double V_max(double v_x[n_x + 1][n_y + 1], double v_y[n_x + 1][n_y + 1])
{
    double maximum = 0.;
    for(int i = 0; i <= n_x; i++)
    {
        for(int j = 0; j <= n_y; j++)
        {
            if(maximum < sqrt(v_x[i][j] * v_x[i][j] + v_y[i][j] * v_y[i][j]))
            {
                maximum = sqrt(v_x[i][j] * v_x[i][j] + v_y[i][j] * v_y[i][j]);
            }
        }
    }
    return maximum;
}

void fill_u(double u_0[n_x+ 1][n_y + 1])
{
    for(int i = 0; i <= n_x; i++)
    {
        for(int j = 0; j <= n_y; j++)
        {
            u_0[i][j] = 1./(2. * M_PI * sigma * sigma) * exp( -(pow((delta * i) - x_A, 2) + pow((delta * j) - y_A, 2))/(2. * sigma * sigma));
        }
    }
}

void advection_diffusion(double u_0[n_x + 1][n_y + 1], double u_1[n_x + 1][n_y + 1], double v_x[n_x + 1][n_y + 1], double v_y[n_x + 1][n_y + 1], double psi[n_x + 1][n_y + 1], FILE*  fp_1, FILE* fp_2, double D, double delta_t)
{
    for(int it = 0; it <= IT_MAX; it++)
    {
        for(int i = 0; i <= n_x; i++)
        {
            for(int j = 0; j <= n_y; j++)
            {
                u_1[i][j] = u_0[i][j];
            }
        }

        for(int k = 1; k <= 20; k++)
        {
            for(int i = 0; i <= n_x; i++)
            {
                for(int j = 1; j <= n_y - 1; j++)
                {
                    if(i < i_1 || i > i_2 || j > j_1)
                    {
                        if(i == 0)
                        {
                            u_1[i][j] = (1. / (1. + ((2. * D * delta_t) / (delta * delta)))) * (u_0[i][j] - (delta_t / 2.) * v_x[i][j] * (((u_0[i + 1][j] - u_0[n_x][j])/(delta * 2.)) + (u_1[i + 1][j] - u_1[n_x][j])/ (2. * delta)) - (delta_t / 2.) * v_y[i][j] *((u_0[i][j + 1] - u_0[i][j - 1])/(2. * delta) + (u_1[i][j + 1] - u_1[i][j - 1])/(2. * delta)) + (delta_t / 2.) * D *((u_0[i + 1][j] + u_0[n_x][j] + u_0[i][j + 1] + u_0[i][j - 1] - 4. * u_0[i][j])/(delta * delta) + (u_1[i + 1][j] + u_1[n_x][j] + u_1[i][j + 1] + u_1[i][j - 1])/(delta * delta)));
                        }
                        else if(i == n_x)
                        {
                            u_1[i][j] = (1. / (1. + ((2. * D * delta_t) / (delta * delta)))) * (u_0[i][j] - (delta_t / 2.) * v_x[i][j] * (((u_0[0][j] - u_0[i - 1][j])/(delta * 2.)) + (u_1[0][j] - u_1[i - 1][j])/(2. * delta)) - (delta_t / 2.) * v_y[i][j] * ((u_0[i][j + 1] - u_0[i][j - 1])/(2. * delta) + (u_1[i][j + 1] - u_1[i][j - 1])/(2. * delta)) + (delta_t /2.) * D * ((u_0[0][j] + u_0[i - 1][j] + u_0[i][j + 1] + u_0[i][j - 1] - 4. * u_0[i][j])/(delta * delta) + (u_1[0][j] + u_1[i - 1][j] + u_1[i][j + 1] + u_1[i][j - 1])/(delta * delta)));
                        }
                        else
                        {
                            u_1[i][j] = (1. / (1. +((2. * D * delta_t ) / (delta * delta)))) * (u_0[i][j] - (delta_t / 2.) * v_x[i][j] * (((u_0[i + 1][j] - u_0[i - 1][j])/(delta * 2.)) + (u_1[i + 1][j] - u_1[i - 1][j])/(2. * delta)) - (delta_t / 2.) * v_y[i][j] * ((u_0[i][j + 1] - u_0[i][j - 1])/(2. * delta) + (u_1[i][j + 1] - u_1[i][j - 1])/(2. * delta)) + (delta_t / 2.) * D * ((u_0[i + 1][j] + u_0[i - 1][j] + u_0[i][j + 1] + u_0[i][j - 1] - 4. * u_0[i][j])/(delta * delta) + (u_1[i + 1][j] + u_1[i - 1][j] + u_1[i][j + 1] + u_1[i][j - 1])/(delta * delta)));
                        }
                    }
                }
            }
        }

        if(it == 4000) //zmiana od 800 do 4000 co 800 (5 map)
        {
            for(int i = 0; i <= n_x; i++)
            {
                for(int j = 0; j <= n_y; j++)
                {
                    fprintf(fp_1, "%d \t %d \t %g \n", i, j, u_1[i][j]);
                }
                fprintf(fp_1, "\n");
            }
        }

        for(int i = 0; i <= n_x; i++)
        {
            for(int j = 0; j <= n_y; j++)
            {
                u_0[i][j] = u_1[i][j];
            }
        }

        double c = 0.;
        double x_sr = 0.;

        for(int i = 0; i <= n_x; i++)
        {
            for(int j = 0; j <= n_y; j++)
            {
                c += delta * delta * u_0[i][j];
                x_sr += delta * delta * delta * i * u_0[i][j];
            }
        }

        fprintf(fp_2, "%g \t %g \t %g \n", delta_t * it, c, x_sr);

    }
}


int main()
{
    FILE* vx = fopen("vx.txt", "w");
    FILE* vy = fopen("vy.txt", "w");
    FILE* integral = fopen("integral1.txt", "w");
    FILE* read = fopen("psi.dat", "r");

    double psi[n_x + 1][n_y + 1];
    double v_x[n_x + 1][n_y + 1];
    double v_y[n_x + 1][n_y + 1];
    double u_0[n_x + 1][n_y + 1];
    double u_1[n_x + 1][n_y + 1];

    double temp;
    int i, j;

    while(fscanf(read, "%d \t %d \t %lf", &i, &j, &temp) == 3)
    {
        psi[i][j] = temp;
    }

    fclose(read);

    V_field(v_x, v_y, psi, vx, vy);
    fclose(vx);
    fclose(vy);

    double delta_t = delta / (4 * V_max(v_x, v_y));

    //FILE* map1_1 = fopen("map1_1.txt", "w");
    //FILE* map1_2 = fopen("map1_2.txt", "w");
    //FILE* map1_3 = fopen("map1_3.txt", "w");
    //FILE* map1_4 = fopen("map1_4.txt", "w");
    FILE* map1_5 = fopen("map1_5.txt", "w");

    fill_u(u_0);

    advection_diffusion(u_0, u_1, v_x, v_y, psi, map1_5, integral, 0., delta_t);

    //fclose(map1_1);
    //fclose(map1_2);
    //fclose(map1_3);
    //fclose(map1_4);
    fclose(map1_5);
    fclose(integral);

    FILE* integral1 = fopen("integral2.txt", "w");

    //FILE* map2_1 = fopen("map2_1.txt", "w");
    //FILE* map2_2 = fopen("map2_2.txt", "w");
    //FILE* map2_3 = fopen("map2_3.txt", "w");
    //FILE* map2_4 = fopen("map2_4.txt", "w");
    FILE* map2_5 = fopen("map2_5.txt", "w");


    fill_u(u_0);

    advection_diffusion(u_0, u_1, v_x, v_y, psi, map2_5, integral1, 0.1, delta_t);

    //fclose(map2_1);
    //fclose(map2_2);
    //fclose(map2_3);
    //fclose(map2_4);
    //fclose(map2_5);
    fclose(integral1);

    return 0;
}
