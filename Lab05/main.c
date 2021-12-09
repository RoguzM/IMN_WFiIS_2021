#include <stdio.h>
#include <math.h>
#include <stdlib.h>

const double delta = 0.2;
const int n_x = 128;
const int n_y = 128;
const double x_max = 25.6;
const double y_max = 25.6;
const double TOL = pow(10, -8);

void multi_grid_relax(int k, double V[n_x+1][n_y+1], FILE* fp_S, FILE* fp_V)
{
    double S_it = 0.,  S_it_1 = 0.;
    int iter = 0;

    do
    {
        for (int i = k; i < n_x - k + 1; i += k)
        {
            for (int j = k; j < n_y - k + 1; j += k)
            {
                V[i][j] = 0.25 * (V[i + k][j] + V[i - k][j] + V[i][j + k] + V[i][j - k]);
            }
        }

        S_it_1 = S_it;
        S_it = 0.;

        for(int i = 0; i < n_x - k + 1; i += k)
        {
            for(int j = 0; j < n_y - k + 1; j += k)
            {
                S_it += (0.5 * k * k * delta * delta) * (pow(((V[i + k][j] - V[i][j])/(2 * k * delta) + (V[i + k][j + k] - V[i][j + k])/(2 * k * delta)), 2) + pow(((V[i][j + k] - V[i][j]) / (2 * k * delta) + (V[i + k][j + k] - V[i + k][j]) / (2 * k * delta)), 2));
            }
        }

        fprintf(fp_S, "%d \t %f \n", iter, S_it);
        iter++;

    }while(fabs((S_it - S_it_1)/S_it_1) > TOL);


    for (int i =0; i < n_x + 1; i += k)
    {
        for (int j = 0; j < n_y + 1; j += k)
        {
            fprintf(fp_V, "%f \t %f \t %f \n", i * delta, j * delta, V[i][j]);
        }
    }

    if (k != 1)
    {
        for (int i = 0; i < n_x - k + 1; i += k)
        {
            for (int j = 0; j < n_y - k + 1; j += k)
            {
                V[i + k / 2][j + k / 2] = (V[i][j] + V[i + k][j] + V[i][j + k] + V[i + k][j + k]) / 4.;

                if (i != n_x - k)
                {
                    V[i + k][j + k / 2] = (V[i + k][j] + V[i + k][j + k]) / 2.;
                }

                if (i != 0)
                {
                    V[i][j + k / 2] = (V[i][j] + V[i][j + k]) / 2.;
                }

                if (j != n_y - k)
                {
                    V[i + k / 2][j + k] = (V[i][j + k] + V[i + k][j + k]) / 2.;
                }

                if (j != 0)
                {
                    V[i + k / 2][j] = (V[i][j] + V[i + k][j]) / 2.;
                }
            }
        }
    }
}

int main()
{
    double V[n_x + 1][n_y + 1];

    for(int i = 0; i < n_x + 1; i++)
    {
        for(int j =0; j < n_y + 1; j++)
        {
            V[i][j] = 0.;
        }
    }

    for(int i = 0; i < n_y + 1; i++)
    {
        V[0][i] = sin(M_PI * i * delta / y_max);
        V[n_x][i] = sin(M_PI * i * delta / y_max);
    }

    for(int i = 0; i < n_x + 1; i++)
    {
        V[i][0] = sin(2. * M_PI * i * delta/ x_max);
        V[i][n_y] = -sin(2. * M_PI * i * delta / x_max);
    }

    FILE* multi_grid_S_16 = fopen("multi_grid_S_16.txt", "w");
    FILE* multi_grid_V_16 = fopen("multi_grid_V_16.txt", "w");
    multi_grid_relax(16, V, multi_grid_S_16, multi_grid_V_16);
    fclose(multi_grid_S_16);
    fclose(multi_grid_V_16);

    FILE* multi_grid_S_8 = fopen("multi_grid_S_8.txt", "w");
    FILE* multi_grid_V_8 = fopen("multi_grid_V_8.txt", "w");
    multi_grid_relax(8, V, multi_grid_S_8, multi_grid_V_8);
    fclose(multi_grid_S_8);
    fclose(multi_grid_V_8);

    FILE* multi_grid_S_4 = fopen("multi_grid_S_4.txt", "w");
    FILE* multi_grid_V_4 = fopen("multi_grid_V_4.txt", "w");
    multi_grid_relax(4, V, multi_grid_S_4, multi_grid_V_4);
    fclose(multi_grid_S_4);
    fclose(multi_grid_V_4);

    FILE* multi_grid_S_2 = fopen("multi_grid_S_2.txt", "w");
    FILE* multi_grid_V_2 = fopen("multi_grid_V_2.txt", "w");
    multi_grid_relax(2, V, multi_grid_S_2, multi_grid_V_2);
    fclose(multi_grid_S_2);
    fclose(multi_grid_V_2);

    FILE* multi_grid_S_1 = fopen("multi_grid_S_1.txt", "w");
    FILE* multi_grid_V_1 = fopen("multi_grid_V_1.txt", "w");
    multi_grid_relax(1, V, multi_grid_S_1, multi_grid_V_1);
    fclose(multi_grid_S_1);
    fclose(multi_grid_V_1);

    return 0;
}

