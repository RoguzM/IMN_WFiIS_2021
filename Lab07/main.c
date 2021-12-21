#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>

const double delta = 0.01;
const double ro = 1.;
const double mi = 1.;
const int n_x = 200;
const int n_y = 90;
const int i_1 = 50;
const int j_1 = 55;
const int IT_MAX = 20000;

double fill_psi(int i, int j, double psi[n_x + 1][n_y + 1], double zeta[n_x + 1][n_y + 1])
{
    return (psi[i + 1][j] + psi[i - 1][j] + psi[i][j + 1] + psi[i][j - 1] - zeta[i][j] * pow(delta, 2))/4.;
}

double fill_zeta(int i, int j, double omega, double psi[n_x + 1][n_y + 1], double zeta[n_x + 1][n_y + 1])
{
    return (zeta[i + 1][j] + zeta[i - 1][j] + zeta[i][j + 1] + zeta[i][j - 1])/4. - omega * ro / (16 * mi) * ((psi[i][j + 1] - psi[i][j - 1]) * (zeta[i + 1][j] - zeta[i - 1][j]) - (psi[i + 1][j] - psi[i - 1][j]) * (zeta[i][j + 1] - zeta[i][j - 1]));
}

double calculate_Q_wy(double Q_we, double y_ny, double y_j1)
{
    return Q_we * (pow(y_ny, 3) - pow(y_j1, 3) - 3 * y_j1 * pow(y_ny, 2) + 3 * pow(y_j1, 2) * y_ny) / pow(y_ny, 3);
}

void boundary_psi(double Q_we, double Q_wy, double psi[n_x + 1][n_y + 1], double y[n_y + 1])
{
    for(int j = j_1; j <= n_y; j++)
    {
        psi[0][j] = Q_we/(2 * mi) * (pow(y[j], 3)/3. - pow(y[j], 2)/2. * (y[j_1] + y[n_y]) + y[j] * y[j_1] * y[n_y]);
    }

    for(int j = 0; j <= n_y; j++)
    {
        psi[n_x][j] = Q_wy/ (2 * mi) * (pow(y[j], 3)/3. - pow(y[j], 2)/2. * y[n_y]) + (Q_we * pow(y[j_1], 2) * (-y[j_1] + 3 * y[n_y]))/ (12 * mi);
    }

    for(int i = 1; i <= n_x - 1; i++)
    {
        psi[i][n_y] = psi[0][n_y];
    }

    for(int i = i_1; i <= n_x - 1; i++)
    {
        psi[i][0] = psi[0][j_1];
    }

    for(int j = 1; j <= j_1; j++)
    {
        psi[i_1][j] = psi[0][j_1];
    }

    for(int i = 1; i <= i_1; i++)
    {
        psi[i][j_1] = psi[0][j_1];
    }
}

void boundary_zeta(double Q_we, double Q_wy, double psi[n_x + 1][n_y + 1], double zeta[n_x + 1][n_y + 1], double y[n_y + 1])
{
    for(int j = j_1; j <= n_y; j++)
    {
        zeta[0][j] = Q_we/ (2 * mi) * (2 * y[j] - y[j_1] - y[n_y]);
    }

    for(int j = 0; j <= n_y; j++)
    {
        zeta[n_x][j] = Q_wy / (2 * mi) * (2 * y[j] - y[n_y]);
    }

    for(int i = 1; i <= n_x - 1; i++)
    {
        zeta[i][n_y] = 2. / (delta * delta) * (psi[i][n_y - 1] - psi[i][n_y]);
    }

    for(int i = i_1 + 1; i <= n_x -1; i++)
    {
        zeta[i][0] = 2. / (delta * delta) * (psi[i][1] - psi[i][0]);
    }

    for(int j = 1; j <= j_1 - 1; j++)
    {
        zeta[i_1][j] = 2. / (delta * delta) * (psi[i_1 + 1][j] - psi[i_1][j]);
    }

    for(int i = 1; i <= i_1; i++)
    {
        zeta[i][j_1] = 2. / (delta * delta) * (psi[i][j_1 + 1] - psi[i][j_1]);
    }

    zeta[i_1][j_1] = (zeta[i_1 - 1][j_1] + zeta[i_1][j_1 - 1])/2.;
}

bool border_check(int i, int j)
{
    return (i == 0) || (i == n_x) || (j == 0) || (j == n_y) || (i <= i_1 && j == j_1) || (i == i_1 && j <= j_1);
}

bool inside_check(int i, int j)
{
    return !(i <= i_1 && j < j_1);
}

void nav_stokes(double Q_we, FILE* fp)
{
    double psi[n_x + 1][n_y + 1];
    double zeta[n_x + 1][n_y + 1];
    double y[n_y + 1];

    for(int i = 0; i <= n_x ; i++)
    {
        for(int j = 0; j <= n_y; j++)
        {
            psi[i][j] = 0;
            zeta[i][j] = 0;
            y[j] = delta * j;
        }
    }

    double Q_wy = calculate_Q_wy(Q_we, y[n_y], y[j_1]);

    int omega;
    bool border, inside;
    double error_control;
    int j_2 = j_1 + 2;

    boundary_psi(Q_we, Q_wy, psi, y);

    for(int it = 1; it < IT_MAX; it++)
    {
        if(it < 2000)
        {
            omega = 0;
        }
        else
        {
            omega = 1;
        }

        for(int i = 1; i <= n_x - 1; i++)
        {
            for(int j = 1; j <= n_y - 1; j++)
            {
                border = border_check(i, j);
                inside = inside_check(i, j);

                if(!border && inside)
                {
                    psi[i][j] = fill_psi(i, j, psi, zeta);
                    zeta[i][j] = fill_zeta(i, j, omega, psi, zeta);
                }
            }
        }

    boundary_zeta(Q_we, Q_wy, psi, zeta, y);
    error_control = 0.;

    for(int i = 1; i <= n_x - 1; i++)
    {
        error_control += (psi[i + 1][j_2] + psi[i - 1][j_2] + psi[i][j_2 + 1] + psi[i][j_2 - 1] - 4 * psi[i][j_2] - delta * delta * zeta[i][j_2]);
    }

    printf("It: %d \t Error_control: %f \n", it, error_control);

    }

    double u;
    double v;

    for(int i = 1; i <= n_x - 1; i++)
    {
        for(int j = 1; j <= n_y - 1; j++)
        {
            border = border_check(i, j);
            inside = inside_check(i, j);

            if(!border && inside)
            {
                u = (psi[i][j + 1] - psi[i][j - 1])/ (2 * delta);
                v = -(psi[i + 1][j] - psi[i - 1][j]) / (2 * delta);
            }
            else
            {
                u = 0;
                v = 0;
            }

            fprintf(fp, "%f \t %f \t %f \t %f \t %f \t %f \n", i * delta, j* delta, psi[i][j], zeta[i][j], u, v);
        }
        fprintf(fp, "\n");
    }

}



int main()
{
    FILE* Q_minus_1000 = fopen("Q_minus_1000.txt", "w");
    nav_stokes(-1000, Q_minus_1000);
    fclose(Q_minus_1000);

    FILE* Q_minus_4000 = fopen("Q_minus_4000.txt", "w");
    nav_stokes(-4000, Q_minus_4000);
    fclose(Q_minus_4000);

    FILE* Q_4000 = fopen("Q_4000.txt", "w");
    nav_stokes(4000, Q_4000);
    fclose(Q_4000);


    return 0;
}
