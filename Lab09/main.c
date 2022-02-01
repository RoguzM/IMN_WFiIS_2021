#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>

 const int n_x = 40;
 const int n_y = 40;
 const int N = (n_x + 1) * (n_y + 1);
 const double delta = 1.;
 const double dt = 1.;
 const int T_A = 40;
 const int T_B = 0;
 const int T_C = 30;
 const int T_D = 0;
 const double k_B = 0.1;
 const double k_D = 0.6;
 const int IT_MAX = 2000;

int reindex(int i, int j)
{
    return i + j * (n_x + 1);
}

void fill(gsl_matrix* A, gsl_matrix* B,  gsl_vector* c, gsl_vector* T)
{
    for(int i = 1; i <= n_x - 1; i++)
    {
        for(int j = 1; j <= n_y - 1; j++)
        {
            gsl_matrix_set(A, reindex(i,j), reindex(i,j) - n_x -1, dt / (2. * delta * delta));
            gsl_matrix_set(A, reindex(i,j), reindex(i,j) - 1, dt / (2. * delta * delta));
            gsl_matrix_set(A, reindex(i,j), reindex(i,j) + 1 , dt / (2. * delta * delta));
            gsl_matrix_set(A, reindex(i,j), reindex(i,j) + n_x + 1, dt / (2. * delta * delta));
            gsl_matrix_set(A, reindex(i,j), reindex(i,j), -2. * dt / (delta * delta) - 1);
            gsl_matrix_set(B, reindex(i,j), reindex(i,j) - n_x - 1, -dt / (2. * delta * delta));
            gsl_matrix_set(B, reindex(i,j), reindex(i,j) - 1, -dt / (2. * delta * delta));
            gsl_matrix_set(B, reindex(i,j), reindex(i,j) + 1, -dt / (2. *delta * delta));
            gsl_matrix_set(B, reindex(i,j), reindex(i,j) + n_x + 1, -dt / (2. * delta * delta));
            gsl_matrix_set(B, reindex(i,j), reindex(i, j), 2. * dt / (delta * delta) - 1);
        }
    }

    for(int i = 0; i <= n_x; i += n_x)
    {
        for(int j = 0; j <= n_y; j++)
        {
            gsl_matrix_set(A, reindex(i, j), reindex(i, j), 1);
            gsl_matrix_set(B, reindex(i, j), reindex(i, j), 1);
            gsl_vector_set(c, reindex(i, j), 0);
        }
    }

    for(int i = 0; i <= n_x - 1; i++)
    {
        gsl_matrix_set(A, reindex(i, n_y), reindex(i, n_y) - n_x - 1, -1./(k_B * delta));
        gsl_matrix_set(A, reindex(i, n_y), reindex(i, n_y), 1 + 1./(k_B * delta));
        gsl_vector_set(c, reindex(i, n_y), T_B);

        for(int j = 0; j < N; j++)
        {
            gsl_matrix_set(B, reindex(i, n_y), j, 0);
        }
    }

    for(int i = 1; i <= n_x - 1; i++)
    {
        gsl_matrix_set(A, reindex(i, 0), reindex(i, 0), 1 + 1. / (k_D * delta));
        gsl_matrix_set(A, reindex(i, 0), reindex(i, 0) + n_x + 1, -1. / (k_D * delta));
        gsl_vector_set(c, reindex(i, 0), T_D);

        for(int j = 0; j < N; j++)
        {
            gsl_matrix_set(B, reindex(i, 0), j, 0);
        }
    }

    for(int j = 0; j <= n_y; j++)
    {
        gsl_vector_set(T, reindex(0, j), T_A);
        gsl_vector_set(T, reindex(n_x, j), T_C);
    }

    for(int i = 1; i <= n_x - 1; i++)
    {
        for(int j = 0; j <= n_y; j++)
        {
            gsl_vector_set(T, reindex(i, j), 0);
        }
    }
}

double differential_T(gsl_vector* T, int l)
{
    return ((gsl_vector_get(T, l + 1) - 2 * gsl_vector_get(T, l) + gsl_vector_get(T, l - 1)) / (delta * delta))
        + ((gsl_vector_get(T, l + n_x + 1) - 2 * gsl_vector_get(T, l) + gsl_vector_get(T, l - n_x - 1)) / (delta * delta));
}


int main()
{
    gsl_matrix* A = gsl_matrix_calloc(N, N);
    gsl_matrix* B = gsl_matrix_calloc(N,N);
    gsl_vector* c = gsl_vector_calloc(N);
    gsl_vector* T = gsl_vector_calloc(N);

    fill(A, B, c, T);

    gsl_permutation *p = gsl_permutation_alloc(N);
    int signum = 0;

    gsl_linalg_LU_decomp(A, p, &signum);

    FILE* fp_T100 = fopen("T100.txt", "w");
    FILE* fp_2T100 = fopen("2T100.txt", "w");
    FILE* fp_T200 = fopen("T200.txt", "w");
    FILE* fp_2T200 = fopen("2T200.txt", "w");
    FILE* fp_T500 = fopen("T500.txt", "w");
    FILE* fp_2T500 = fopen("2T500.txt", "w");
    FILE* fp_T1000 = fopen("T1000.txt", "w");
    FILE* fp_2T1000 = fopen("2T1000.txt", "w");
    FILE* fp_T2000 = fopen("T2000.txt", "w");
    FILE* fp_2T2000 = fopen("2T2000.txt", "w");

    gsl_vector* d = gsl_vector_calloc(N);

    for(int it = 0; it <= IT_MAX; it++)
    {
        gsl_blas_dgemv(CblasNoTrans, 1, B, T, 0, d);
        gsl_blas_daxpy(1, c, d);
        gsl_linalg_LU_solve(A, p, d, T);

        if(it == 100 || it == 200 || it == 500 || it == 1000 || it == 2000)
        {
            for(int i = 1; i <= n_x - 1; i++)
            {
                for(int j = 1; j <= n_y - 1; j++)
                {
                    if(it == 100)
                    {
                        fprintf(fp_T100, "%.12lf \t %.12lf \t %.12lf \n", i * delta, j * delta, gsl_vector_get(T, reindex(i, j)));
                        fprintf(fp_2T100, "%.12lf \t %.12lf \t %.12lf \n", i * delta, j * delta, differential_T(T, reindex(i, j)));
                    }
                    if(it == 200)
                    {
                        fprintf(fp_T200, "%.12lf \t %.12lf \t %.12lf \n", i * delta, j * delta, gsl_vector_get(T, reindex(i, j)));
                        fprintf(fp_2T200, "%.12lf \t %.12lf \t %.12lf \n", i * delta, j * delta, differential_T(T, reindex(i, j)));
                    }
                    if(it == 500)
                    {
                        fprintf(fp_T500, "%.12lf \t %.12lf \t %.12lf \n", i * delta, j * delta, gsl_vector_get(T, reindex(i, j)));
                        fprintf(fp_2T500, "%.12lf \t %.12lf \t %.12lf \n", i * delta, j * delta, differential_T(T, reindex(i, j)));
                    }
                    if(it == 1000)
                    {
                        fprintf(fp_T1000, "%.12lf \t %.12lf \t %.12lf \n", i * delta, j * delta, gsl_vector_get(T, reindex(i, j)));
                        fprintf(fp_2T1000, "%.12lf \t %.12lf \t %.12lf \n", i * delta, j * delta, differential_T(T, reindex(i, j)));
                    }
                    if(it == 2000)
                    {
                        fprintf(fp_T2000, "%.12lf \t %.12lf \t %.12lf \n", i * delta, j * delta, gsl_vector_get(T, reindex(i, j)));
                        fprintf(fp_2T2000, "%.12lf \t %.12lf \t %.12lf \n", i * delta, j * delta, differential_T(T, reindex(i, j)));
                    }
                }
                if(it == 100)
                    {
                        fprintf(fp_T100, "\n");
                        fprintf(fp_2T100, "\n");
                    }
                    if(it == 200)
                    {
                        fprintf(fp_T200, "\n");
                        fprintf(fp_2T200, "\n");
                    }
                    if(it == 500)
                    {
                        fprintf(fp_T500, "\n");
                        fprintf(fp_2T500, "\n");
                    }
                    if(it == 1000)
                    {
                        fprintf(fp_T1000, "\n");
                        fprintf(fp_2T1000, "\n");
                    }
                    if(it == 2000)
                    {
                        fprintf(fp_T2000, "\n");
                        fprintf(fp_2T2000, "\n");
                    }
            }
        }
    }

    fclose(fp_T100);
    fclose(fp_2T100);
    fclose(fp_T200);
    fclose(fp_2T200);
    fclose(fp_T500);
    fclose(fp_2T500);
    fclose(fp_T1000);
    fclose(fp_2T1000);
    fclose(fp_T2000);
    fclose(fp_2T2000);

    gsl_permutation_free(p);
    gsl_matrix_free(A);
    gsl_matrix_free(B);
    gsl_vector_free(c);
    gsl_vector_free(T);
    gsl_vector_free(d);

    return 0;
}
