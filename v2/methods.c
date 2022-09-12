#include "methods.h"

// --------------------------- Funções Internas ---------------------------
/*!
    \brief Função que encontra o valor máximo da coluna
    \param matrix Ponteiro para a Matriz
    \param line linha atual da execução da fatoração LU
    \return linha do pivo de maior valor
*/
int find_max (matrix_t *matrix, int line) {
    int count, max_line;
    double iterator, value;
    iterator = fabs (matrix->coef[line][line]);
    max_line = line;

    for (count = line + 1; count < matrix->n; count++) {
        value = fabs (matrix->coef[count][line]);
        if (value > iterator) {
            max_line = count;
            iterator = value;
        }
    }
    return max_line;
}

/*!
    \brief função que calcula o determiante e retorna o seu valor
    \param matrix Matriz que vai ser calculada o determinante
    \returns valor do determinante
*/
double calc_determinant (matrix_t *matrix) {
    int count, size;
    int unroll_limit;
    double determinant = 1.0;
    size = matrix->n;
    unroll_limit = size - (size % UNROLL_SIZE);

    for (count = 0; count < unroll_limit; count += UNROLL_SIZE) 
        determinant = determinant * matrix->coef[count][count] * matrix->coef[count + 1][count + 1] * matrix->coef[count + 2][count + 2] * matrix->coef[count + 3][count + 3];
    for (count = unroll_limit; count < size; count++) 
        determinant = determinant * matrix->coef[count][count];

    return determinant;
}

/*!
    \brief Função que calcula a matrix de residuos do resultado
    \param residue_matrix Ponteiro para a Matriz de Residuos (Retorno)
    \param matrix Ponteiro para Matriz enviada pelo usuário
    \param inv_matrix Ponteiro para a Matriz inversa calculada
    \return 0 em caso de sucesso, ALLOC_ERROR em caso de falha de alocação
*/
void calc_residue (matrix_t *restrict residue_matrix, matrix_t *restrict matrix, matrix_t *restrict inv_matrix) {
    int istart = 0, iend = 0, jstart = 0, jend = 0, kstart = 0, kend = 0;
    int ii, jj, kk, i, j, k;
    int size, loops, unroll_limit;
    __m256d avx_mul;
    __m128d sum_aux_1, sum_aux_2, final_sum;
    size = matrix->n;
    generate_identity_matrix (residue_matrix);

    loops = size / BLOCK_SIZE;
    unroll_limit = size - (size % UNROLL_SIZE);
    for (ii = 0; ii < loops; ii++) {
        istart = ii * BLOCK_SIZE;
        iend = istart + BLOCK_SIZE;
        for (jj = 0; jj < loops; jj++) {
            jstart = jj * BLOCK_SIZE;
            jend = jstart + BLOCK_SIZE;
            for (kk = 0; kk < loops; kk++) {
                kstart = kk * BLOCK_SIZE;
                kend = kstart + BLOCK_SIZE;
                for (i = istart; i < iend; i++) {
                    for (j = jstart; j < jend; j += UNROLL_SIZE) {
                        for (k = kstart; k < kend; k++) {
                            avx_mul = _mm256_mul_pd (_mm256_loadu_pd (&matrix->coef[i][k]), _mm256_loadu_pd (&inv_matrix->coef[k][j]));
                            sum_aux_1 = _mm256_extractf128_pd (avx_mul, 0);
                            sum_aux_2 = _mm256_extractf128_pd (avx_mul, 1);
                            final_sum = _mm_add_pd (sum_aux_1, sum_aux_2);
                            residue_matrix->coef[i][j] -= _mm_cvtsd_f64 (final_sum);
                        }
                    }
                }
            }
            for (i = istart; i < iend; i++) {
                for (j = jstart; j < jend; j += UNROLL_SIZE) {
                    for (k = kend; k < size; k++) {
                        avx_mul = _mm256_mul_pd (_mm256_loadu_pd (&matrix->coef[i][k]), _mm256_loadu_pd (&inv_matrix->coef[k][j]));
                        sum_aux_1 = _mm256_extractf128_pd (avx_mul, 0);
                        sum_aux_2 = _mm256_extractf128_pd (avx_mul, 1);
                        final_sum = _mm_add_pd (sum_aux_1, sum_aux_2);
                        residue_matrix->coef[i][j] -= _mm_cvtsd_f64 (final_sum);
                    }
                }
            }
        }
        for (i = iend; i < size; i++) {
            for (j = jend; j < unroll_limit; j += 4) {
                for (k = kstart; k < kend; k++) {
                    avx_mul = _mm256_mul_pd (_mm256_loadu_pd (&matrix->coef[i][k]), _mm256_loadu_pd (&inv_matrix->coef[k][j]));
                    sum_aux_1 = _mm256_extractf128_pd (avx_mul, 0);
                    sum_aux_2 = _mm256_extractf128_pd (avx_mul, 1);
                    final_sum = _mm_add_pd (sum_aux_1, sum_aux_2);
                    residue_matrix->coef[i][j] -= _mm_cvtsd_f64 (final_sum);
                }
            }
            for (j = unroll_limit; j < size; j++)
                for (k = kstart; k < kend; k++) 
                    residue_matrix->coef[i][j] -= matrix->coef[i][k] * inv_matrix->coef[k][j];
        }
    }
    for (i = iend; i < size; i++) {
         for (j = jend; j < unroll_limit; j += 4) {
            for (k = kend; k < size; k++) {
                residue_matrix->coef[i][j] -= matrix->coef[i][k] * inv_matrix->coef[k][j];
                residue_matrix->coef[i][j + 1] -= matrix->coef[i][k] * inv_matrix->coef[k][j + 1];
                residue_matrix->coef[i][j + 2] -= matrix->coef[i][k] * inv_matrix->coef[k][j + 2];
                residue_matrix->coef[i][j + 3] -= matrix->coef[i][k] * inv_matrix->coef[k][j + 3];
            }
        }
        for (j = unroll_limit; j < size; j++)
            for (k = kend; k < size; k++)
                residue_matrix->coef[i][j] -= matrix->coef[i][k] * inv_matrix->coef[k][j];
    }
}

/*!
    \brief Função que calcula a norma L2 do resíduo
    \param residue_matrix Ponteiro para a Matriz de resíduos
    \param residue_time tempo para calcular a norma do resíduo
    \return o valor da norma
*/
double calc_norma (matrix_t *residue_matrix, double *residue_time) {
    int line, col, size;
    int unroll_limit;
    double norma = 0.0;

    size = residue_matrix->n;
    unroll_limit = size - (size % UNROLL_SIZE);

    *residue_time = timestamp ();

    for (line = 0; line < unroll_limit; line += UNROLL_SIZE) {
        for (col = 0; col < size; col++) {
            norma += residue_matrix->coef[line][col] * residue_matrix->coef[line][col];
            norma += residue_matrix->coef[line + 1][col] * residue_matrix->coef[line + 1][col];
            norma += residue_matrix->coef[line + 2][col] * residue_matrix->coef[line + 2][col];
            norma += residue_matrix->coef[line + 3][col] * residue_matrix->coef[line + 3][col];
        }
    }
    for (line = unroll_limit; line < size; line++) 
        for (col = 0; col < size; col++) 
            norma += residue_matrix->coef[line][col] * residue_matrix->coef[line][col];

    *residue_time = timestamp () - *residue_time;
    return sqrt (norma);
}

// --------------------------- Funções Externas ---------------------------
int lu_factorization (matrix_t *restrict u_matrix, matrix_t *restrict l_matrix, pivot_steps_t *pivot_steps, double *lu_time) {
    int line, i_line, col, pivot, size;
    double m;

    size = u_matrix->n;
    memset (l_matrix->coef[0], 0.0, sizeof (double) * size * size);

    *lu_time = timestamp ();

    for (line = 0; line < size; line++) {
        pivot = find_max (u_matrix, line);
        if (line != pivot) {
            swap_line (u_matrix, line, pivot);
            swap_line (l_matrix, line, pivot);
            add_pivot_step (pivot_steps, line, pivot);
        }

        l_matrix->coef[line][line] = 1.0;
        for (i_line = line + 1; i_line < size; i_line++) {
            m = u_matrix->coef[i_line][line] / u_matrix->coef[line][line];
            // if (isnan (m) || isinf (m)) 
            //     return NAN_INF_ERROR;
            
            u_matrix->coef[i_line][line] = 0.0;
            l_matrix->coef[i_line][line] = m;
            for (col = line + 1; col < size; col++) 
                u_matrix->coef[i_line][col] -= u_matrix->coef[line][col] * m;
        }
    }
    *lu_time = timestamp () - *lu_time;

    if (fabs (calc_determinant (u_matrix) - 0) < DBL_EPSILON)
        return INVERTIBLE_ERROR;

    return EXIT_SUCCESS;
}

int calc_inverse_matrix (matrix_t *restrict inv_matrix, matrix_t *restrict l_matrix, matrix_t *restrict u_matrix, pivot_steps_t *steps) {
    matrix_t *i_matrix, *solution;
    double *temp_sol, *aux_avx;
    int count, sol_count; 
    int size, line, col;
    int unroll_limit, retro_limit, line_size;
    __m256d sol_avx;

    size = inv_matrix->n;
    unroll_limit = size - (size % UNROLL_SIZE);

    if (! (i_matrix = alloc_matrix (size)) || ! (solution = alloc_matrix (size)))
        return ALLOC_ERROR;
    if (! (temp_sol = malloc (sizeof (double) * size))) 
        return ALLOC_ERROR;
    memset (temp_sol, 0.0, sizeof (double) * size);
    memset (solution->coef[0], 0.0, sizeof (double) * size * size);
    generate_identity_matrix (i_matrix);

    LIKWID_MARKER_START ("retrosubs");
    apply_pivot_steps (i_matrix, steps);
    for (count = 0; count < size; count++) {
        // Retrosubstituição inversa
        for (line = 0; line < size; line++) {
            temp_sol[line] = i_matrix->coef[line][count];
            for (col = line - 1; col >= 0; col--)
                temp_sol[line] -= l_matrix->coef[line][col] * temp_sol[col];
            temp_sol[line] /= l_matrix->coef[line][line];
            // if (isnan (solution[line]) || isinf (solution[line]))
            //     return NAN_INF_ERROR;
        }

        // adição da solução temporaria da execução do refinamento
        for (sol_count = 0; sol_count < unroll_limit; sol_count += UNROLL_SIZE) {
            sol_avx = _mm256_add_pd (_mm256_loadu_pd (&solution->coef[count][sol_count]), _mm256_loadu_pd (&temp_sol[sol_count]));
            aux_avx = (double *)&sol_avx;
            solution->coef[count][sol_count] = aux_avx[sol_count];
            solution->coef[count][sol_count + 1] = aux_avx[sol_count + 1];
            solution->coef[count][sol_count + 2] = aux_avx[sol_count + 2];
            solution->coef[count][sol_count + 3] = aux_avx[sol_count + 3];
        }
        for (sol_count = unroll_limit; sol_count < size; sol_count++)
            solution->coef[count][sol_count] += temp_sol[sol_count];

        // Retrosubstituição
        for (line = size - 1; line >= 0; line--) {
            line_size = size - line - 1;
            retro_limit = size - (linesize % 4);
            
            inv_matrix->coef[line][count] = solution->coef[count][line];
            for (col = line + 1; col < retro_limit; col += 4) {
                inv_matrix->coef[line][count] -= u_matrix->coef[line][col] * inv_matrix->coef[col][count];
                inv_matrix->coef[line][count] -= u_matrix->coef[line][col + 1] * inv_matrix->coef[col + 1][count];
                inv_matrix->coef[line][count] -= u_matrix->coef[line][col + 2] * inv_matrix->coef[col + 2][count];
                inv_matrix->coef[line][count] -= u_matrix->coef[line][col + 3] * inv_matrix->coef[col + 3][count];
            }
            for (col = retro_limit; col < size; col++) 
                inv_matrix->coef[line][count] -= u_matrix->coef[line][col] * inv_matrix->coef[col][count];
            inv_matrix->coef[line][count] /= u_matrix->coef[line][line];
            // if (isnan (inv_matrix->coef[line][count]) || isinf (inv_matrix->coef[line][count]))
            //     return NAN_INF_ERROR;
        }
    }
    LIKWID_MARKER_STOP ("retrosubs");
    
    free_matrix (solution);
    free_matrix (i_matrix);
    free (temp_sol);
    return EXIT_SUCCESS;
}

int matrix_refinement (matrix_t *restrict inv_matrix, matrix_t *restrict matrix, matrix_t *restrict l_matrix, matrix_t *restrict u_matrix, pivot_steps_t *steps, 
                       double *restrict norma_vet, int iterations, double *restrict iter_time, double *restrict residue_time) {
    matrix_t *residue_matrix, *solution;
    double act_iter_time, act_residue_time;
    double *temp_sol, *aux_avx;
    int size, count, ls_count, sol_count;
    int line, col;
    int unroll_limit, retro_limit, line_size;
     __m256d sol_avx;

    *iter_time = 0;
    *residue_time = 0;

    size = matrix->n;
    unroll_limit = size - (size % UNROLL_SIZE);

    // Alocação e inicialização dos vetores/Matrizes internas
    if (! (solution = alloc_matrix (size)) || ! (residue_matrix = alloc_matrix (size)))
        return ALLOC_ERROR;
    if (! (temp_sol = malloc (sizeof (double) * size))) 
        return ALLOC_ERROR;
    memset (temp_sol, 0.0, sizeof (double) * size);
    memset (solution->coef[0], 0.0, sizeof (double) * size * size);

    // passos do refinamento
    for (count = 0; count < iterations; count++) {
        LIKWID_MARKER_START ("calc-residue");
        calc_residue (residue_matrix, matrix, inv_matrix);
        LIKWID_MARKER_STOP ("calc-residue");
        
        apply_pivot_steps (residue_matrix, steps);

        act_iter_time = timestamp ();
        for (ls_count = 0; ls_count < size; ls_count++) {
            // Retrosubstituição inversa
            for (line = 0; line < size; line++) {
                temp_sol[line] = residue_matrix->coef[line][ls_count];
                for (col = line - 1; col >= 0; col--)
                    temp_sol[line] -= l_matrix->coef[line][col] * temp_sol[col];
                temp_sol[line] /= l_matrix->coef[line][line];
                // if (isnan (solution[line]) || isinf (solution[line]))
                //     return NAN_INF_ERROR;
            }

            // adição da solução temporaria da execução do refinamento
            for (sol_count = 0; sol_count < unroll_limit; sol_count += UNROLL_SIZE) {
                sol_avx = _mm256_add_pd (_mm256_loadu_pd (&solution->coef[ls_count][sol_count]), _mm256_loadu_pd (&temp_sol[sol_count]));
                aux_avx = (double *)&sol_avx;
                solution->coef[ls_count][sol_count] = aux_avx[sol_count];
                solution->coef[ls_count][sol_count + 1] = aux_avx[sol_count + 1];
                solution->coef[ls_count][sol_count + 2] = aux_avx[sol_count + 2];
                solution->coef[ls_count][sol_count + 3] = aux_avx[sol_count + 3];
            }
            for (sol_count = unroll_limit; sol_count < size; sol_count++)
                solution->coef[ls_count][sol_count] += temp_sol[sol_count];

            // Retrosubstituição
            for (line = size - 1; line >= 0; line--) {
                line_size = size - line - 1;
                retro_limit = size - (linesize % 4);
                
                inv_matrix->coef[line][ls_count] = solution->coef[ls_count][line];
                for (col = line + 1; col < retro_limit; col += 4) {
                    inv_matrix->coef[line][ls_count] -= u_matrix->coef[line][col] * inv_matrix->coef[col][ls_count];
                    inv_matrix->coef[line][ls_count] -= u_matrix->coef[line][col + 1] * inv_matrix->coef[col + 1][ls_count];
                    inv_matrix->coef[line][ls_count] -= u_matrix->coef[line][col + 2] * inv_matrix->coef[col + 2][ls_count];
                    inv_matrix->coef[line][ls_count] -= u_matrix->coef[line][col + 3] * inv_matrix->coef[col + 3][ls_count];
                }
                for (col = retro_limit; col < size; col++) 
                    inv_matrix->coef[line][ls_count] -= u_matrix->coef[line][col] * inv_matrix->coef[col][ls_count];
                inv_matrix->coef[line][ls_count] /= u_matrix->coef[line][line];
                // if (isnan (inv_matrix->coef[line][ls_count]) || isinf (inv_matrix->coef[line][ls_count]))
                //     return NAN_INF_ERROR;
            }
        }
        act_iter_time = timestamp () - act_iter_time;

        norma_vet[count] = calc_norma (residue_matrix, &act_residue_time);

        *residue_time += act_residue_time;
        *iter_time += act_iter_time;
    }

    *iter_time = *iter_time / count;
    *residue_time = *residue_time / count;

    free_matrix (solution);
    free_matrix (residue_matrix);
    free (temp_sol);
    return EXIT_SUCCESS;
}