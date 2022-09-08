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
    double iterator;
    iterator = ABS (matrix->coef[line][line]);
    max_line = line;

    for (count = line + 1; count < matrix->n; count++) {
        max_line = iterator >= ABS (matrix->coef[count][line]) ? max_line : count;
        iterator = iterator >= ABS (matrix->coef[count][line]) ? iterator : ABS (matrix->coef[count][line]);
    }
    return max_line;
}

/*!
    \brief função que calcula o determiante e retorna o seu valor
    \param matrix Matriz que vai ser calculada o determinante
    \returns valor do determinante
*/
double calc_determinant (matrix_t *matrix) {
    int count;
    double determinant = 1.0;

    for (count = 0; count < matrix->n; count++)
        determinant = determinant * matrix->coef[count][count];

    return determinant;
}

/*!
    \brief Função que calcula o valor do sistema linear de uma matriz triangular inferior
    \param l_matrix Ponteiro para a Matriz L
    \param b_matrix Ponteiro para a Matriz de termos independentes
    \param solution vetor com a solução da retrosubstituição
    \param ls_b coluna da Matriz identidade (B do sistema linear)
    \return 0 em caso de sucesso. NAN_INF_ERROR em caso de falha.
*/
int inv_retro_subs (matrix_t *restrict l_matrix, matrix_t *restrict b_matrix, double *solution, int ls_b) {    
    int line, col, size;

    size = l_matrix->n;

    for (line = 0; line < size; line++) {
        solution[line] = b_matrix->coef[line][ls_b];
        for (col = line - 1; col >= 0; col--) 
            solution[line] -= l_matrix->coef[line][col] * solution[col];
        solution[line] /= l_matrix->coef[line][line];
        if (isnan (solution[line]) || isinf (solution[line]))
            return NAN_INF_ERROR;
    }
    return EXIT_SUCCESS;
}

/*!
    \brief Função que calcula o valor do sistema linear de uma matriz triangular superior
    \param u_matrix Ponteiro para a Matriz U
    \param inv_matrix Ponteiro para Matriz Inversa
    \param solution matrix com a solução da retrosubstituição anterior
    \param ls_b coluna da Matriz identidade (B do sistema linear)
    \return 0 em caso de sucesso. NAN_INF_ERROR em caso de falha.
*/
int retro_subs (matrix_t *restrict u_matrix, matrix_t *restrict inv_matrix, matrix_t *restrict solution, int ls_b) {    
    int line, col, size;

    size = u_matrix->n;

    for (line = size - 1; line >= 0; line--) {
        inv_matrix->coef[line][ls_b] = solution->coef[line][ls_b];
        for (col = line + 1; col < size; col++) 
            inv_matrix->coef[line][ls_b] -= u_matrix->coef[line][col] * inv_matrix->coef[col][ls_b];
        inv_matrix->coef[line][ls_b] /= u_matrix->coef[line][line];
        if (isnan (inv_matrix->coef[line][ls_b]) || isinf (inv_matrix->coef[line][ls_b]))
            return NAN_INF_ERROR;
    }
    return EXIT_SUCCESS;
}

/*!
    \brief Função que calcula a matrix de residuos do resultado
    \param residue_matrix Ponteiro para a Matriz de Residuos (Retorno)
    \param matrix Ponteiro para Matriz enviada pelo usuário
    \param inv_matrix Ponteiro para a Matriz inversa calculada
    \return 0 em caso de sucesso, ALLOC_ERROR em caso de falha de alocação
*/
void calc_residue (matrix_t *restrict residue_matrix, matrix_t *restrict matrix, matrix_t *restrict inv_matrix) {
    int line, col, pivot_col, size;
    double temp = 0;

    size = matrix->n;
    generate_identity_matrix (residue_matrix);

    for (pivot_col = 0; pivot_col < size ; pivot_col++) {
        for (line = 0; line < size; line++) {
            for (col = 0; col < size; col++) 
                temp += matrix->coef[line][col] * inv_matrix->coef[col][pivot_col];
            residue_matrix->coef[line][pivot_col] -= temp;
            temp = 0;
        }
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
    double norma = 0.0;

    size = residue_matrix->n;

    *residue_time = timestamp ();
    for (line = 0; line < size; line++)
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
            if (isnan (m) || isinf (m)) 
                return NAN_INF_ERROR;
            
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
    double *temp_sol;
    int count, sol_count; 
    int size;

    size = inv_matrix->n;

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
        if (inv_retro_subs (l_matrix, i_matrix, temp_sol, count) == NAN_INF_ERROR)
            return NAN_INF_ERROR;

        // adição da solução temporaria da execução do refinamento
        for (sol_count = 0; sol_count < size; sol_count++)
            solution->coef[sol_count][count] += temp_sol[sol_count];

        if (retro_subs (u_matrix, inv_matrix, solution, count) == NAN_INF_ERROR)
            return NAN_INF_ERROR;
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
    double *temp_sol;
    int size, count, ls_count, sol_count;

    *iter_time = 0;
    *residue_time = 0;

    size = matrix->n;

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
            if (inv_retro_subs (l_matrix, residue_matrix, temp_sol, ls_count) == NAN_INF_ERROR)
                return NAN_INF_ERROR;

            // adição da solução temporaria da execução do refinamento
            for (sol_count = 0; sol_count < size; sol_count++)
                solution->coef[sol_count][ls_count] += temp_sol[sol_count];

            if (retro_subs (u_matrix, inv_matrix, solution, ls_count) == NAN_INF_ERROR)
                return NAN_INF_ERROR;
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