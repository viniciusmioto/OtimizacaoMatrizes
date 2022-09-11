#include "steps.h"

pivot_steps_t *alloc_steps (int n) {
    pivot_steps_t *pivot_steps;

    if ((pivot_steps = malloc (sizeof (pivot_steps_t)))) {
        pivot_steps->lines = calloc (n, sizeof (int));
        pivot_steps->pivots = calloc (n, sizeof (int));
        if (! pivot_steps->lines || !pivot_steps->pivots)
            return NULL;

        memset (pivot_steps->lines, 0, n * sizeof (int));
        memset (pivot_steps->pivots, 0, n * sizeof (int));
        pivot_steps->qtd = 0;
        return pivot_steps;
    } else return NULL;
}

void add_pivot_step (pivot_steps_t *steps, int line, int pivot) {
    steps->lines[steps->qtd] = line;
    steps->pivots[steps->qtd] = pivot;
    steps->qtd++;
}

void apply_pivot_steps (matrix_t *matrix, pivot_steps_t *steps) {
    int count;
    for (count = 0; count < steps->qtd; count++) 
        swap_line (matrix, steps->lines[count], steps->pivots[count]);
}

void swap_line (matrix_t *matrix, int line, int pivot) {
    double m_aux;
    int count, size;
    int unroll_limit;
    size = matrix->n;
    unroll_limit = size - (size % UNROLL_SIZE);

    for (count = 0; count < unroll_limit; count += UNROLL_SIZE) {
        m_aux = matrix->coef[line][count];
        matrix->coef[line][count] = matrix->coef[pivot][count];
        matrix->coef[pivot][count] = m_aux;

        m_aux = matrix->coef[line][count + 1];
        matrix->coef[line][count + 1] = matrix->coef[pivot][count + 1];
        matrix->coef[pivot][count + 1] = m_aux;

        m_aux = matrix->coef[line][count + 2];
        matrix->coef[line][count + 2] = matrix->coef[pivot][count + 2];
        matrix->coef[pivot][count + 2] = m_aux;

        m_aux = matrix->coef[line][count + 3];
        matrix->coef[line][count + 3] = matrix->coef[pivot][count + 3];
        matrix->coef[pivot][count + 3] = m_aux;
    }

    for (count = unroll_limit; count < size; count++) {
        m_aux = matrix->coef[line][count];
        matrix->coef[line][count] = matrix->coef[pivot][count];
        matrix->coef[pivot][count] = m_aux;
    }
}

void free_steps (pivot_steps_t *steps) {
    free (steps->lines);
    free (steps->pivots);
    free (steps);
}