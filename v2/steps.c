#include "steps.h"

pivot_steps_t *alloc_steps (int n) {
    pivot_steps_t *pivot_steps;
    double log_2;
    int internal_n;

    if ((pivot_steps = malloc (sizeof (pivot_steps_t)))) {
        internal_n = n;
        log_2 = log10 (n) / log10 (2);
        if (ceil (log_2) == floor (log_2)) internal_n++;

        pivot_steps->lines = calloc (internal_n, sizeof (int));
        pivot_steps->pivots = calloc (internal_n, sizeof (int));
        if (! pivot_steps->lines || !pivot_steps->pivots)
            return NULL;

        memset (pivot_steps->lines, 0, internal_n * sizeof (int));
        memset (pivot_steps->pivots, 0, internal_n * sizeof (int));
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
    int count;

    for (count = 0; count < matrix->n; count++) {
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