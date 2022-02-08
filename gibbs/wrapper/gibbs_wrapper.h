
#include <stdint.h>
#include <complex.h>
#include <stdbool.h>

void gibbs_gibbs_energy(int64_t *species_len, char *species, double *T, double *P, int64_t *err_len, void *err, double *G);
void gibbs_gibbs_energy_err(int64_t *err_len, void *err_cp, char *err);
void gibbs_load_spronsbl(int64_t *path_len, char *path);