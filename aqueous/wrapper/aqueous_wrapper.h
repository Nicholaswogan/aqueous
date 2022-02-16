
#include <stdint.h>
#include <complex.h>
#include <stdbool.h>

void aqueous_gibbs_energy(int64_t *species_len, char *species, double *T, double *P, int64_t *err_len, void *err, double *G);
void aqueous_err(int64_t *err_len, void *err_cp, char *err);
void aqueous_load_spronsbl(int64_t *path_len, char *path);

void aqueous_alloc_aqueoussolution(void *ptr);
void aqueous_dealloc_aqueoussolution(void *ptr);
void aqueous_aqueoussolution_init(void *ptr, int64_t *species_dim, char* species, int64_t *err_len, void *err);
void aqueous_aqueoussolution_equilibrate(void *ptr, int64_t *m_dim, double *m, double *T, double *P, int64_t *err_len, void *err);
void aqueous_aqueoussolution_xtol_get(void *ptr, double *val);
void aqueous_aqueoussolution_xtol_set(void *ptr, double *val);
void aqueous_aqueoussolution_conserv_tol_get(void *ptr, double *val);
void aqueous_aqueoussolution_conserv_tol_set(void *ptr, double *val);
void aqueous_aqueoussolution_g_init_get(void *ptr, double *val);
void aqueous_aqueoussolution_g_opt_get(void *ptr, double *val);
void aqueous_aqueoussolution_algorithm_get(void *ptr, char *val);
void aqueous_aqueoussolution_algorithm_set(void *ptr, int64_t *val_len, char *val);