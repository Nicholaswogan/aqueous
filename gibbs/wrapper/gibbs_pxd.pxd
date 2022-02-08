
from libc.stdint cimport int8_t, int16_t, int32_t, int64_t
from libcpp cimport bool, complex

cdef extern from "gibbs_wrapper.h":
  cdef void gibbs_gibbs_energy(int64_t *species_len, char *species, double *T, double *P, int64_t *err_len, void *err, double *G);
  cdef void gibbs_gibbs_energy_err(int64_t *err_len, void *err_cp, char *err)
  cdef void gibbs_load_spronsbl(int64_t *path_len, char *path);