
from numpy cimport ndarray, int8_t, int16_t, int32_t, int64_t
from libcpp cimport bool, complex
import numpy as np
cimport gibbs_pxd as pxd

import os


cdef load_spronsbl():
  cdef int64_t path_len
  cdef char *path_c
  cdef bytes path = (os.path.dirname(os.path.realpath(__file__))+'/data/').encode()
  path_c = path
  path_len = len(path)
  pxd.gibbs_load_spronsbl(&path_len, path_c)
  
load_spronsbl()

cpdef api double gibbs_energy(str species, double T, double P) except? -1.0:
  cdef int64_t species_len
  cdef char *species_c
  cdef int64_t err_len
  cdef void *err
  cdef double G
  cdef ndarray err_
  
  cdef bytes species_b = species.encode()
  species_c = species_b
  species_len = len(species)
  pxd.gibbs_gibbs_energy(&species_len, species_c, &T, &P, &err_len, &err, &G)
  if err:
    err_ = np.zeros((),dtype=np.dtype(('S', err_len)))
    pxd.gibbs_gibbs_energy_err(&err_len, err, <char *> err_.data)
    raise Exception(err_.item().decode())
  return G
  



  
