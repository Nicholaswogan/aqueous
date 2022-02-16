
from numpy cimport ndarray, int8_t, int16_t, int32_t, int64_t
from libcpp cimport bool, complex
import numpy as np
cimport gibbs_pxd as pxd
import os

DEF STR_LEN = 20

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
    pxd.gibbs_err(&err_len, err, <char *> err_.data)
    raise Exception(err_.item().decode())
  return G
  
cdef class AqueousSolution:
  cdef void *_ptr
  
  def __cinit__(self, species_l = None):
    pxd.gibbs_alloc_aqueoussolution(&self._ptr)
    cdef int64_t species_dim = len(species_l)
    cdef ndarray species = np.array(species_l,dtype=np.dtype(('S', 20)),order='F')
    cdef int64_t err_len;
    cdef void *err
    cdef ndarray err_
    
    pxd.gibbs_aqueoussolution_init(&self._ptr, &species_dim, <char *> species.data, &err_len, &err)
    
    if err:
      err_ = np.zeros((),dtype=np.dtype(('S', err_len)))
      pxd.gibbs_err(&err_len, err, <char *> err_.data)
      raise Exception(err_.item().decode())
    
  def __dealloc__(self):
    pxd.gibbs_dealloc_aqueoussolution(&self._ptr)
    
  property xtol:
    def __get__(self):
      cdef double val;
      pxd.gibbs_aqueoussolution_xtol_get(&self._ptr, &val)
      return val
    def __set__(self, double val):
      pxd.gibbs_aqueoussolution_xtol_set(&self._ptr, &val)
      
  property conserv_tol:
    def __get__(self):
      cdef double val;
      pxd.gibbs_aqueoussolution_conserv_tol_get(&self._ptr, &val)
      return val
    def __set__(self, double val):
      pxd.gibbs_aqueoussolution_conserv_tol_set(&self._ptr, &val)
      
  property G_init:
    def __get__(self):
      cdef double val;
      pxd.gibbs_aqueoussolution_g_init_get(&self._ptr, &val)
      return val
      
  property G_opt:
    def __get__(self):
      cdef double val;
      pxd.gibbs_aqueoussolution_g_opt_get(&self._ptr, &val)
      return val
      
  property algorithm:
    def __get__(self):
      cdef char val[STR_LEN];
      pxd.gibbs_aqueoussolution_algorithm_get(&self._ptr, val)
      return val.decode()
    def __set__(self, str val):
      cdef int64_t val_len
      cdef char *val_c
      cdef bytes val_b = val.encode()
      val_c = val_b
      val_len = len(val)
      pxd.gibbs_aqueoussolution_algorithm_set(&self._ptr, &val_len, val_c)  
      
  def equilibrate(self, ndarray[double, ndim=1] m, double T, double P):
    cdef int64_t m_dim = m.shape[0]
    cdef int64_t err_len;
    cdef void *err
    cdef ndarray err_
    
    pxd.gibbs_aqueoussolution_equilibrate(&self._ptr, &m_dim, <double *> m.data, &T, &P, &err_len, &err);
    if err:
      err_ = np.zeros((),dtype=np.dtype(('S', err_len)))
      pxd.gibbs_err(&err_len, err, <char *> err_.data)
      raise Exception(err_.item().decode())
  



  
