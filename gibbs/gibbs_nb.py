
import numba as nb
import ctypes as ct
import _gibbs

addr = nb.extending.get_cython_function_address("_gibbs","gibbs_energy")
gibbs_energy_functype = ct.CFUNCTYPE(ct.c_double, ct.py_object, ct.c_double, ct.c_double)
gibbs_energy_for_numba = gibbs_energy_functype(addr)

@nb.njit
def gibbs_energy(species, T, P):
    return gibbs_energy_for_numba(species, T, P)


