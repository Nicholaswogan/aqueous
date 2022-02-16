
import numba as nb
import ctypes as ct
import _aqueous

addr = nb.extending.get_cython_function_address("_aqueous","gibbs_energy")
gibbs_energy_functype = ct.CFUNCTYPE(ct.c_double, ct.py_object, ct.c_double, ct.c_double)

gibbs_energy = gibbs_energy_functype(addr)



