#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import triqstools.tools as tools
import triqstools.pytools as pytools
from triqstools.pytools import print_mpi
from pytriqs.gf import Gf, BlockGf, iOmega_n, Idx
from pytriqs.gf import MeshBrillouinZone, MeshImFreq, MeshImTime, MeshProduct, MeshCyclicLattice, MeshPoint
from pytriqs.lattice import BravaisLattice, BrillouinZone
from pytriqs.archive import HDFArchive, hdf_archive_schemes
from pytriqs.utility import mpi
from pytriqs.lattice.tight_binding import TightBinding, dos
from triqs_ctint import SolverCore
import numpy as np
from math import pi
from math import sqrt
from scipy.optimize import curve_fit, approx_fprime
from scipy import polyfit, polyval, polyder


# In[ ]:


beta = 2.
mu   = 0.


# In[ ]:


n_iw = 4
n_iW = 4
n_k  = 4
n_q  = 4


# In[ ]:


# define the Bravais lattice: a square lattice in 2d
BL = BravaisLattice(units = [(1, 0), (0, 1)])

# define the tight-binding model, i.e., the hopping parameters
t = -1.0               # nearest neighbor hopping
tp = -0.00*t           # next-nearest neighbor hopping

# hopping[displacement on the lattice] = [[t11, t12, t13....], [t21, t22, t23....], ..., [...., tnn]] 
# where n=Number_Orbitals
hop= {  (1, 0)  :  [[ t]],       
        (-1, 0) :  [[ t]],     
        (0, 1)  :  [[ t]], 
        (0, -1) :  [[ t]], 
        (1, 1)  :  [[ tp]], 
        (-1, -1):  [[ tp]], 
        (1, -1) :  [[ tp]], 
        (-1, 1) :  [[ tp]]}

TB = TightBinding(BL, hop)


# In[ ]:


iw_mesh = MeshImFreq(beta, 'Fermion', n_iw)
iW_mesh = MeshImFreq(beta, 'Boson', n_iW)
iw_iW_mesh = MeshProduct(iw_mesh, iW_mesh)
k_mesh = MeshBrillouinZone(BrillouinZone(BL), n_k)
q_mesh = MeshBrillouinZone(BrillouinZone(BL), n_q)
k_iw_mesh = MeshProduct(k_mesh, iw_mesh)
q_iW_mesh = MeshProduct(q_mesh, iW_mesh)


# In[ ]:


G0 = pytools.make_G0(Gf(mesh=k_iw_mesh, target_shape=[]), TB, mu)


# In[ ]:


print_mpi("calculating total self-energy")
Sigma = tools.self_energy_weak_coupling(G0, q_mesh, iw_mesh, 2.0)


# In[ ]:


print_mpi("calculating bosonic fluctuation diagnostics for momentum")
Sigma_q = tools.self_energy_weak_coupling_bosonic_momentum(G0, q_mesh, iw_mesh, 2.0)


# In[ ]:


print_mpi("calculating bosonic fluctuation diagnostics for frequency")
Sigma_iW = tools.self_energy_weak_coupling_bosonic_frequency(G0, q_mesh, iW_mesh, 2.0)


# In[ ]:


print_mpi("calculating fermionic fluctuation diagnostics for momentum")
Sigma_kp = tools.self_energy_weak_coupling_fermionic_momentum(G0, q_mesh, iw_mesh, 2.0)


# In[ ]:


mpi.barrier()
print_mpi("saving results")
if mpi.is_master_node():
    with HDFArchive('weak_coupling_fluct_diag.h5', 'a') as A:
        A['beta']     = beta
        A['n_k']      = n_k
        A['k_mesh']   = k_mesh
        A['n_q']      = n_q
        A['q_mesh']   = q_mesh
        A['Sigma']    = Sigma
        A['Sigma_q']  = Sigma_q
        A['Sigma_kp'] = Sigma_kp
        A['Sigma_iW'] = Sigma_iW

