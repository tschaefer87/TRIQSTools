from __future__ import print_function
import pytriqs.utility.mpi
import pytriqs.lattice.tight_binding
import pytriqs.lattice.lattice_tools
from pytriqs.gf import Gf, BlockGf, Block2Gf, Idx, MeshProduct, map_block, MeshCyclicLattice
import numpy
import numpy.linalg
import itertools
import scipy.optimize
import scipy.special
import triqstools.tools
import sys

def print_mpi(*args, **kwargs):
    """ Printing function, which prints a string on the MPI master node only.
    
    """
    if pytriqs.utility.mpi.is_master_node():
        print(*args, file=sys.stdout, **kwargs)
        sys.stdout.flush()
        
def print_err_mpi(*args, **kwargs):
    """ Printing function, which prints a string on the MPI master node only.
    
    """
    if pytriqs.utility.mpi.is_master_node():
        print(*args, file=sys.stderr, **kwargs)
        sys.stderr.flush()

def max_norm(G):
    """ Calculates the maximum (infinity) norm of G
    
    Parameters
    ----------
    G : Gf, BlockGf
        Object, for which the max_norm shall be calculated
    
    Returns
    -------
    float
        infinity norm of the data array of G
    
    Raises
    ------
    NotImplementedError
        If G is neither Gf nor BlockGf
    """
    if isinstance(G, Gf):
        return numpy.amax(numpy.abs(G.data))
    elif isinstance(G, BlockGf):
        return max(max_norm(G_bl) for bl, G_bl in G)
    else:
        raise NotImplementedError, " max_norm not implemented for this type "

def find_chem_pot(filling, G, mu, accuracy=1E-12, max_it=10000):
    """ Determines the chemical potential, which leads to a certain filling calculated by the density of G, by means of a bisectioning algorithm.
    
    Parameters
    ----------
    filling : float
        The filling to be achieved by shifting the chemical potential. Half-filling is given by 0.5
    G : Gf
        The Green function for which mu has to be adopted
    accuracy : float
        Threshold how accurate the density should be determined
    max_it : int
        Maximum number of iterations
        
    Returns
    -------
    float
        chemical potential, which gives a density of filling with the Green function G
    
    
    """
    return mu + scipy.optimize.bisect(lambda mu_diff: (trilex.tools.calc_filling(G,mu_diff) - filling).real,-10.0,10.0,xtol=accuracy,rtol=accuracy,maxiter=max_it)

# FIXME: implementation of this routine should be simpler --> tight-binding functionality of TRIQS
# FIXME: TB should also know its dim
# FIXME: eps(TB,k)
def eps(BL, TB, k):
    """ Calculates the dispersion at a specific reciprocal k-point.

    Parameters
    ----------
    BL : BravaisLattice
        Bravais lattice for which the dispersion shall be determined
    TB : TightBinding
        Tight binding model for which the dispersion shall be determined
    k : MeshPoint
        Brillouin zone mesh point

    Returns
    -------
    float
        Dispersion energy at k

    """
    return pytriqs.lattice.tight_binding.energies_on_bz_path(TB, \
            [k[i]*1./(2.*numpy.pi) for i in range(0,BL.dim)], \
            [k[i]*1./(2.*numpy.pi) for i in range(0,BL.dim)], \
            1)[0,0]

def make_G0(G, TB, mu=0.):
    """ Calculates and returns a (non-interacting) lattice Green function.
    
    Parameters
    ----------
    G : Gf, BlockGf
        The structure of this Gf is used for the construction of G0
    TB : TightBinding
        Used for determining the dispersion via the function eps
    mu : float
        chemical potential
    
    Returns
    -------
    Gf
        (non-interacting) lattice Green function
    """
    G0 = G.copy()
    G0.zero()

    k_mesh, iw_mesh = G0.mesh[0], G0.mesh[1]

    iw_vec = numpy.array([iw.value for iw in iw_mesh])
    k_vec = numpy.array([k.value for k in k_mesh])
    e_k_vec = pytriqs.lattice.hopping_stack(TB, k_vec.T / 2. / numpy.pi).transpose(2, 0, 1)[::,0,0]

    G0.data[:] = 1.0 / (iw_vec[None,::] - e_k_vec[::,None] + mu)
    return G0

def make_vq(V_Hubbard, q_mesh, BL, cutoff=0.05):
    """ Calculates the Ewald summation (PhD Thesis Thomas Ayral, Eq. (6.5)) for treating non-local Coulomb interactions
    
    Parameters
    ----------
    V_Hubbard : 
        prefactor of the non-local part of the Coulomb interaction V(q)

    q_mesh :
        Brioullin zone mesh on which to evaluate V(q)

    BL :
        Bravais lattice accompanying q_mesh

    cutoff :
        The cutoff of the Ewald summation
        
    Returns
    -------
    v
        Gf v(q_mesh), i.e. V(q)
    
    """
    v = Gf(mesh=q_mesh, target_shape=[])
    n_q = numpy.sqrt(len(q_mesh))
    eta = 0.1 * numpy.sqrt(n_q) + 0.9 * n_q
    
    periodization_matrix = numpy.matrix([[2*n_q+1, 0, 0], [0, 2*n_q+1, 0], [0, 0, 1]], numpy.int32)
    r_mesh = MeshCyclicLattice(BL, periodization_matrix)
    
    v.zero()
    for q in q_mesh:
        norm_q = numpy.linalg.norm(q.value)
        if norm_q <= cutoff:
            norm_q = cutoff
        for r in r_mesh:
            r_value = r.value - ([n_q, n_q, 0] if len(r.value) == 3 else [n_q, n_q] if len(r_value) == 2 else [n_q])
            norm_r = numpy.linalg.norm(r_value)
            if norm_r > 1E-3:
                v[q] += V_Hubbard * (scipy.special.erfc(norm_r / eta) * numpy.exp(1j * numpy.dot(q.value, r_value)) / norm_r).real
        v[q] += V_Hubbard * (2. * numpy.pi * scipy.special.erfc(0.5 * norm_q * eta) / norm_q - 2. / (eta * numpy.sqrt(numpy.pi))).real
        
    return v


def get_charge_channel(chi):
    """ Calculates the linear combination suitable for going from spin-indices to the spin-diagonalized (physical) charge channel. This routine assumes SU(2) symmetry.
    
    Parameters
    ----------
    chi : Block2Gf
        Green function in spin basis
        
    Returns
    -------
    Gf
        Green function in the charge channel
    
    """
    return (chi['up','up'] + chi['dn','dn'] \
          + chi['up','dn'] + chi['dn','up'])

def get_spin_channel(chi):
    """ Calculates the linear combination suitable for going from spin-indices to the spin-diagonalized (physical) spin channel. This routine assumes SU(2) symmetry.
    
    Parameters
    ----------
    chi : Block2Gf
        Green function in spin basis
        
    Returns
    -------
    Gf
        Green function in the spin channel
    """
    return (chi['up','up'] + chi['dn','dn'] \
          - chi['up','dn'] - chi['dn','up'])
 
def make_G_ij_from_G_R(G_R):
    """ Computes the Green function in an absolute site representation from the one of a relative site representation.
    
        In the absolute site representation, the sites are given by the target space structure, in the relative site representation by an additional cyclic mesh.
        
    Parameters
    ----------
    G_R : Gf, BlockGf, Block2Gf
        Green function in relative site representation
        
    Returns
    -------
    Gf, BlockGf, Block2Gf
        Green function in absolute site representation
        
    Raises
    ------
    Exception
        If the type of G_R is not Gf, BlockGf or Block2Gf
    """
    if isinstance(G_R,Gf): 
        assert G_R.target_shape==()
        
        R_points = list(G_R.mesh[0])
        n_c = len(R_points)

        if G_R.mesh.rank == 2:
            new_mesh = G_R.mesh[1]
        else:
            new_mesh = MeshProduct(*G_R.mesh[1:])

        G_ij= Gf(name=G_R.name.replace("_R","_ij"), mesh=new_mesh, target_shape=[n_c,n_c], is_real=(G_R.data.dtype==numpy.dtype('float64')))
        G_ij.zero()
        for i in range(0,n_c):
            for j in range(0,n_c):
                # FIXME: can we get rid of cast (long -> int), R = R_points[i] - R_points[j]
                R = [ int(R_points[i][k]-R_points[j][k]) for k in range(3) ]
                G_ij[i,j] = G_R[Idx(*R),:]
        return G_ij
                    
    elif isinstance(G_R, (BlockGf, Block2Gf)):
        return map_block(make_G_ij_from_G_R, G_R)
    
    else:
        raise Exception("type of GF unknown")
        

def make_G_R_from_G_ij(G_ij, R_mesh):
    """ Computes the Green function in a relative site representation from the one of an absolute site representation.
    
        In the absolute site representation, the sites are given by the target space structure, in the relative site representation by an additional cyclic mesh.
        
    Parameters
    ----------
    G_ij : Gf, BlockGf, Block2Gf
        Green function in absolute site representation
    R_mesh : MeshCyclicLattice
        Mesh for the construction of the relative site Green function
        
    Returns
    -------
    Gf, BlockGf, Block2Gf
        Green function in relative site representation
        
    Raises
    ------
    Exception
        If the type of G_R is not Gf, BlockGf or Block2Gf
    """
    R_points = list(R_mesh)
    n_c = len(R_points)
    if isinstance(G_ij, Gf):
        assert len(R_mesh) == G_ij.target_shape[0]
        if not isinstance(G_ij.mesh, MeshProduct):
            G_R = Gf(name=G_ij.name.replace("_ij","_R"), mesh=MeshProduct(R_mesh, G_ij.mesh), target_shape=[])
            G_R.zero()
            for i in range(0,n_c):
                for j in range(0,n_c):
                    # FIXME: can we get rid of cast (long -> int), R = R_points[i] - R_points[j]
                    R = [ int(R_points[i][k]-R_points[j][k]) for k in range(3) ]
                    if G_ij.target_rank == 2:
                        G_R[Idx(*R),:] = G_R[Idx(*R),:] + G_ij[i,j]
                    elif G_ij.target_rank == 4:
                        G_R[Idx(*R),:] = G_R[Idx(*R),:] + G_ij[i,i,j,j]
                    else:
                        raise Exception("type of GF not implemented")
            return G_R / n_c
        elif G_ij.mesh.rank == 2:
            assert G_ij.target_rank == 4
            G_R = Gf(name=G_ij.name.replace("_ij","_R"), mesh=MeshProduct(R_mesh, R_mesh, G_ij.mesh[0], G_ij.mesh[1]), target_shape=[])
            G_R.zero()
            for i in range(0,n_c):
                for j in range(0,n_c):
                    # FIXME: can we get rid of cast (long -> int), R = R_points[i] - R_points[j]
                    R1 = [ int(R_points[i][l]-R_points[j][l]) for l in range(3) ]
                    for k in range(0,n_c):
                        # FIXME: can we get rid of cast (long -> int), R = R_points[i] - R_points[j]
                        R2 = [ int(R_points[k][l]-R_points[j][l]) for l in range(3) ]
                        G_R[Idx(*R1),Idx(*R2),:,:] = G_R[Idx(*R1),Idx(*R2),:,:] + G_ij[i,j,k,k]                
            return G_R / n_c
        
        else:
            raise Exception("mesh not implemented")
            
        
    elif isinstance(G_ij, (BlockGf, Block2Gf)):
        return map_block(lambda G: make_G_R_from_G_ij(G, R_mesh), G_ij)

    elif isinstance(G_ij, numpy.ndarray):
        assert len(R_mesh) == G_ij.shape[0]
        G_R = Gf(name='density_R', mesh=R_mesh, target_shape=[])
        for i in range(0,n_c):
            for j in range(0,n_c):
            # FIXME: can we get rid of cast (long -> int), R = R_points[i] - R_points[j]
                R = [ int(R_points[i][k]-R_points[j][k]) for k in range(3) ]
                G_R[Idx(*R)] = G_R[Idx(*R)] + G_ij[i,j]
        return G_R / n_c
    
         
    else:
        raise Exception("type of GF unknown: ", type(G_ij))
