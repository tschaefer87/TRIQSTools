# Generated automatically using the command :
# c++2py ../../c++/triqstools/tools.hpp -p --members_read_only -N triqstools -a triqstools -m tools -o tools -C pytriqs --cxxflags="-std=c++17 "
from cpp2py.wrap_generator import *

# The module
module = module_(full_name = "tools", doc = r"", app_name = "triqstools")

# Imports
module.add_imports(*['pytriqs.gf', 'pytriqs.lattice'])

# Add here all includes
module.add_include("triqstools/tools.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <cpp2py/converters/complex.hpp>
#include <triqs/cpp2py_converters/arrays.hpp>
#include <triqs/cpp2py_converters/gf.hpp>

using namespace triqstools;
""")


module.add_function ("triqstools::g_r_t triqstools::make_gf_from_fourier_k (triqstools::g_k_cvt G)", doc = r"""""")

module.add_function ("triqstools::g_k_t triqstools::make_gf_from_fourier_r (triqstools::g_r_cvt G)", doc = r"""""")

module.add_function ("triqstools::g_r_iw_t triqstools::make_gf_from_fourier_k (triqstools::g_k_iw_cvt G)", doc = r"""""")

module.add_function ("triqstools::g_k_iw_t triqstools::make_gf_from_fourier_r (triqstools::g_r_iw_cvt G)", doc = r"""""")

module.add_function ("triqstools::chi2_r_iW_t triqstools::make_gf_from_fourier_k (triqstools::chi2_q_iW_cvt W)", doc = r"""""")

module.add_function ("triqstools::chi2_q_iW_t triqstools::make_gf_from_fourier_r (triqstools::chi2_r_iW_cvt W)", doc = r"""""")

module.add_function ("triqstools::g_R_R_iw_iW_t triqstools::make_gf_from_fourier_k (triqstools::g_K_Q_iw_iW_cvt G)", doc = r"""""")

module.add_function ("triqstools::g_K_Q_iw_iW_t triqstools::make_gf_from_fourier_r (triqstools::g_R_R_iw_iW_cvt G)", doc = r"""""")

module.add_function ("triqstools::chi3_R_iw_t triqstools::make_gf_from_fourier_k (triqstools::chi3_Q_iw_cvt G)", doc = r"""""")

module.add_function ("triqstools::chi3_Q_iw_t triqstools::make_gf_from_fourier_r (triqstools::chi3_R_iw_cvt G)", doc = r"""""")

module.add_function ("triqstools::g_k_iw_t triqstools::dyson_k_iw (triqstools::g_k_iw_cvt G0, triqstools::g_k_iw_cvt Sigma, double mu = 0.)", doc = r"""Dyson equation for a one-particle lattice Green function

 :math:`G({\mathbf k},i\omega)=\frac{1}{i\omega+\mu-\epsilon_{\mathbf k}-\Sigma({\mathbf k},i\omega)}`

Parameters
----------
G0
     non-interacting Green function :math:`G_{0}(\mathbf{k},i\omega)`

Sigma
     self-energy :math:`\Sigma(\mathbf{k},i\omega)`

mu
     additional constant real shift (e.g., for chemical potential)

Returns
-------
out
     interacting one-particle fermionic lattice Green function :math:`G(\mathbf{k},i\omega)`""")

module.add_function ("triqstools::g_k_iw_t triqstools::dyson_k_iw (triqstools::g_k_iw_cvt G0, triqstools::g_iw_cvt Sigma, double mu = 0.)", doc = r"""Dyson equation for a one-particle lattice Green function with local self-energy

 :math:`G({\mathbf k},i\omega)=\frac{1}{i\omega+\mu-\epsilon_{\mathbf k}-\Sigma(i\omega)}`

Parameters
----------
G0
     non-interacting Green function :math:`G_{0}(\mathbf{k},i\omega)`

Sigma
     local self-energy :math:`\Sigma(i\omega)`

mu
     additional constant real shift (e.g., for chemical potential)

Returns
-------
out
     interacting one-particle lattice Green function :math:`G(\mathbf{k},i\omega)`""")

module.add_function ("triqstools::chi2_q_iW_t triqstools::dyson_k_iw (triqstools::g_k_iw_cvt G0, triqstools::chi2_iW_cvt Sigma, double mu = 0.)", doc = r"""""")

module.add_function ("triqstools::chi2_q_iW_t triqstools::dyson_q_iW (array<double,1> U, triqstools::chi2_q_iW_cvt P)", doc = r"""Dyson equation for a one-particle lattice block Green function

 :math:`W^{\eta}({\mathbf q},i\Omega)=\frac{U^{\eta}}{1-U^{\eta}P^{\eta}({\mathbf q},i\Omega)}`

Parameters
----------
U
     bare propagator :math:`U^{\eta}`

P
     polarization :math:`P^{\eta}(\mathbf{q},i\Omega)`

Returns
-------
out
     interacting one-particle lattice block Green function :math:`W^{\eta}(i\Omega)`""")

module.add_function ("triqstools::chi2_q_iW_t triqstools::dyson_q_iW (triqstools::chi2_q_cvt U, triqstools::chi2_q_iW_cvt P)", doc = r"""""")

module.add_function ("triqstools::chi2_iW_t triqstools::dyson_iW (array<double,1> U, triqstools::chi2_iW_cvt P)", doc = r"""local Dyson equation for a one-particle lattice block Green function

 :math:`W^{\eta}(i\Omega)=\frac{U^{\eta}}{1-U^{\eta}P^{\eta}(i\Omega)}`

Parameters
----------
U
     bare propagator :math:`U^{\eta}`

P
     local polarization :math:`P^{\eta}(i\Omega)`

Returns
-------
out
     interacting one-particle lattice block Green function :math:`W^{\eta}(\mathbf{q},i\Omega)`""")

module.add_function ("triqstools::chi2_q_iW_t triqstools::dyson_iW (triqstools::chi2_q_cvt U, triqstools::chi2_iW_cvt P)", doc = r"""""")

module.add_function ("triqstools::g_iw_t triqstools::inverse_dyson (triqstools::g_iw_cvt G, triqstools::g_iw_cvt Sigma)", doc = r"""inverse local Dyson equation

 :math:`G_{0}(i\omega) = \left(G^{-1}(i\omega) + \Sigma(i\omega)\right)^{-1}`

Parameters
----------
G
     interacting one-particle Green function :math:`G(i\omega)`

Sigma
     local self-energy :math:`\Sigma(i\omega)`

Returns
-------
out
     non-interacting Green function :math:`G_{0}(i\omega)`""")

module.add_function ("triqstools::g_k_iw_t triqstools::inverse_dyson (triqstools::g_k_iw_cvt G, triqstools::g_k_iw_cvt Sigma)", doc = r"""inverse lattice Dyson equation

 :math:`G_{0}(\mathbf{k},i\omega) = \left(G^{-1}(\mathbf{k},i\omega) + \Sigma(\mathbf{k},i\omega)\right)^{-1}`

Parameters
----------
G
     interacting one-particle Green function :math:`G(\mathbf{k},i\omega)`

Sigma
     self-energy :math:`\Sigma(\mathbf{k},i\omega)`

Returns
-------
out
     non-interacting Green function :math:`G_{0}(\mathbf{k},i\omega)`""")

module.add_function ("triqstools::chi2_iW_t triqstools::inverse_dyson (triqstools::chi2_iW_cvt W, triqstools::chi2_iW_cvt P)", doc = r"""inverse local Dyson equation

 :math:`W_{0}(i\omega) = \left(W^{-1}(i\Omega) + P(i\Omega)\right)^{-1}`

Parameters
----------
G
     interacting one-particle block Green function :math:`W^{\eta}(i\Omega)`

Sigma
     local self-energy :math:`P^{\eta}(i\Omega)`

Returns
-------
out
     non-interacting Green function :math:`W_{0}^{\eta}(i\Omega)`""")

module.add_function ("triqstools::chi2_q_iW_t triqstools::inverse_dyson (triqstools::chi2_q_iW_cvt W, triqstools::chi2_q_iW_cvt P)", doc = r"""inverse Dyson equation

 :math:`W_{0}({\mathbf q},i\omega) = \left(W^{-1}({\mathbf q},i\Omega) + P({\mathbf q},i\Omega)\right)^{-1}`

Parameters
----------
G
     interacting one-particle block Green function :math:`W^{\eta}({\mathbf q},i\Omega)`

Sigma
     local self-energy :math:`P^{\eta}({\mathbf q},i\Omega)`

Returns
-------
out
     non-interacting Green function :math:`W_{0}^{\eta}({\mathbf q},i\Omega)`""")

module.add_function ("triqstools::g_iw_t triqstools::make_local_gf (triqstools::g_k_iw_cvt G)", doc = r"""localizes the lattice Green function by a normalized sum over the k-mesh

Parameters
----------
G
     lattice Green function :math:`G(\mathbf{k},i\omega)`

Returns
-------
out
     local Green function :math:`G(i\omega)`""")

module.add_function ("std::complex<double> triqstools::make_local_gf (triqstools::g_iw_cvt G)", doc = r"""""")

module.add_function ("triqstools::g_w_t triqstools::make_local_gf (triqstools::g_k_w_cvt G)", doc = r"""""")

module.add_function ("triqstools::g_K_iw_t triqstools::sum_over_k_mesh (triqstools::g_K_k_iw_cvt G)", doc = r"""(normalized) sum over the second of two Brillouin zone meshes of a theta-weighted Green function

Parameters
----------
G
     theta-weighted Green function :math:`G(\mathbf{K},\mathbf{k},i\omega)`

Returns
-------
out
     cluster Green function :math:`G(\mathbf{K},i\omega)`""")

module.add_function ("std::complex<double> triqstools::sum_k (triqstools::g_k_cvt G)", doc = r"""""")

module.add_function ("triqstools::chi2_iW_t triqstools::make_local_gf (triqstools::chi2_q_iW_cvt W)", doc = r"""localizes the lattice block Green function by a normalized sum over the q-mesh

Parameters
----------
G
     lattice block Green function :math:`W^{\eta}(\mathbf{q},i\Omega)`

Returns
-------
out
     local block Green function :math:`W^{\eta}(i\Omega)`""")

module.add_function ("triqstools::chi2_Q_iW_t triqstools::sum_over_k_mesh (triqstools::chi2_Q_q_iW_cvt W)", doc = r"""(normalized) sum over the second of two Brillouin zone meshes of a theta-weighted block Green function

Parameters
----------
W
     theta-weighted block Green function :math:`W^{\eta}(\mathbf{Q},\mathbf{q},i\Omega)`

Returns
-------
out
     cluster block Green function :math:`W^{\eta}(\mathbf{Q},i\Omega)`""")

module.add_function ("triqstools::chi3_iw_t triqstools::make_Lambda_reg (triqstools::chi3_iw_cvt chi3_tilde, triqstools::g_iw_cvt G_imp, triqstools::chi2_iW_cvt U_Weiss, triqstools::chi2_iW_cvt chi2, triqstools::chi2_iW_cvt l, double n)", doc = r"""extracts the regular part of the electron-boson coupling vertex :math:`\Lambda_{reg}^{\eta}(i\omega,i\Omega)`
 from the equal-time correlation function :math:`\tilde{\chi}^{\eta}(i\omega,i\Omega)`
 T. Ayral and O. Parcollet,  Phys. Rev. B 93, 235124 (2016), Eq. (76) combined with Eq. (66)

Parameters
----------
chi3_tilde
     the equal-time correlation function :math:`\tilde{\chi}^{\eta}(i\omega,i\Omega)`

G_imp
     the local Green function (e.g. impurity model Green function) :math:`G(i\omega)`

U_Weiss
     the bosonic Weiss field :math:`U^{\eta}(i\Omega)`

chi2
     the local connected density-density correlation function :math:`\chi^{\eta}(i\Omega)`

l
     part of the vertex, so that :math:`\Lambda^{\eta}(i\omega,i\Omega)=\Lambda^{\eta}_{reg}(i\omega,i\Omega)+l^{\eta}(i\Omega)`

Returns
-------
out
     regular part of the electron-boson coupling vertex :math:`Lambda^{\eta}_{reg}(i\omega,i\Omega)`""")

module.add_function ("triqstools::chi3_Q_iw_t triqstools::make_Lambda_reg (triqstools::chi3_Q_iw_cvt chi3_tilde, triqstools::g_k_iw_cvt G_imp, triqstools::chi2_q_iW_cvt U_Weiss, triqstools::chi2_q_iW_cvt chi2, triqstools::chi2_q_iW_cvt l, triqstools::g_k_cvt n_K)", doc = r"""extracts the regular part of the electron-boson coupling vertex :math:`\Lambda_{reg}^{\eta}(\mathbf{K},\mathbf{Q},i\omega,i\Omega)`
 from the equal-time correlation function :math:`\tilde{\chi}^{\eta}(\mathbf{K},\mathbf{Q},i\omega,i\Omega)`
 T. Ayral and O. Parcollet,  Phys. Rev. Lett. 119, 166401 (2017)

Parameters
----------
chi3_tilde
     the equal-time correlation function :math:`\tilde{\chi}^{\eta}(\mathbf{K},\mathbf{Q},i\omega,i\Omega)`

G_imp
     the cluster impurity Green function (e.g. impurity model Green function) :math:`G(\mathbf{K},i\omega)`

U_Weiss
     the bosonic Weiss field :math:`U^{\eta}(\mathbf{Q},i\Omega)`

chi2
     the cluster connected density-density correlation function :math:`\chi^{\eta}(\mathbf{Q},i\Omega)`

l
     part of the vertex, so that :math:`\Lambda^{\eta}(\mathbf{K},\mathbf{Q},i\omega,i\Omega)=\Lambda^{\eta}_{reg}(\mathbf{K},\mathbf{Q},i\omega,i\Omega)+l^{\eta}(\mathbf{Q},i\Omega)`

Returns
-------
out
     regular part of the electron-boson coupling vertex :math:`Lambda^{\eta}_{reg}(\mathbf{K},\mathbf{Q},i\omega,i\Omega)`""")

module.add_function ("triqstools::chi2_q_iW_t triqstools::make_chi_lambda (triqstools::chi2_q_iW_cvt chi, array<double,1> lambda)", doc = r"""""")

module.add_function ("triqstools::chi2_q_iW_t triqstools::make_chi_from_W (array<double,1> U, triqstools::chi2_q_iW_cvt W)", doc = r"""""")

module.add_function ("triqstools::chi2_q_iW_t triqstools::make_chi_from_W (triqstools::chi2_q_cvt U, triqstools::chi2_q_iW_cvt W)", doc = r"""""")

module.add_function ("triqstools::chi2_q_iW_t triqstools::make_P_from_W (array<double,1> U, triqstools::chi2_q_iW_cvt W)", doc = r"""""")

module.add_function ("triqstools::chi2_q_iW_t triqstools::make_W_from_chi (array<double,1> U, triqstools::chi2_q_iW_cvt chi)", doc = r"""""")

module.add_function ("std::complex<double> triqstools::density (triqstools::g_k_iw_cvt G)", doc = r"""calculates the density for a given lattice Green function

Parameters
----------
G
     Green function :math:`G(\mathbf{k},i\omega)` for which the density is determined. The necessary tail will automatically be determined.

Returns
-------
out
     density, half-filling corresponds to 0.5""")

module.add_function ("std::complex<double> triqstools::density (triqstools::g_iw_cvt G)", doc = r"""""")

module.add_function ("std::complex<double> triqstools::sum_k_iw (triqstools::g_k_iw_cvt G)", doc = r"""(normalized) sum over both, the Brillouin zone and Matsubara mesh of a given lattice Green function.
 Be careful: no tail is used.

Parameters
----------
G
     Green function :math:`G(\mathbf{k},i\omega)`.

Returns
-------
out
     (normalized) sum over both meshes without tail""")

module.add_function ("std::complex<double> triqstools::calc_filling (triqstools::g_k_iw_cvt G, double mu_diff)", doc = r"""calculates the density for a given lattice Green function including a shift in the chemical potential

Parameters
----------
G
     Green function :math:`G(\mathbf{k},i\omega)` for which the density is determined. The necessary tail will automatically be determined.

mu_diff
     additional shift in the chemical potential

Returns
-------
out
     density, half-filling corresponds to 0.5""")

module.add_function ("triqstools::g_k_iw_t triqstools::copy_gf_to_smaller_intervall (triqstools::g_k_iw_cvt G_big, triqstools::g_k_iw_cvt G_small)", doc = r"""copies Green function into a Gf of smaller mesh size

Parameters
----------
G_big
     Green function :math:`G(\mathbf{k},i\omega)` to be copied

G_small
     Green function, which structure serves as a template

Returns
-------
out
     Green function with data from G_big, but structure of G_small""")

module.add_function ("triqstools::g_r_iw_t triqstools::copy_gf_to_smaller_intervall (triqstools::g_r_iw_cvt G_big, triqstools::g_r_iw_cvt G_small)", doc = r"""copies Green function into a Gf of smaller mesh size

Parameters
----------
G_big
     Green function :math:`G(\mathbf{r},i\omega)` to be copied

G_small
     Green function, which structure serves as a template

Returns
-------
out
     Green function with data from G_big, but structure of G_small""")

module.add_function ("triqstools::g_r_tau_t triqstools::copy_gf_to_smaller_intervall (triqstools::g_r_tau_cvt G_big, triqstools::g_r_tau_cvt G_small)", doc = r"""copies Green function into a Gf of smaller mesh size

Parameters
----------
G_big
     Green function :math:`G(\mathbf{r},\tau)` to be copied

G_small
     Green function, which structure serves as a template

Returns
-------
out
     Green function with data from G_big, but structure of G_small""")

module.add_function ("triqstools::g_iw_iW_t triqstools::copy_gf_to_smaller_intervall (triqstools::g_iw_iW_cvt G_big, triqstools::g_iw_iW_cvt G_small)", doc = r"""copies Green function into a Gf of smaller mesh size

Parameters
----------
G_big
     Green function :math:`G(i\omega,i\Omega)` to be copied

G_small
     Green function, which structure serves as a template

Returns
-------
out
     Green function with data from G_big, but structure of G_small""")

module.add_function ("triqstools::g_iw_t triqstools::copy_gf_to_smaller_intervall (triqstools::g_iw_cvt G_big, triqstools::g_iw_cvt G_small)", doc = r"""copies Green function into a Gf of smaller mesh size

Parameters
----------
G_big
     Green function :math:`G(i\omega)` to be copied

G_small
     Green function, which structure serves as a template

Returns
-------
out
     Green function with data from G_big, but structure of G_small""")

module.add_function ("triqstools::chi2_q_iW_t triqstools::copy_gf_to_smaller_intervall (triqstools::chi2_q_iW_cvt W_big, triqstools::chi2_q_iW_cvt W_small)", doc = r"""copies block Green function into a Gf of smaller mesh size

Parameters
----------
W_big
     Green function :math:`W^{\eta}(\mathbf{q},i\Omega)` to be copied

W_small
     Green function, which structure serves as a template

Returns
-------
out
     Green function with data from W_big, but structure of W_small""")

module.add_function ("triqstools::g_iw_t triqstools::evaluate_gf_on_mesh_same_beta (triqstools::g_iw_cvt gf, triqstools::iw_mesh_t mesh)", doc = r"""""")

module.add_function ("triqstools::g_iw_t triqstools::evaluate_gf_on_mesh (triqstools::g_iw_cvt gf, triqstools::iw_mesh_t mesh)", doc = r"""""")

module.add_function ("triqstools::g_k_iw_t triqstools::evaluate_gf_on_mesh (triqstools::g_k_iw_cvt gf, triqstools::k_iw_mesh_t mesh)", doc = r"""""")

module.add_function ("triqstools::g_r_iw_t triqstools::evaluate_gf_on_mesh (triqstools::g_r_iw_cvt gf, triqstools::r_iw_mesh_t mesh)", doc = r"""""")

module.add_function ("triqstools::g_iw_iW_t triqstools::evaluate_gf_on_mesh (triqstools::g_iw_iW_cvt gf, triqstools::iw_iW_mesh_t mesh)", doc = r"""""")

module.add_function ("triqstools::g_iw_t triqstools::make_gf_hermitian (triqstools::g_iw_cvt G)", doc = r"""""")

module.add_function ("triqstools::chi2_iW_t triqstools::make_gf_hermitian (triqstools::chi2_iW_cvt W)", doc = r"""""")

module.add_function ("triqstools::g_k_iw_t triqstools::make_gf_hermitian (triqstools::g_k_iw_cvt G)", doc = r"""""")

module.add_function ("triqstools::chi2_q_iW_t triqstools::make_gf_hermitian (triqstools::chi2_q_iW_cvt W)", doc = r"""""")

module.add_function ("bool triqstools::is_gf_hermitian (triqstools::g_k_iw_cvt G, double tolerance = 1E-8)", doc = r"""""")

module.add_function ("triqstools::g_k_iw_t triqstools::self_energy_weak_coupling (triqstools::g_k_iw_cvt G, triqstools::q_mesh_t q_mesh, triqstools::iW_mesh_t iW_mesh, double coupling)", doc = r"""""")

module.add_function ("triqstools::g_q_k_iw_t triqstools::self_energy_weak_coupling_bosonic_momentum (triqstools::g_k_iw_cvt G, triqstools::q_mesh_t q_mesh, triqstools::iW_mesh_t iW_mesh, double coupling)", doc = r"""""")

module.add_function ("triqstools::g_k_k_iw_t triqstools::self_energy_weak_coupling_fermionic_momentum (triqstools::g_k_iw_cvt G, triqstools::q_mesh_t q_mesh, triqstools::iW_mesh_t iW_mesh, double coupling)", doc = r"""""")



module.generate_code()