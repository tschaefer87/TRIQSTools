//
// Created by thomas on 13.11.17.
//

#pragma once
#include "./types.hpp"
#include <triqs/gfs/transform/fourier.hpp>

namespace triqstools {
  //FIXME: wrap template instantiation
  inline g_r_t make_gf_from_fourier_k(g_k_cvt G) { return triqs::gfs::make_gf_from_fourier<0>(G); }
  inline g_k_t make_gf_from_fourier_r(g_r_cvt G) { return triqs::gfs::make_gf_from_fourier<0>(G); }

  inline g_r_iw_t make_gf_from_fourier_k(g_k_iw_cvt G) { return triqs::gfs::make_gf_from_fourier<0>(G); }
  inline g_k_iw_t make_gf_from_fourier_r(g_r_iw_cvt G) { return triqs::gfs::make_gf_from_fourier<0>(G); }

  inline chi2_r_iW_t make_gf_from_fourier_k(chi2_q_iW_cvt W) {
    return map_block_gf(static_cast<g_r_iw_t (*)(g_k_iw_cvt)>(&make_gf_from_fourier_k), W);
  }
  inline chi2_q_iW_t make_gf_from_fourier_r(chi2_r_iW_cvt W) {
    return map_block_gf(static_cast<g_k_iw_t (*)(g_r_iw_cvt)>(&make_gf_from_fourier_r), W);
  }

  inline g_R_R_iw_iW_t make_gf_from_fourier_k(g_K_Q_iw_iW_cvt G) { return triqs::gfs::make_gf_from_fourier<0, 1>(G); }
  inline g_K_Q_iw_iW_t make_gf_from_fourier_r(g_R_R_iw_iW_cvt G) { return triqs::gfs::make_gf_from_fourier<0, 1>(G); }

  inline chi3_R_iw_t make_gf_from_fourier_k(chi3_Q_iw_cvt G) {
    return map_block_gf(static_cast<g_R_R_iw_iW_t (*)(g_K_Q_iw_iW_cvt)>(&make_gf_from_fourier_k), G);
  }
  inline chi3_Q_iw_t make_gf_from_fourier_r(chi3_R_iw_cvt G) {
    return map_block_gf(static_cast<g_K_Q_iw_iW_t (*)(g_R_R_iw_iW_cvt)>(&make_gf_from_fourier_r), G);
  }

  /**
   * Dyson equation for a one-particle lattice Green function
   *
   * :math:`G({\mathbf k},i\omega)=\frac{1}{i\omega+\mu-\epsilon_{\mathbf k}-\Sigma({\mathbf k},i\omega)}`
   *
   * @param G0 non-interacting Green function :math:`G_{0}(\mathbf{k},i\omega)`
   * @param Sigma self-energy :math:`\Sigma(\mathbf{k},i\omega)`
   * @param mu additional constant real shift (e.g., for chemical potential)
   * @return interacting one-particle fermionic lattice Green function :math:`G(\mathbf{k},i\omega)`
   * @remark
   */
  g_k_iw_t dyson_k_iw(g_k_iw_cvt G0, g_k_iw_cvt Sigma, double mu = 0.);

  /**
   * Dyson equation for a one-particle lattice Green function with local self-energy
   *
   * :math:`G({\mathbf k},i\omega)=\frac{1}{i\omega+\mu-\epsilon_{\mathbf k}-\Sigma(i\omega)}`
   *
   * @param G0 non-interacting Green function :math:`G_{0}(\mathbf{k},i\omega)`
   * @param Sigma local self-energy :math:`\Sigma(i\omega)`
   * @param mu additional constant real shift (e.g., for chemical potential)
   * @return interacting one-particle lattice Green function :math:`G(\mathbf{k},i\omega)`
   * @remark
   */
  g_k_iw_t dyson_k_iw(g_k_iw_cvt G0, g_iw_cvt Sigma, double mu = 0.);

  chi2_q_iW_t dyson_k_iw(g_k_iw_cvt G0, chi2_iW_cvt Sigma, double mu = 0.);
  /**
   * Dyson equation for a one-particle lattice block Green function
   *
   * :math:`W^{\eta}({\mathbf q},i\Omega)=\frac{U^{\eta}}{1-U^{\eta}P^{\eta}({\mathbf q},i\Omega)}`
   *
   * @param U bare propagator :math:`U^{\eta}`
   * @param P polarization :math:`P^{\eta}(\mathbf{q},i\Omega)`
   * @return interacting one-particle lattice block Green function :math:`W^{\eta}(i\Omega)`
   * @remark
   */
  chi2_q_iW_t dyson_q_iW(array<double, 1> U, chi2_q_iW_cvt P);

  chi2_q_iW_t dyson_q_iW(chi2_q_cvt U, chi2_q_iW_cvt P);

  /**
   * local Dyson equation for a one-particle lattice block Green function
   *
   * :math:`W^{\eta}(i\Omega)=\frac{U^{\eta}}{1-U^{\eta}P^{\eta}(i\Omega)}`
   *
   * @param U bare propagator :math:`U^{\eta}`
   * @param P local polarization :math:`P^{\eta}(i\Omega)`
   * @return interacting one-particle lattice block Green function :math:`W^{\eta}(\mathbf{q},i\Omega)`
   * @remark
   */
  chi2_iW_t dyson_iW(array<double, 1> U, chi2_iW_cvt P);

  chi2_q_iW_t dyson_iW(chi2_q_cvt U, chi2_iW_cvt P);

  /**
   * inverse local Dyson equation
   *
   * :math:`G_{0}(i\omega) = \left(G^{-1}(i\omega) + \Sigma(i\omega)\right)^{-1}`
   *
   * @param G interacting one-particle Green function :math:`G(i\omega)`
   * @param Sigma local self-energy :math:`\Sigma(i\omega)`
   * @return non-interacting Green function :math:`G_{0}(i\omega)`
   * @remark
   */
  g_iw_t inverse_dyson(g_iw_cvt G, g_iw_cvt Sigma);

  /**
   * inverse lattice Dyson equation
   *
   * :math:`G_{0}(\mathbf{k},i\omega) = \left(G^{-1}(\mathbf{k},i\omega) + \Sigma(\mathbf{k},i\omega)\right)^{-1}`
   *
   * @param G interacting one-particle Green function :math:`G(\mathbf{k},i\omega)`
   * @param Sigma self-energy :math:`\Sigma(\mathbf{k},i\omega)`
   * @return non-interacting Green function :math:`G_{0}(\mathbf{k},i\omega)`
   * @remark
   */
  g_k_iw_t inverse_dyson(g_k_iw_cvt G, g_k_iw_cvt Sigma);

  /**
   * inverse local Dyson equation
   *
   * :math:`W_{0}(i\omega) = \left(W^{-1}(i\Omega) + P(i\Omega)\right)^{-1}`
   *
   * @param G interacting one-particle block Green function :math:`W^{\eta}(i\Omega)`
   * @param Sigma local self-energy :math:`P^{\eta}(i\Omega)`
   * @return non-interacting Green function :math:`W_{0}^{\eta}(i\Omega)`
   * @remark
   */
  chi2_iW_t inverse_dyson(chi2_iW_cvt W, chi2_iW_cvt P);

  /**
   * inverse Dyson equation
   *
   * :math:`W_{0}({\mathbf q},i\omega) = \left(W^{-1}({\mathbf q},i\Omega) + P({\mathbf q},i\Omega)\right)^{-1}`
   *
   * @param G interacting one-particle block Green function :math:`W^{\eta}({\mathbf q},i\Omega)`
   * @param Sigma local self-energy :math:`P^{\eta}({\mathbf q},i\Omega)`
   * @return non-interacting Green function :math:`W_{0}^{\eta}({\mathbf q},i\Omega)`
   * @remark
   */
  chi2_q_iW_t inverse_dyson(chi2_q_iW_cvt W, chi2_q_iW_cvt P);

  /**
   * localizes the lattice Green function by a normalized sum over the k-mesh
   *
   * @param G lattice Green function :math:`G(\mathbf{k},i\omega)`
   * @return local Green function :math:`G(i\omega)`
   * @remark
   */
  g_iw_t make_local_gf(g_k_iw_cvt G);

  std::complex<double> make_local_gf(g_iw_cvt G);

  g_w_t make_local_gf(g_k_w_cvt G);

  /**
   * (normalized) sum over the second of two Brillouin zone meshes of a theta-weighted Green function
   *
   * @param G theta-weighted Green function :math:`G(\mathbf{K},\mathbf{k},i\omega)`
   * @return cluster Green function :math:`G(\mathbf{K},i\omega)`
   * @remark
   */
  g_K_iw_t sum_over_k_mesh(g_K_k_iw_cvt G);

  std::complex<double> sum_k(g_k_cvt G);

  /**
   * localizes the lattice block Green function by a normalized sum over the q-mesh
   *
   * @param G lattice block Green function :math:`W^{\eta}(\mathbf{q},i\Omega)`
   * @return local block Green function :math:`W^{\eta}(i\Omega)`
   * @remark
   */
  chi2_iW_t make_local_gf(chi2_q_iW_cvt W);

  /**
   * (normalized) sum over the second of two Brillouin zone meshes of a theta-weighted block Green function
   *
   * @param W theta-weighted block Green function :math:`W^{\eta}(\mathbf{Q},\mathbf{q},i\Omega)`
   * @return cluster block Green function :math:`W^{\eta}(\mathbf{Q},i\Omega)`
   * @remark
   */
  chi2_Q_iW_t sum_over_k_mesh(chi2_Q_q_iW_cvt W);

  /**
   * extracts the regular part of the electron-boson coupling vertex :math:`\Lambda_{reg}^{\eta}(i\omega,i\Omega)`
   * from the equal-time correlation function :math:`\tilde{\chi}^{\eta}(i\omega,i\Omega)`
   * T. Ayral and O. Parcollet,  Phys. Rev. B 93, 235124 (2016), Eq. (76) combined with Eq. (66)
   *
   * @param chi3_tilde the equal-time correlation function :math:`\tilde{\chi}^{\eta}(i\omega,i\Omega)`
   * @param G_imp the local Green function (e.g. impurity model Green function) :math:`G(i\omega)`
   * @param U_Weiss the bosonic Weiss field :math:`U^{\eta}(i\Omega)`
   * @param chi2 the local connected density-density correlation function :math:`\chi^{\eta}(i\Omega)`
   * @param l part of the vertex, so that :math:`\Lambda^{\eta}(i\omega,i\Omega)=\Lambda^{\eta}_{reg}(i\omega,i\Omega)+l^{\eta}(i\Omega)`
   * @return regular part of the electron-boson coupling vertex :math:`Lambda^{\eta}_{reg}(i\omega,i\Omega)`
   * @remark
   */
  chi3_iw_t make_Lambda_reg(chi3_iw_cvt chi3_tilde, g_iw_cvt G_imp, chi2_iW_cvt U_Weiss, chi2_iW_cvt chi2, chi2_iW_cvt l, double n);

  /**
   * extracts the regular part of the electron-boson coupling vertex :math:`\Lambda_{reg}^{\eta}(\mathbf{K},\mathbf{Q},i\omega,i\Omega)`
   * from the equal-time correlation function :math:`\tilde{\chi}^{\eta}(\mathbf{K},\mathbf{Q},i\omega,i\Omega)`
   * T. Ayral and O. Parcollet,  Phys. Rev. Lett. 119, 166401 (2017)
   *
   * @param chi3_tilde the equal-time correlation function :math:`\tilde{\chi}^{\eta}(\mathbf{K},\mathbf{Q},i\omega,i\Omega)`
   * @param G_imp the cluster impurity Green function (e.g. impurity model Green function) :math:`G(\mathbf{K},i\omega)`
   * @param U_Weiss the bosonic Weiss field :math:`U^{\eta}(\mathbf{Q},i\Omega)`
   * @param chi2 the cluster connected density-density correlation function :math:`\chi^{\eta}(\mathbf{Q},i\Omega)`
   * @param l part of the vertex, so that :math:`\Lambda^{\eta}(\mathbf{K},\mathbf{Q},i\omega,i\Omega)=\Lambda^{\eta}_{reg}(\mathbf{K},\mathbf{Q},i\omega,i\Omega)+l^{\eta}(\mathbf{Q},i\Omega)`
   * @return regular part of the electron-boson coupling vertex :math:`Lambda^{\eta}_{reg}(\mathbf{K},\mathbf{Q},i\omega,i\Omega)`
   * @remark
   */
  chi3_Q_iw_t make_Lambda_reg(chi3_Q_iw_cvt chi3_tilde, g_k_iw_cvt G_imp, chi2_q_iW_cvt U_Weiss, chi2_q_iW_cvt chi2, chi2_q_iW_cvt l, g_k_cvt n_K);

  chi2_q_iW_t make_chi_lambda(chi2_q_iW_cvt chi, array<double, 1> lambda);

  chi2_q_iW_t make_chi_from_W(array<double, 1> U, chi2_q_iW_cvt W);

  chi2_q_iW_t make_chi_from_W(chi2_q_cvt U, chi2_q_iW_cvt W);

  chi2_q_iW_t make_P_from_W(array<double, 1> U, chi2_q_iW_cvt W);

  chi2_q_iW_t make_W_from_chi(array<double, 1> U, chi2_q_iW_cvt chi);

  g_q_iW_t make_chi_ornstein_zernike(q_iW_mesh_t q_iW_mesh, double a, double xi, double gamma, double Qx, double Qy);
  /**
   * calculates the density for a given lattice Green function
   *
   * @param G Green function :math:`G(\mathbf{k},i\omega)` for which the density is determined. The necessary tail will automatically be determined.
   * @return density, half-filling corresponds to 0.5
   * @remark
   */
  std::complex<double> density(g_k_iw_cvt G);

  std::complex<double> density(g_iw_cvt G);
  /**
   * (normalized) sum over both, the Brillouin zone and Matsubara mesh of a given lattice Green function.
   * Be careful: no tail is used.
   *
   * @param G Green function :math:`G(\mathbf{k},i\omega)`.
   * @return (normalized) sum over both meshes without tail
   * @remark
   */
  std::complex<double> sum_k_iw(g_k_iw_cvt G);

  /**
   * calculates the density for a given lattice Green function including a shift in the chemical potential
   *
   * @param G Green function :math:`G(\mathbf{k},i\omega)` for which the density is determined. The necessary tail will automatically be determined.
   * @param mu_diff additional shift in the chemical potential
   * @return density, half-filling corresponds to 0.5
   * @remark
   */
  std::complex<double> calc_filling(g_k_iw_cvt G, double mu_diff);

  /**
   * copies Green function into a Gf of smaller mesh size
   *
   * @param G_big Green function :math:`G(\mathbf{k},i\omega)` to be copied
   * @param G_small Green function, which structure serves as a template
   * @return Green function with data from G_big, but structure of G_small
   * @remark
   */
  g_k_iw_t copy_gf_to_smaller_intervall(g_k_iw_cvt G_big, g_k_iw_cvt G_small);

  /**
   * copies Green function into a Gf of smaller mesh size
   *
   * @param G_big Green function :math:`G(\mathbf{r},i\omega)` to be copied
   * @param G_small Green function, which structure serves as a template
   * @return Green function with data from G_big, but structure of G_small
   * @remark
   */
  g_r_iw_t copy_gf_to_smaller_intervall(g_r_iw_cvt G_big, g_r_iw_cvt G_small);

  /**
   * copies Green function into a Gf of smaller mesh size
   *
   * @param G_big Green function :math:`G(\mathbf{r},\tau)` to be copied
   * @param G_small Green function, which structure serves as a template
   * @return Green function with data from G_big, but structure of G_small
   * @remark
   */
  g_r_tau_t copy_gf_to_smaller_intervall(g_r_tau_cvt G_big, g_r_tau_cvt G_small);

  /**
  * copies Green function into a Gf of smaller mesh size
  *
  * @param G_big Green function :math:`G(i\omega,i\Omega)` to be copied
  * @param G_small Green function, which structure serves as a template
  * @return Green function with data from G_big, but structure of G_small
  * @remark
  */
  g_iw_iW_t copy_gf_to_smaller_intervall(g_iw_iW_cvt G_big, g_iw_iW_cvt G_small);

  /**
  * copies Green function into a Gf of smaller mesh size
  *
  * @param G_big Green function :math:`G(i\omega)` to be copied
  * @param G_small Green function, which structure serves as a template
  * @return Green function with data from G_big, but structure of G_small
  * @remark
  */
  g_iw_t copy_gf_to_smaller_intervall(g_iw_cvt G_big, g_iw_cvt G_small);

  /**
  * copies block Green function into a Gf of smaller mesh size
  *
  * @param W_big Green function :math:`W^{\eta}(\mathbf{q},i\Omega)` to be copied
  * @param W_small Green function, which structure serves as a template
  * @return Green function with data from W_big, but structure of W_small
  * @remark
  */
  chi2_q_iW_t copy_gf_to_smaller_intervall(chi2_q_iW_cvt W_big, chi2_q_iW_cvt W_small);

  g_iw_t evaluate_gf_on_mesh_same_beta(g_iw_cvt gf, iw_mesh_t mesh);

  g_iw_t evaluate_gf_on_mesh(g_iw_cvt gf, iw_mesh_t mesh);

  g_k_iw_t evaluate_gf_on_mesh(g_k_iw_cvt gf, k_iw_mesh_t mesh);

  g_r_iw_t evaluate_gf_on_mesh(g_r_iw_cvt gf, r_iw_mesh_t mesh);

  g_iw_iW_t evaluate_gf_on_mesh(g_iw_iW_cvt gf, iw_iW_mesh_t mesh);

  g_iw_t make_gf_hermitian(g_iw_cvt G);

  chi2_iW_t make_gf_hermitian(chi2_iW_cvt W);

  g_k_iw_t make_gf_hermitian(g_k_iw_cvt G);

  chi2_q_iW_t make_gf_hermitian(chi2_q_iW_cvt W);

  bool is_gf_hermitian(g_k_iw_cvt G, double tolerance = 1E-8);

  g_k_iw_t make_bubble_from_G(g_k_iw_t G);

  g_k_iw_t self_energy_ornstein_zernike(g_k_iw_cvt G, q_iW_mesh_t q_iW_mesh, double a, double xi, double Qx, double Qy);

  g_k_iw_t self_energy_ornstein_zernike_iW(g_k_iw_cvt G, q_iW_mesh_t q_iW_mesh, double a, double xi, double Qx, double Qy);

  g_k_iw_t self_energy_chi(g_k_iw_cvt G, g_q_iW_cvt chi, k_iw_mesh_t k_iw_mesh, double coupling);

  g_k_iw_t self_energy_chi_restricted(g_k_iw_cvt G, g_q_iW_cvt chi, k_iw_mesh_t k_iw_mesh, double coupling, double xi);

  g_k_iw_t bubble2(g_k_iw_cvt chi, g_k_iw_cvt g0);

  g_k_iw_t self_energy_weak_coupling(g_k_iw_cvt G, q_mesh_t q_mesh, iW_mesh_t iW_mesh, double coupling);

  g_iW_k_iw_t self_energy_weak_coupling_bosonic_frequency(g_k_iw_cvt G, q_mesh_t q_mesh, iW_mesh_t iW_mesh, double coupling);

  g_q_k_iw_t self_energy_weak_coupling_bosonic_momentum(g_k_iw_cvt G, q_mesh_t q_mesh, iW_mesh_t iW_mesh, double coupling);

  g_k_k_iw_t self_energy_weak_coupling_fermionic_momentum(g_k_iw_cvt G, q_mesh_t q_mesh, iW_mesh_t iW_mesh, double coupling);
} // namespace triqstools
