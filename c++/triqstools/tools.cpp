//
// Created by thomas on 13.11.17.
//
#include <triqs/utility/itertools.hpp>
#include <mpi/mpi.hpp>
#include <array>
#include <triqs/gfs.hpp>
#include <triqs/lattice/tight_binding.hpp>
#include "./types.hpp"
#include "./tools.hpp"

namespace triqstools {
  using std::pow;
  using std::sin;
  using triqs::utility::enumerate;
  using triqs::utility::product;
  using triqs::utility::zip;

  g_k_iw_t dyson_k_iw(g_k_iw_cvt G0, g_k_iw_cvt Sigma, double mu) {
    auto G = g_k_iw_t{G0};
    G[k_, iw_] << 1.0 / (1.0 / G0[k_, iw_] + mu - Sigma[k_, iw_]);
    return G;
  }

  g_k_iw_t dyson_k_iw(g_k_iw_cvt G0, g_iw_cvt Sigma, double mu) {
    auto G = g_k_iw_t{G0};
    G[k_, iw_] << 1.0 / (1.0 / G0[k_, iw_] + mu - Sigma[iw_]);
    return G;
  }

  chi2_q_iW_t dyson_k_iw(g_k_iw_cvt G0, chi2_iW_cvt Sigma, double mu) {
    auto G_bl = g_k_iw_t{G0};
    auto G    = make_block_gf({"up", "dn"}, {G_bl, G_bl});
    G[eta_][k_, iw_] << 1.0 / (1.0 / G0[k_, iw_] + mu - Sigma[eta_][iw_]);
    return G;
  }

  chi2_q_iW_t dyson_q_iW(array<double, 1> U, chi2_q_iW_cvt P) {
    auto W = chi2_q_iW_t{P};
    W[eta_][q_, iW_] << U(eta_) / (1.0 - U(eta_) * P[eta_][q_, iW_]);
    return W;
  }

  chi2_q_iW_t dyson_q_iW(chi2_q_cvt U, chi2_q_iW_cvt P) {
    auto W = chi2_q_iW_t{P};
    W[eta_][q_, iW_] << U[eta_][q_] / (1.0 - U[eta_][q_] * P[eta_][q_, iW_]);
    return W;
  }

  chi2_iW_t dyson_iW(array<double, 1> U, chi2_iW_cvt P) {
    auto W = chi2_iW_t{P};
    W[eta_][iW_] << U(eta_) / (1.0 - U(eta_) * P[eta_][iW_]);
    return W;
  }

  chi2_q_iW_t dyson_iW(chi2_q_cvt U, chi2_iW_cvt P) {
    auto W_bl = g_q_iW_t{{U[0].mesh(), P[0].mesh()}};
    auto W    = make_block_gf({"ch", "sp"}, {W_bl, W_bl});
    W[eta_][q_, iW_] << U[eta_][q_] / (1.0 - U[eta_][q_] * P[eta_][iW_]);
    return W;
  }

  g_iw_t inverse_dyson(g_iw_cvt G, g_iw_cvt Sigma) {
    auto G0 = g_iw_t{G};
    G0[iw_] << 1.0 / (1.0 / G[iw_] + Sigma[iw_]);
    return G0;
  }

  g_k_iw_t inverse_dyson(g_k_iw_cvt G, g_k_iw_cvt Sigma) {
    auto G0 = g_k_iw_t{G};
    G0[k_, iw_] << 1.0 / (1.0 / G[k_, iw_] + Sigma[k_, iw_]);
    return G0;
  }

  chi2_iW_t inverse_dyson(chi2_iW_cvt W, chi2_iW_cvt P) {
    auto W0 = chi2_iW_t{W};
    W0[eta_][iW_] << 1.0 / (1.0 / W[eta_][iW_] + P[eta_][iW_]);
    return W0;
  }

  chi2_q_iW_t inverse_dyson(chi2_q_iW_cvt W, chi2_q_iW_cvt P) {
    auto W0 = chi2_q_iW_t{W};
    W0[eta_][q_, iW_] << 1.0 / (1.0 / W[eta_][q_, iW_] + P[eta_][q_, iW_]);
    return W0;
  }

  g_iw_t make_local_gf(g_k_iw_cvt G) {
    auto const &[k_mesh, iw_mesh] = G.mesh();
    auto G_loc                    = g_iw_t{iw_mesh};
    for (auto const &iw : iw_mesh) G_loc[iw] = sum(G[k_, iw], k_ = k_mesh) / k_mesh.size();
    return G_loc;
  }

  std::complex<double> make_local_gf(g_iw_cvt G) { return sum(G[iw_], iw_ = G.mesh()) / G.mesh().domain().beta; }

  g_w_t make_local_gf(g_k_w_cvt G) {
    auto const &[k_mesh, w_mesh] = G.mesh();
    auto G_loc                   = g_w_t{w_mesh};
    for (auto const w : w_mesh) G_loc[w] = sum(G[k_, w], k_ = k_mesh) / k_mesh.size();
    return G_loc;
  }

  g_K_iw_t sum_over_k_mesh(g_K_k_iw_cvt G) {
    auto [K_mesh, k_mesh, iw_mesh] = G.mesh();

    auto G_loc = g_K_iw_t{{K_mesh, iw_mesh}};
    for (auto const &[K, iw] : G_loc.mesh()) G_loc[K, iw] = sum(G[K, k_, iw], k_ = k_mesh) / k_mesh.size();

    return G_loc;
  }

  chi2_iW_t make_local_gf(chi2_q_iW_cvt W) { return map_block_gf(static_cast<g_iw_t (*)(g_k_iw_cvt)>(&make_local_gf), W); }

  chi2_Q_iW_t sum_over_k_mesh(chi2_Q_q_iW_cvt W) { return map_block_gf(static_cast<g_K_iw_t (*)(g_K_k_iw_cvt)>(&sum_over_k_mesh), W); }

  chi3_iw_t make_Lambda_reg(chi3_iw_cvt chi3_tilde, g_iw_cvt G_imp, chi2_iW_cvt U_Weiss, chi2_iW_cvt chi2, chi2_iW_cvt l, double n) {
    auto const &iw_mesh = std::get<0>(chi3_tilde[0].mesh());
    double beta         = iw_mesh.domain().beta;

    auto Lambda_reg = chi3_iw_t{chi3_tilde};
    Lambda_reg[eta_][iw_, iW_] << -(chi3_tilde[eta_][iw_, iW_] - 2.0 * beta * G_imp[iw_] * kronecker(iW_) * kronecker(eta_, 0) * n)
             / (G_imp[iw_] * G_imp(iw_ + iW_) * (1. - U_Weiss[eta_][iW_] * chi2[eta_][iW_]))
          - l[eta_][iW_];

    return Lambda_reg;
  }

  chi3_Q_iw_t make_Lambda_reg(chi3_Q_iw_cvt chi3_tilde, g_k_iw_cvt G_imp, chi2_q_iW_cvt U_Weiss, chi2_q_iW_cvt chi2, chi2_q_iW_cvt l, g_k_cvt n_K) {
    all_t _;
    auto const &[K_mesh, Q_mesh, iw_mesh, iW_mesh] = chi3_tilde[0].mesh();
    double beta                                    = iw_mesh.domain().beta;

    auto Lambda_reg = chi3_Q_iw_t{chi3_tilde};
    // FIXME: are channels available also in eta-loop (e.g. if we want to distinguish charge channel in kronecker)?
    // FIXME: kronecker(Q,Idx(0,0))
    for (int eta = 0; eta < 2; ++eta) {
      Lambda_reg[eta]() = 0.;
      for (auto const &[K, Q] : product(K_mesh, Q_mesh)) {
        for (auto const &[iw, iW] : product(iw_mesh, iW_mesh)) {
          Lambda_reg[eta][K, Q, iw, iW] = -(chi3_tilde[eta][K, Q, iw, iW]
                                            - 2.0 * beta * G_imp[K, iw] * kronecker(iW) * Q_mesh.size() * kronecker(Q[0], 0) * kronecker(Q[1], 0)
                                               * kronecker(eta, 0) * triqstools::sum_k(n_K))
                / (G_imp[K, iw] * G_imp[K + Q, _](iw + iW) * (1. - U_Weiss[eta][Q, iW] * chi2[eta][Q, iW]))
             - l[eta][Q, iW];
        }
      }
    }

    return Lambda_reg;
  }

  chi2_q_iW_t make_chi_from_W(array<double, 1> U, chi2_q_iW_cvt W) {
    auto chi = chi2_q_iW_t{W};

    chi[eta_][q_, iW_] << -1. / U(eta_) * W[eta_][q_, iW_] * 1. / U(eta_) + 1. / U(eta_);

    return chi;
  }

  chi2_q_iW_t make_W_from_chi(array<double, 1> U, chi2_q_iW_cvt chi) {
    auto W = chi2_q_iW_t{chi};

    W[eta_][q_, iW_] << -U(eta_) * chi[eta_][q_, iW_] * U(eta_) + U(eta_);

    return W;
  }

  chi2_q_iW_t make_P_from_W(array<double, 1> U, chi2_q_iW_cvt W) {
    auto P = chi2_q_iW_t{W};

    P[eta_][q_, iW_] << 1. / U(eta_) - 1. / W[eta_][q_, iW_];

    return W;
  }

  chi2_q_iW_t make_chi_lambda(chi2_q_iW_cvt chi, array<double, 1> lambda) {
    auto chi_lambda = chi2_q_iW_t{chi};

    chi_lambda[eta_][q_, iW_] << 1. / (1. / chi[eta_][q_, iW_] + lambda(eta_));

    return chi_lambda;
  }

  chi2_q_iW_t make_chi_from_W(chi2_q_cvt U, chi2_q_iW_cvt W) {
    auto chi = chi2_q_iW_t{W};

    chi[eta_][q_, iW_] << -1. / U[eta_][q_] * W[eta_][q_, iW_] * 1. / U[eta_][q_] + 1. / U[eta_][q_];

    return chi;
  }

  g_iW_iw_iw_t make_F_updn_from_chi_updn(g_iW_iw_iw_cvt chi_updn, g_iw_cvt G) {
    auto const &iW_iw_iwp_mesh = chi_updn.mesh();

    auto F_updn = g_iW_iw_iw_t{chi_updn};
    F_updn()    = 0.;

    for (auto const &[iW, iw, iwp] : iW_iw_iwp_mesh) { F_updn[iW, iw, iwp] = -chi_updn[iW, iw, iwp] / (G(iw + iW) * G(iw) * G(iwp + iW) * G(iwp)); }

    return F_updn;
  }

  g_iW_iw_iw_mat_t make_chi_from_G2c(g_iW_iw_iw_mat_cvt G2c, g_iw_mat_cvt G) {
    double beta = G.domain().beta;

    auto chi = g_iW_iw_iw_mat_t{G2c};

    // Calculate generalized susceptibility in the ph channel from G2c and G
    //chi[iW_, iw_, iwp_][i_, j_, k_, l_] << G2c[iW_, iw_, iwp_][i_, j_, k_, l_] - beta * kronecker(iw_, iwp_) * G[iw_][l_, i_] * G(iwp_ + iW_)[j_, k_];

    return chi;
  }

  std::complex<double> sum_k(g_k_cvt G) {
    std::complex<double> summed = 0.;

    for (auto k : G.mesh()) summed += G[k];
    return summed / G.mesh().size();
  }

  std::complex<double> sum_k_iw(g_k_iw_cvt G) {
    auto const &[k_mesh, iw_mesh] = G.mesh();
    const double beta             = iw_mesh.domain().beta;

    std::complex<double> summed = 0.;

    for (auto [k, iw] : G.mesh()) summed += G[k, iw];
    return summed / (k_mesh.size() * beta);
  }

  std::complex<double> density(g_k_iw_cvt G) {
    g_iw_t G_loc = make_local_gf(G);
    return triqstools::density(G_loc);
  }

  std::complex<double> density(g_iw_cvt G) {
    auto known_moments = make_zero_tail(G, 2);
    known_moments(1)   = 1;
    auto [tail, err]   = triqs::gfs::fit_hermitian_tail(G, known_moments);
    return triqs::gfs::density(G, tail);
  }

  std::complex<double> calc_filling(g_k_iw_cvt G, double mu_diff) {
    if (std::abs(mu_diff) < 1E-15) { return triqstools::density(G); }

    auto G_adapted_mu = g_k_iw_t{G};
    for (auto [k, iw] : G_adapted_mu.mesh()) { G_adapted_mu[k, iw] = 1. / (1. / G_adapted_mu[k, iw] + mu_diff); }

    return triqstools::density(G_adapted_mu);
  }

  g_k_iw_t copy_gf_to_smaller_intervall(g_k_iw_cvt G_big, g_k_iw_cvt G_small) {
    auto G = g_k_iw_t{G_small};

    for (auto [k, iw] : G.mesh()) { G[k, iw] = G_big(k, iw); }
    return G;
  }

  g_r_iw_t copy_gf_to_smaller_intervall(g_r_iw_cvt G_big, g_r_iw_cvt G_small) {
    auto G = g_r_iw_t{G_small};

    for (auto [r, iw] : G.mesh()) { G[r, iw] = G_big(r, iw); }
    return G;
  }

  g_r_tau_t copy_gf_to_smaller_intervall(g_r_tau_cvt G_big, g_r_tau_cvt G_small) {
    auto G = g_r_tau_t{G_small};

    for (auto [r, tau] : G.mesh()) { G[r, tau] = G_big(r, tau); }
    return G;
  }

  g_iw_iW_t copy_gf_to_smaller_intervall(g_iw_iW_cvt G_big, g_iw_iW_cvt G_small) {
    auto G = g_iw_iW_t{G_small};

    for (auto [iw, iW] : G.mesh()) { G[iw, iW] = G_big(iw, iW); }
    return G;
  }

  g_iw_t copy_gf_to_smaller_intervall(g_iw_cvt G_big, g_iw_cvt G_small) {
    auto G = g_iw_t{G_small};

    for (auto iw : G.mesh()) { G[iw] = G_big(iw); }
    return G;
  }

  chi2_q_iW_t copy_gf_to_smaller_intervall(chi2_q_iW_cvt W_big, chi2_q_iW_cvt W_small) {
    auto W = chi2_q_iW_t{W_small};

    for (auto const &[W_eta, W_big_eta, W_small_eta] : zip(W, W_big, W_small)) { W_eta = copy_gf_to_smaller_intervall(W_big_eta, W_small_eta); }

    return W;
  }

  g_iw_t evaluate_gf_on_mesh_same_beta(g_iw_cvt gf, iw_mesh_t mesh) {
    auto G = g_iw_t{mesh};

    for (auto iw : mesh) { G[iw] = gf(iw); }

    return G;
  }

  g_iw_t evaluate_gf_on_mesh(g_iw_cvt gf, iw_mesh_t iw_mesh) {
    // assume that meshes have same statistics
    assert(gf.mesh().domain().statistic == iw_mesh.domain().statistic);
    const double beta_ratio = iw_mesh.domain().beta / gf.mesh().domain().beta;

    // if meshes are at same temperature: just evaluate
    if (abs(1. - beta_ratio) < 1E-14) { return evaluate_gf_on_mesh_same_beta(gf, iw_mesh); }

    auto gf_in          = g_iw_t{gf};
    auto tau_mesh       = make_adjoint_mesh(iw_mesh);
    auto tau_small_mesh = make_adjoint_mesh(gf_in.mesh());
    auto G_tau          = g_tau_t{tau_mesh};
    auto known_moments  = triqs::gfs::make_zero_tail(gf_in, 1);

    // fit and subtract constant offset
    gf_in.mesh().set_tail_fit_parameters(0.2, 30, 8);
    auto [tail, err] = triqs::gfs::fit_tail(gf_in);
    if (err > 1E-6 && !mpi::communicator().rank()) { std::cerr << "WARNING: offset fitting with error > 1E-6" << std::endl; }
    auto offset = tail(0).real();
    gf_in[iw_] << gf_in[iw_] - offset;

    // Fourier transform to imaginary time, evaluate, transform back
    auto gf_tau = make_gf_from_fourier(gf_in, tau_small_mesh, known_moments);
    for (auto tau : tau_mesh) { G_tau[tau] = gf_tau(tau / beta_ratio); }
    auto G = make_gf_from_fourier(G_tau, iw_mesh, triqs::gfs::make_zero_tail(G_tau));

    // add again constant offset
    G[iw_] << G[iw_] + offset;

    return G;
  }

  g_k_iw_t evaluate_gf_on_mesh(g_k_iw_cvt gf, k_iw_mesh_t mesh) {
    all_t _;
    auto G         = g_k_iw_t{mesh};
    auto G_large_k = g_k_iw_t{{std::get<0>(mesh), std::get<1>(gf.mesh())}};

    auto comm   = mpi::communicator();
    G_large_k() = 0.;
    for (auto iw : std::get<1>(gf.mesh())) {
      auto gf_slice = gf[_, iw];
      for (auto k : std::get<0>(mesh)) {
        if (k.linear_index() % comm.size() == comm.rank()) { G_large_k[k, iw] = gf_slice(k); }
      }
    }
    G_large_k = mpi::all_reduce(G_large_k, comm);

    G() = 0.;
    for (auto k : std::get<0>(mesh)) {
      if (k.linear_index() % comm.size() == comm.rank()) { G[k, _] = evaluate_gf_on_mesh(G_large_k[k, _], std::get<1>(mesh)); }
    }
    G = mpi::all_reduce(G, comm);

    return G;
  }

  g_r_iw_t evaluate_gf_on_mesh(g_r_iw_cvt gf, r_iw_mesh_t mesh) {
    all_t _;
    auto G         = g_r_iw_t{mesh};
    auto G_large_r = g_r_iw_t{{std::get<0>(mesh), std::get<1>(gf.mesh())}};

    auto comm   = mpi::communicator();
    G_large_r() = 0.;
    for (auto iw : std::get<1>(gf.mesh())) {
      auto gf_slice = gf[_, iw];
      for (auto r : std::get<0>(mesh)) {
        if (r.linear_index() % comm.size() == comm.rank()) { G_large_r[r, iw] = gf_slice(r); }
      }
    }
    G_large_r = mpi::all_reduce(G_large_r, comm);

    G() = 0.;
    for (auto r : std::get<0>(mesh)) {
      if (r.linear_index() % comm.size() == comm.rank()) { G[r, _] = evaluate_gf_on_mesh(G_large_r[r, _], std::get<1>(mesh)); }
    }
    G = mpi::all_reduce(G, comm);

    return G;
  }

  g_iw_iW_t evaluate_gf_on_mesh(g_iw_iW_cvt gf, iw_iW_mesh_t mesh) {
    all_t _;
    auto G          = g_iw_iW_t{mesh};
    auto G_large_iW = g_iw_iW_t{{std::get<0>(gf.mesh()), std::get<1>(mesh)}};

    auto comm    = mpi::communicator();
    G_large_iW() = 0.;
    for (auto iw : std::get<0>(gf.mesh())) {
      if (iw.linear_index() % comm.size() == comm.rank()) { G_large_iW[iw, _] = evaluate_gf_on_mesh(gf[iw, _], std::get<1>(mesh)); }
    }
    G_large_iW = mpi::all_reduce(G_large_iW, comm);

    G() = 0.;
    for (auto iW : std::get<1>(mesh)) {
      if (iW.linear_index() % comm.size() == comm.rank()) { G[_, iW] = evaluate_gf_on_mesh(G_large_iW[_, iW], std::get<0>(mesh)); }
    }
    G = mpi::all_reduce(G, comm);

    return G;
  }

  g_iw_t make_gf_hermitian(g_iw_cvt G) { return triqs::gfs::make_hermitian(G); }

  chi2_iW_t make_gf_hermitian(chi2_iW_cvt W) { return map_block_gf(static_cast<g_iw_t (*)(g_iw_cvt)>(&make_gf_hermitian), W); }

  g_k_iw_t make_gf_hermitian(g_k_iw_cvt G) {
    all_t _;
    auto G_hermitian = g_k_iw_t{G};

    auto comm     = mpi::communicator();
    G_hermitian() = 0.;
    for (auto k : std::get<0>(G.mesh())) {
      if (k.linear_index() % comm.size() == comm.rank()) { G_hermitian[k, _] = triqs::gfs::make_hermitian(G[k, _]); }
    }
    G_hermitian = mpi::all_reduce(G_hermitian, comm);

    return G_hermitian;
  }

  chi2_q_iW_t make_gf_hermitian(chi2_q_iW_cvt W) { return map_block_gf(static_cast<g_k_iw_t (*)(g_k_iw_cvt)>(&make_gf_hermitian), W); }

  bool is_gf_hermitian(g_k_iw_cvt G, double tolerance) {
    all_t _;
    bool hermitian = true;

    for (auto k : std::get<0>(G.mesh())) { hermitian = hermitian && triqs::gfs::is_gf_hermitian(G[k, _], tolerance); }

    return hermitian;
  }

  g_q_iW_t make_chi_ornstein_zernike(q_iW_mesh_t q_iW_mesh, double a, double xi, double gamma, double z, double Qx, double Qy) {
    auto chi = g_q_iW_t{q_iW_mesh};

    for (auto const &[q, iW] : q_iW_mesh) {
      chi[q, iW] = a
         / (4. * pow(sin((q[0] - Qx) / 2.), 2) + 4. * pow(sin((q[1] - Qy) / 2.), 2) + pow(std::abs(std::sqrt(iW * iW)), 2. / z) / gamma
            + pow(xi, -2));
    }
    return chi;
  }

  g_iw_t self_energy_from_dyson_schwinger(g_iW_iw_iw_cvt F_updn, g_iw_cvt G, double U) {
    auto const &iW_iw_iwp_mesh = F_updn.mesh();
    const double beta          = std::get<0>(iW_iw_iwp_mesh).domain().beta;

    auto Sigma = g_iw_t{G};
    Sigma()    = 0.;

    for (auto const &[iW, iw, iwp] : iW_iw_iwp_mesh) { Sigma[iw] += F_updn[iW, iw, iwp] * G(iwp) * G(iwp + iW) * G(iw + iW); }

    return -Sigma * U / (beta * beta);
  }

  g_k_iw_t self_energy_ornstein_zernike(g_k_iw_cvt G, q_iW_mesh_t q_iW_mesh, double a, double xi, double Qx, double Qy) {
    auto const &[k_mesh, iw_mesh] = G.mesh();
    auto const &[q_mesh, iW_mesh] = q_iW_mesh;
    const double beta             = iw_mesh.domain().beta;

    auto Sigma = g_k_iw_t{G};
    Sigma()    = 0.;

    auto comm = mpi::communicator();
    for (auto k : k_mesh) {
      if (k.linear_index() % comm.size() == comm.rank()) {
        for (auto const &iw : iw_mesh) {
          for (auto const &[q, iW] : q_iW_mesh) {
            Sigma[k, iw] += G(k + q, iw) * a
               / (4. * pow(sin((q[0] - Qx) / 2.), 2) + 4. * pow(sin((q[1] - Qy) / 2.), 2) + std::abs(std::sqrt(iW * iW)) + pow(xi, -2));
          }
          Sigma[k, iw] /= (beta * q_mesh.size());
        }
      }
    }
    Sigma = mpi::all_reduce(Sigma, comm);
    return Sigma;
  }

  g_k_iw_t self_energy_ornstein_zernike_iW(g_k_iw_cvt G, q_iW_mesh_t q_iW_mesh, double a, double xi, double Qx, double Qy) {
    auto const &[k_mesh, iw_mesh] = G.mesh();
    auto const &[q_mesh, iW_mesh] = q_iW_mesh;
    const double beta             = iw_mesh.domain().beta;

    auto Sigma = g_k_iw_t{G};
    Sigma()    = 0.;

    auto comm = mpi::communicator();
    for (auto k : k_mesh) {
      if (k.linear_index() % comm.size() == comm.rank()) {
        for (auto const &iw : iw_mesh) {
          for (auto const &[q, iW] : q_iW_mesh) {
            Sigma[k, iw] += G(k + q, iw) * a / (4. * pow(sin((q[0] - Qx) / 2.), 2) + 4. * pow(sin((q[1] - Qy) / 2.), 2) + iW + pow(xi, -2));
          }
          Sigma[k, iw] /= (beta * q_mesh.size());
        }
      }
    }
    Sigma = mpi::all_reduce(Sigma, comm);
    return Sigma;
  }
  g_k_iw_t self_energy_chi(g_k_iw_cvt G, g_q_iW_cvt chi, k_iw_mesh_t k_iw_mesh, double coupling) {
    auto const &k_mesh            = std::get<0>(k_iw_mesh);
    auto const &iw_mesh           = std::get<1>(k_iw_mesh);
    auto const &[q_mesh, iW_mesh] = chi.mesh();
    const double beta             = iw_mesh.domain().beta;

    auto Sigma = g_k_iw_t{{k_iw_mesh}};

    auto comm = mpi::communicator();
    for (auto k : k_mesh) {
      if (k.linear_index() % comm.size() == comm.rank()) {
        for (auto iw : iw_mesh) {
          Sigma[k, iw] = coupling * coupling * sum(sum(G(k + q_, iw + iW_) * chi(q_, iW_), q_ = q_mesh), iW_ = iW_mesh) / (beta * q_mesh.size());
        }
      }
    }
    Sigma = mpi::all_reduce(Sigma, comm);
    return Sigma;
  }

  g_k_iw_t self_energy_chi_restricted(g_k_iw_cvt G, g_q_iW_cvt chi, k_iw_mesh_t k_iw_mesh, double coupling, double xi) {
    auto const &k_mesh            = std::get<0>(k_iw_mesh);
    auto const &iw_mesh           = std::get<1>(k_iw_mesh);
    auto const &[q_mesh, iW_mesh] = chi.mesh();
    const double beta             = iw_mesh.domain().beta;

    auto Sigma = g_k_iw_t{{k_iw_mesh}};
    const double pi=std::atan(1.0)*4.0;

    auto comm = mpi::communicator();
    for (auto k : k_mesh) {
      if (k.linear_index() % comm.size() == comm.rank()) {
        for (auto iw : iw_mesh) {
          for (auto q : q_mesh) {
            if (std::abs(q[0] - pi) < 2. / std::abs(xi) && std::abs(q[1] - pi) < 2. / std::abs(xi)) {
              Sigma[k, iw] = coupling * coupling * sum(G(k + q, iw + iW_) * chi(q, iW_), iW_ = iW_mesh) / (beta * q_mesh.size());
            }
          }
        }
      }
    }
    Sigma = mpi::all_reduce(Sigma, comm);
    return Sigma;
  }

  g_k_iw_t self_energy_weak_coupling(g_k_iw_cvt G, q_mesh_t q_mesh, iW_mesh_t iW_mesh, double coupling) {
    auto const &[k_mesh, iw_mesh] = G.mesh();
    const double beta             = iw_mesh.domain().beta;

    auto Sigma = g_k_iw_t{G};
    Sigma()    = 0.;

    auto comm = mpi::communicator();
    for (auto k : k_mesh) {
      if (k.linear_index() % comm.size() == comm.rank()) {
        for (auto iw : iw_mesh) {
          Sigma[k, iw] = -2. * coupling * coupling
             * sum(sum(sum(sum(G(k + q_, iw + iW_) * G(kp_ + q_, iwp_ + iW_), q_ = q_mesh), iW_ = iW_mesh) * G[kp_, iwp_], iwp_ = iw_mesh),
                   kp_ = k_mesh)
             / (beta * beta * k_mesh.size() * q_mesh.size());
        }
      }
    }
    Sigma = mpi::all_reduce(Sigma, comm);
    return Sigma;
  }

  g_q_k_iw_t self_energy_weak_coupling_bosonic_momentum(g_k_iw_cvt G, q_mesh_t q_mesh, iW_mesh_t iW_mesh, double coupling) {
    auto const &[k_mesh, iw_mesh] = G.mesh();
    const double beta             = iw_mesh.domain().beta;

    auto Sigma = g_q_k_iw_t{{q_mesh, k_mesh, iw_mesh}};
    Sigma()    = 0.;

    auto comm = mpi::communicator();
    for (auto k : k_mesh) {
      if (k.linear_index() % comm.size() == comm.rank()) {
        for (auto iw : iw_mesh) {
          for (auto q : q_mesh) {
            Sigma[q, k, iw] = -2. * coupling * coupling
               * sum(sum(sum(G(k + q, iw + iW_) * G(kp_ + q, iwp_ + iW_), iW_ = iW_mesh) * G[kp_, iwp_], kp_ = k_mesh), iwp_ = iw_mesh)
               / (beta * beta * k_mesh.size());
          }
        }
      }
    }
    Sigma = mpi::all_reduce(Sigma, comm);
    return Sigma;
  }

  g_iW_k_iw_t self_energy_weak_coupling_bosonic_frequency(g_k_iw_cvt G, q_mesh_t q_mesh, iW_mesh_t iW_mesh, double coupling) {
    auto const &[k_mesh, iw_mesh] = G.mesh();
    const double beta             = iw_mesh.domain().beta;

    auto Sigma = g_iW_k_iw_t{{iW_mesh, k_mesh, iw_mesh}};
    Sigma()    = 0.;

    auto comm = mpi::communicator();
    for (auto k : k_mesh) {
      if (k.linear_index() % comm.size() == comm.rank()) {
        for (auto iw : iw_mesh) {
          for (auto iW : iW_mesh) {
            Sigma[iW, k, iw] = -2. * coupling * coupling
               * sum(sum(sum(G(k + q_, iw + iW) * G(kp_ + q_, iwp_ + iW), q_ = q_mesh) * G[kp_, iwp_], kp_ = k_mesh), iwp_ = iw_mesh)
               / (beta * q_mesh.size() * k_mesh.size());
          }
        }
      }
    }
    Sigma = mpi::all_reduce(Sigma, comm);
    return Sigma;
  }

  g_k_k_iw_t self_energy_weak_coupling_fermionic_momentum(g_k_iw_cvt G, q_mesh_t q_mesh, iW_mesh_t iW_mesh, double coupling) {
    auto const &[k_mesh, iw_mesh] = G.mesh();
    const double beta             = iw_mesh.domain().beta;

    auto Sigma = g_k_k_iw_t{{q_mesh, k_mesh, iw_mesh}};
    Sigma()    = 0.;

    auto comm = mpi::communicator();
    for (auto k : k_mesh) {
      if (k.linear_index() % comm.size() == comm.rank()) {
        for (auto iw : iw_mesh) {
          for (auto kp : k_mesh) {
            Sigma[kp, k, iw] = -2. * coupling * coupling
               * sum(sum(sum(G(k + q_, iw + iW_) * G(kp + q_, iwp_ + iW_), q_ = q_mesh), iW_ = iW_mesh) * G[kp, iwp_], iwp_ = iw_mesh)
               / (beta * beta * q_mesh.size());
          }
        }
      }
    }
    Sigma = mpi::all_reduce(Sigma, comm);
    return Sigma;
  }

  g_k_iw_t make_bubble_from_G(g_k_iw_t G) {
    auto G_r_tau = make_gf_from_fourier<0, 1>(G);

    auto [r_mesh, tau_mesh] = G_r_tau.mesh();
    double beta             = tau_mesh.domain().beta;

    auto tau_mesh_bosonic = gf_mesh<imtime>{beta, Boson, tau_mesh.size()};

    auto chi0 = g_r_tau_t{{r_mesh, tau_mesh_bosonic}};

    chi0[r_, tau_] << -G_r_tau(-r_, beta - tau_) * G_r_tau(r_, tau_);

    return make_gf_from_fourier<0, 1>(chi0);
  }

  g_k_iw_t bubble2(g_k_iw_cvt chi, g_k_iw_cvt g0) {

    // Fourier Transformation of k, \omega to obtain g(r,t)
    auto chirt = make_gf_from_fourier<0, 1>(chi);
    auto g0rt  = make_gf_from_fourier<0, 1>(g0);

    auto sigma = g0rt;
    sigma()    = 0;

    // The mesh of gtr is a cartesian product mt x mr. We decompose it.
    auto [mr, mt] = g0rt.mesh();

    // we fill sigma : sigma(r, tau) = chi(r, tau) * g0(r, tau)
    // I am not sure this is correct .......
    for (auto const &r : mr)
      for (auto const &t : mt) sigma[r, t] = chirt(r, t) * g0rt(r, t);

    // Fourier transform back to k, \omega space and return
    return make_gf_from_fourier<0, 1>(sigma);
  }

} // namespace triqstools
