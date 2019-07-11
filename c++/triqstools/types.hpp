/***************************************************************************
 *
 *    file:          types.hpp
 *    contents:      definitions of types
 *
 ***************************************************************************/

#pragma once
#include <triqs/gfs.hpp>
#include <vector>
#include <utility>

namespace triqstools {
  using namespace triqs::gfs;
  using namespace triqs::arrays;

  // Type to describe the index- and block-structure of a Green function
  using gf_struct_t = std::vector<std::pair<std::string /*block name*/, std::vector<std::string> /*indices in block*/>>;

  //======= Aliases for meshes ========================================
  using iw_mesh_t = gf_mesh<imfreq>;
  using iW_mesh_t = iw_mesh_t;

  using k_mesh_t = gf_mesh<brillouin_zone>;

  using k_iw_mesh_t = gf_mesh<cartesian_product<brillouin_zone, imfreq>>;
  using q_iW_mesh_t = k_iw_mesh_t;

  using r_iw_mesh_t = gf_mesh<cartesian_product<cyclic_lattice, imfreq>>;
  using r_iW_mesh_t = r_iw_mesh_t;

  using iw_iW_mesh_t = gf_mesh<cartesian_product<imfreq, imfreq>>;

  //======= Aliases for the different Green function container's
  // Container type of static single-particle fermionic propagator G(k)
  using g_k_t   = gf<brillouin_zone, scalar_valued>;
  using g_k_cvt = g_k_t::const_view_type;

  // Container type of static single-particle fermionic propagator G(r)
  using g_r_t   = gf<cyclic_lattice, scalar_valued>;
  using g_r_cvt = g_r_t::const_view_type;

  // Container type of local single-particle fermionic propagator G(iw)
  using g_iw_t   = gf<imfreq, scalar_valued>;
  using g_iw_cvt = g_iw_t::const_view_type;

  // Container type of local single-particle fermionic propagator G(w)
  using g_w_t   = gf<refreq, scalar_valued>;
  using g_w_cvt = g_w_t::const_view_type;

  // Container type of local single-particle fermionic propagator G(tau)
  using g_tau_t   = gf<imtime, scalar_valued>;
  using g_tau_cvt = g_tau_t::const_view_type;

  // Container type of local single-particle bosonic propagator chi(iW)
  using g_iW_t   = g_iw_t;
  using g_iW_cvt = g_iw_cvt;

  // Container type of local boson-fermion vertex Lambda(iw,iW)
  using g_iw_iW_t   = gf<cartesian_product<imfreq, imfreq>, scalar_valued>;
  using g_iw_iW_cvt = g_iw_iW_t::const_view_type;

  // Container type of single-particle fermionic propagator G(k,iw)
  using g_k_iw_t   = gf<cartesian_product<brillouin_zone, imfreq>, scalar_valued>;
  using g_k_iw_cvt = g_k_iw_t::const_view_type;

  // Container type of single-particle fermionic propagator G(k,w)
  using g_k_w_t   = gf<cartesian_product<brillouin_zone, refreq>, scalar_valued>;
  using g_k_w_cvt = g_k_w_t::const_view_type;

  // Container type of single-particle fermionic propagator G(k,iw)
  using g_K_iw_t   = g_k_iw_t;
  using g_K_iw_cvt = g_k_iw_cvt;

  // Container type of single-particle fermionic propagator G(K,r,tau)
  using g_K_r_tau_t   = gf<cartesian_product<brillouin_zone, cyclic_lattice, imtime>, scalar_valued>;
  using g_K_r_tau_cvt = g_K_r_tau_t::const_view_type;

  // Container type of single-particle fermionic propagator G(K,r,iw)
  using g_K_r_iw_t   = gf<cartesian_product<brillouin_zone, cyclic_lattice, imfreq>, scalar_valued>;
  using g_K_r_iw_cvt = g_K_r_iw_t::const_view_type;

  // Container type of single-particle fermionic theta-weighted propagator G(K,k,iw)
  using g_K_k_iw_t   = gf<cartesian_product<brillouin_zone, brillouin_zone, imfreq>, scalar_valued>;
  using g_K_k_iw_cvt = g_K_k_iw_t::const_view_type;

  // Container type of boson-electron vertex Lambda(K,Q,iw,iW)
  using g_K_Q_iw_iW_t   = gf<cartesian_product<brillouin_zone, brillouin_zone, imfreq, imfreq>, scalar_valued>;
  using g_K_Q_iw_iW_cvt = g_K_Q_iw_iW_t::const_view_type;

  // Container type of boson-electron vertex Lambda(R,R,iw,iW)
  using g_R_R_iw_iW_t   = gf<cartesian_product<cyclic_lattice, cyclic_lattice, imfreq, imfreq>, scalar_valued>;
  using g_R_R_iw_iW_cvt = g_R_R_iw_iW_t::const_view_type;

  // Container type of single-particle fermionic propagator G(r,tau)
  using g_r_tau_t   = gf<cartesian_product<cyclic_lattice, imtime>, scalar_valued>;
  using g_r_tau_cvt = g_r_tau_t::const_view_type;

  // Container type of single-particle bosonic propagator G(q,iW)
  using g_q_iW_t   = g_k_iw_t;
  using g_q_iW_cvt = g_k_iw_cvt;

  // Container type of single-particle fermionic propagator G(r,tau)
  using g_r_tau_t   = gf<cartesian_product<cyclic_lattice, imtime>, scalar_valued>;
  using g_r_tau_cvt = g_r_tau_t::const_view_type;

  // Container type of single-particle fermionic propagator G(r,iw)
  using g_r_iw_t   = gf<cartesian_product<cyclic_lattice, imfreq>, scalar_valued>;
  using g_r_iw_cvt = g_r_iw_t::const_view_type;

  // Container type of single-particle bosonic propagator G(r,iW)
  using g_r_iW_t   = g_r_iw_t;
  using g_r_iW_cvt = g_r_iw_cvt;

  // Container type of local single-particle bosonic propagator W[channel](iW)
  using chi2_iW_t   = block_gf<imfreq, scalar_valued>;
  using chi2_iW_cvt = chi2_iW_t::const_view_type;

  // Container type of single-particle static bosonic propagator W[channel](q)
  using chi2_q_t   = block_gf<brillouin_zone, scalar_valued>;
  using chi2_q_cvt = chi2_q_t::const_view_type;

  // Container type of single-particle bosonic propagator W[channel](q,iW)
  using chi2_q_iW_t   = block_gf<cartesian_product<brillouin_zone, imfreq>, scalar_valued>;
  using chi2_q_iW_cvt = chi2_q_iW_t::const_view_type;

  // Container type of single-particle bosonic propagator W[channel](q,iW)
  using chi2_Q_iW_t   = chi2_q_iW_t;
  using chi2_Q_iW_cvt = chi2_q_iW_cvt;

  // Container type of single-particle bosonic theta-weighted propagator W[channel](Q,q,iW)
  using chi2_Q_q_iW_t   = block_gf<cartesian_product<brillouin_zone, brillouin_zone, imfreq>, scalar_valued>;
  using chi2_Q_q_iW_cvt = chi2_Q_q_iW_t::const_view_type;

  // Container type of single-particle bosonic propagator W[channel](k,iW)
  using chi2_r_iW_t   = block_gf<cartesian_product<cyclic_lattice, imfreq>, scalar_valued>;
  using chi2_r_iW_cvt = chi2_r_iW_t::const_view_type;

  // Container type of local fermion-boson vertex Lambda[channel](iw,iW)
  using chi3_iw_t   = block_gf<cartesian_product<imfreq, imfreq>, scalar_valued>;
  using chi3_iw_cvt = chi3_iw_t::const_view_type;

  // Container type of fermion-boson vertex Lambda[channel](k,q,iw,iW)
  using chi3_Q_iw_t   = block_gf<cartesian_product<brillouin_zone, brillouin_zone, imfreq, imfreq>, scalar_valued>;
  using chi3_Q_iw_cvt = chi3_Q_iw_t::const_view_type;

  // Container type of fermion-boson vertex Lambda[channel](k,q,iw,iW)
  using chi3_R_iw_t   = block_gf<cartesian_product<cyclic_lattice, cyclic_lattice, imfreq, imfreq>, scalar_valued>;
  using chi3_R_iw_cvt = chi3_R_iw_t::const_view_type;

  // Declare some placeholders for the rest of the code. Use anonymous namespace for proper linkage
  // in this code, all variables with trailing _ are placeholders by convention.
  namespace {
    triqs::clef::placeholder<0> iw_;
    triqs::clef::placeholder<1> k_;
    triqs::clef::placeholder<2> eta_;
    triqs::clef::placeholder<3> q_;
    triqs::clef::placeholder<4> iW_;
    triqs::clef::placeholder<5> r_;
    triqs::clef::placeholder<6> tau_;
    triqs::clef::placeholder<7> K_;
    triqs::clef::placeholder<8> Q_;
    triqs::clef::placeholder_prime<0> iwp_;
    triqs::clef::placeholder_prime<1> iWp_;
    triqs::clef::placeholder_prime<2> kp_;
  } // anonymous namespace
} // namespace trilex
