// CPU core of the FPWPM/BWPM in 2D and 3D (stage B4c).
//
// Mirrors bwpm/bidirectional.py (2D) and bwpm/volume3d.py (3D); 2D is the case nx = 1.
// Specification: docs/spezifikation.md, sections 4 and 5; table 1 of wende2023 (i = 1 … 4).
// Planes p = 0 … nz, plane p belongs to slab p − 1 (plane 0 to `index_before`); fields (ny, nx, 2)
// in C order with components (E_x, E_y); centred grids; k in 1/mm, lengths in mm.
// Optional weak index fluctuations δn (nz, ny, nx) scatter with the vector Born approximation in
// both passes: forward increment and back scattering R1, increment R2 and forward scattering R3 of
// the backward field (stages B7b, B8b; docs/spezifikation.md §10, §11; bwpm/fluctuation.py).
#pragma once

#include <complex>
#include <cstdint>
#include <map>
#include <tuple>
#include <vector>

#include "bwpm/interface.hpp"

namespace bwpm {

struct Grid {
  std::size_t ny, nx;
  double dy, dx, dz, k0;
};

struct RunResult {
  std::vector<cplx> transmitted, reflected;          // (ny, nx, 2) after the last iteration
  std::vector<std::vector<cplx>> transmitted_iter;   // per iteration (energy residual)
  std::vector<std::vector<cplx>> reflected_iter;
  std::map<std::size_t, std::vector<cplx>> forward;  // requested planes
  std::map<std::size_t, std::vector<cplx>> backward;
  std::size_t stored_values = 0;                     // stored reflection source vectors
};

class Solver {
 public:
  // labels: (nz, ny, nx) material labels into `indices`; mask: (ny, nx) or empty (no mask);
  // fluctuation: δn (nz, ny, nx) or empty, within a slab only on pixels of one index (the host).
  Solver(Grid grid, std::vector<cplx> indices, std::vector<std::int32_t> labels, std::size_t nz,
         cplx index_before, std::vector<double> mask, int threads,
         std::vector<cplx> fluctuation = {});

  // iterations = 0: forward FPWPM only; keep: planes to return.
  RunResult run(const std::vector<cplx>& field_in, int iterations,
                const std::vector<std::size_t>& keep);

 private:
  enum Kind { FwdT = 0, FwdR = 1, BwdT = 2, BwdR = 3 };
  using Key = std::tuple<int, std::int32_t, std::int32_t>;

  const std::int32_t* plane_labels(std::size_t p) const;
  cplx index_of(std::int32_t label) const;
  const std::vector<cplx>& op(Kind kind, std::int32_t r, std::int32_t s);
  void fft(std::vector<cplx>& data, bool forward) const;
  // out = Σ_pairs I_r,s F⁻¹{M F{in}}; reflections skip equal indices; mask on transmissions.
  void apply(Kind kind, const std::vector<cplx>& in, std::size_t p_r, std::size_t p_s,
             std::vector<cplx>& out);
  // FFT over (ny, nx) of `ncomp` interleaved components.
  void fft_components(std::vector<cplx>& data, std::size_t ncomp, bool forward) const;
  // Born sources F{V·E_m} of slab `slab` at the two Gauss points (§11) for the incident field
  // `in` of `direction` (+1 at the left slab edge, −1 at the right edge); (2, ny, nx, 3).
  std::vector<cplx> born_sources(const std::vector<cplx>& in, int direction, std::size_t slab);
  // Field radiated by the sources in `direction`, at the right slab edge, real space (ny, nx, 2).
  std::vector<cplx> born_emission(const std::vector<cplx>& sources, int direction,
                                  std::size_t slab);

  Grid g_;
  std::vector<cplx> indices_;
  std::vector<std::int32_t> labels_;
  std::vector<std::int32_t> before_;
  std::size_t nz_;
  cplx index_before_;
  std::vector<double> mask_;
  int threads_;
  std::vector<double> kx_, ky_;  // (ny, nx) frequencies in numpy order
  std::vector<cplx> fluct_;       // δn (nz, ny, nx) or empty
  std::vector<char> has_host_;    // per slab: δn ≠ 0 somewhere
  std::vector<cplx> host_;        // per slab: host index ñ₀
  std::map<Key, std::vector<cplx>> cache_;
};

}  // namespace bwpm
