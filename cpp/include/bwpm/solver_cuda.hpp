// CUDA core of the FPWPM/BWPM in 2D and 3D (stage B4d, decisions E7 D4–D6, D8).
//
// Same algorithm, inputs and outputs as bwpm::Solver (cpp/include/bwpm/solver.hpp): planes
// p = 0 … nz, table 1 of wende2023, sparse reflection sources. Differences: FFTs with cuFFT
// (double complex), T and R with propagation factor computed inside the kernel instead of a
// cache (D6), labels uploaded per plane. Fields stay on the device during a run.
#pragma once

#include <cstdint>
#include <vector>

#include "bwpm/solver.hpp"

namespace bwpm {

// Number of usable CUDA devices (0 if none or no driver).
int cuda_device_count();

class CudaSolver {
 public:
  CudaSolver(Grid grid, std::vector<cplx> indices, std::vector<std::int32_t> labels,
             std::size_t nz, cplx index_before, std::vector<double> mask, int device);
  ~CudaSolver();
  CudaSolver(const CudaSolver&) = delete;
  CudaSolver& operator=(const CudaSolver&) = delete;

  RunResult run(const std::vector<cplx>& field_in, int iterations,
                const std::vector<std::size_t>& keep);

 private:
  struct Impl;
  Impl* impl_;
};

}  // namespace bwpm
