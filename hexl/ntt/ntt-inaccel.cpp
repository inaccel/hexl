// Copyright (C) 2018-2022 InAccel
// SPDX-License-Identifier: Apache-2.0

#include <inaccel/coral>

#include "hexl/ntt/ntt.hpp"
#include "hexl/util/check.hpp"

namespace intel {
namespace hexl {

void ForwardTransformToBitReverseInaccel(
    uint64_t* result, const uint64_t* operand, uint64_t n,
    const uint64_t *modulus, const uint64_t* root_of_unity_powers,
    const uint64_t* precon_root_of_unity_powers, uint64_t input_mod_factor,
    uint64_t output_mod_factor, int batch) {
  HEXL_CHECK(NTT::CheckArguments(n, modulus), "");
  HEXL_CHECK(n == 16834, "Require n == 16834");
  HEXL_CHECK(modulus != nullptr, "modulus == nullptr");
  HEXL_CHECK_BOUNDS(operand, batch * n, batch * modulus * input_mod_factor,
                "operand exceeds bound " << batch * modulus * input_mod_factor);
  HEXL_CHECK(root_of_unity_powers != nullptr,
             "root_of_unity_powers == nullptr");
  HEXL_CHECK(precon_root_of_unity_powers != nullptr,
             "precon_root_of_unity_powers == nullptr");
  HEXL_CHECK(
      input_mod_factor == 1 || input_mod_factor == 2 || input_mod_factor == 4,
      "input_mod_factor must be 1, 2, or 4; got " << input_mod_factor);
  HEXL_UNUSED(input_mod_factor);
  HEXL_CHECK(output_mod_factor == 1,
             "output_mod_factor must be 1 for FPGA; got " << output_mod_factor);
  HEXL_UNUSED(output_mod_factor);

  inaccel::request ntt("hexl.ntt.ForwardTransformToBitReverseRadix2");
  ntt.arg(batch)
      .arg_array<uint64_t>(operand, operand + batch * n)
      .arg_array<uint64_t>(modulus, modulus + 1)
      .arg_array<uint64_t>(root_of_unity_powers, root_of_unity_powers + n)
      .arg_array<uint64_t>(precon_root_of_unity_powers,
          precon_root_of_unity_powers + n)
      .arg_array<uint64_t>(result, result + batch * n);
  inaccel::submit(ntt).get();
}

void InverseTransformFromBitReverseInaccel(
    uint64_t* result, const uint64_t* operand, uint64_t n, uint64_t *modulus,
    const uint64_t* inv_root_of_unity_powers,
    const uint64_t* precon_inv_root_of_unity_powers, uint64_t input_mod_factor,
    uint64_t output_mod_factor, uint64_t *inv_n, uint64_t *inv_n_w, int batch) {
  HEXL_CHECK(NTT::CheckArguments(n, modulus), "");
  HEXL_CHECK(n == 16834, "Require n == 16834");
  HEXL_CHECK(modulus != nullptr, "modulus == nullptr");
  HEXL_CHECK(inv_root_of_unity_powers != nullptr,
             "inv_root_of_unity_powers == nullptr");
  HEXL_CHECK(precon_inv_root_of_unity_powers != nullptr,
             "precon_inv_root_of_unity_powers == nullptr");
  HEXL_CHECK(operand != nullptr, "operand == nullptr");
  HEXL_CHECK(input_mod_factor == 1 || input_mod_factor == 2,
             "input_mod_factor must be 1 or 2; got " << input_mod_factor);
  HEXL_UNUSED(input_mod_factor);
  HEXL_CHECK(output_mod_factor == 1 || output_mod_factor == 2,
             "output_mod_factor must be 1 or 2; got " << output_mod_factor);
  HEXL_UNUSED(output_mod_factor);
  HEXL_CHECK(inv_n != nullptr, "inv_n == nullptr");
  HEXL_CHECK(inv_n_w != nullptr, "inv_n_w == nullptr");

  inaccel::request ntt("hexl.ntt.InverseTransformFromBitReverseRadix2");
  ntt.arg(batch)
      .arg_array<uint64_t>(operand, operand + batch * n)
      .arg_array<uint64_t>(modulus, modulus + 1)
      .arg_array<uint64_t>(inv_n, inv_n + 1)
      .arg_array<uint64_t>(inv_n_w, inv_n_w + 1)
      .arg_array<uint64_t>(inv_root_of_unity_powers,
          inv_root_of_unity_powers + n)
      .arg_array<uint64_t>(precon_inv_root_of_unity_powers,
          precon_inv_root_of_unity_powers + n)
      .arg_array<uint64_t>(result, result + batch * n);
  inaccel::submit(ntt).get();
}

}  // namespace hexl
}  // namespace intel
