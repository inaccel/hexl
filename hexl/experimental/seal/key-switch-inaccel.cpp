#include "hexl/number-theory/number-theory.hpp"
#include "hexl/util/check.hpp"
#include <inaccel/coral>
#include <iostream>

#define BIT_MASK(BITS) ((1UL << BITS) - 1)
#define MAX_DEGREE 16384

typedef struct {
    uint64_t key1 : 52;
    uint64_t key2 : 52;
    uint64_t key3 : 52;
    uint64_t key4 : 52;
    uint64_t key5 : 48;
} __attribute__((packed)) DyadmultKeys1_t;

/// @brief
/// Struct DyadmultKeys2_t
/// @param[in] key1-6 stores the bits of compressed switch key data
typedef struct {
    uint64_t key1 : 4;
    uint64_t key2 : 52;
    uint64_t key3 : 52;
    uint64_t key4 : 52;
    uint64_t key5 : 52;
    uint64_t key6 : 44;
} __attribute__((packed)) DyadmultKeys2_t;

/// @brief
/// Struct DyadmultKeys3_t
/// @param[in] key1-5 stores the bits of compressed switch key data
typedef struct {
    uint64_t key1 : 8;
    uint64_t key2 : 52;
    uint64_t key3 : 52;
    uint64_t key4 : 52;
    uint64_t key5 : 52;
    uint64_t NOT_USED : 40;
} __attribute__((packed)) DyadmultKeys3_t;

typedef struct {
    uint64_t data[8][4];
} moduli_t;


namespace intel {
namespace hexl {
namespace internal {


void ComputeRootOfUnityPowers(uint64_t m_q, uint64_t m_degree,
                              uint64_t m_degree_bits, uint64_t m_w,
                              uint64_t* inv_root_of_unity_powers,
                              uint64_t* precon64_inv_root_of_unity_powers,
                              uint64_t* root_of_unity_powers,
                              uint64_t* precon64_root_of_unity_powers) {
    uint64_t inv_root_of_unity_powers_pre[MAX_DEGREE];

    // 64-bit preconditioning
    root_of_unity_powers[0] = 1;
    inv_root_of_unity_powers_pre[0] = 1;
    uint64_t idx = 0;
    uint64_t prev_idx = idx;

    for (size_t i = 1; i < m_degree; i++) {
        idx = ReverseBits(i, m_degree_bits);
        root_of_unity_powers[idx] =
            MultiplyMod(root_of_unity_powers[prev_idx], m_w, m_q);
        inv_root_of_unity_powers_pre[idx] =
            InverseMod(root_of_unity_powers[idx], m_q);

        prev_idx = idx;
    }

    precon64_root_of_unity_powers[0] = 0;
    for (size_t i = 1; i < m_degree; i++) {
        precon64_root_of_unity_powers[i] =
            MultiplyFactor(root_of_unity_powers[i], 64, m_q).BarrettFactor();
    }

    idx = 0;

    for (size_t m = (m_degree >> 1); m > 0; m >>= 1) {
        for (size_t i = 0; i < m; i++) {
            inv_root_of_unity_powers[idx] = inv_root_of_unity_powers_pre[m + i];
            idx++;
        }
    }

    inv_root_of_unity_powers[m_degree - 1] = 0;

    for (uint64_t i = 0; i < m_degree; i++) {
        precon64_inv_root_of_unity_powers[i] =
            MultiplyFactor(inv_root_of_unity_powers[i], 64, m_q)
                .BarrettFactor();
    }
}


void loadTwiddleFactors(uint64_t n, uint64_t key_modulus_size,
                          const uint64_t *moduli,
                          uint64_t *root_of_unity_powers_ptr) {
    for (uint64_t i = 0; i < key_modulus_size; i++) {
        ComputeRootOfUnityPowers(
            moduli[i], n, Log2(n),
            MinimalPrimitiveRoot(2 * n, moduli[i]),
            root_of_unity_powers_ptr + i * n * 4,
            root_of_unity_powers_ptr + i * n * 4 + n,
            root_of_unity_powers_ptr + i * n * 4 + n * 2,
            root_of_unity_powers_ptr + i * n * 4 + n * 3);
    }
}


uint64_t precomputeModulusK(uint64_t modulus) {
    uint64_t k = 0;
    for (uint64_t i = 64; i > 0; i--) {
        if ((1UL << i) >= modulus) {
            k = i;
        }
    }
    return k;
}


void buildModulusMeta(uint64_t key_modulus_size, const uint64_t *moduli,
                        const uint64_t* modswitch_factors,
                        moduli_t *modulus_meta) {
    for (uint64_t i = 0; i < key_modulus_size; i++) {
        (*modulus_meta).data[i][0] = moduli[i];
        (*modulus_meta).data[i][1] =
            MultiplyFactor(1, 64, moduli[i]).BarrettFactor();
        uint64_t modulus = moduli[i];
        uint64_t twice_modulus = 2 * modulus;
        uint64_t four_times_modulus = 4 * modulus;
        uint64_t arg2 = modswitch_factors[i];
        const int InputModFactor = 8;
        arg2 = ReduceMod<InputModFactor>(arg2, modulus, &twice_modulus,
                                         &four_times_modulus);
        (*modulus_meta).data[i][2] = arg2;
        uint64_t k = precomputeModulusK(moduli[i]);
        __extension__ __int128 a = 1;
        uint64_t r = (uint64_t) ((a << (2 * k)) / moduli[i]);
        (*modulus_meta).data[i][3] = (r << 8) | k;
    }
}


void buildInvnMeta(uint64_t n, uint64_t key_modulus_size,
                     const uint64_t *moduli,
                     const uint64_t* root_of_unity_powers_ptr, moduli_t *invn) {
    for (uint64_t i = 0; i < key_modulus_size; i++) {
        uint64_t inv_n = InverseMod(n, moduli[i]);
        uint64_t W_op =
            root_of_unity_powers_ptr[i * n * 4 + n - 1];
        uint64_t inv_nw = MultiplyMod(inv_n, W_op, moduli[i]);
        uint64_t y_barrett_n = DivideUInt128UInt64Lo(inv_n, 0, moduli[i]);
        uint64_t y_barrett_nw =
            DivideUInt128UInt64Lo(inv_nw, 0, moduli[i]);
        (*invn).data[i][0] = inv_n;
        unsigned long k = precomputeModulusK(moduli[i]);
        __extension__ __int128 a = 1;
        uint64_t r = (uint64_t) ((a << (2 * k)) / moduli[i]);
        (*invn).data[i][1] = (r << 8) | k;
        (*invn).data[i][2] = y_barrett_n;
        (*invn).data[i][3] = y_barrett_nw;
    }
}


void loadKeys(uint64_t n, uint64_t decomp_modulus_size,
               uint64_t key_modulus_size, const uint64_t** k_switch_keys,
               DyadmultKeys1_t* key_vector1, DyadmultKeys2_t* key_vector2,
               DyadmultKeys3_t* key_vector3) {
    size_t kv_idx = 0;
    for (uint64_t k = 0; k < decomp_modulus_size; k++) {
        for (uint64_t j = 0; j < n; j++, kv_idx++) {
            for (uint64_t i = 0; i < key_modulus_size; i++) {
                uint64_t key1 = k_switch_keys[k][i * n + j];
                uint64_t key2 =
                    k_switch_keys[k][(i + key_modulus_size) * n + j];
                if (i == 0) {
                    key_vector1[kv_idx].key1 = key1 & BIT_MASK(52);
                    key_vector1[kv_idx].key2 = key2 & BIT_MASK(52);
                } else if (i == 1) {
                    key_vector1[kv_idx].key3 = key1 & BIT_MASK(52);
                    key_vector1[kv_idx].key4 = key2 & BIT_MASK(52);
                } else if (i == 2) {
                    key_vector1[kv_idx].key5 = key1 & BIT_MASK(48);
                    key_vector2[kv_idx].key1 = (key1 >> 48) & BIT_MASK(4);
                    key_vector2[kv_idx].key2 = key2 & BIT_MASK(52);
                } else if (i == 3) {
                    key_vector2[kv_idx].key3 = key1 & BIT_MASK(52);
                    key_vector2[kv_idx].key4 = key2 & BIT_MASK(52);
                } else if (i == 4) {
                    key_vector2[kv_idx].key5 = key1 & BIT_MASK(52);
                    key_vector2[kv_idx].key6 = key2 & BIT_MASK(44);
                    key_vector3[kv_idx].key1 = (key2 >> 44) & BIT_MASK(8);
                } else if (i == 5) {
                    key_vector3[kv_idx].key2 = key1 & BIT_MASK(52);
                    key_vector3[kv_idx].key3 = key2 & BIT_MASK(52);
                } else if (i == 6) {
                    key_vector3[kv_idx].key4 = key1 & BIT_MASK(52);
                    key_vector3[kv_idx].key5 = key2 & BIT_MASK(52);
                    key_vector3[kv_idx].NOT_USED = 0;
                } else {
                        HEXL_CHECK(0, "Not supported keys");
                }
            }
        }
    }
}


void readOutput(uint64_t n, uint64_t decomp_modulus_size,
                 const uint64_t *moduli, const uint64_t *fpga_result,
                 uint64_t* result, uint64_t batch_size) {
    size_t output_size = n * decomp_modulus_size * 2;
    for (size_t off = 0; off < batch_size * output_size; off += output_size) {
        for (size_t i = 0; i < decomp_modulus_size; i++) {
            uint64_t modulus = moduli[i];
            for (size_t j = 0; j < n; j++) {
                size_t k = i * n + j;
                result[off + k] += fpga_result[off + 2 * k];
                if (result[off + k] >= modulus) result[off + k] -= modulus;

                size_t idx = k + n * decomp_modulus_size;
                result[off + idx] += fpga_result[off + 2 * k + 1];
                if (result[off + idx] >= modulus) result[off + idx] -= modulus;
            }
        }
    }
}

void KeySwitchInAccel(uint64_t* result, const uint64_t* t_target_iter_ptr,
               uint64_t n, uint64_t decomp_modulus_size,
               uint64_t key_modulus_size, uint64_t rns_modulus_size,
               uint64_t key_component_count, const uint64_t* moduli,
               const uint64_t** k_switch_keys,
               const uint64_t* modswitch_factors,
               const uint64_t* root_of_unity_powers_ptr,
               uint64_t batch_size) {
    HEXL_UNUSED(rns_modulus_size);
    HEXL_CHECK(key_component_count == 2, "key_component_count must be 2");

    int wmem = 1;
    uint64_t nil = 0;
    size_t root_of_unity_powers_ptr_size = n * key_modulus_size * 4;
    size_t key_size = decomp_modulus_size * n;
    size_t t_target_iter_ptr_size = batch_size * n * decomp_modulus_size;
    size_t fpga_result_size =
        batch_size * n * decomp_modulus_size * key_component_count;

    inaccel::allocator<uint64_t> alloc_uint64_t;
    inaccel::allocator<DyadmultKeys1_t> alloc_DyadmultKeys1_t;
    inaccel::allocator<DyadmultKeys2_t> alloc_DyadmultKeys2_t;
    inaccel::allocator<DyadmultKeys3_t> alloc_DyadmultKeys3_t;

    moduli_t modulus_meta = {}, invn = {};
    uint64_t *fpga_result = alloc_uint64_t.allocate(fpga_result_size);
    uint64_t *fpga_root_of_unity_powers_ptr =
        alloc_uint64_t.allocate(root_of_unity_powers_ptr_size);
    uint64_t *fpga_t_target_iter_ptr =
        alloc_uint64_t.allocate(t_target_iter_ptr_size);
    DyadmultKeys1_t *key1 = alloc_DyadmultKeys1_t.allocate(key_size);
    DyadmultKeys2_t *key2 = alloc_DyadmultKeys2_t.allocate(key_size);
    DyadmultKeys3_t *key3 = alloc_DyadmultKeys3_t.allocate(key_size);

    memcpy(fpga_t_target_iter_ptr, t_target_iter_ptr,
            t_target_iter_ptr_size * sizeof(uint64_t));

    if (!root_of_unity_powers_ptr) {
        loadTwiddleFactors(n, key_modulus_size, moduli,
            fpga_root_of_unity_powers_ptr);
    } else {
        memcpy(fpga_root_of_unity_powers_ptr, root_of_unity_powers_ptr,
            root_of_unity_powers_ptr_size * sizeof(uint64_t));
    }

    buildModulusMeta(key_modulus_size, moduli, modswitch_factors,
        &modulus_meta);
    buildInvnMeta(n, key_modulus_size, moduli, fpga_root_of_unity_powers_ptr,
        &invn);

    loadKeys(n, decomp_modulus_size, key_modulus_size, k_switch_keys, key1,
        key2, key3);

   inaccel::request keyswitch("hexl.experimental.seal.KeySwitch");
   keyswitch.arg((int) batch_size)
           .arg( batch_size)
           .arg((int) n)
           .arg(n)
           .arg(modulus_meta)
           .arg(invn)
           .arg(decomp_modulus_size)
           .arg_array<uint64_t>(fpga_root_of_unity_powers_ptr, fpga_root_of_unity_powers_ptr + root_of_unity_powers_ptr_size)
           .arg(root_of_unity_powers_ptr_size)
           .arg_array<DyadmultKeys1_t>(key1, key1 + key_size)
           .arg_array<DyadmultKeys2_t>(key2, key2 + key_size)
           .arg_array<DyadmultKeys3_t>(key3, key3 + key_size)
           .arg(key_size)
           .arg_array<uint64_t>(fpga_t_target_iter_ptr, fpga_t_target_iter_ptr + t_target_iter_ptr_size)
           .arg(t_target_iter_ptr_size)
           .arg_array<uint64_t>(fpga_result, fpga_result + fpga_result_size)
           .arg(fpga_result_size)
           .arg(nil)
           .arg(wmem);

    std::cerr << keyswitch << std::endl;

    inaccel::submit(keyswitch).get();

    readOutput(n, decomp_modulus_size, moduli, fpga_result, result, batch_size);

    alloc_DyadmultKeys3_t.deallocate(key3, key_size);
    alloc_DyadmultKeys2_t.deallocate(key2, key_size);
    alloc_DyadmultKeys1_t.deallocate(key1, key_size);
    alloc_uint64_t.deallocate(fpga_t_target_iter_ptr, t_target_iter_ptr_size);
    alloc_uint64_t.deallocate(fpga_root_of_unity_powers_ptr,
                              root_of_unity_powers_ptr_size);
    alloc_uint64_t.deallocate(fpga_result, fpga_result_size);
}

}  // namespace internal
}  // namespace hexl
}  // namespace intel
