#pragma once

#include "bigint.hpp"
#include <vector>
#include <algorithm>
#include <type_traits> // if constexpr 用

// =================================================================
// 【多項式演算エンジン】 (テンプレート版)
// =================================================================
template <typename IntType>
struct Poly {
    std::vector<IntType> coef; // coef[i] が X^i の係数

    Poly() {}
    Poly(const std::vector<IntType>& c) : coef(c) {}

    int degree() const {
        return coef.empty() ? -1 : (int)coef.size() - 1;
    }

    void clean() {
        while (!coef.empty() && coef.back() == IntType(0)) {
            coef.pop_back();
        }
    }

    Poly add(const Poly& other, const IntType& N) const {
        Poly res;
        int max_deg = std::max(degree(), other.degree());
        res.coef.resize(max_deg + 1, IntType(0));
        for (int i = 0; i <= max_deg; ++i) {
            IntType a = (i <= degree()) ? coef[i] : IntType(0);
            IntType b = (i <= other.degree()) ? other.coef[i] : IntType(0);
            res.coef[i] = mod(a + b, N);
        }
        res.clean();
        return res;
    }

    Poly sub(const Poly& other, const IntType& N) const {
        Poly res;
        int max_deg = std::max(degree(), other.degree());
        res.coef.resize(max_deg + 1, IntType(0));
        for (int i = 0; i <= max_deg; ++i) {
            IntType a = (i <= degree()) ? coef[i] : IntType(0);
            IntType b = (i <= other.degree()) ? other.coef[i] : IntType(0);
            res.coef[i] = mod(a - b, N);
        }
        res.clean();
        return res;
    }
};

// =================================================================
// 補助関数：分割統治法による O(N log N) パッキング
// 戻り値は巨大になるため、常に cpp_int を返す
// =================================================================
template <typename IntType>
cpp_int pack_poly(const std::vector<IntType>& coef, int left, int right, int K) {
    if (left == right) {
        if constexpr (std::is_same_v<IntType, cpp_int>) {
            return coef[left];
        } else {
            return cpp_int(coef[left].str()); // uint384_t等からのブリッジ
        }
    }
    int mid = left + (right - left) / 2;
    cpp_int P_left = pack_poly(coef, left, mid, K);
    cpp_int P_right = pack_poly(coef, mid + 1, right, K);
    
    int shift = (mid - left + 1) * K;
    return P_left + (P_right << shift);
}

// =================================================================
// 補助関数：分割統治法による O(N log N) アンパッキング
// =================================================================
template <typename IntType>
void unpack_poly(const cpp_int& packed, int left, int right, int K, const cpp_int& mask, const IntType& N, std::vector<IntType>& result) {
    if (left == right) {
        cpp_int rem = packed & mask;
        if constexpr (std::is_same_v<IntType, cpp_int>) {
            result[left] = mod(rem, N);
        } else {
            cpp_int N_cpp(N.str());
            rem = mod(rem, N_cpp);
            result[left] = IntType(rem.str()); // uint384_t等へのブリッジ
        }
        return;
    }
    int mid = left + (right - left) / 2;
    int shift = (mid - left + 1) * K;
    
    cpp_int lower_mask = (cpp_int(1) << shift) - 1;
    cpp_int P_left = packed & lower_mask;
    cpp_int P_right = packed >> shift;
    
    unpack_poly(P_left, left, mid, K, mask, N, result);
    unpack_poly(P_right, mid + 1, right, K, mask, N, result);
}

// =================================================================
// 超高速化版: クロネッカー代入を用いた高速多項式乗算 (mod N)
// =================================================================
template <typename IntType>
Poly<IntType> kronecker_multiply(const Poly<IntType>& A, const Poly<IntType>& B, const IntType& N) {
    if (A.degree() < 0 || B.degree() < 0) return Poly<IntType>();

    int n = A.degree();
    int m = B.degree();
    int min_len = std::min(n + 1, m + 1);

    int bits_N = N.bit_length(); // ★共通インターフェースを使用
    int log2_min = 0;
    int temp_min = min_len;
    while (temp_min > 0) { log2_min++; temp_min >>= 1; }

    int K = 2 * bits_N + log2_min + 2;

    cpp_int packed_A = pack_poly(A.coef, 0, n, K);
    cpp_int packed_B = pack_poly(B.coef, 0, m, K);

    // GMPの超巨大整数FFT乗算
    cpp_int packed_C = packed_A * packed_B;

    Poly<IntType> C;
    C.coef.resize(n + m + 1, IntType(0));
    cpp_int mask = (cpp_int(1) << K) - 1;
    
    unpack_poly(packed_C, 0, n + m, K, mask, N, C.coef);

    C.clean();
    return C;
}

// =================================================================
// 分割統治法による剰余木 (Subproduct Tree) の構築
// =================================================================
template <typename IntType>
Poly<IntType> build_poly_tree(const std::vector<IntType>& roots, int left, int right, const IntType& N) {
    if (left == right) {
        return Poly<IntType>({mod(-roots[left], N), IntType(1)});
    }
    
    int mid = left + (right - left) / 2;
    Poly<IntType> left_poly = build_poly_tree(roots, left, mid, N);
    Poly<IntType> right_poly = build_poly_tree(roots, mid + 1, right, N);
    
    return kronecker_multiply(left_poly, right_poly, N);
}

// =================================================================
// モニック多項式による多項式剰余算 (A mod B)
// ★超高速化: inner-loop内のメモリ確保を完全に防ぐ最適化
// =================================================================
template <typename IntType>
Poly<IntType> poly_mod_monic(Poly<IntType> A, const Poly<IntType>& B, const IntType& N) {
    int degB = B.degree();
    if (degB < 0) return A;
    
    if constexpr (std::is_same_v<IntType, cpp_int>) {
        // cpp_int の場合は malloc を回避するため mpz_t を直接回す
        mpz_t temp, factor_mpz, mod_N;
        mpz_init(temp);
        mpz_init(factor_mpz);
        mpz_init_set(mod_N, N.get_mpz_t());
        
        while (A.degree() >= degB) {
            int d = A.degree() - degB;
            mpz_set(factor_mpz, A.coef.back().get_mpz_t());
            
            for (int i = 0; i <= degB; ++i) {
                // A.coef[i+d] -= B.coef[i] * factor (mod N)
                mpz_mul(temp, B.coef[i].get_mpz_t(), factor_mpz);
                mpz_sub(A.coef[i + d].get_mpz_t(), A.coef[i + d].get_mpz_t(), temp);
                // mpz_fdiv_r は常に正の剰余を返す
                mpz_fdiv_r(A.coef[i + d].get_mpz_t(), A.coef[i + d].get_mpz_t(), mod_N);
            }
            A.clean();
        }
        mpz_clear(temp);
        mpz_clear(factor_mpz);
        mpz_clear(mod_N);
        return A;
    } else {
        // uint384_t などの固定長整数の場合はそのまま演算子を使う (元々mallocが発生しないため)
        while (A.degree() >= degB) {
            int d = A.degree() - degB;
            IntType factor = A.coef.back();
            for (int i = 0; i <= degB; ++i) {
                A.coef[i + d] -= B.coef[i] * factor;
                A.coef[i + d] = mod(A.coef[i + d], N);
            }
            A.clean();
        }
        return A;
    }
}

// =================================================================
// 評価対象の点群から、多項式除算用の木を構築
// =================================================================
template <typename IntType>
void build_eval_tree(int node, const std::vector<IntType>& points, int left, int right, const IntType& N, std::vector<Poly<IntType>>& tree) {
    if (left == right) {
        tree[node] = Poly<IntType>({mod(-points[left], N), IntType(1)});
        return;
    }
    int mid = left + (right - left) / 2;
    build_eval_tree(2 * node, points, left, mid, N, tree);
    build_eval_tree(2 * node + 1, points, mid + 1, right, N, tree);
    
    tree[node] = kronecker_multiply(tree[2 * node], tree[2 * node + 1], N);
}

// =================================================================
// 木構造の上から下へ多項式の割り算を繰り返し、最終的な代入値を抽出する
// =================================================================
template <typename IntType>
void evaluate_down(int node, Poly<IntType> F, int left, int right, const IntType& N, const std::vector<Poly<IntType>>& tree, std::vector<IntType>& results) {
    F = poly_mod_monic(F, tree[node], N);
    
    if (left == right) {
        results.push_back(F.degree() >= 0 ? F.coef[0] : IntType(0));
        return;
    }
    
    int mid = left + (right - left) / 2;
    evaluate_down(2 * node, F, left, mid, N, tree, results);
    evaluate_down(2 * node + 1, F, mid + 1, right, N, tree, results);
}