#pragma once

#include "bigint.hpp"
#include <vector>

// --- テンプレート化: モントゴメリ曲線用のX, Z座標のみの構造体 ---
template <typename IntType>
struct PointXZ {
    IntType X, Z;
};

// =================================================================
// 【将来実装】384ビット専用のCIOSモントゴメリコンテキスト
// =================================================================
struct MontgomeryContext384 {
    MontgomeryContext384(const uint384_t& modulus) {}
    // 次回ここに xADD, xDBL などの実装を追加します
};

// =================================================================
// 従来の GMPベース モントゴメリ乗算コンテキスト
// =================================================================
struct MontgomeryContextGMP {
    cpp_int n;
    cpp_int r;     
    int r_bits;    
    cpp_int mask;  
    cpp_int n_prime; 
    
    // 事前確保バッファ（malloc/freeを完全に防ぐ）
    mutable mpz_t tmp_t, tmp_m, tmp_u, tmp_temp;
    mutable mpz_t t1, t2, t3, t4, t5, t6; 
    
    MontgomeryContextGMP(const cpp_int& modulus) {
        mpz_init(tmp_t); mpz_init(tmp_m); mpz_init(tmp_u); mpz_init(tmp_temp);
        mpz_init(t1); mpz_init(t2); mpz_init(t3); mpz_init(t4); mpz_init(t5); mpz_init(t6);
        
        n = modulus;
        int bits = 0;
        cpp_int temp = n;
        while (temp > 0) { bits++; temp >>= 1; }
        r_bits = ((bits + 63) / 64) * 64;
        r = cpp_int(1) << r_bits;
        mask = (cpp_int(1) << r_bits) - 1;
        
        cpp_int t = 0, newt = 1;
        cpp_int r_val = r, newr = n;
        while (newr != 0) {
            cpp_int quotient = r_val / newr;
            cpp_int temp_t = t - quotient * newt;
            t = newt; newt = temp_t;
            cpp_int temp_r = r_val - quotient * newr;
            r_val = newr; newr = temp_r;
        }
        if (t < 0) t += r;
        n_prime = (r - t) % r;
    }

    ~MontgomeryContextGMP() {
        mpz_clear(tmp_t); mpz_clear(tmp_m); mpz_clear(tmp_u); mpz_clear(tmp_temp);
        mpz_clear(t1); mpz_clear(t2); mpz_clear(t3); mpz_clear(t4); mpz_clear(t5); mpz_clear(t6);
    }

    cpp_int to_mont(const cpp_int& a) const { return (a << r_bits) % n; }
    cpp_int from_mont(const cpp_int& a) const { return mont_mul(a, 1); }

    cpp_int mont_mul(const cpp_int& a, const cpp_int& b) const {
        cpp_int res;
        mont_mul_raw(res.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
        return res;
    }

    void mont_mul_raw(mpz_t res, const mpz_t a, const mpz_t b) const {
        mpz_mul(tmp_t, a, b);
        mpz_mul(tmp_temp, tmp_t, n_prime.get_mpz_t());
        mpz_and(tmp_m, tmp_temp, mask.get_mpz_t());
        mpz_mul(tmp_temp, tmp_m, n.get_mpz_t());
        mpz_add(tmp_temp, tmp_t, tmp_temp);
        mpz_fdiv_q_2exp(tmp_u, tmp_temp, r_bits);
        if (mpz_cmp(tmp_u, n.get_mpz_t()) >= 0) {
            mpz_sub(tmp_u, tmp_u, n.get_mpz_t());
        }
        mpz_set(res, tmp_u);
    }

    void add_mod_raw(mpz_t res, const mpz_t a, const mpz_t b) const {
        mpz_add(res, a, b);
        if (mpz_cmp(res, n.get_mpz_t()) >= 0) mpz_sub(res, res, n.get_mpz_t());
    }

    void sub_mod_raw(mpz_t res, const mpz_t a, const mpz_t b) const {
        mpz_sub(res, a, b);
        if (mpz_sgn(res) < 0) mpz_add(res, res, n.get_mpz_t());
    }

    void xADD(PointXZ<cpp_int>& out, const PointXZ<cpp_int>& P, const PointXZ<cpp_int>& Q, const PointXZ<cpp_int>& P_minus_Q) const {
        sub_mod_raw(t1, P.X.get_mpz_t(), P.Z.get_mpz_t());
        add_mod_raw(t2, Q.X.get_mpz_t(), Q.Z.get_mpz_t());
        mont_mul_raw(t3, t1, t2);

        add_mod_raw(t1, P.X.get_mpz_t(), P.Z.get_mpz_t());
        sub_mod_raw(t2, Q.X.get_mpz_t(), Q.Z.get_mpz_t());
        mont_mul_raw(t4, t1, t2);

        add_mod_raw(t5, t3, t4);
        sub_mod_raw(t6, t3, t4);

        mont_mul_raw(t1, t5, t5);
        mont_mul_raw(out.X.get_mpz_t(), P_minus_Q.Z.get_mpz_t(), t1);

        mont_mul_raw(t2, t6, t6);
        mont_mul_raw(out.Z.get_mpz_t(), P_minus_Q.X.get_mpz_t(), t2);
    }

    void xDBL(PointXZ<cpp_int>& out, const PointXZ<cpp_int>& P, const cpp_int& A24_mont) const {
        add_mod_raw(t1, P.X.get_mpz_t(), P.Z.get_mpz_t());
        sub_mod_raw(t2, P.X.get_mpz_t(), P.Z.get_mpz_t());

        mont_mul_raw(t3, t1, t1);
        mont_mul_raw(t4, t2, t2);
        sub_mod_raw(t5, t3, t4);

        mont_mul_raw(out.X.get_mpz_t(), t3, t4);

        mont_mul_raw(t6, A24_mont.get_mpz_t(), t5);
        add_mod_raw(t6, t4, t6);

        mont_mul_raw(out.Z.get_mpz_t(), t5, t6);
    }
};

// --- テンプレート化: モントゴメリ・ラダーによる超高速スカラー倍 ---
// IntType と ContextType に依存し、GMP固有の処理を持たないクリーンな実装
template <typename IntType, typename ContextType>
PointXZ<IntType> mont_ladder(const PointXZ<IntType>& P, const IntType& k, const IntType& A24_mont, const ContextType& ctx, const IntType& N) {
    if (k == 0) return {ctx.to_mont(1), 0};
    
    PointXZ<IntType> R0 = P;
    PointXZ<IntType> R1;
    ctx.xDBL(R1, P, A24_mont);
    
    // ★ 追加した共通インターフェースを使用
    size_t bits = k.bit_length();
    
    for (int i = (int)bits - 2; i >= 0; --i) {
        bool bit = k.test_bit(i); // ★ 追加した共通インターフェースを使用
        if (!bit) {
            ctx.xADD(R1, R0, R1, P);    
            ctx.xDBL(R0, R0, A24_mont); 
        } else {
            ctx.xADD(R0, R0, R1, P);    
            ctx.xDBL(R1, R1, A24_mont); 
        }
    }
    return R0;
}