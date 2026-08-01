#pragma once

#include <string>
#include <cstdlib>
#include <gmp.h>
#include <numeric>

// =================================================================
// 【将来実装】384ビット専用の固定長整数
// =================================================================
struct uint384_t {
    uint64_t limbs[6] = {0};

    uint384_t() {}
    uint384_t(int v) { limbs[0] = v; }
    uint384_t(const std::string& s) { /* 次回実装 */ }
    
    std::string str() const { return "not implemented"; }
    
    bool operator==(int b) const { return limbs[0] == (uint64_t)b && limbs[1] == 0 && limbs[2] == 0 && limbs[3] == 0 && limbs[4] == 0 && limbs[5] == 0; }
    bool operator<(int b) const { return false; /* 次回実装 */ }
    bool operator>(int b) const { return false; /* 次回実装 */ }

    // テンプレート互換用メソッド
    size_t bit_length() const { return 0; /* 次回実装 */ }
    bool test_bit(int i) const { return false; /* 次回実装 */ }
};

// =================================================================
// 【mini-gmp C++ Wrapper】
// 既存の cpp_int の振る舞いを再現し、数式演算子をそのまま使えるようにする
// =================================================================
class cpp_int {
private:
    mpz_t val;
public:
    cpp_int() { mpz_init(val); }
    cpp_int(int v) { mpz_init_set_si(val, v); }
    cpp_int(long long v) { 
        std::string s = std::to_string(v);
        mpz_init_set_str(val, s.c_str(), 10); 
    }
    cpp_int(const char* s) { mpz_init_set_str(val, s, 10); }
    cpp_int(const std::string& s) { mpz_init_set_str(val, s.c_str(), 10); }
    
    cpp_int(const cpp_int& other) { mpz_init_set(val, other.val); }
    cpp_int(cpp_int&& other) noexcept { 
        val[0] = other.val[0]; 
        mpz_init(other.val); 
    }
    
    ~cpp_int() { mpz_clear(val); }

    cpp_int& operator=(const cpp_int& other) {
        if (this != &other) mpz_set(val, other.val);
        return *this;
    }
    cpp_int& operator=(int v) { mpz_set_si(val, v); return *this; }

    std::string str() const {
        char* s = mpz_get_str(NULL, 10, val);
        std::string res(s);
        free(s);
        return res;
    }

    explicit operator int() const { return (int)mpz_get_si(val); }
    explicit operator long long() const { 
        char* s = mpz_get_str(NULL, 10, val);
        long long res = std::stoll(s);
        free(s);
        return res;
    }

    const mpz_t& get_mpz_t() const { return val; }
    mpz_t& get_mpz_t() { return val; }

    // ★追加: テンプレートアルゴリズム（mont_ladder等）から呼ばれる共通インターフェース
    size_t bit_length() const { return mpz_sizeinbase(val, 2); }
    bool test_bit(int i) const { return mpz_tstbit(val, i); }

    // 算術演算子 (cpp_int同士)
    cpp_int operator+(const cpp_int& b) const { cpp_int r; mpz_add(r.val, val, b.val); return r; }
    cpp_int operator-(const cpp_int& b) const { cpp_int r; mpz_sub(r.val, val, b.val); return r; }
    cpp_int operator*(const cpp_int& b) const { cpp_int r; mpz_mul(r.val, val, b.val); return r; }
    cpp_int operator/(const cpp_int& b) const { cpp_int r; mpz_tdiv_q(r.val, val, b.val); return r; }
    cpp_int operator%(const cpp_int& b) const { cpp_int r; mpz_tdiv_r(r.val, val, b.val); return r; }
    cpp_int operator-() const { cpp_int r; mpz_neg(r.val, val); return r; } 
    
    // 算術演算子 (intとの計算)
    cpp_int operator+(int b) const { return *this + cpp_int(b); }
    cpp_int operator-(int b) const { return *this - cpp_int(b); }
    cpp_int operator*(int b) const { return *this * cpp_int(b); }
    cpp_int operator/(int b) const { return *this / cpp_int(b); }
    cpp_int operator%(int b) const { return *this % cpp_int(b); }

    friend cpp_int operator+(int a, const cpp_int& b) { return cpp_int(a) + b; }
    friend cpp_int operator-(int a, const cpp_int& b) { return cpp_int(a) - b; }
    friend cpp_int operator*(int a, const cpp_int& b) { return cpp_int(a) * b; }
    friend cpp_int operator/(int a, const cpp_int& b) { return cpp_int(a) / b; }
    friend cpp_int operator%(int a, const cpp_int& b) { return cpp_int(a) % b; }

    cpp_int operator<<(int b) const { cpp_int r; mpz_mul_2exp(r.val, val, b); return r; }
    cpp_int operator>>(int b) const { cpp_int r; mpz_fdiv_q_2exp(r.val, val, b); return r; }
    cpp_int operator&(const cpp_int& b) const { cpp_int r; mpz_and(r.val, val, b.val); return r; }
    cpp_int operator|(const cpp_int& b) const { cpp_int r; mpz_ior(r.val, val, b.val); return r; }

    cpp_int& operator+=(const cpp_int& b) { mpz_add(val, val, b.val); return *this; }
    cpp_int& operator-=(const cpp_int& b) { mpz_sub(val, val, b.val); return *this; }
    cpp_int& operator*=(const cpp_int& b) { mpz_mul(val, val, b.val); return *this; }
    cpp_int& operator/=(const cpp_int& b) { mpz_tdiv_q(val, val, b.val); return *this; }
    cpp_int& operator%=(const cpp_int& b) { mpz_tdiv_r(val, val, b.val); return *this; }
    cpp_int& operator>>=(int b) { mpz_fdiv_q_2exp(val, val, b); return *this; }
    cpp_int& operator<<=(int b) { mpz_mul_2exp(val, val, b); return *this; }

    cpp_int& operator+=(int b) { return *this += cpp_int(b); }
    cpp_int& operator-=(int b) { return *this -= cpp_int(b); }
    cpp_int& operator*=(int b) { return *this *= cpp_int(b); }
    cpp_int& operator/=(int b) { return *this /= cpp_int(b); }
    cpp_int& operator%=(int b) { return *this %= cpp_int(b); }

    bool operator==(const cpp_int& b) const { return mpz_cmp(val, b.val) == 0; }
    bool operator!=(const cpp_int& b) const { return mpz_cmp(val, b.val) != 0; }
    bool operator<(const cpp_int& b) const { return mpz_cmp(val, b.val) < 0; }
    bool operator>(const cpp_int& b) const { return mpz_cmp(val, b.val) > 0; }
    bool operator<=(const cpp_int& b) const { return mpz_cmp(val, b.val) <= 0; }
    bool operator>=(const cpp_int& b) const { return mpz_cmp(val, b.val) >= 0; }

    bool operator==(int b) const { return mpz_cmp_si(val, b) == 0; }
    bool operator!=(int b) const { return mpz_cmp_si(val, b) != 0; }
    bool operator<(int b) const { return mpz_cmp_si(val, b) < 0; }
    bool operator>(int b) const { return mpz_cmp_si(val, b) > 0; }
};

// 安全なモジュロ演算（負の数に対応）
inline cpp_int mod(cpp_int a, cpp_int n) {
    a %= n;
    if (a < 0) a += n;
    return a;
}

// Boostの gcd をエミュレートするグローバル関数
inline cpp_int gcd(const cpp_int& a, const cpp_int& b) {
    cpp_int r;
    mpz_gcd(r.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
    return r;
}

// 将来の uint384_t 用ダミー
inline uint384_t mod(uint384_t a, uint384_t n) { return uint384_t(); }
inline uint384_t gcd(const uint384_t& a, const uint384_t& b) { return uint384_t(); }