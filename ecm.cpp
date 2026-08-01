#include <emscripten/bind.h>
#include <emscripten/val.h>
#include <string>
#include <cstdlib>
#include <vector>
#include <map>
#include <gmp.h>
#include <numeric> // std::gcd用に追加

using namespace emscripten;

// =================================================================
// 【mini-gmp C++ Wrapper】
// 既存の cpp_int の振る舞いを完全に再現し、数式演算子をそのまま使えるようにする
// =================================================================
class cpp_int {
private:
    mpz_t val;
public:
    // コンストラクタ
    cpp_int() { mpz_init(val); }
    cpp_int(int v) { mpz_init_set_si(val, v); }
    cpp_int(long long v) { 
        std::string s = std::to_string(v);
        mpz_init_set_str(val, s.c_str(), 10); 
    }
    cpp_int(const char* s) { mpz_init_set_str(val, s, 10); }
    cpp_int(const std::string& s) { mpz_init_set_str(val, s.c_str(), 10); }
    
    // コピーとムーブ
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

    // ★追加: 型変換演算子 ( (int)(k % 4) などのキャスト用 )
    explicit operator int() const { return (int)mpz_get_si(val); }
    explicit operator long long() const { 
        char* s = mpz_get_str(NULL, 10, val);
        long long res = std::stoll(s);
        free(s);
        return res;
    }

    // 生の mpz_t へのアクセス（内部用）
    const mpz_t& get_mpz_t() const { return val; }
    mpz_t& get_mpz_t() { return val; }

    // 算術演算子 (cpp_int同士)
    cpp_int operator+(const cpp_int& b) const { cpp_int r; mpz_add(r.val, val, b.val); return r; }
    cpp_int operator-(const cpp_int& b) const { cpp_int r; mpz_sub(r.val, val, b.val); return r; }
    cpp_int operator*(const cpp_int& b) const { cpp_int r; mpz_mul(r.val, val, b.val); return r; }
    cpp_int operator/(const cpp_int& b) const { cpp_int r; mpz_tdiv_q(r.val, val, b.val); return r; }
    cpp_int operator%(const cpp_int& b) const { cpp_int r; mpz_tdiv_r(r.val, val, b.val); return r; }
    // 算術演算子 (cpp_int同士)
    cpp_int operator-() const { cpp_int r; mpz_neg(r.val, val); return r; } // ★この1行を追加
    
    // 算術演算子 (intとの計算)
    cpp_int operator+(int b) const { return *this + cpp_int(b); }
    cpp_int operator-(int b) const { return *this - cpp_int(b); }
    cpp_int operator*(int b) const { return *this * cpp_int(b); }
    cpp_int operator/(int b) const { return *this / cpp_int(b); }
    cpp_int operator%(int b) const { return *this % cpp_int(b); }

    // friend関数 (intが左辺に来る場合)
    friend cpp_int operator+(int a, const cpp_int& b) { return cpp_int(a) + b; }
    friend cpp_int operator-(int a, const cpp_int& b) { return cpp_int(a) - b; }
    friend cpp_int operator*(int a, const cpp_int& b) { return cpp_int(a) * b; }
    friend cpp_int operator/(int a, const cpp_int& b) { return cpp_int(a) / b; }
    friend cpp_int operator%(int a, const cpp_int& b) { return cpp_int(a) % b; }

    // ビットシフト・論理演算
    cpp_int operator<<(int b) const { cpp_int r; mpz_mul_2exp(r.val, val, b); return r; }
    cpp_int operator>>(int b) const { cpp_int r; mpz_fdiv_q_2exp(r.val, val, b); return r; }
    cpp_int operator&(const cpp_int& b) const { cpp_int r; mpz_and(r.val, val, b.val); return r; }
    cpp_int operator|(const cpp_int& b) const { cpp_int r; mpz_ior(r.val, val, b.val); return r; }

    // 代入演算子 (cpp_intとの計算)
    cpp_int& operator+=(const cpp_int& b) { mpz_add(val, val, b.val); return *this; }
    cpp_int& operator-=(const cpp_int& b) { mpz_sub(val, val, b.val); return *this; }
    cpp_int& operator*=(const cpp_int& b) { mpz_mul(val, val, b.val); return *this; }
    // ★追加: / と % の代入
    cpp_int& operator/=(const cpp_int& b) { mpz_tdiv_q(val, val, b.val); return *this; }
    cpp_int& operator%=(const cpp_int& b) { mpz_tdiv_r(val, val, b.val); return *this; }
    cpp_int& operator>>=(int b) { mpz_fdiv_q_2exp(val, val, b); return *this; }
    cpp_int& operator<<=(int b) { mpz_mul_2exp(val, val, b); return *this; }

    // ★追加: 代入演算子 (intとの計算。 k /= 2 等に必須)
    cpp_int& operator+=(int b) { return *this += cpp_int(b); }
    cpp_int& operator-=(int b) { return *this -= cpp_int(b); }
    cpp_int& operator*=(int b) { return *this *= cpp_int(b); }
    cpp_int& operator/=(int b) { return *this /= cpp_int(b); }
    cpp_int& operator%=(int b) { return *this %= cpp_int(b); }

    // 比較演算子
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
// =================================================================

// =================================================================
// 【将来実装】384ビット専用の固定長整数とCIOSモントゴメリコンテキスト
// =================================================================
struct uint384_t {
    // 次回ここに実装を追加します
    uint384_t() {}
    uint384_t(int v) {}
    uint384_t(const std::string& s) {}
    std::string str() const { return "not implemented"; }
    bool operator==(int b) const { return false; }
    // ...
};

struct MontgomeryContext384 {
    MontgomeryContext384(const uint384_t& modulus) {}
    // 次回ここに実装を追加します
};

// Boostの gcd をエミュレートするグローバル関数
inline cpp_int gcd(const cpp_int& a, const cpp_int& b) {
    cpp_int r;
    mpz_gcd(r.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
    return r;
}
// =======================================================


// --- モントゴメリ曲線用のX, Z座標のみの構造体 ---
struct PointXZ {
    cpp_int X, Z;
};

// --- モントゴメリ乗算用クラス・ヘルパー (ゼロ・アロケーション最適化版) ---
struct MontgomeryContext {
    cpp_int n;
    cpp_int r;     
    int r_bits;    
    cpp_int mask;  
    cpp_int n_prime; 
    
    // 事前確保バッファ（malloc/freeを完全に防ぐ）
    mutable mpz_t tmp_t, tmp_m, tmp_u, tmp_temp;
    mutable mpz_t t1, t2, t3, t4, t5, t6; // ★追加: xADD/xDBL 専用の使い回しバッファ
    
    MontgomeryContext(const cpp_int& modulus) {
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

    ~MontgomeryContext() {
        mpz_clear(tmp_t); mpz_clear(tmp_m); mpz_clear(tmp_u); mpz_clear(tmp_temp);
        mpz_clear(t1); mpz_clear(t2); mpz_clear(t3); mpz_clear(t4); mpz_clear(t5); mpz_clear(t6);
    }

    cpp_int to_mont(const cpp_int& a) const { return (a << r_bits) % n; }
    cpp_int from_mont(const cpp_int& a) const { return mont_mul(a, 1); }

    // 外部（C++コード）から呼ばれる用のラッパー
    cpp_int mont_mul(const cpp_int& a, const cpp_int& b) const {
        cpp_int res;
        mont_mul_raw(res.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
        return res;
    }

    // =================================================================
    // ★ 内部計算用の超高速 raw 関数群 (割り算レス & インプレース)
    // =================================================================
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

    // =================================================================
    // ★ 究極の xADD / xDBL (結果を参照渡しして alloc を完全に回避)
    // =================================================================
    void xADD(PointXZ& out, const PointXZ& P, const PointXZ& Q, const PointXZ& P_minus_Q) const {
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

    void xDBL(PointXZ& out, const PointXZ& P, const cpp_int& A24_mont) const {
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

// --- 追加: 計算状態を記録する構造体（射影座標版） ---
struct ECMState {
    cpp_int found_factor;
    long long current_prime;
    cpp_int current_multiplier;
    cpp_int final_X, final_Y, final_Z; // 射影座標
    int phase; // ★追加: どちらのフェーズで見つかったか記録
};

ECMState global_ecm_state;

// 安全なモジュロ演算（負の数に対応）
cpp_int mod(cpp_int a, cpp_int n) {
    a %= n;
    if (a < 0) a += n;
    return a;
}

// --- 追加: グローバルな素数キャッシュ ---
std::vector<int> cached_primes;
int max_sieved_B1 = 1;

// 定数ではなく、動的に書き換えるグローバル変数に変更
int current_phase2_D = 0;
std::vector<int> current_phase2_coprimes;

// B2 のサイズに応じて最適な D を選び、Baby-step のリストを自動生成する
void setup_phase2_params(int B2) {
    int target_D = 210;
    if (B2 > 50000000) target_D = 510510;       // 2 * 3 * 5 * 7 * 11 * 13 * 17
    else if (B2 > 2000000) target_D = 30030;    // 2 * 3 * 5 * 7 * 11 * 13
    else if (B2 > 50000) target_D = 2310;       // 2 * 3 * 5 * 7 * 11
    
    // すでに計算済みならスキップ
    if (current_phase2_D == target_D) return;
    
    current_phase2_D = target_D;
    current_phase2_coprimes.clear();
    
    // D/2 以下の、D と互いに素な整数をすべてリストアップ（これがそのまま素数＆合成数の探索候補になる）
    for (int i = 1; i <= target_D / 2; i += 2) {
        if (std::gcd(i, target_D) == 1) {
            current_phase2_coprimes.push_back(i);
        }
    }
}

// 必要なB1までエラトステネスの篩で素数テーブルを構築・拡張する
void ensure_primes_up_to(int B1) {
    if (B1 <= max_sieved_B1) return; // 既に計算済みなら何もしない
    
    std::vector<bool> is_p(B1 + 1, true);
    is_p[0] = is_p[1] = false;
    for (int p = 2; p * p <= B1; p++) {
        if (is_p[p]) {
            for (int i = p * p; i <= B1; i += p) {
                is_p[i] = false;
            }
        }
    }
    
    cached_primes.clear();
    for (int p = 2; p <= B1; p++) {
        if (is_p[p]) cached_primes.push_back(p);
    }
    max_sieved_B1 = B1;
}

// 素数判定用の繰り返し二乗法 (変更なし)
cpp_int powm(cpp_int base, cpp_int exp, cpp_int mod) {
    cpp_int res = 1;
    base %= mod;
    while (exp > 0) {
        if (exp % 2 == 1) res = (res * base) % mod;
        base = (base * base) % mod;
        exp /= 2;
    }
    return res;
}

// 巨大な多倍長整数用の乱数生成関数 (変更なし)
cpp_int generate_random_base(cpp_int limit) {
    if (limit <= 0) return 0;
    cpp_int res = 0;
    cpp_int temp = limit;
    int chunks = 0;
    while(temp > 0) { temp >>= 30; chunks++; }
    for(int i = 0; i < chunks; i++) res = (res << 30) | (std::rand() & 0x3FFFFFFF);
    return res % limit;
}

// 本格版ミラー・ラビン素数判定 (変更なし)
bool is_prime_cpp(const std::string& n_str) {
    cpp_int n(n_str);
    if (n < 2) return false;
    int small_primes[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97};
    for (int p : small_primes) {
        if (n == p) return true;
        if (n % p == 0) return false;
    }
    cpp_int d = n - 1;
    int s = 0;
    while (d % 2 == 0) { d /= 2; s++; }
    
    int k = 20; 
    for (int i = 0; i < k; i++) {
        cpp_int limit = n - 3;
        cpp_int a = 2 + generate_random_base(limit);
        cpp_int x = powm(a, d, n);
        if (x == 1 || x == n - 1) continue;
        bool composite = true;
        for (int r = 1; r < s; r++) {
            x = (x * x) % n;
            if (x == n - 1) { composite = false; break; }
        }
        if (composite) return false;
    }
    return true;
}

// 拡張ユークリッドの互除法による逆元計算
cpp_int modInverse(cpp_int a, cpp_int m) {
    cpp_int m0 = m, y = 0, x = 1;
    if (m == 1) return 0;
    cpp_int q, t;
    while (a > 1) {
        if (m0 == 0) break;
        q = a / m0;
        t = m0;
        m0 = a % m0;
        a = t;
        t = y;
        y = x - q * y;
        x = t;
    }
    if (a > 1) { 
        // 【素因数発見】グローバル状態に保存
        global_ecm_state.found_factor = a; 
        return 1; 
    }
    if (x < 0) x += m;
    return x;
}

// --- モントゴメリのバッチ逆元 (一括逆元計算) ---
// 複数のZ座標の逆元を、1回の modInverse でまとめて計算する超高速化手法
bool batch_inverse(std::vector<cpp_int>& vec, const cpp_int& N) {
    int n = vec.size();
    if (n == 0) return true;
    std::vector<cpp_int> prod(n);
    cpp_int acc = 1;
    for (int i = 0; i < n; ++i) {
        prod[i] = acc;
        acc = (acc * vec[i]) % N;
    }
    cpp_int acc_inv = modInverse(acc, N);
    
    // 逆元計算中に素因数を引き当てた場合
    if (global_ecm_state.found_factor > 0) return false;
    
    for (int i = n - 1; i >= 0; --i) {
        cpp_int current = vec[i];
        vec[i] = (acc_inv * prod[i]) % N;
        if (vec[i] < 0) vec[i] += N;
        acc_inv = (acc_inv * current) % N;
    }
    return true;
}

// --- (追加) ポラード・ロー法 ---
// 小さな素因数を高速に探索する。見つかったらその文字列表現を、見つからなければ空文字を返す。
std::string pollard_rho(const std::string& n_str, int max_steps) {
    cpp_int n(n_str);
    if (n % 2 == 0) return "2";
    
    // x = (x^2 + c) mod n という関数を使う
    cpp_int x = 2;
    cpp_int y = 2;
    cpp_int d = 1;
    cpp_int c = cpp_int(std::rand()) % (n - 1) + 1;

    int step = 0;
    while (d == 1 && step < max_steps) {
        // 亀(x)は1歩進む
        x = (x * x + c) % n;
        // 兎(y)は2歩進む
        y = (y * y + c) % n;
        y = (y * y + c) % n;
        
        cpp_int diff = x > y ? x - y : y - x;
        d = gcd(diff, n);
        
        step++;
    }

    if (d > 1 && d < n) {
        return d.str(); // 素因数発見！
    }
    return ""; // 見つからなかった（諦めてECMに任せる）
}


// モントゴメリ・ラダーによる超高速スカラー倍 (Zero-Allocation 版)
PointXZ mont_ladder(const PointXZ& P, cpp_int k, const cpp_int& A24_mont, const MontgomeryContext& ctx, const cpp_int& N) {
    if (k == 0) return {ctx.to_mont(1), 0};
    
    PointXZ R0 = P;
    PointXZ R1;
    ctx.xDBL(R1, P, A24_mont);
    
    // ★ C++の演算子を使わず、GMPAPIで直接ビット数とビット値を抽出 (malloc回避)
    size_t bits = mpz_sizeinbase(k.get_mpz_t(), 2);
    
    for (int i = (int)bits - 2; i >= 0; --i) {
        bool bit = mpz_tstbit(k.get_mpz_t(), i);
        if (!bit) {
            ctx.xADD(R1, R0, R1, P);    // R1 を上書き
            ctx.xDBL(R0, R0, A24_mont); // R0 を上書き
        } else {
            ctx.xADD(R0, R0, R1, P);    // R0 を上書き
            ctx.xDBL(R1, R1, A24_mont); // R1 を上書き
        }
    }
    return R0;
}

// =================================================================
// 【多項式演算エンジン (FFTフェーズ2用)】
// =================================================================
struct Poly {
    std::vector<cpp_int> coef; // coef[i] が X^i の係数

    Poly() {}
    Poly(const std::vector<cpp_int>& c) : coef(c) {}

    int degree() const {
        return coef.empty() ? -1 : (int)coef.size() - 1;
    }

    // 最高次の0を削除して正規化する
    void clean() {
        while (!coef.empty() && coef.back() == 0) {
            coef.pop_back();
        }
    }

    // 多項式の加算 (mod N)
    Poly add(const Poly& other, const cpp_int& N) const {
        Poly res;
        int max_deg = std::max(degree(), other.degree());
        res.coef.resize(max_deg + 1, 0);
        for (int i = 0; i <= max_deg; ++i) {
            cpp_int a = (i <= degree()) ? coef[i] : 0;
            cpp_int b = (i <= other.degree()) ? other.coef[i] : 0;
            res.coef[i] = mod(a + b, N);
        }
        res.clean();
        return res;
    }

    // 多項式の減算 (mod N)
    Poly sub(const Poly& other, const cpp_int& N) const {
        Poly res;
        int max_deg = std::max(degree(), other.degree());
        res.coef.resize(max_deg + 1, 0);
        for (int i = 0; i <= max_deg; ++i) {
            cpp_int a = (i <= degree()) ? coef[i] : 0;
            cpp_int b = (i <= other.degree()) ? other.coef[i] : 0;
            res.coef[i] = mod(a - b, N);
        }
        res.clean();
        return res;
    }
};

// =================================================================
// 補助関数：分割統治法による O(N log N) パッキング
// =================================================================
cpp_int pack_poly(const std::vector<cpp_int>& coef, int left, int right, int K) {
    if (left == right) {
        return coef[left];
    }
    int mid = left + (right - left) / 2;
    cpp_int P_left = pack_poly(coef, left, mid, K);
    cpp_int P_right = pack_poly(coef, mid + 1, right, K);
    
    // 右半分をシフトして左半分と足す
    int shift = (mid - left + 1) * K;
    return P_left + (P_right << shift);
}

// =================================================================
// 補助関数：分割統治法による O(N log N) アンパッキング
// =================================================================
void unpack_poly(const cpp_int& packed, int left, int right, int K, const cpp_int& mask, const cpp_int& N, std::vector<cpp_int>& result) {
    if (left == right) {
        result[left] = mod(packed & mask, N);
        return;
    }
    int mid = left + (right - left) / 2;
    int shift = (mid - left + 1) * K;
    
    // 下位ビットを抽出するためのマスク
    cpp_int lower_mask = (cpp_int(1) << shift) - 1;
    
    // 巨大な数を「左半分」と「右半分」に分割
    cpp_int P_left = packed & lower_mask;
    cpp_int P_right = packed >> shift;
    
    unpack_poly(P_left, left, mid, K, mask, N, result);
    unpack_poly(P_right, mid + 1, right, K, mask, N, result);
}

// =================================================================
// 超高速化版: クロネッカー代入を用いた高速多項式乗算 (mod N)
// =================================================================
Poly kronecker_multiply(const Poly& A, const Poly& B, const cpp_int& N) {
    if (A.degree() < 0 || B.degree() < 0) return Poly();

    int n = A.degree();
    int m = B.degree();
    int min_len = std::min(n + 1, m + 1);

    int bits_N = 0;
    cpp_int temp_N = N;
    while (temp_N > 0) { bits_N++; temp_N >>= 1; }

    int log2_min = 0;
    int temp_min = min_len;
    while (temp_min > 0) { log2_min++; temp_min >>= 1; }

    int K = 2 * bits_N + log2_min + 2;

    // ★ 改善点 1: 分割統治法で高速パッキング
    cpp_int packed_A = pack_poly(A.coef, 0, n, K);
    cpp_int packed_B = pack_poly(B.coef, 0, m, K);

    // GMPによる本物のFFT巨大整数乗算 (ここが一番重くなるのが正常)
    cpp_int packed_C = packed_A * packed_B;

    // ★ 改善点 2: 分割統治法で高速アンパッキング
    Poly C;
    C.coef.resize(n + m + 1, 0);
    cpp_int mask = (cpp_int(1) << K) - 1;
    
    unpack_poly(packed_C, 0, n + m, K, mask, N, C.coef);

    C.clean();
    return C;
}

// 分割統治法による剰余木 (Subproduct Tree) の構築
// 根のリスト [x_1, x_2, ... x_n] から F(X) = Π (X - x_i) を O(n log^2 n) で計算
Poly build_poly_tree(const std::vector<cpp_int>& roots, int left, int right, const cpp_int& N) {
    if (left == right) {
        // (X - root) を作成: coef[0] = -root, coef[1] = 1
        return Poly({mod(-roots[left], N), 1});
    }
    
    int mid = left + (right - left) / 2;
    Poly left_poly = build_poly_tree(roots, left, mid, N);
    Poly right_poly = build_poly_tree(roots, mid + 1, right, N);
    
    return kronecker_multiply(left_poly, right_poly, N);
}

// =================================================================
// 【高速多点評価 (Fast Multipoint Evaluation)】
// =================================================================

// モニック多項式（最高次の係数が1）による多項式剰余算 (A mod B)
Poly poly_mod_monic(Poly A, const Poly& B, const cpp_int& N) {
    int degB = B.degree();
    if (degB < 0) return A;
    
    // 多項式の筆算（テンポラリ変数を排除し、直接配列を書き換える）
    while (A.degree() >= degB) {
        int d = A.degree() - degB;
        cpp_int factor = A.coef.back();
        for (int i = 0; i <= degB; ++i) {
            // mod() 関数の呼び出しを避け、直接インプレースで減算と剰余を実行
            A.coef[i + d] -= B.coef[i] * factor;
            A.coef[i + d] %= N;
            if (A.coef[i + d] < 0) A.coef[i + d] += N;
        }
        A.clean();
    }
    return A;
}

// 評価対象の点群から、多項式除算用の木（Subproduct Tree）を構築
void build_eval_tree(int node, const std::vector<cpp_int>& points, int left, int right, const cpp_int& N, std::vector<Poly>& tree) {
    if (left == right) {
        tree[node] = Poly({mod(-points[left], N), 1});
        return;
    }
    int mid = left + (right - left) / 2;
    build_eval_tree(2 * node, points, left, mid, N, tree);
    build_eval_tree(2 * node + 1, points, mid + 1, right, N, tree);
    
    tree[node] = kronecker_multiply(tree[2 * node], tree[2 * node + 1], N);
}

// 木構造の上から下へ多項式の割り算を繰り返し、最終的な代入値を抽出する
void evaluate_down(int node, Poly F, int left, int right, const cpp_int& N, const std::vector<Poly>& tree, std::vector<cpp_int>& results) {
    // まず自分自身のノードの多項式で割った余りを取る
    F = poly_mod_monic(F, tree[node], N);
    
    // 葉に到達したら、Fは定数（評価値）になっている
    if (left == right) {
        results.push_back(F.degree() >= 0 ? F.coef[0] : 0);
        return;
    }
    
    int mid = left + (right - left) / 2;
    evaluate_down(2 * node, F, left, mid, N, tree, results);
    evaluate_down(2 * node + 1, F, mid + 1, right, N, tree, results);
}
// =================================================================

// ECMを1サイクル実行する関数
// --- バッチ処理対応の ECM 実行関数 ---
// 指定された数 (num_curves) だけ、C++内部で高速にガチャを回し続けます
val run_ecm_batch(std::string n_str, int B1, int seed, int num_curves) {
    std::srand(seed);
    cpp_int N(n_str);
    val result = val::object();

    if (N % 2 == 0 || N % 3 == 0) {
        cpp_int factor = (N % 2 == 0) ? 2 : 3;
        result.set("status", "found");
        result.set("factor", factor.str());
        result.set("phase", 1);
        result.set("prime", factor.str());
        result.set("k", "1");
        result.set("Z", factor.str());
        result.set("A", "0"); result.set("B", "0");
        result.set("P_x", "0"); result.set("P_y", "0");
        result.set("curves_tried", 1);
        return result;
    }

    ensure_primes_up_to(B1);
    MontgomeryContext ctx(N); 
    
    for (int iter = 0; iter < num_curves; ++iter) {
        global_ecm_state.found_factor = 0;
        global_ecm_state.phase = 1;
        
        // ==========================================
        // 【最適化】 Suyama's Parametrization
        // ==========================================
        cpp_int sigma = std::rand() % 1000000 + 6; // σ >= 6 を選択
        cpp_int u = mod(sigma * sigma - 5, N);
        cpp_int v = mod(4 * sigma, N);
        
        cpp_int u3 = mod(u * u * u, N);
        cpp_int v3 = mod(v * v * v, N);
        
        // 初期点 P = (X0, Z0) = (u^3, v^3)
        PointXZ P_XZ = {ctx.to_mont(u3), ctx.to_mont(v3)};
        
        // A24 = (A+2)/4 = (v-u)^3 * (3u+v) / (16 * u^3 * v)
        cpp_int v_minus_u = mod(v - u, N);
        cpp_int term1 = mod(v_minus_u * v_minus_u * v_minus_u, N);
        cpp_int term2 = mod(3 * u + v, N);
        cpp_int num = mod(term1 * term2, N);
        
        cpp_int den = mod(16 * u3 * v, N);
        cpp_int den_inv = modInverse(den, N);
        
        // ※ もしここで modInverse が N と共通因数を見つけた場合、
        // global_ecm_state.found_factor に記録されるためそのまま処理を進めてOK
        if (global_ecm_state.found_factor > 0) break;
        
        cpp_int A24 = mod(num * den_inv, N);
        cpp_int A24_mont = ctx.to_mont(A24);
        
        // ==========================================
        // 【新生フェーズ 1】 モントゴメリ・ラダー
        // ==========================================
        for (int p : cached_primes) {
            if (p > B1) break;
            cpp_int q = p, q_max = p;
            while (q_max * p <= B1) q_max *= p;
            
            P_XZ = mont_ladder(P_XZ, q_max, A24_mont, ctx, N);
            
            cpp_int Z_normal = ctx.from_mont(P_XZ.Z);
            cpp_int d = gcd(Z_normal, N);
            if (d > 1 && d < N) {
                global_ecm_state.found_factor = d;
                global_ecm_state.current_prime = p;
                global_ecm_state.current_multiplier = q_max;
                global_ecm_state.final_Z = Z_normal;
                break;
            } else if (d == N) {
                global_ecm_state.found_factor = N;
                break;
            }
        }

        // ==========================================
        // 【新生フェーズ 2】 多項式多点評価 (FFTベース)
        // ==========================================
        if (global_ecm_state.found_factor == 0) {
            // ★ FFTの恩恵を最大化するため、B2 を B1 の 500倍 に設定（GMP-ECM標準クラス）
            int B2 = B1 * 500; 
            setup_phase2_params(B2);
            
            // 1. Baby-step のアフィンX座標をバッチ逆元で集める
            std::vector<PointXZ> baby_steps;
            std::vector<cpp_int> baby_z;
            for (int j : current_phase2_coprimes) { // ★ 動的リストを参照
                PointXZ step = mont_ladder(P_XZ, j, A24_mont, ctx, N);
                baby_steps.push_back(step);
                baby_z.push_back(step.Z);
            }
            
            // ★ ここが消えていました：Z座標を一括で逆元に変換
            if (!batch_inverse(baby_z, N)) continue;
            
            // ★ ここが消えていました：X * Z^-1 でアフィンX座標 (baby_x) を計算
            std::vector<cpp_int> baby_x;
            for (size_t i = 0; i < baby_steps.size(); ++i) {
                baby_x.push_back(mod(baby_steps[i].X * baby_z[i], N));
            }

            // 2. Baby-stepの点から 多項式 F(X) = Π(X - x_j) を構築
            Poly F = build_poly_tree(baby_x, 0, baby_x.size() - 1, N);

            // 3. Giant-step のアフィンX座標をバッチ逆元で集める
            PointXZ DQ = mont_ladder(P_XZ, current_phase2_D, A24_mont, ctx, N); // ★ D を current_phase2_D に変更
            int i_min = B1 / current_phase2_D;
            int i_max = B2 / current_phase2_D;
            
            int prev_idx = (i_min == 0) ? 1 : i_min - 1;
            PointXZ T_prev = mont_ladder(DQ, prev_idx, A24_mont, ctx, N);
            PointXZ T_curr = mont_ladder(DQ, i_min, A24_mont, ctx, N);
            
            std::vector<PointXZ> giant_steps;
            std::vector<cpp_int> giant_z;
            for (int i = i_min; i <= i_max; ++i) {
                giant_steps.push_back(T_curr);
                giant_z.push_back(T_curr.Z);
                PointXZ T_next;
                ctx.xADD(T_next, T_curr, DQ, T_prev); // ★ ここを変更
                T_prev = T_curr;
                T_curr = T_next;
            }
            
            // Z座標を一括で逆元に変換
            if (!batch_inverse(giant_z, N)) continue;
            
            // アフィンX座標 (giant_x) を計算
            std::vector<cpp_int> giant_x;
            for (size_t i = 0; i < giant_steps.size(); ++i) {
                giant_x.push_back(mod(giant_steps[i].X * giant_z[i], N));
            }

            // 4. Giant-stepの点から評価ツリーを構築し、F(X) を一気に代入して割り算
            std::vector<Poly> eval_tree;
            eval_tree.resize(4 * giant_x.size()); 
            build_eval_tree(1, giant_x, 0, giant_x.size() - 1, N, eval_tree);
            
            std::vector<cpp_int> eval_results;
            eval_results.reserve(giant_x.size()); 
            evaluate_down(1, F, 0, giant_x.size() - 1, N, eval_tree, eval_results);

            // 5. すべての評価値 F(x_giant) の積を取り、定期的に GCD を計算
            cpp_int accum = 1;
            int batch_count = 0;
            for (const cpp_int& res : eval_results) {
                if (res == 0) continue;
                accum = (accum * res) % N;
                
                if (++batch_count % 128 == 0) {
                    cpp_int d = gcd(accum, N);
                    if (d > 1 && d < N) {
                        global_ecm_state.found_factor = d;
                        global_ecm_state.phase = 2; 
                        global_ecm_state.final_Z = accum;
                        break;
                    }
                }
            }
            
            if (global_ecm_state.found_factor == 0 && batch_count % 128 != 0) {
                cpp_int d = gcd(accum, N);
                if (d > 1 && d < N) {
                    global_ecm_state.found_factor = d;
                    global_ecm_state.phase = 2; 
                    global_ecm_state.final_Z = accum;
                } else if (d == N) {
                    global_ecm_state.found_factor = N;
                }
            }
        }

        // ==========================================
        // 結果の返却
        // ==========================================
        if (global_ecm_state.found_factor > 0 && global_ecm_state.found_factor < N) {
            // ★ UI表示用に、モントゴメリ曲線のパラメータ A とアフィン座標 x を逆算して復元する
            // A24 = (A + 2) / 4 より、A = 4 * A24 - 2
            cpp_int A = mod(4 * A24 - 2, N);
            // X = u^3, Z = v^3 より、x = X / Z = u^3 * (v^3)^-1
            cpp_int x = mod(u3 * modInverse(v3, N), N);

            result.set("status", "found");
            result.set("factor", global_ecm_state.found_factor.str());
            result.set("phase", global_ecm_state.phase);
            if (global_ecm_state.phase == 1) {
                result.set("prime", std::to_string(global_ecm_state.current_prime));
                result.set("k", global_ecm_state.current_multiplier.str());
            }
            result.set("Z", global_ecm_state.final_Z.str());
            result.set("A", A.str());
            result.set("B", "0");       // モントゴメリ曲線移行によりY関連は非表示
            result.set("P_x", x.str());
            result.set("P_y", "0");     // 同上
            result.set("curves_tried", iter + 1);
            return result;
        }
    }
    
    result.set("status", "continue");
    result.set("curves_tried", num_curves);
    return result;
}

EMSCRIPTEN_BINDINGS(ecm_module) {
    function("is_prime_cpp", &is_prime_cpp);
    function("run_ecm_batch", &run_ecm_batch); // バッチ処理版に変更
    function("pollard_rho", &pollard_rho);
}

