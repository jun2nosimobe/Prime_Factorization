#include <emscripten/bind.h>
#include <emscripten/val.h>
#include <string>
#include <cstdlib>
#include <vector>
#include <map>

// Boostの代わりにmini-gmpをインクルード
#include "mini-gmp.h"

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

// Boostの gcd をエミュレートするグローバル関数
inline cpp_int gcd(const cpp_int& a, const cpp_int& b) {
    cpp_int r;
    mpz_gcd(r.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
    return r;
}
// =======================================================

// --- モントゴメリ乗算用クラス・ヘルパー ---
struct MontgomeryContext {
    cpp_int n;
    cpp_int r;     
    int r_bits;    
    cpp_int mask;  // ★必須: 下位ビット抽出用マスク
    cpp_int n_prime; 
    
    MontgomeryContext(const cpp_int& modulus) {
        n = modulus;
        int bits = 0;
        cpp_int temp = n;
        while (temp > 0) { bits++; temp >>= 1; }
        r_bits = ((bits + 63) / 64) * 64;
        r = cpp_int(1) << r_bits;
        
        // ★修正: 正しいマスク (2^r_bits - 1) を生成
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

    cpp_int to_mont(const cpp_int& a) const {
        return (a << r_bits) % n;
    }

    cpp_int from_mont(const cpp_int& a) const {
        return mont_mul(a, 1);
    }

    cpp_int mont_mul(const cpp_int& a, const cpp_int& b) const {
        cpp_int t = a * b;
        // ★修正: -1 が抜けていた誤った式を削除し、正しく mask を使用する
        cpp_int m = cpp_int((t * n_prime) & mask);
        cpp_int u = cpp_int(t + m * n) >> r_bits; 
        if (u >= n) {
            u -= n;
        }
        return u;
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




// --- モントゴメリ曲線用のX, Z座標のみの構造体 ---
struct PointXZ {
    cpp_int X, Z;
};

// モントゴメリ・ラダー: 差分加算 (xADD)
PointXZ xADD(const PointXZ& P, const PointXZ& Q, const PointXZ& P_minus_Q, const MontgomeryContext& ctx, const cpp_int& N) {
    cpp_int U = mod(P.X - P.Z, N);
    cpp_int V = mod(Q.X + Q.Z, N);
    cpp_int U_V = ctx.mont_mul(U, V);
    
    cpp_int W = mod(P.X + P.Z, N);
    cpp_int Y = mod(Q.X - Q.Z, N);
    cpp_int W_Y = ctx.mont_mul(W, Y);
    
    cpp_int sum = mod(U_V + W_Y, N);
    cpp_int diff = mod(U_V - W_Y, N);
    
    cpp_int X_plus = ctx.mont_mul(P_minus_Q.Z, ctx.mont_mul(sum, sum));
    cpp_int Z_plus = ctx.mont_mul(P_minus_Q.X, ctx.mont_mul(diff, diff));
    
    return {X_plus, Z_plus};
}

// モントゴメリ・ラダー: 2倍算 (xDBL)
PointXZ xDBL(const PointXZ& P, const cpp_int& A24_mont, const MontgomeryContext& ctx, const cpp_int& N) {
    cpp_int U = mod(P.X + P.Z, N);
    cpp_int V = mod(P.X - P.Z, N);
    
    cpp_int U2 = ctx.mont_mul(U, U);
    cpp_int V2 = ctx.mont_mul(V, V);
    
    cpp_int diff = mod(U2 - V2, N);
    
    cpp_int X_out = ctx.mont_mul(U2, V2);
    cpp_int term = mod(V2 + ctx.mont_mul(A24_mont, diff), N);
    cpp_int Z_out = ctx.mont_mul(diff, term);
    
    return {X_out, Z_out};
}

// モントゴメリ・ラダーによる超高速スカラー倍
PointXZ mont_ladder(const PointXZ& P, cpp_int k, const cpp_int& A24_mont, const MontgomeryContext& ctx, const cpp_int& N) {
    if (k == 0) return {ctx.to_mont(1), 0};
    
    PointXZ R0 = P;
    PointXZ R1 = xDBL(P, A24_mont, ctx, N);
    
    int bits = 0;
    cpp_int temp = k;
    while(temp > 0) { bits++; temp >>= 1; }
    
    // 最上位ビット(常に1)の次は bits - 2 からスタート
    for (int i = bits - 2; i >= 0; --i) {
        bool bit = ((k >> i) & 1) == 1;
        if (!bit) {
            R1 = xADD(R0, R1, P, ctx, N);
            R0 = xDBL(R0, A24_mont, ctx, N);
        } else {
            R0 = xADD(R0, R1, P, ctx, N);
            R1 = xDBL(R1, A24_mont, ctx, N);
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

// クロネッカー代入を用いた高速多項式乗算 (mod N)
Poly kronecker_multiply(const Poly& A, const Poly& B, const cpp_int& N) {
    if (A.degree() < 0 || B.degree() < 0) return Poly();

    int n = A.degree();
    int m = B.degree();
    int min_len = std::min(n + 1, m + 1);

    // Nのビット数を計算
    int bits_N = 0;
    cpp_int temp_N = N;
    while (temp_N > 0) { bits_N++; temp_N >>= 1; }

    // min_len のビット数を計算
    int log2_min = 0;
    int temp_min = min_len;
    while (temp_min > 0) { log2_min++; temp_min >>= 1; }

    // 係数がオーバーフローしないための安全なシフト幅 K
    int K = 2 * bits_N + log2_min + 2;

    // パッキング
    cpp_int packed_A = 0;
    for (int i = 0; i <= n; ++i) {
        packed_A += (A.coef[i] << (i * K));
    }

    cpp_int packed_B = 0;
    for (int i = 0; i <= m; ++i) {
        packed_B += (B.coef[i] << (i * K));
    }

    // GMPの巨大整数乗算
    cpp_int packed_C = packed_A * packed_B;

    // アンパックと mod N
    Poly C;
    C.coef.resize(n + m + 1, 0);
    cpp_int mask = (cpp_int(1) << K) - 1;

    for (int i = 0; i <= n + m; ++i) {
        cpp_int coeff = (packed_C >> (i * K)) & mask;
        C.coef[i] = mod(coeff, N);
    }
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
    
    // 多項式の筆算（長除算）
    while (A.degree() >= degB) {
        int d = A.degree() - degB;
        cpp_int factor = A.coef.back();
        for (int i = 0; i <= degB; ++i) {
            cpp_int term = mod(B.coef[i] * factor, N);
            A.coef[i + d] = mod(A.coef[i + d] - term, N);
        }
        A.clean();
    }
    return A;
}

// 評価対象の点群から、多項式除算用の木（Subproduct Tree）を構築
void build_eval_tree(int node, const std::vector<cpp_int>& points, int left, int right, const cpp_int& N, std::vector<Poly>& tree) {
    if ((int)tree.size() <= node) tree.resize(node + 1);
    
    if (left == right) {
        tree[node] = Poly({mod(-points[left], N), 1}); // (X - p)
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
            int B2 = B1 * 50; 
            int D = 210;      
            int coprimes[] = {1, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 
                              53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103};
            
            // 1. Baby-step のアフィンX座標を集める
            std::vector<cpp_int> baby_x;
            for (int j : coprimes) {
                PointXZ step = mont_ladder(P_XZ, j, A24_mont, ctx, N);
                // 射影座標からアフィン座標 X/Z への変換
                cpp_int x_affine = mod(step.X * modInverse(step.Z, N), N);
                if (global_ecm_state.found_factor > 0) break;
                baby_x.push_back(x_affine);
            }
            
            // もしアフィン変換時の modInverse で素因数が見つかれば、ループの最後(結果返却)へ飛ぶ
            if (global_ecm_state.found_factor > 0) continue;

            // 2. Baby-stepの点から 多項式 F(X) = Π(X - x_j) を構築
            Poly F = build_poly_tree(baby_x, 0, baby_x.size() - 1, N);

            // 3. Giant-step のアフィンX座標を集める (差分加算の漸化式を使用)
            PointXZ DQ = mont_ladder(P_XZ, D, A24_mont, ctx, N);
            int i_min = B1 / D;
            int i_max = B2 / D;
            
            int prev_idx = (i_min == 0) ? 1 : i_min - 1;
            PointXZ T_prev = mont_ladder(DQ, prev_idx, A24_mont, ctx, N);
            PointXZ T_curr = mont_ladder(DQ, i_min, A24_mont, ctx, N);
            
            std::vector<cpp_int> giant_x;
            for (int i = i_min; i <= i_max; ++i) {
                cpp_int x_affine = mod(T_curr.X * modInverse(T_curr.Z, N), N);
                if (global_ecm_state.found_factor > 0) break;
                giant_x.push_back(x_affine);
                
                PointXZ T_next = xADD(T_curr, DQ, T_prev, ctx, N);
                T_prev = T_curr;
                T_curr = T_next;
            }
            if (global_ecm_state.found_factor > 0) continue;

            // 4. Giant-stepの点から評価ツリーを構築し、F(X) を一気に代入して割り算
            std::vector<Poly> eval_tree;
            build_eval_tree(1, giant_x, 0, giant_x.size() - 1, N, eval_tree);
            
            std::vector<cpp_int> eval_results;
            evaluate_down(1, F, 0, giant_x.size() - 1, N, eval_tree, eval_results);

            // 5. すべての評価値 F(x_giant) の積を取り、NとのGCDを計算
            cpp_int accum = 1;
            for (const cpp_int& res : eval_results) {
                if (res == 0) continue;
                accum = mod(accum * res, N);
            }
            
            cpp_int d = gcd(accum, N);
            if (d > 1 && d < N) {
                global_ecm_state.found_factor = d;
                global_ecm_state.phase = 2; 
                global_ecm_state.final_Z = accum; // 代替として記録
            } else if (d == N) {
                global_ecm_state.found_factor = N;
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

