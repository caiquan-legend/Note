#### min25筛

min25筛用来求积性函数前缀和 $\sum_{i = 1}^nf(i)$, 时间复杂度约为：$O(\frac {n^\frac 34}{logn})$，$10^{10}$数量级1s



##### 要求

积性函数 $f(x)$, 如果能够快速知道 $f(p^k)$ 处的值，且 $f(p)$ 可以表示为低阶多项式（项数较少或能表示成若干项使得每一项都是完全积性函数），那么min_25筛可以在 $O(\frac{n^ \frac34}{logn})$ 时间复杂度内求出 $\sum_{i = 1}^nf(i)$

由于其思想与埃氏筛类似，所以亦被称作扩展埃氏筛。



### 内容

#### 第一步

###### 推导

目标：$g(n) = \sum_{p \le n}f(p)$

不妨设 $f(p)$ 是完全积性函数，如果不是则可以尝试拆成若干项完全积性函数，分别求然后相加。

首先要线性筛求出 $\sqrt n$以内的质数。

$g(n)$ 难以直接求解，考虑 $dp$。设$g(n, j) = \sum_{i = 1}^n f(i)[i \in \mathbb P \or i的最小质因子 > p_j]$，其中 $p_j$ 表示第 $j$ 个素数，那么我们需要的就是 $g(n, k), k为最小的满足 p_k \ge \sqrt n$, 考虑从 $j - 1$ 到 $j$, 那么最小质因子为 $p_j$ 的合数会被筛掉，那么它们的贡献需要减去。

则有状态转移 ：
$$
g(n, j) = g(n, j - 1) - f(p_j)\Big(g(\lfloor \frac n{p_j} \rfloor, j - 1) - g(p_{j - 1}, j - 1) \Big)
$$
系数 $f(p_j)$ 表示由 $f(p)$ 是完全积性函数，所以可以把它分解为$f(n) = f(p)f(\frac np)$。

###### 个人证明：

$$
\begin{flalign}
g(n, j) = & \ \ g(n, j - 1) - \sum_{i = 1}^n F(i)[i的最小质因子是p_j] & \\
= & \ \ g(n, j - 1) - \sum_{i = 1}^n F(p_j \cdot \frac i{p_j})[i的最小质因子是p_j] & \\
= & \ \ g(n, j - 1) - F(p_j)\sum_{i = 1}^n F(\frac i{p_j})[i的最小质因子是p_j] & \\
= & \ \ g(n, j - 1) - F(p_j)\sum_{i = 1}^{\lfloor \frac n{p_j} \rfloor} F(i)[i的最小质因子 > p_{j - 1}] & \\
= & \ \ g(n, j - 1) - F(p_j)(\sum_{i = 1}^{\lfloor \frac n{p_j} \rfloor}F(i)[i \in \mathbb P \or i的最小质因子 > p_{j - 1}] - \sum_{i = 1}^{\lfloor \frac n{p_j} \rfloor}F(i)[i \in \mathbb P \and i的最小质因子 \le p_{j - 1}]) \ \ (容斥原理)& \\
= & \ \ g(n, j - 1) - F(p_j)\Big(g(n, j - 1) - g(p_{j - 1}, j - 1)\Big)
\end{flalign}
$$





由于有公式 $\lfloor \frac{\lfloor \frac ab \rfloor}c \rfloor = \lfloor \frac a{bc} \rfloor$, 因此容易发现上述式子只会用到形如 $\frac nx, x \le n$ 的点处的 $dp$ 值，即第一项的状态数是 $O(\sqrt n)$, （实际实现时的时候注意状态数是 $2 \sqrt n$。我们预处理出 $O(\sqrt n)$ 个数，将它们离散化顺带求出 $g(x, 0)$，然后 $dp$ 即可。

###### 注意到上述转移式还有 $g(p_{j - 1}, j -1)$ 一项，我们只处理了所有形如 $\lfloor \frac nx \rfloor$的数，实际上并不会漏掉某些质数 $p$。注意到能用来转移的 $p$ 一定满足 $p \le \sqrt n$, 我们只需证明 $\forall k \le \sqrt n, \exists x, k= \lfloor \frac nx \rfloor$.

###### 证明:

1. 若 $k = \lfloor \sqrt n \rfloor$, 设 $n = k^2 + d$, 则由于 $k^2 \le n < (k+1)^2$, 故 $d \in [0, 2k]$.
   1. $ d \in [0, k)$, 那么 $\lfloor \frac nk \rfloor = k + \lfloor \frac dk \rfloor = k$
   2. $d \in [k, 2k]$, 那么 $\lfloor \frac n{k + 1} \rfloor = \lfloor \frac {k^2 + k + d - k}{k + 1} \rfloor = k + \lfloor \frac{d - k}{k + 1} \rfloor = k$
2. 若 $k < \lfloor \sqrt n \rfloor$, 即 $k \le \sqrt n - 1$, 假设存在$i$, $\lfloor \frac n{i + 1} < k < \lfloor \frac ni \rfloor$, 此时 $k$ 恰好在两个连续的 $\lfloor \frac nx \rfloor$ 之间，即不可被表出。则 $ \frac n{i + 1} < k$, 故$n < k(i + 1)$, 从而$k < \lfloor \frac ni \rfloor < \lfloor \frac {k(i + 1)}i \rfloor = k + \lfloor \frac ki \rfloor$。另一方面，$\frac n{i + 1} < k < \sqrt n$, 所以 $i + 1 > \sqrt n$, 于是 $i > \sqrt n - 1 \ge k$, 因此$\lfloor \frac ki \rfloor = 0$，于是得到 $k < k$，假设不成立，原命题成立。 





### 

#### 第二步

###### 推导

目标：求$S(n) = \sum_{i \le n} f(i)$, 与第一步类似，设 $S(n, j) = \sum_{i = 1}^nf(i)[i的最小质因子 > p_j]$, 但 $f$ 不需要再拆分成单项式，直接是原函数即可（不需要依赖完全积性，只需要积性即可，但要能快速计算 $f(p^k)$的值）



##### 方法一		一般情况下复杂度$O(1^{1 - \varepsilon})$, $\varepsilon$是一个无穷小量，但在 $n \le 10^{13}$ 时，复杂度为 $O(\frac {n^\frac 34}{logn})$

考虑把贡献拆成质数和合数的，合数枚举最小质因子以及次数，于是有转移：
$$
S(n, j) = g(n) - g(p_j) + \sum_{k = j + 1}^{\sqrt n}\sum_{e = 1}^{p_k^e \le n}f(p_k^e) \Big( S(\Big\lfloor \frac n{p_k^e} \Big\rfloor, k) + [e \neq 1] \Big)
$$
最后一项$[e \neq 1]$的意思是，对于 $e = 1$ 的情况，$S$ 没有计算 `1` 的贡献，因为此时 $p_k \times 1$ 是质数，其贡献在之前计算过。

对于 $e > 1$ 的情况， $p_k^1 \times 1$ 是合数，贡献没有计算需要补上。直接暴力递归计算，且不需要记忆化。

###### 个人证明：

$$
\begin{flalign}
S(n, j) = & \ \ g(n) - (p_j) + \sum_{i = 1}^n{f(i)[i的最小质因子 > p_j] - \sum_{k = j + 1}^{p_k^2 \le n}f(p_k)} & \\
= & \ \ g(n) - g(p_j) + \sum_{i = 1}^n{f(i = p_{j+1}^{r_1}p_{j+2}^{r_2}\cdots p_{k}^{r_k})[i的最小质因子 > p_j] - \sum_{k = j + 1}^{p_k^2 \le n}f(p_k)} & \\ 
= & \ \ g(n) - g(p_j) + \sum_{k = j + 1}^{p_k^2 \le n}\sum_{e = 1}^{p_k^e \le n}\sum_{i = 1}^n{f(p_k^e)\big(1 + f(\frac{i}{p_k^e})[i的最小质因子 > p_j]\big)} - \sum_{k = j + 1}^{p_k^2 \le n}f(p_k) \ \ (需要加上i=p^k时的贡献) & \\
= & \ \ g(n) - g(p_j) + \sum_{k = j + 1}^{p_k^2 \le n}\sum_{e = 1}^{p_k^e \le n}f(p_k^e)\sum_{i = 1}^n\big(1 + f(\frac{i}{p_k^e})[i的最小质因子 > p_j]\big) - \sum_{k = j + 1}^{p_k^2 \le n}f(p_k) & \\
= & \ \ g(n) - g(p_j) + \sum_{k = j + 1}^{p_k^2 \le n}\sum_{e = 1}^{p_k^e \le n}f(p_k^e)\sum_{i = 1}^{n / p_k^e}\big(1 + f(i)[i的最小质因子 > p_k]\big) - \sum_{k = j + 1}^{p_k^2 \le n}f(p_k) & \\
= & \ \ g(n) - g(p_j) + \sum_{k = j + 1}^{p_k^2 \le n}\sum_{e = 1}^{p_k^e \le n}f(p_k^e)\Big(S(\Big \lfloor \frac{n}{p_k^e} \Big \rfloor, k) + 1\Big) - \sum_{k = j + 1}^{p_k^2 \le n}f(p_k) & \\ 
= & \ \ g(n) - g(p_j) + \sum_{k = j + 1}^{p_k^2 \le n}\sum_{e = 1}^{p_k^e \le n}f(p_k^e)\Big(S(\Big \lfloor \frac{n}{p_k^e} \Big \rfloor, k) + 1 - [e = 1]\Big) & \\ 
= & \ \ g(n) - g(p_j) + \sum_{k = j + 1}^{p_k^2 \le n}\sum_{e = 1}^{p_k^e \le n}f(p_k^e)\Big(S(\Big \lfloor \frac{n}{p_k^e} \Big \rfloor, k) + [e \neq 1]\Big) & \\
\end{flalign}
$$





##### 方法二		复杂度$O(\frac{n^ \frac34}{logn})$, 常数略大

也是将贡献拆成质数和合数，只是采用类似于第一步的递推方式：
$$
S(n, j) = f(p_{j + 1}) + S(n, j + 1) + \sum_{p_{j + 1} \le \sqrt n，e = 1}^{p_{j + 1}^e \le n}f(p_{j + 1}^e) \Big( S(\Big\lfloor \frac n{p_{j + 1}^e} \Big \rfloor, j + 1) + [e \neq 1] \Big)
$$
直接转移的复杂度有误，我们要严格控制只有 $ \le \sqrt n$ 的质数才转移，对于 $p_{j + 1} > \sqrt n$ 的状态 $S(n, j)$，不显式计算，用到的时候特判为$g(n) - g(p_j)$。所以我们要更新状态 $S(n, j)$ 时，一定有 $p_{j + 1} \le \sqrt n$, 此时需要用到 $S(n, j + 1)$, 即便要特判为 $g(n) - g(p_{j + 1})$, 由步骤一的论述可知 $g(p_{j + 1})$ 已经算出，所以不会有问题。需要注意 $n = 2, 3$ 的状态不会被任何 $p$ 更新到，所以也要特判。

###### 改进

既然采用了类似第一步的递推方法，干脆直接套用第一步的状态设计。

设 $S(n, j) = \sum_{i = 1}^nf(i)[i是质数或i的最小质因子 > p_j]$, 那么参考第一步的方程，有转移：
$$
S(n, j) = S(n, j + 1) + \sum_{e = 1}^{p_{j + 1} \le \sqrt n，p_{j + 1}^e \le n}f(p_{j + 1}^e) \Big(S(\Big\lfloor \frac n{p_{j + 1}^e} \Big\rfloor, j + 1) - g(min\{\Big\lfloor \frac n{p_{j + 1}^e} \Big\rfloor p_{j + 1}\}) + [e \neq 1] \Big)
$$
最后面 $-g(\cdots)$ 是因为 $S(\lfloor \frac n{p_{j + 1}^e} \rfloor， j + 1)$ 中包含了 $\le p_{j + 1}$ 的质数的贡献，与转移式的意义不符需要减去。我们可以进一步简化上式：
$$
S(n, j) = S(n, j + 1) + \sum_{p_{j + 1} \le \sqrt n，e = 1}^{p_{j + 1}^{e + 1} \le n}f(p_{j + 1}^e) \Big(S(\Big\lfloor \frac n{p_{j + 1}^e} \Big\rfloor, j + 1) - g(p_{j + 1}) \Big) + f(p_{j + 1}^{e + 1})
$$
这是因为当 $p_{j + 1}^{e + 1} > n$ 的时候 $s(\cdots) - g(p_{j + 1})$ 这一项不会产生贡献。

所以只需要设置 $S(n, +\infty)$ 初值为 $g(n)$, 即可进行和第一步几乎一样的 DP过程了.

此外，这种方法不仅计算出了 $S(n)$, 也顺带计算出了所有 $S(\lfloor n/x \rfloor)$.





###

#### 实现细节

第一步中要使空间为$O(\sqrt n)$，考虑根号分治表示下标。即对 $x < \sqrt n$ 用 $x$ 映射到下标，对于 $x > \sqrt n$ 用 $n / x$ 映射下标

[P5325 【模板】Min_25 筛 - 洛谷](https://www.luogu.com.cn/problem/P5325)

第一步将原函数拆成$g_1(p) = p, g_2(p) = p^2$ 两个完全积性函数计算

##### 方法一

```cpp
constexpr int N = 2e5 + 10, mod = 1e9 + 7, i2 = 5e8 + 4, i6 = 166666668;
i64 s[N], g1[N], g2[N], v[N];
int p[N], vis[N], cnt;
int id1[N], id2[N], m;
i64 n;
int get(i64 x) { return x < N ? id1[x] : id2[n / x]; }
i64 S1(i64 x) { return x %= mod, x * (x + 1) % mod * i2 % mod; }
i64 S2(i64 x) { return x %= mod, x * (x + 1) % mod * (2 * x + 1) % mod * i6 % mod; }
i64 sq(i64 x) { return x %= mod, x * x % mod; }
i64 F(i64 x) { return x %= mod, (sq(x) - x + mod) % mod; }

void init(int n) {
    for(int i = 2; i <= n; ++ i) {
        if(!vis[i]) p[ ++ cnt] = i;
        for(int j = 1; p[j] <= n / i; ++ j) {
            vis[p[j] * i] = 1;
            if(i % p[j] == 0) {
                break;
            }
        }
    }
}

i64 S(i64 x, i64 y) {
    if(p[y] > x) return 0;
    i64 res = (g2[get(x)] - g1[get(x)] - g2[get(p[y])] + g1[get(p[y])] + 2 * mod) % mod;
    for(int i = y + 1; i <= cnt && p[i] <= x / p[i]; ++ i) {
        i64 w = p[i];
        for(int j = 1; w <= x / p[i]; ++ j, w = w * p[i]) {
            res = (res + F(w) * S(x / w, i) % mod + F(w * p[i])) % mod;
        }
    }
    return res;
}

void solve() {
    init(sqrt(n) + 1);
    for(i64 l = 1, r; l <= n; l = r + 1) {
        r = n / (n / l), v[ ++ m] = n / l;
        if(v[m] < N) id1[v[m]] = m;
        else id2[n / v[m]] = m;
        g1[m] = (S1(v[m]) - 1 + mod) % mod, g2[m] = (S2(v[m]) - 1 + mod) % mod;
    }
    
    for(int j = 1; j <= cnt; ++ j) {
        for(int i = 1; i <= m && p[j] <= v[i] / p[j]; ++ i) {
            g1[i] = (g1[i] - p[j] * (g1[get(v[i] / p[j])] - g1[get(p[j - 1])]) % mod + mod) % mod;
            g2[i] = (g2[i] - sq(p[j]) * (g2[get(v[i] / p[j])] - g2[get(p[j - 1])]) % mod + mod) % mod;
        }
    }
    cout << (S(n, 0) + 1) % mod << endl;
}
```



##### 方法二

```cpp
constexpr int N = 2e5 + 10, mod = 1e9 + 7, i2 = 5e8 + 4, i6 = 166666668;
i64 s[N], g1[N], g2[N], v[N];
int p[N], vis[N], cnt;
int id1[N], id2[N], m;
i64 n;
int get(i64 x) { return x < N ? id1[x] : id2[n / x]; }
i64 S1(i64 x) { return x %= mod, x * (x + 1) % mod * i2 % mod; }
i64 S2(i64 x) { return x %= mod, x * (x + 1) % mod * (2 * x + 1) % mod * i6 % mod; }
i64 sq(i64 x) { return x %= mod, x * x % mod; }
i64 F(i64 x) { return x %= mod, (sq(x) - x + mod) % mod; }

void init(int n) {
    for(int i = 2; i <= n; ++ i) {
        if(!vis[i]) p[ ++ cnt] = i;
        for(int j = 1; p[j] <= n / i; ++ j) {
            vis[p[j] * i] = 1;
            if(i % p[j] == 0) {
                break;
            }
        }
    }
}

void solve() {
    init(sqrt(n) + 1);
    for(i64 l = 1, r; l <= n; l = r + 1) {
        r = n / (n / l), v[ ++ m] = n / l;
        if(v[m] < N) id1[v[m]] = m;
        else id2[n / v[m]] = m;
        g1[m] = (S1(v[m]) - 1 + mod) % mod, g2[m] = (S2(v[m]) - 1 + mod) % mod;
    }
    
    for(int j = 1; j <= cnt; ++ j) {
        for(int i = 1; i <= m && p[j] <= v[i] / p[j]; ++ i) {
            g1[i] = (g1[i] - p[j] * (g1[get(v[i] / p[j])] - g1[get(p[j - 1])]) % mod + mod) % mod;
            g2[i] = (g2[i] - sq(p[j]) * (g2[get(v[i] / p[j])] - g2[get(p[j - 1])]) % mod + mod) % mod;
        }
    }
    
    for(int i = 1; i <= m; ++ i) {
        s[i] = g1[i] = (g2[i] - g1[i] + mod) % mod;
    }
    for(int j = cnt; j; -- j) {
        for(int i = 1; i <= m && p[j] <= v[i] / p[j]; ++ i) {
            i64 w = p[j];
            for(int k = 1; w <= v[i] / p[j]; ++ k, w *= p[j]) {
                s[i] = (s[i] + F(w) * (s[get(v[i] / w)] - g1[get(p[j])] + mod) % mod + F(w * p[j])) % mod;
            }
        }
    }
    cout << (s[get(n)] + 1) % mod << endl;
}
```

