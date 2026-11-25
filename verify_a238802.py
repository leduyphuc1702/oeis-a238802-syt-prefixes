"""Verification for the OEIS A238802 note (independent of Jaggard's theorem).

1. Brute-force SYT enumeration for n<=13 -> triangle T(n,k), compared with
   the OEIS data and with B(n,k)=sum_j C(n-k,j) I_{n-k-j}.
2. Brute-force over involutions: w(1)>...>w(k) count equals B(n,k).
3. Symbolic (generic k) expansion of B(n,k)/I_n in powers of n^{-1/2},
   from the asymptotic expansion of rho_n = I_n/I_{n-1}.
4. High-precision numeric check at n up to 20000.
"""
import itertools
from math import comb, factorial
import sympy as sp
import mpmath as mp

OEIS = [1,0,1,0,1,1,0,2,1,1,0,5,3,1,1,0,13,8,3,1,1,0,38,24,9,3,1,1,0,116,74,
        28,9,3,1,1,0,382,246,93,29,9,3,1,1,0,1310,848,321,98,29,9,3,1,1,0,
        4748,3088,1168,350,99,29,9,3,1,1,0,17848,11644,4404,1302,356,99,29,9,3,1,1]

I = [1, 1]
for n in range(2, 20001):
    I.append(I[-1] + (n - 1) * I[-2])


def B(n, k):
    return sum(comb(n - k, j) * I[n - k - j] for j in range(min(k, n - k) + 1))


# --- 1. SYT enumeration -------------------------------------------------
def syt_first_column_runs(N):
    """For each n<=N, histogram of k = maximal run 1..k in column 1."""
    hist = [[0] * (n + 1) for n in range(N + 1)]
    pass

    def rec(shape, n, run, alive):
        # shape: tuple of row lengths; entry n has just been placed
        hist[n][run] += 1
        if n == N:
            return
        for i in range(len(shape) + 1):
            row_len = shape[i] if i < len(shape) else 0
            if i == 0 or shape[i - 1] > row_len:
                new = list(shape)
                if i < len(shape):
                    new[i] += 1
                else:
                    new.append(1)
                in_col1 = row_len == 0
                nalive = alive and in_col1
                rec(tuple(new), n + 1, run + 1 if nalive else run, nalive)

    rec((), 0, 0, True)
    return hist


N = 13
hist = syt_first_column_runs(N)
flat = [x for row in hist for x in row]
assert flat[:len(OEIS)] == OEIS, "SYT enumeration disagrees with OEIS"
for n in range(1, N + 1):
    for k in range(0, n + 1):
        b = B(n, k) - (B(n, k + 1) if k + 1 <= n else 0)
        assert hist[n][k] == b, (n, k, hist[n][k], b)
print("1. SYT enumeration n<=13 matches OEIS triangle and B(n,k)-B(n,k+1)")

# --- 2. involution brute force ------------------------------------------
for n in range(1, 10):
    invs = [p for p in itertools.permutations(range(n)) if all(p[p[i]] == i for i in range(n))]
    assert len(invs) == I[n]
    for k in range(1, n + 1):
        c = sum(all(p[i] > p[i + 1] for i in range(k - 1)) for p in invs)
        assert c == B(n, k), (n, k)
print("2. involutions with w(1)>...>w(k) counted by B(n,k), n<=9")

# --- 3. symbolic expansion ----------------------------------------------
e = sp.symbols('e', positive=True)          # e = n^{-1/2}
ORDER = 6                                   # keep terms up to e^ORDER
cs = sp.symbols('c1:%d' % (ORDER + 2))
# ansatz rho_n = n^{1/2} (1 + sum c_i e^i)


def rho_shift(i, cvals):
    """rho_{n-i} / n^{1/2} as a series in e, with n = e^{-2}."""
    # (n-i)^{1/2} = n^{1/2} (1 - i e^2)^{1/2};  (n-i)^{-1/2} = e (1-i e^2)^{-1/2}
    f = sp.sqrt(1 - i * e**2)
    g = e / sp.sqrt(1 - i * e**2)            # (n-i)^{-1/2}
    ser = f * (1 + sum(c * g**(m + 1) for m, c in enumerate(cvals)))
    return sp.series(ser, e, 0, ORDER + 1).removeO()


# solve rho_n rho_{n-1} = rho_{n-1} + (n-1) for c's
expr = sp.expand(rho_shift(0, cs) * rho_shift(1, cs) * e**-2) - \
    sp.expand(rho_shift(1, cs) * e**-1) - (e**-2 - 1)
expr = sp.expand(sp.series(sp.expand(expr * e**2), e, 0, ORDER + 2).removeO())
sol = {}
for p in range(1, ORDER + 1):
    eq = sp.expand(expr.coeff(e, p).subs(sol))
    var = cs[p - 1]
    s = sp.solve(eq, var)
    assert len(s) == 1, (p, s)
    sol[var] = s[0]
cv = [sol[c] for c in cs[:ORDER]]
print("3. rho_n = sqrt(n) (1 + %s ...)" % sp.expand(sum(c * e**(m + 1) for m, c in enumerate(cv))))

i, M, kk, d = sp.symbols('i M k d', integer=True, nonnegative=True)
# log(1/rho_{n-i}) + (1/2) log n  as polynomial in i per power of e
logterm = sp.expand(sp.series(-sp.log(rho_shift(i, cv)), e, 0, 5).removeO())
# sum over i = 0..M-1 of logterm, symbolic in M
S_log = 0
for p in range(0, 5):
    coeff = sp.expand(logterm.coeff(e, p))
    S_log += sp.summation(coeff, (i, 0, M - 1)) * e**p
S_log = sp.expand(S_log)


def ratio_term(dd):
    """C(n-k,j) I_{n-k-j} / I_n for j=k-dd, times n^{dd/2}, series in e (generic k)."""
    j = kk - dd
    Mv = kk + j
    # prod_{i<M} 1/rho_{n-i} = n^{-M/2} exp(S_log)
    prod_rho = sp.exp(S_log.subs(M, Mv))
    # C(n-k,j) = n^j/j! prod_{l<j} (1 - (k+l) e^2)
    logC = sp.summation(sp.expand(sp.series(sp.log(1 - (kk + i) * e**2), e, 0, 5).removeO()), (i, 0, j - 1))
    # n^j n^{-M/2} = n^{-dd/2}... wait: j - M/2 = (k-dd) - (2k-dd)/2 = -dd/2
    full = sp.exp(sp.expand(S_log.subs(M, Mv) + logC))
    return sp.series(full, e, 0, 4).removeO()   # times 1/j!


total = 0
for dd in range(0, 4):
    j = kk - dd
    rt = sp.expand(ratio_term(dd))
    print("   d=%d: (k-%d)! n^(%d/2) C(n-k,k-%d) I_(n-2k+%d)/I_n = 1 + %s" % (
        dd, dd, dd, dd, dd, " + ".join("(%s) n^-%d/2" % (sp.factor(rt.coeff(e, p)), p)
                                       for p in range(1, 4 - dd))))
    total += e**dd * rt / sp.factorial(j)
total = sp.series(sp.expand(total), e, 0, 4).removeO()
norm = sp.simplify(sp.expand(total * sp.factorial(kk)))
coeffs = [sp.factor(sp.simplify(sp.combsimp(sp.expand_func(norm.coeff(e, p))))) for p in range(4)]
print("   k! B(n,k)/I_n = 1 + (%s) n^-1/2 + (%s) n^-1 + (%s) n^-3/2 + O(n^-2)" % tuple(coeffs[1:]))
target = kk * (kk - 1) * (kk - 2) / 3
for kv in range(3, 13):
    vals = [sp.nsimplify(sp.simplify(c.subs(kk, kv))) for c in coeffs]
    assert vals[0] == 1 and vals[1] == 0 and vals[2] == 0, (kv, vals)
    assert vals[3] == target.subs(kk, kv), (kv, vals)
print("   checked: coefficient of n^-3/2 equals k(k-1)(k-2)/3 for k=3..12 (generic formula above)")

# --- 4. numerics ----------------------------------------------------------
mp.mp.dps = 60
print("4. n^2 * (n^1.5 (T/I_n (k+1)!/k - 1) - (k^2-1)(k-3)/3) should stay bounded:")
for n in (1000, 4000, 10000, 20000):
    row = []
    for k in (1, 2, 3, 4, 8, 10, 12):
        T = B(n, k) - B(n, k + 1)
        D = mp.mpf(n)**mp.mpf(1.5) * (mp.mpf(factorial(k + 1) * T) / (k * I[n]) - 1)
        row.append(mp.nstr((D - mp.mpf((k * k - 1) * (k - 3)) / 3) * mp.sqrt(n), 6))
    print("  n=%d  sqrt(n)*(D - target), k=1,2,3,4,8,10,12:" % n, row)
