from sage.all import (
    GF,
    PolynomialRing,
    PowerSeriesRing,
    divisors,
    Integer,
    set_random_seed,
    random_prime,
    is_prime,
    next_prime,
)
import random
import functools
import collections

@functools.cache
def ntdivs(n):
    """Returns a tuple of the non-trivial divisors of integer n."""
    return tuple(divisors(n)[1:-1])

def decomp_dense_s(f, s):
    """
    Complete fast decomposition: f(x) = g(h(x))
    f: polynomial
    s: target degree of h
    """
    n = f.degree()
    r = n // s
    PR = f.parent()
    F = PR.base_ring()

    # 1. Fast h search via Power Series
    f_monic = f / f.leading_coefficient()
    S_ring = PowerSeriesRing(F, 'y', default_prec=s)

    # Reverse f to treat as power series at infinity
    P_ser = S_ring(f_monic.list()[::-1][:s])
    h_rev = P_ser.nth_root(r).list()
    pad = s - len(h_rev)
    h_rev.extend(F.zero() for _ in range(pad))
    h = PR(h_rev[::-1]).shift(1) # Convert back to polynomial

    # 2. Fast g search via D&C
    g = find_g(f, h)

    # 3. Verification (Crucial for decomposition)
    if g is not None and g(h) == f:
        return g, h
    return None

def find_g(f, h):
    """
    Finds g such that f(x) = g(h(x)) using a Divide and Conquer approach.
    This is much faster than iterative subtraction for high degrees.
    """
    PR = f.parent()

    if f.is_constant():
        return f
    elif h.is_constant():
        return None

    r, rem = divmod(f.degree(), h.degree())

    if rem:
        return None
    elif r == 1:
        # f = g1*h + g0. Since f = g(h), and deg(g)=1,
        # g1 is leading_coeff(f)/leading_coeff(h)
        g1 = f.leading_coefficient() / h.leading_coefficient()
        g0 = (f - g1 * h)[0]
        return PR([g0, g1])

    # Divide and Conquer step
    m = r // 2
    h_m = h**m

    # f = q * h_m + r
    q, r_poly = f.quo_rem(h_m)

    # Recurse
    # g = g_upper * x^m + g_lower
    g_lower = find_g(r_poly, h)
    g_upper = find_g(q, h)

    if g_lower is None or g_upper is None:
        return None
    else:
        # Shift g_upper by m and add
        return g_lower + g_upper.shift(m)

def decomp_dense(f):
    """Computes a decomposition of f = g(h), if it exists."""
    d = f.degree()
    PR = f.parent()
    F = PR.base_ring()

    for s in ntdivs(d):
        res = decomp_dense_s(f, s)
        if res is not None:
            return res

def vss(f, R):
    """Computes the size of the value set of f in the given (finite) ring."""
    return len(set(f(a) for a in R))

def vss_est(f, R, k):
    """Estimates vss(f,GF(p)) by taking k samples and using Chao1."""
    sampled = set()
    while len(sampled) < k:
        sampled.add(R.random_element())
    counts = collections.Counter()
    for a in sampled:
        counts[f(a)] += 1
    uniq = len(counts)
    singles = sum(1 for c in counts.values() if c == 1)
    doubles = sum(1 for c in counts.values() if c == 2)
    collisions = sum(c*(c-1)//2 for c in counts.values())
    chao1 = uniq + singles*singles / (2 * doubles)
    gt = uniq / (1 - singles / len(sampled))
    bday = k*(k-1) / (2 * collisions)
    return round(chao1,1), round(gt,1), round(bday,1)

def rand_poly(PR, d, t):
    """Generates a random monic polynomial with degree t and sparsity t."""
    F = PR.base_ring()
    coeffs = {d: F(1)}
    for exp in random.sample(range(d), min(t-1,d)):
        while True:
            c = F.random_element()
            if c: break
        coeffs[exp] = c
    return PR(coeffs)

def rand_comp(PR, d, t):
    """Generates two random polynomials, each with sparsity t, to give a composed polynomial of degree d."""
    try:
        r = random.choice(ntdivs(d))
    except IndexError:
        raise ValueError(f"Can't compose to prime degree: {d}") from None
    s = d // r
    g = rand_poly(PR, r, min(r+1,t))
    h = rand_poly(PR, s, min(s+1,t))
    return g(h), g, h

def test_dense_decomp():
    """Performs a bunch of random decomposable / indecomposable tests."""
    global f1, g1, h1
    for _ in range(1000):
        F = GF(random_prime(1000,lbound=500))
        PR = PolynomialRing(F,'x')
        d = random.randrange(4, 500)
        f0 = rand_poly(PR, d, random.randrange(3,d))
        dd = decomp_dense(f0)
        if dd is not None:
            g0,h0 = dd
            assert f0 == g0(h0)
            print('--------------------------------------------')
            print("found a random composition:")
            print('f', f0)
            print('g', g0)
            print('h', h0)
            print('--------------------------------------------')
        while True:
            d = random.randrange(4,500)
            if not is_prime(d): break
        f1, g1, h1 = rand_comp(PR, d, random.randrange(3,d))
        g2,h2 = decomp_dense(f1)
        assert f1 == g2(h2)


if __name__ == '__main__':
    seed = random.randrange(10**10)
    print("USING SEED:", seed)
    random.seed(seed)
    set_random_seed(seed)
    #test_dense_decomp()

    p = next_prime(10000)
    F = GF(p)
    PR = PolynomialRing(F,'x')
    d = 5000
    t = 10
    k = 500
    print("indecomp:")
    for _ in range(10):
        while True:
            f = rand_poly(PR, d, t)
            if decomp_dense(f) is None: break
        print(vss(f, F), *vss_est(f, F, k), sep='\t')
    print("decomp:")
    for _ in range(10):
        f,g,h = rand_comp(PR, d, t)
        print(vss(f, F), *vss_est(f, F, k), sep='\t')
