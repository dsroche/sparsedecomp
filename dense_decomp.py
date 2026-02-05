from sage.all import (
    GF,
    PolynomialRing,
    PowerSeriesRing,
    divisors,
    Integer,
    set_random_seed,
    random_prime,
    is_prime,
)
import random
import functools

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
    d = f.degree()
    PR = f.parent()
    F = PR.base_ring()

    for s in ntdivs(d):
        res = decomp_dense_s(f, s)
        if res is not None:
            return res

        """
        r = d // s
        f_monic = f / f.leading_coefficient()

        # Explicitly define S and its generator y
        S = PowerSeriesRing(F, 'y', default_prec=s+1)
        y = S.gen()

        # PowerSeries constructor from list
        coeffs_rev = f_monic.list()[::-1][:s+1]
        P_ser = S(coeffs_rev)

        try:
            h_ser = P_ser.nth_root(r)
            # Reconstruct h as a polynomial
            h = PR(h_ser.list()[::-1])

            g = find_g(f, h)
            if g(h) == f:
                return g, h
        except (ValueError, ArithmeticError):
            continue
            """

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
    for _ in range(1000):
        F = GF(random_prime(1000,lbound=100))
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
    #seed = 1234
    print("USING SEED:", seed)
    random.seed(seed)
    set_random_seed(seed)
    test_dense_decomp()


    F = GF(101)
    PR = PolynomialRing(F,'x')
    x = PR.gen()
    f0 = rand_poly(PR,100,5)
    f1,g1,h1 = rand_comp(PR,100,5)
    #run_experiment(100,10,101)
