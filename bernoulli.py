from sympy import var, simplify, S, bernoulli, binomial, zoo, Rational, factor, factor_list, numer, denom
var('p')

def Vt(x, t=var('p')):
    """ Computes the t-adic valuation of a polynomial x"""
    if x == 0:
        return 99999
    for L in factor_list(x)[1]:
        if L[0] == t:
            return L[1]
    return 0

def val(x, t=var('p')):
    """ Computes the t-adic valuation of a rational function x"""
    if x == 0:
        return 99999
    return Vt(numer(factor(x)), t)-Vt(denom(factor(x)), t)

def v(p, x):
    """Compute p-adic valuation of x.
    Returns 99999 if x is 0"""
    if p <= 1:
        raise ValueError('Invalid prime for computing valuation.')
        return None
    if x == 0:
        return 99999
    if isinstance(x, int):
        x = Rational(x)
    if denom(x) % p == 0:
        return v(p, x * p) - 1
    if numer(x) % p == 0:
        return v(p, x / p) + 1
    return 0

def B_test(LL1,LL2,N):
    assert len(LL2) == len(LL1) + 1
    L1 = [S(x) for x in LL1]
    L2 = [S(x) for x in LL2]
    bb = True
    for P in [11,13,17]:
        x = sum(L2[i+1].subs(p,P)*bernoulli(L1[i].subs(p,P)) for i in range(len(L1)))
        x -= L2[0].subs(p,P)
        if v(P,x) < N:
            bb = False
            break
    return bb

def B(LL1,LL2,N,alert=False):
    assert len(LL2) == len(LL1) + 1
    L1 = []
    L2 = [S(LL2[0])]
    for i in range(len(LL1)):
        a = factor_list(LL1[i])
        if a[1] == []:
            L2[0] -= bernoulli(a[0])*LL2[i+1]
        else:
            L1.append(S(LL1[i]))
            L2.append(S(LL2[i+1]))
    y = [f.subs(p,1) for f in L1]
    MIN = 99999
    for f in L2:
        MIN = min(MIN, val(f))
    SET = []
    for f in L1:
        q = f.subs(p,1)
        if q not in SET:
            SET.append(q)
    T = L2[0]
    for i in range(len(L1)):
        if y[i] == 0:
            T -= (1-1/p)*L2[i+1]
        elif y[i] >= 2 and y[i]%2 == 0:
            T -= (1-p**(y[i]-1))*bernoulli(y[i])/y[i]*L1[i]*L2[i+1]
    if val(T) < N:
        if alert:
            print "First test: %i"%val(T)
            print T
        return False
    for k in SET:
        if k <= 0 and k%2 == 0:
            for m in range(1,N-MIN+1):
                T = S(0)
                for i in range(len(L1)):
                    if y[i] == k:
                        T += L2[i+1]*L1[i]**m
                if val(T) < N+1-m:
                    if alert:
                        print "Second test: %i, %i, %i"%(k,m,val(T))
                        print T
                    return False
    for k in SET:
        if k >= 2 and k%2 == 0:
            for m in range(2,N-MIN+1):
                T = S(0)
                for i in range(len(L1)):
                    if y[i] == k:
                        T += L2[i+1]*(L1[i]**m - k**(m-1)*L1[i])
                if val(T) < N+1-m:
                    if alert:
                        print "Third test: %i, %i, %i"%(k,m,val(T))
                        print T
                    return False
    return True

def tests():
    """Should do nothing if program is working correctly"""
    assert B([p-1,2*p-2],[p-1,2*p,-p],2)
    assert not B([p-1,2*p-2],[p-1,2*p,-p],3)
    for k in range(2,7):
        for b in range(2,8,2):
            assert B([k*(p-1)+b,p-1+b],[-(k-1)*(1-p**(b-1))/b*bernoulli(b),1/(k*(p-1)+b),-k/(p-1+b)],2)
            assert not B([k*(p-1)+b,p-1+b],[-(k-1)*(1-p**(b-1))/b*bernoulli(b),1/(k*(p-1)+b),-k/(p-1+b)],3)
    for k in range(3,10):
        for b in range(2,10,2):
            assert B([k*(p-1)+b,2*(p-1)+b,p-1+b,b],[0,1/(k*(p-1)+b),-binomial(k,2)/(2*(p-1)+b),k*(k-2)/(p-1+b),-binomial(k-1,2)*(1-p**(b-1))/b],3)
            assert not B([k*(p-1)+b,2*(p-1)+b,p-1+b,b],[0,1/(k*(p-1)+b),-binomial(k,2)/(2*(p-1)+b),k*(k-2)/(p-1+b),-binomial(k-1,2)*(1-p**(b-1))/b],4)
    print "Success"