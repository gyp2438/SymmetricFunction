from sage.combinat.q_analogues import *

q_bin=q_binomial

def q_int(n):
	return add(q**i for i in range(n))

def q_cat(a,b=None):
    if b==None:
        return q_cat(a,a+1)
    return 1/q_int(a+b)*q_binomial(a+b,a)


Phivar=var(['varphi_%s'%i for i in range(100)])
Qvar=var(['q_%s'%i for i in range(100)])
Factvar=var(['f_%s'%i for i in range(100)])

@cached_function
def sk_sj(k,j):
    return q**(k*j)+add(c*q**mu[0] for mu,c in s[k](s[j]).restrict_partition_lengths(2))

def to_q_fact(pol):
    expr=to_qn(pol)
    expr=expr.substitute({Qvar[i]:Factvar[i]/Factvar[i-1] for i in range(2,100)})
    return expr.substitute({Factvar[1]:1})


def to_qn(pol):
    if pol==0 or pol==1:
        return pol
    else:
        expr=qcyclo_decomp(pol)
        if expr.parent()==factor(1).parent():
            return expr
        else:
            return expr.substitute({Phivar[i]:phi_to_qn(i) for i in range(1,100)})

@cached_function
def inv_cyclo(d):
    return [n for n in range(1,1260) if d==euler_phi(n)]

@cached_function
def qcyclo_decomp(pol):
    if pol==0:
        return 0
    else:
        fpol=factor(pol)
        return fpol.unit()*factor(mul(to_qcyclo(c[0])**c[1] for c in list(fpol)))
    
    
@cached_function
def to_qcyclo(pol):
    pol=SR(pol)
    d=pol.degree(SR(q))
    for n in inv_cyclo(d):
        if cyclotomic_polynomial(n,SR(q))==pol:
            return Phivar[n]
    return pol

def phi_to_qn(k):
    if Integer(k).is_prime():
        return Qvar[k]
    else:
        return Qvar[k]/mul(phi_to_qn(d) for d in divisors(k) if not d==k and not d==1)  


def phi_en_qn(n):
    return mul(Qvar[d]**moebius(n/d) for d in divisors(n))


def CC(k,j): return SR(k).binomial(j, hold=True)

def tobinom(pol):
    pol=numerator(pol)
    return add(pol.coefficient({q:k})*binqn(k) for k in range(pol.degree()+1))

def ToBin(pol):
    pol=numerator(pol)
    return add(a*CC('k',b.degree()) for a,b in pol)

def binqn(n):
    return add(stirling_number2(n,k)*factorial(k)*q**k for k in range(n+1))


def qt_cyclo(n):
    if n==1:
        return 1
    return t**euler_phi(n)*cyclotomic_polynomial(n,SR(q)).substitute({q:q/t})

@cached_function
def s_cyclo(n):
    if n==1:
        return 1
    return InSchur(qt_cyclo(n))

def qt_binomial(n,k):
    return t**((n-k)*k)*q_binomial(n,k).substitute({q:q/t})
