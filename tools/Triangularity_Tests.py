@cached_function
def TriangularPartitions(n,Method='Safarul'):
    if n<4:
        return sorted([Partition(list(tau)) for tau in Partitions(n)])
    else:
        return sorted(Set([Partition(list(mu.add_cell(c[0])))
                       for mu in TriangularPartitions(n-1)
                       for c in mu.addable_cells()
                       if Is_Triangular(mu.add_cell(c[0]),Method=Method)]))

@cached_function
def Is_Triangular(tau,Method='Safarul'):
    tau=Partition(tau)
    if tau.size()==0:
        return True
    if Method=='Sim':
        return t_min(tau)<t_max(tau)
    elif Method=='Sturm_Factor':
        return part_to_word(tau).is_sturmian_factor()
    elif Method=='ConcaveConvex':
        return Is_Convex(tau) and Is_Concave(tau)
    elif Method=='Safarul':
        return Safarul(tau)
    else:
        return False


@cached_function
def part_to_word(mu):
    var('a','b')
    if mu.length()==0:
        return Word([])
    else:
        L=[0]+list(reversed(mu))
        DL=[L[i]-L[i-1] for i in range(1,len(L))]
        return mul(Word([a for i in range(DL[k])]+[b]) for k in range(len(DL)))

@cached_function
def Is_Concave(tau):
    tau=Partition(tau)
    if tau.size()<=3:
        return True
    else:
        tau=Partition(tau)
        cell_tau=Set([tuple(c) for c in tau.cells()])
        rho=Partition([max(tau)+1 for i in range(tau.length()+1)])
        cell_rho=Set([tuple(c) for c in rho.cells()])
        P_tau=Polyhedron(vertices = cell_rho.difference(cell_tau))
        cell_P_tau=Set([tuple(c) for c in P_tau.integral_points()])
        return cell_tau==cell_rho.difference(cell_P_tau)

@cached_function
def Is_Convex(mu):
    mu=Partition(mu)
    cell_mu=Set([tuple(c) for c in mu.cells()])
    P_mu=Polyhedron(vertices = cell_mu)
    cell_P_mu=Set([tuple(c) for c in P_mu.integral_points()])
    return cell_mu==cell_P_mu

def t_min(tau):
    tau=Partition(tau)
    if tau.is_empty():
        return 0
    else:
        return QQ(Max([low_t(c,tau) for c in Cellules(tau)]))

def t_max(tau):
    tau=Partition(tau)
    if tau.is_empty():
        return 1
    else:
        return QQ(Min([top_t(c,tau) for c in Cellules(tau)]))

def top_t(c,alpha):
    leg=alpha.leg_length(c[1]-1,c[0]-1)
    arm=alpha.arm_length(c[1]-1,c[0]-1)
    return (leg+1)/(arm+leg+1)

def low_t(c,alpha):
    leg=alpha.leg_length(c[1]-1,c[0]-1)
    arm=alpha.arm_length(c[1]-1,c[0]-1)
    return leg/(arm+leg+1)

def Cellules(mu):
    return [(b+1,a+1) for a in range(mu.length()) for b in range(mu[a])]


@cached_function
def Safarul(tau): #checks triangularity using the second sturmian characterization described in the source
    tau = Partition(tau)
    if tau.length() > tau[0]:
        tau = tau.conjugate()
    p = tau.length()
    S = S_tau(tau)
    m = min(S)    
    if S not in [{m},{m,m+1}] or tau[-1]>m+1:
        return false    
    s = chi_2(tau)
    w = chi_3(tau)
    a = Word(w)
    if w!=[1]*p and (w[0]==1 and s!=0 and tau[-1] == s+1) or (w[0]==0 and ((w==[0]*p and tau[-1] == s) or (w != [0]*p and tau[-1] <= s))):
        if Partition([tau[-1]+(p-i)*s+sum(w[j] for j in range(1,p-i+1)) for i in range(1,p+1)]) == tau:                                 
            return (a.is_balanced() or (a[0]==0 and a[1::].is_balanced()))
    else:
        return false

def S_tau(tau): #set of steps of tau
    tau = Partition(tau)
    p = tau.length()
    if p == 1:
        return Set([tau[0]])
    return Set([tau[p-i]-tau[p-i+1] for i in range(2,p+1)])

def chi_2(tau): #function chi_2 defined in the source
    tau = Partition(tau)
    m = min(S_tau(tau))
    if S_tau(tau).cardinality() == 2 or tau[-1] >= m:
        return m
    else:
        return m-1

def chi_3(tau): #function chi_3 defined in the source
    tau = Partition(tau)
    p = tau.length()
    s = chi_2(tau)
    if tau[-1] == s+1:
        w = [1]
    else:
        w = [0]
    for i in range(2,p+1):
        if tau[p-i]-tau[p-i+1] == s+1:
            w.append(1)
        else:
            w.append(0)
    return w