import numpy as np
import time
from numpy.random import randint,seed
from fractions import Fraction
from itertools import combinations, product
from  scipy.special import comb
from graph import Graph


def PolyArea(x,y):
    return Fraction(1,2)*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))

def PolyArea3(x,y):
    S2 = abs(x[0]*y[1] + x[1]*y[2] + x[2]*y[0] - x[1]*y[0] - x[2]*y[1] - x[0]*y[2])
    assert(S2>=0)
    return Fraction(1,2)*S2


def TreugArea(Tr):
    S = PolyArea3(Tr[:,0],Tr[:,1])
    assert(S>=0)
    return S

def random_points(N, **kwargs):
    size = kwargs.get('size', 2**30)
    a = randint(-size,size+1, size= [N,2])
    b = randint(-size,size+1, size= [N,2])
    b[b==0] = 1
    points = np.empty((N,2), dtype = object)
    for i in range(N):
        for j in range(2):
            points[i,j] = Fraction(int(a[i,j]),int(b[i,j]))
    return points   

def erased_code(Points):
    N = Points.shape[0]
    S_rel=[]
    for tr_pair in combinations(combinations(range(N), 3), 2):
        Tr  = [ Points[tr_pair[0],:], Points[tr_pair[1],:]] # triangle coordinates
        S = [TreugArea(tr) for tr in Tr]
        Sr = max(S)/min(S)
        S_rel.append(Sr)
    return np.sort(S_rel)

def allS(Points):
    # get all triangles areas
    N=Points.shape[0]
    S = [TreugArea(Points[tr,:]) for tr in combinations(range(N), 3)]
    return np.sort(S)

def find_set_areas(CodeE, way, N):
    '''
    Determines original areas by erased code
    CodeE - Erased Code
    way = 0 or 1 - branch
    N - number of points
    '''
    CodeE=np.sort(CodeE)
    Smax = CodeE[-1]
    M=CodeE.shape[0]
    CodeE_rev = Smax / CodeE[::-1]
    assert(np.all(CodeE_rev[:-1] < CodeE_rev[1:])) # check that sorted too
    ##% find pairs (TODO: separate function)
    i=j=0
    pairs0=[]
    while i < M and j < M:
        if CodeE[i] == CodeE_rev[j]:
            if CodeE[i]*CodeE[M-1-j] == Smax:
                ind_pair = np.sort([i, M-1-j])
                pairs0.append(ind_pair)
                #print(ind_pair)
            i+=1
            j+=1
        else:
            if CodeE[i] < CodeE_rev[j]:
                i+=1
            else:
                j+=1
    pairs1 = np.vstack(pairs0)
    pairs = np.unique(pairs1, axis=0);
    assert(pairs.shape[0] == comb(N,3)-2); #% - 2, because Smin,Smax excluded
    #% then test for Sbase - either CodeEnorm(i) or CodeEnorm(j)
    Sbase=CodeE[pairs[0][way]] #way must be 0 or 1
    assert(Sbase>1)
    Sset = [1, Sbase, Smax]
    for p in pairs[1:]:
        i , j =  p
        S_i = max(CodeE[i]/Sbase, Sbase/CodeE[i])
        S_j = max(CodeE[j]/Sbase, Sbase/CodeE[j])
        b_i = np.any(S_i == CodeE)
        b_j = np.any(S_j == CodeE)
        if b_i+b_j != 1:#if fail - not general position triangles
            return None # recovery failed
        if b_i:
            Sset.append(CodeE[i])
        else:
            Sset.append(CodeE[j])
    assert(len(Sset) == comb(N,3))
    return np.sort(Sset)

def find_triples(Sset, N):
    '''
    function [Triples, TripleTypes] = find_triples(Sset, Npt)
    N - number of points
    '''
    n = len(Sset); # % number of triangles, should be C_N^3
    assert n == comb(N,3)
    Triples =[]; TripleTypes=[]
    Sset = np.sort(Sset)
    assert Sset[0] == 1
    for main in range(1,n): # 2:n % main triangle
        Smain = Sset[main]
        for s0,s1 in combinations(range(n), 2):
            if s0==main or s1==main: 
                continue
            #% type 1 Main = S1+S2+1, tyhpe 2 Main = S1+S2-1
            S12 = Sset[s0] + Sset[s1]
            if Smain == S12+1  or Smain ==  S12-1:
                Triples.append([main, s0, s1])
                TripleTypes.append(1 if Smain == S12+1 else 2)
    #% each triple correspond to point (exclude central triangle)
    assert len(Triples) == N-3, 'Not general position triangles!'
    assert len(TripleTypes) == N-3
    return Triples, TripleTypes

def have_common_side(i1,i2, Sset):
    '''
    Checks if Sset(i1) and Sset(i2) have a common side
    First sign always "+", others - 7 combinations
    '''
    r1 = Fraction(1,1)
    SignsCombine = np.array([
        [r1,  r1,  r1, -r1],
        [r1,  r1, -r1,  r1],
        [r1, -r1,  r1,  r1], # % 1"-"
        [r1,  r1, -r1, -r1],
        [r1, -r1, -r1,  r1],
        [r1, -r1,  r1, -r1], #% 2"-"
        [r1, -r1, -r1, -r1]# % 3"-"
        ], dtype=object)
    nTr = len(Sset)
    ind = [i for i in range(1,nTr) if (i!=i1 and i!=i2)]
    S = np.array([ Sset[i1], Sset[i2], 0, 0 ])
    for pair in combinations(ind, 2):
        S[2] = Sset[pair[0]]
        S[3] = Sset[pair[1]]
        if any(S@SignsCombine.T==0):
            return True
    return False



def have_common_side_opt(i1,i2, Sset):
    '''
    Checks if Sset(i1) and Sset(i2) have a common side
    First sign always "+", others - 7 combinations
    TODO: precompute all sums and diffs in Sset
    use "from bisect import bisect_left"
    '''
    ind = [i for i in range(1,len(Sset)) if (i!=i1 and i!=i2)]
    S = np.array([ Sset[i1], Sset[i2], 0, 0 ])
    S12sum, S12diff  = Sset[i1] + Sset[i2], Sset[i1] - Sset[i2]
    for pair in combinations(ind, 2):
        S[2], S[3] = Sset[pair[0]], Sset[pair[1]]
        S23sum, S23diff = S[2] + S[3],  S[2] - S[3] 
        if S12sum in [S23diff, -S23diff, S23sum]:
            return True
        if S12diff in [S23diff, -S23diff,S23sum, -S23sum]:
            return True
    return False



def factorize_adj_class(Sset, Triples):
    n_Triples = len(Triples)
    pos_max = len(Sset)
    conn = np.zeros((pos_max,pos_max))
    #% set -1 for pos if two triangles in the same triple
    for i in range(n_Triples): 
        for j1,j2 in product(range(3), range(3)):
            if j1==j2: continue
            conn[Triples[i][j1], Triples[i][j2]] = -1

    TriplesInd  = np.unique(np.array(Triples).ravel())

    for i, j  in product(TriplesInd, TriplesInd):
        if conn[i,j]!=-1 and have_common_side_opt(i, j, Sset):
            conn[i,j] = conn[j,i] = 1
    conn_red1 = conn[:,TriplesInd]
    conn_red = conn_red1[TriplesInd, :]
    n_triples = len(TriplesInd)
    G = Graph(n_triples)
    for i,j in combinations(range(n_triples), 2):
        if conn_red[i,j]==1:
            G.addEdge(i,j)
    cc = G.connectedComponents()
    ClassIndExt = np.zeros((pos_max,))
    for i, comp in enumerate(cc):
        ClassIndExt[TriplesInd[comp]] = i+1 #starting from 1
    #%% check that there are 3 types in each triple
    for tr in Triples:
        tri_bin = ClassIndExt[tr]
        #print(tri_bin)
        assert  all(np.sort(tri_bin) == [1,2,3]), 'Not 1,2,3 class in triple'
    return ClassIndExt

def  type_class2sign(triangle_type, triangle_class):
    if triangle_class==0:
        return (0, 0)
    elif triangle_class==1:
        return (-1,1) if triangle_type == 1 else (1,-1)
    elif triangle_class==2:
        return (-1,-1) if triangle_type == 1 else (1,1)
    else:
        return (1,-1) if triangle_type == 1 else (-1,1)

def getSignFromClass(Triples, TripleTypes, ClassInd):
    n_triples = len(Triples)
    Triples2D = np.array(Triples)
    assert len(TripleTypes)==n_triples
    n_triangle = len(ClassInd)
    assert np.max(Triples2D)<=n_triangle
    ClsInd =  ClassInd[Triples2D[:,0]]
    res = [type_class2sign(x,y) for x,y in zip(TripleTypes, ClsInd)]
    return res

def run_experiment(N):
    seed(42)
    Points = random_points(N, size=2**24)
    S = allS(Points)

    S = S / S[0] # normalize it
    CodeE = erased_code(Points)
    print(S.shape, CodeE.shape)
    for way in [1,0]:
        Sset = find_set_areas(CodeE, way, N)
        if np.all(S==Sset):
            print('Recover original')
        else:
            print('Recover Ghost')
            continue
        Triples, TripleTypes =  find_triples(Sset, N)
        ClassInd = factorize_adj_class(Sset, Triples)    
        #print(ClassInd)
        ClassMap = {x:y for x,y in zip(Sset, ClassInd)}#% TODO: check syntax!
        xy_sign = getSignFromClass(Triples, TripleTypes, ClassInd)
        #% Let the initial triangle has coordinates (0,0)-(0,1)-(1,0)
        #% then S'=S/2 (as it has area=0.5)
        #% Suppose S1 - prilegaet k (0,0)-(0,1) - triple of class 1.
        #% then geometric set is a line such that 1/2*|y|*1 = 1/2S1 =>
        #% y=S1*y_sign(i), ditto for x = S2*x_sign(i)
        XY_recover = np.zeros((N,2) , dtype=object)
        XY_recover[1,:] = [1,0] 
        XY_recover[2,:] = [0,1]
        Triples2D = np.array(Triples)
        for i in range(3,N):
            #% recover i-th point from (i-3)th triple
            i_t = i-3
            triple = Triples2D[i_t,:]
            x_s, y_s =xy_sign[i_t]
            #y_s = y_sign(i_t);
            assert x_s!=0 and y_s !=0
            for tr in triple:
                cls_tr = ClassMap[Sset[tr]]
                if cls_tr==1:
                    XY_recover[i,1] =  Sset[tr]*y_s
                elif cls_tr==2:
                    XY_recover[i,0] = - Sset[tr]*x_s #;%% TODO: why minus here???
                elif cls_tr==3:
                    pass
                else:
                    raise f'Wrong class ID = {cls_tr}'
        CodeErec = erased_code(XY_recover)
        if np.all(CodeErec == CodeE):
            print('Recover success')
        else:
            print('Recover fail')



if __name__ == "__main__":
    N=8
    t0 = time.time()
    run_experiment(N)
    t1 = time.time()
    print(N, t1-t0)
