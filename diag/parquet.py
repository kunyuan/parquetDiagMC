#!/usr/bin/env python
import numpy as np
# import unionfind
import random
# from nullspace import rank, nullspace
from numpy.linalg import matrix_rank
import itertools 
import sys,os
import unittest
from nullspace import rank, nullspace

# permutation meaning:
# l_in ---> Per[0]
# r_in ---> Per[0]
# l_out=-1, r_out=-2

Order=2
I, S, C, R, L, U, Gamma={}, {}, {}, {}, {}, {}, {}

for i in range(1,Order+1):
    I[i]=[]

for i in range(Order):
    I[1].append((i*2,i*2+1,-1,-2))

S[1], C[1], R[1], L[1], U[1]=[], [], [], [], []
Gamma[1]=I[1]

LOUT=-1
ROUT=-2
LIN=0
RIN=1

VerDict={}
VerMap={}
VerIndex=0

def GetVer(diag):
    return (diag[0], diag[1], diag.index(-1)+diag[0]-2, diag.index(-2)+diag[0]-2)

for e in I[1]:
    VerDict[(1, 0, GetVer(e))]=[]
    VerMap[(1,0,GetVer(e))]=(1,0,VerIndex)
    VerIndex+=1

def FindIndex(diag, order):
    if diag in I[order]:
        return 1
    if diag in S[order]:
        return 2
    if diag in C[order]:
        return 3
    if diag in R[order]:
        return 4
    if diag in L[order]:
        return 5
    if diag in U[order]:
        return 6

def Iterate(Type, order):
    for l in range(1, order):
        r=order-l
        Operation(Type, l, r)

def Operation(Type, orderL, orderR):
    global VerIndex

    if Type is S:
        Left=I[orderL]+C[orderL]+R[orderL]+L[orderL]+U[orderL]
    elif Type is C or Type is R:
        Left=I[orderL]+S[orderL]+U[orderL]
    elif Type is L or Type is U:
        Left=I[orderL]+S[orderL]+C[orderL]+R[orderL]+L[orderL]
    # Gamma=I[orderR]+S[orderR]+C[orderR]+R[orderR]+L[orderR]+U[orderR]
    Gam=Gamma[orderR]
    order=orderL+orderR
    if not Type.has_key(order):
        Type[order]=[]

    typeindex=0
    for ldiag in Left:
        for gamma in Gam:
            if max(ldiag)+1==gamma[0]:
                new=list(ldiag)
                new_l_out=new.index(LOUT)
                new_r_out=new.index(ROUT)
                rhalf=list(gamma[2:])
                rhalf_l_out=rhalf.index(LOUT)
                rhalf_r_out=rhalf.index(ROUT)
                if Type is S:
                    new[new_l_out]=gamma[LIN]
                    new[new_r_out]=gamma[RIN]
                    new+=rhalf
                    typeindex=2
                elif Type is C:
                    new[RIN]=gamma[RIN] #r_in of new is r_in of the gamma
                    new[new_r_out]=gamma[LIN] #r_out of ldiag is redirected to gamma l_in
                    rhalf[rhalf_l_out]=ldiag[RIN] #l_out of gamma is to ldiag r_in
                    new+=rhalf
                    typeindex=3
                elif Type is R:
                    new[RIN]=gamma[RIN]
                    new[new_r_out]=gamma[LIN]
                    rhalf[rhalf_r_out]=ldiag[RIN]
                    rhalf[rhalf_l_out]=ROUT
                    new+=rhalf
                    typeindex=4
                elif Type is L:
                    new[RIN]=gamma[RIN]
                    new[new_l_out]=gamma[LIN]
                    rhalf[rhalf_l_out]=ldiag[RIN]
                    new[new_r_out]=LOUT
                    new+=rhalf
                    typeindex=5
                elif Type is U:
                    new[LIN]=ldiag[LIN]
                    new[RIN]=gamma[RIN]
                    new[new_l_out]=gamma[LIN]
                    new[new_r_out]=ROUT
                    rhalf[rhalf_l_out]=LOUT
                    rhalf[rhalf_r_out]=ldiag[RIN]
                    new+=rhalf
                    typeindex=6
                    
                Type[order].append(tuple(new))

                NewVer=(typeindex, order-1, GetVer(new))
                if not VerDict.has_key(NewVer):
                    VerDict[NewVer]=[]
                    VerMap[NewVer]=(typeindex, order-1, VerIndex)
                    VerIndex+=1

                LVer=(FindIndex(ldiag, orderL), orderL-1, GetVer(ldiag))
                RVer=(FindIndex(gamma, orderR), orderR-1, GetVer(gamma))
                if (LVer, RVer) not in VerDict[NewVer]:
                    VerDict[NewVer].append((LVer, RVer))

                print LVer, RVer, "tooooooooooo", NewVer
                # Index+=1
    # print Type[order]

def SetGamma(order):
    Gamma[order]=I[order]+S[order]+C[order]+R[order]+L[order]+U[order]
    Test(Gamma[order])
    # print "Gamma at Order {0}: {1}".format(order, Gamma[order])

    # print "The order ", order
    Equal=[]
    for e in L[order]:
        # print "IN, OUT", e[LIN], e[RIN], e.index(LOUT), e.index(ROUT)
        Equal.append((e[LIN], e[RIN], e.index(LOUT), e.index(ROUT)))
    # print "Size", len(Equal), len(set(Equal))

def GetFullGamma(order):
    Exchange=[]
    for d in Gamma[order]:
        diag=list(d)
        l_out=diag.index(LOUT)
        r_out=diag.index(ROUT)
        diag[l_out], diag[r_out]=diag[r_out], diag[l_out]
        Exchange.append(tuple(diag))
    Total=Gamma[order]+Exchange
    Test(Total)
    return Total

def GetSelfEnergy(order, TotalGamma):
    SelfEnergy=[]
    l_ver=2*(order-1)
    r_ver=2*(order-1)+1
    for diag in TotalGamma:
        new1=list(diag)
        l_out=new1.index(LOUT)
        r_out=new1.index(ROUT)
        new1[l_out]=l_ver
        new1[r_out]=r_ver
        new1.append(-1) #add vertex 2*order
        new1.append(new1[RIN]) #add vertex 2*order+1
        del new1[1]
        SelfEnergy.append(tuple(new1))

        new2=list(diag)
        l_out=new2.index(LOUT)
        r_out=new2.index(ROUT)
        new2[l_out]=l_ver
        new2[r_out]=r_ver
        new2.append(new2[RIN]) #add vertex 2*order
        new2.append(-1) #add vertex 2*order+1
        del new2[1]

        SelfEnergy.append(tuple(new2))
    return SelfEnergy

def swap(array, i, j):
    array = list(array)
    array[i], array[j] = array[j], array[i]
    return tuple(array)

def swap_interaction(permutation, m, n, k, l):
    permutation = list(permutation)
    mp,np,kp,lp=(permutation.index(e) for e in (m,n,k,l))
    permutation[mp]=k
    permutation[kp]=m
    permutation[np]=l
    permutation[lp]=n
    permutation[m],permutation[k]=permutation[k],permutation[m]
    permutation[n],permutation[l]=permutation[l],permutation[n]
    return tuple(permutation)

def swap_LR(permutation, i, j):
    # print permutation, i, j
    permutation = list(permutation)
    ip,jp=permutation.index(i),permutation.index(j)
    permutation[ip]=j
    permutation[jp]=i
    permutation[i],permutation[j]=permutation[j],permutation[i]
    # print "after", permutation
    return tuple(permutation)


def Test(Collection):
    if len(Collection)!=len(set(Collection)):
        print "Elements are not unique!\n", Collection
        sys.exit(0)

def TestConnection(BuildTree, VerTable, Index2Ver):
    order=VerTable.shape[0]
    VerNum=BuildTree.shape[0]
    BranchNum=BuildTree.shape[1]
    WeightTable=np.zeros(VerNum)

    #initialize all first order vertex with 1
    OperNum=0
    for i in range(VerNum):
        index=VerTable[0,i]
        if index<0:
            break
        WeightTable[i]=1
        OperNum+=1

    for o in range(1, order):
        # print "Order {0}".format(o+1)
        Weight=0
        for i in range(VerNum):
            index=VerTable[o, i]
            if index<0:
                break
            # print "{0} with {1}".format(index, Index2Ver[index])
            WeightTable[index]=0.0
            for b in range(BranchNum):
                lver=BuildTree[index, b, 1] 
                rver=BuildTree[index, b, 2] 
                if lver<0 or rver<0:
                    break
                # print "left {0}, right {1}: {2}x{3}".format(lver, rver, WeightTable[lver], WeightTable[rver])
                WeightTable[index]+=WeightTable[lver]*WeightTable[rver]
                OperNum+=1
            Weight+=WeightTable[index]
        print "Order {0} total Weight: {1}\n".format(o+1, Weight)

    print "Operation Num: {0}".format(OperNum)
    # print WeightTable

def TestSelfEnergyDiag(SelfEnergy):
    Diag={}
    for d in SelfEnergy: 
        dd=list(d)
        dd[d.index(-1)]=d[0]
        dd=dd[1:]
        Diag[tuple(dd)]=None

    PermutationDict=dict(Diag)
    for permutation in Diag.keys():
        del PermutationDict[permutation]
        order = len(permutation)/2
        measure_in=permutation[0]

        Deformation = [permutation]
        for idx in range(1, Order):
            if idx==measure_in:
                continue
            for i in range(len(Deformation)):
                for j in range(idx):
                    if j==measure_in:
                        continue
                    Deformation.append(swap_interaction(Deformation[i], idx*2, idx*2+1, j*2, j*2+1))

        for idx in range(Order):
            if idx==measure_in:
                continue
            for i in range(len(Deformation)):
                Deformation.append(swap_LR(Deformation[i], idx*2, idx*2+1))

        Deformation = set(Deformation)
        for p in Deformation:
            if p in PermutationDict:
                print "Duplicate diagram is found! {0}, {1}".format(permutation, p)
                del PermutationDict[p]
                sys.exit(0)

        kG, kW=AssignMomentums(permutation)

        for i in range(0,len(kG)):
            if abs(kG[i])<1e-10:
                print "Tadpole detected!", permutation
                sys.exit(0)

            for j in range(i+1,len(kG)):
                if abs(kG[i]-kG[j])<1e-10:
                    print "Same k on G for {0}: {1} on {2}; {3} on {4}".format(permutation, kG[i],i,kG[j],j)
                    # print "Same k on W for {0}: {1}; 1, {2}".format(p, kG[i],kG[j])
                    sys.exit(0)

        #check reducibility
    return

def AssignMomentums(permutation):
    N=len(permutation)/2
    vectors=FindIndependentK(permutation)
    freedoms=vectors.shape[1]
    karray=np.array([random.random() for _ in range(freedoms)])
    kVector=np.dot(vectors, karray)
    # kVector=vectors[:,0]
    return kVector[:2*N], kVector[2*N:]

def FindIndependentK(permutation):
    # kList=[(random.randint(0, Nmax), random.randint(0,Nmax)) for i in range(len(InteractionPairs)+1)]
    N=len(permutation)/2
    Matrix=np.zeros((2*N,3*N))
    for i in range(2*N):
        interaction=int(i/2)+2*N
        sign=i%2
        Matrix[i,interaction]=-(-1)**sign
        Matrix[i, i]=-1
        Matrix[i, permutation.index(i)]=1
    # print Matrix
    vectors = nullspace(Matrix)
    # print len(vectors)
    # print vectors
    freedoms=vectors.shape[1]
    if freedoms!=N+1:
        print "Warning! Rank is wrong for {0} with \n{1}".format(permutation, vectors)
    return vectors

def Compress(VerDict):
    VerDictSimple={}
    VerNum=len(VerDict.keys())
    BranchNum=max([len(e) for e in VerDict.values()])
    print BranchNum
    Index2Ver=np.zeros((VerNum, 4), dtype=int) #l_in, r_in, l_out, r_out
    BuildTree=np.zeros((VerNum, BranchNum, 7), dtype=int) #VerIndex, BranchIndex, Type/LeftVer/RightVer/lver_l_out/lver_r_out/rver_l_in/rver_r_in
    BuildTree-=1
    VerTable=np.zeros((Order, VerNum), dtype=int) #for a given order, what are the vertex
    VerTable-=1

    # for k in sorted(VerDict.keys(), key = lambda x: x[1]):
    for k in VerMap.keys():
        index=VerMap[k][-1]
        Index2Ver[index]=k[-1]

    for k in VerDict.keys():
        VerDictSimple[VerMap[k]]=[]
        index=VerMap[k][-1]
        bran=0
        for elem in VerDict[k]:
            left=elem[0]
            right=elem[1]
            VerDictSimple[VerMap[k]].append((VerMap[left], VerMap[right]))
            # print left
            BuildTree[index, bran, :]=[k[1], VerMap[left][-1], VerMap[right][-1], left[-1][-2], left[-1][-1], right[-1][0], right[-1][1]] 
            # print [k[1], VerMap[left][-1], VerMap[right][-1], left[-1][-2], left[-1][-1], right[-1][0], right[-1][1]], left, right
            bran+=1

    for o in range(Order):
        ver=[e for e in VerMap.values() if e[1]==o]
        for i in range(len(ver)):
            VerTable[o, i]=ver[i][-1]

    return VerDictSimple, Index2Ver, BuildTree, VerTable


if __name__=="__main__":
    ObjectList=[S, C, R, L, U]
    for order in range(2, Order+1):
        for O in ObjectList:
            Iterate(O, order)
        SetGamma(order)
        Test(Gamma[order])

    TotalGamma=GetFullGamma(Order)
    SelfEnergy=GetSelfEnergy(Order+1, Gamma[Order])
    for i in range(1,Order+1):
        print "Gamma diagram number at order {0}:{1}".format(i, len(Gamma[i])) 
    Test(SelfEnergy)
    TestSelfEnergyDiag(SelfEnergy)
    # print "Self Energy:", SelfEnergy
    print "Total Self Energy diagram number", len(SelfEnergy)

    for k in VerDict.keys():
        print "{0} : {1}".format(k, VerDict[k])
        pass
    print "Total Vertex ", sum([len(e) for e in VerDict.values()])

    VerDictSimple, Index2Ver, BuildTree, VerTable=Compress(VerDict)

    for k in VerDictSimple.keys():
        # print "{0} : {1}".format(k, VerDictSimple[k])
        pass
    print "Total terms ", sum([len(e) for e in VerDictSimple.values()])
    # print VerMap
    print len(VerMap.keys())
    print Index2Ver
    print BuildTree[4,...]
    print VerTable
    np.savez("table.npz", buildtree=BuildTree, vertable=VerTable)
    TestConnection(BuildTree, VerTable, Index2Ver)

    

