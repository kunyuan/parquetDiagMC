#!/usr/bin/env python
import numpy as np
# import unionfind
import random
# from nullspace import rank, nullspace
from numpy.linalg import matrix_rank
import itertools 
import sys,os

# permutation meaning:
# l_in ---> Per[0]
# r_in ---> Per[0]
# l_out=-1, r_out=-2

Order=2
I, S, C, R, L, U, Gamma={}, {}, {}, {}, {}, {}, {}
I[1]=1
I[2]=0
I[3]=0
I[4]=32
I[5]=944
I[6]=22608
S[1], C[1], R[1], L[1], U[1]=0, 0, 0, 0, 0 
Gamma[1]=I[1]

LOUT=-1
ROUT=-2
LIN=0
RIN=1

def Iterate(Type, order):
    for l in range(1, order):
        r=order-l
        Operation(Type, l, r)

def Operation(Type, orderL, orderR):
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
        Type[order]=0

    Type[order]+=Left*Gam
    print Type[order]

def SetGamma(order):
    Gamma[order]=I[order]+S[order]+C[order]+R[order]+L[order]+U[order]
    Test(Gamma[order])
    print "Gamma at Order {0}: {1}".format(order, Gamma[order])

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
    SelfEnergy=TotalGamma*2
    return SelfEnergy


def Test(Collection):
    return
    # if len(Collection)!=len(set(Collection)):
        # print "Elements are not unique!\n", Collection
        # sys.exit(0)

if __name__=="__main__":
    ObjectList=[S, C, R, L, U]
    for order in range(2, Order+1):
        for O in ObjectList:
            Iterate(O, order)
        SetGamma(order)
        Test(Gamma[order])

    # TotalGamma=GetFullGamma(Order)
    SelfEnergy=GetSelfEnergy(Order+1, Gamma[Order])
    Test(SelfEnergy)
    print "Self Energy:", SelfEnergy
    print "Total Self Energy diagram number", SelfEnergy

