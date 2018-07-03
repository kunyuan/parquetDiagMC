#!/usr/bin/python
import pprint
import gzip,os,sys,time
print sys.version
import hickle as hkl
from numpy import *
#all numpy symbols have to be imported as * in order to read "array([...])" in .txt file with LoadDict function

set_printoptions(threshold=nan) #make sure numpy will print all elements, so that SaveDict and LoadDict will work even for very large array

def SaveDict(filename, mode, root):
    if filename[-4:]!=".txt":
        filename+=".txt"
    with open(filename, mode) as f:
        f.write(pprint.pformat(root))

def LoadDict(filename):
    if filename[-4:]!=".txt":
        filename+=".txt"
    with open(filename, "r") as f:
        return eval(f.read())

def SaveBigDict(filename, root):
    if filename[-4:]!=".hkl":
        filename+=".hkl"
    if "GammaWStatis" in root.keys():
        gammaw=root["GammaWStatis"]
        gammaw["WeightAccu"]=array(gammaw["WeightAccu"], dtype=complex64)
    hkl.dump(root, "_"+filename, mode='w', compression='gzip')
    os.rename("_"+filename, filename)

def LoadBigDict(filename):
    if filename[-4:]!=".hkl":
        filename+=".hkl"
    data=hkl.load(filename)
    if "WWGammaW" in data.keys():
        gammaw=data["WWGammaW"]
        if "SmoothT" in gammaw.keys():
            gammaw["SmoothT"]=array(gammaw["SmoothT"], dtype=complex)
    return data
