import numpy as np
from isap_great3_tools import *
from isap_great3_sparse2d import *
import argparse
import subprocess
import pyfits


class test:
   def __init__(self):
      self.bef = self.boostEnumFindObject



   @property
   def boostEnumFindObject(self, typeValue,enumType):
      """! 
         Get the key associated to a boost python type in a 
         boost python enumeration
         @param typeValue: a boost python type 
         @param enumType: a boost python enumeration
         @return key in the boost enumeration
      """
      Nobj=len(enumType.names.keys())
      return [enumType.values.keys()[k] for k in range(Nobj) \
                                   if enumType.values.values()[k] == typeValue][0]

pylab

image=numpy.zeros((192,96), numpy.float32)+1.
for kx in range(192):
   for ky in range(96):
      image[kx,ky]=kx+ky*192


Data=Iflt(image.shape[0],image.shape[1])
Data.data[:]=image


#image[:]=numpy.random.randn(192,96)
#Data.data[:]=image


#Test Parser
command="-t24 -n4 -T1 -L -v"
#Get transform or filter size
nUndecFilter=len(type_undec_filter.names.keys())
nTransform=len(type_transform.names.keys())
nLift=len(type_lift.names.keys())

parser = argparse.ArgumentParser(prog='mr_transform',description="Process mr_transform options.", add_help=False)
parser.add_argument('-U',dest='UFilter',default=bboostEnumFindObject(type_undec_filter.U_B3SPLINE_2,type_undec_filter),type=int,help="Type of Undecimated Filter. Default is SplineB3-Id:  H=[1,4,6,4,1]/16, Ht=H, G=Id-H*H, Gt=Id",choices=range(nUndecFilter))
parser.add_argument('-u',dest='NbrUndec',default="-1",type=int,help="Number of Undecimated Scales",choices=range(nUndecFilter))
parser.add_argument('-L',dest='L2Norm',action='store_true',default=False,help="Use a L2 normalization. Default is L1")
parser.add_argument('-v',dest='Verbose',action='store_true',default=False,help="Verbose. Default is no.")
parser.add_argument('-t',dest='TypeTransform',default=bboostEnumFindObject(type_transform.TO_PAVE_BSPLINE,type_transform),type=int,help="Type of MultiResolution Transform. Default is B-Spline Wavelet Transform with a trous algorithm.",choices=range(nTransform))
parser.add_argument('-l',dest='TypeLifting',default=boostEnumFindObject(type_lift.TL_INT_HAAR,type_lift),type=int,help="Type of Lifting Transform. Default is Lifting scheme: integer Haar WT.",choices=range(nLift))
parser.add_argument('-T',dest='TypeFilter',default=bboostEnumFindObject(type_sb_filter.F_MALLAT_7_9,type_sb_filter),type=int,help="Type of Filter. default is Biorthogonal 7/9 filters.",choices=range(nLift))
parser.add_argument('-n',dest='NbrPlan',default=4,type=int,help="Number of scales used in the multiresolution transform. Default is 4.",choices=range(2,10))
parser.add_argument('-c',dest='NbrIter',default=0,type=int,help="Iterative transformation. Iter = number of iterations. This option is valid only if the chosen multiresolution transform is pyramidal (6,7,8,9,10,11). The reconstruction is not exact and we need few iterations. Generally, we take 3.",choices=range(2,20))

test=parser.parse_args(command.split())




shapeImage=image.shape
FAS= FilterAnaSynt()
FAS.alloc(type_sb_filter.F_MALLAT_7_9)
MR_Data =MultiResol()
MR_Data.Border=type_border.I_CONT
MR_Data.alloc(image.shape[0], image.shape[1],4,type_transform.TO_MALLAT,FAS,sb_type_norm.NORM_L2,0,type_undec_filter.U_B3SPLINE_2)
MR_Data.transform(Data)
imshow(MR_Data.band(2).data)
print MR_Data.Set_Transform

FAS=FilterAnaSynt()
MR_Data.alloc(image.shape[0], image.shape[1],4, type_transform.TO_PAVE_BSPLINE,FAS,sb_type_norm.NORM_L2,0,type_undec_filter.U_B3SPLINE_2)
MR_Data.transform(Data)
imshow(MR_Data.band(2).data)
MR_Data.size_scale_nl(2)


shapeImage=image.shape
FAS= FilterAnaSynt()
FAS.alloc(type_sb_filter.F_MALLAT_7_9)
MR_Data.alloc(image.shape[0], image.shape[1], 2, type_transform.TO_DIADIC_MALLAT,FAS,sb_type_norm.NORM_L2,0,type_undec_filter.U_B3SPLINE_2)
MR_Data.transform(Data)
print MR_Data.band(0).data.shape
print MR_Data.Set_Transform

naxes=numpy.zeros(3, numpy.float32)
naxes[0] = MR_Data.nbr_band()
naxes[1] = MR_Data.size_ima_nl()
naxes[2] = MR_Data.size_ima_nc()

wavecoeffs=numpy.zeros(naxes,numpy.float32)
for kband in range(MR_Data.nbr_band()):
   wavecoeffs[kband,:,:]=MR_Data.band(kband).data[:]

imshow(wavecoeffs[0,:,:])

MR_Data =MultiResol()
MR_Data.Border=type_border.I_CONT
MR_Data.Verbose= Bool.values[1]
FAS= FilterAnaSynt()
MR_Data.alloc(Data.nl(), Data.nc(), 4, type_transform.TO_SEMI_PYR,FAS,sb_type_norm.NORM_L2,-1,type_undec_filter.U_B3SPLINE_2)
MR_Data.transform(Data)
print MR_Data.Set_Transform
print MR_Data.nbr_scale()
print MR_Data.band(0).data.shape
imshow(MR_Data.band(0).data)
naxes=numpy.zeros(3, numpy.float32)
naxes[0] = MR_Data.nbr_band()
naxes[1] = MR_Data.size_ima_nl()
naxes[2] = MR_Data.size_ima_nc()
wavecoeffs=numpy.zeros(naxes,numpy.float32)
for kband in range(MR_Data.nbr_band()):
   wavecoeffs[kband,0:MR_Data.size_band_nl(kband),0:MR_Data.size_band_nc(kband)]=MR_Data.band(kband).data[:]

imshow(wavecoeffs[0,:,:])
pyfits.writeto("test_DATA.fits",Data.data[:],clobber=True)
process = subprocess.Popen("mr_transform -t19 -n4 test_DATA.fits test_wave_DATA.fits", shell=True)
mr=pyfits.getdata("test_wave_DATA.fits.mr").astype("float32")
imshow(mr[0,:,:]-MR_Data.band(0).data)
imshow(wavecoeffs[0,:,:]-MR_Data.band(0).data)
print mr[0,0,0],wavecoeffs[0,0,0],MR_Data.band(0).data[0,0]
TM_TO_PYR

io_write_ima_float("test_DATA_intern.fits",Data)
process = subprocess.Popen("mr_transform -t19 -n4 test_DATA_intern.fits test_wave_DATA_intern.fits", shell=True)
mr_intern=pyfits.getdata("test_wave_DATA_intern.fits.mr").astype("float32")

#####
#TEST BSPLINE
MR_Data =MultiResol()
MR_Data.Border=type_border.I_CONT
MR_Data.Verbose= Bool.values[1]
FAS= FilterAnaSynt()
FAS.alloc(type_sb_filter.F_MALLAT_7_9)
MR_Data.alloc(Data.nl(), Data.nc(), 4, type_transform.TO_MALLAT,FAS,sb_type_norm.NORM_L2,-1,type_undec_filter.U_B3SPLINE_2)
MR_Data.transform(Data)
print MR_Data.Set_Transform
print MR_Data.nbr_scale()
print MR_Data.band(0).data.shape
imshow(MR_Data.band(0).data)
Ima=Iflt(MR_Data.size_ima_nl(), MR_Data.size_ima_nc())
ortho_trans_to_ima(MR_Data,Ima)
wavecoeffs=Ima.data[:]

imshow(wavecoeffs)
shapeImage=image.shape
MRTransObj=dict(nLines=shapeImage[0],nColumns=shapeImage[1],NbrPlan=4,Border=MR_Data.Border,Transform=type_transform.TO_MALLAT, SB_Filter=type_sb_filter.F_MALLAT_7_9 ,U_Filter= type_undec_filter.U_B3SPLINE_2, LiftingTrans=type_lift.TL_INT_HAAR, Norm= sb_type_norm.NORM_L2)
MRTransObj.update(dict(WaveCoeffs=wavecoeffs))

pyfits.writeto("test_DATA.fits",Data.data[:],clobber=True)
process = subprocess.Popen("mr_transform -t2 -n4 test_DATA.fits test_wave_DATA.fits", shell=True)
mr=pyfits.getdata("test_wave_DATA.fits.mr").astype("float32")
imshow(mr[0,:,:]-MR_Data.band(0).data)
print mr[0,0,0],wavecoeffs[0,0,0],MR_Data.band(0).data[0,0]
TM_TO_PYR

io_write_ima_float("test_DATA_intern.fits",Data)
process = subprocess.Popen("mr_transform -t19 -n4 test_DATA_intern.fits test_wave_DATA_intern.fits", shell=True)
mr_intern=pyfits.getdata("test_wave_DATA_intern.fits.mr").astype("float32")




#-----------------------------------------------------------------------
import numpy
import argparse
import subprocess
import pyfits
from IsapPyWrapper.sparse2d.mr_transform import *
from isap_great3_tools import *
from isap_great3_sparse2d import *


image=numpy.zeros((192,96), numpy.float32)+1.
for kx in range(192):
   for ky in range(96):
      image[kx,ky]=kx+ky*192




Data=Iflt(image.shape[0],image.shape[1])
Data.data[:]=image[:]
image[:]=numpy.random.randn(192,96)
Data.data[:]=image[:]

command="-t14 -n4 -v -T1 -L"
PyMR= MRTrans(command)
coef = PyMR.transform(image)
print coef.MRSignature()
PyMR.parseTransCommand("-t2 -n4 -v -u4")
print coef.MRSignature()
coef2=PyMR.transform(image)
PyMR.checkInputObject(coef)

PyMR.setAllCoefs(coef)
t3=PyMR.recons()

image[:]=numpy.random.randn(192,96)
Data.data[:]=image[:]
command="-t29 -v -l6 "
PyMR= MRTrans(command)
coef4 = PyMR.transform(image)
pyfits.writeto("test_DATA.fits",image,clobber=True)
process = subprocess.Popen("mr_transform "+command+" test_DATA.fits test_wave_DATA.fits", shell=True)
mr=pyfits.getdata("test_wave_DATA.fits.mr").astype("float32")

imshow(coef3.Coefs-coef.Coefs)
MR_Data =MultiResol()
MR_Data.Border=type_border.I_CONT
MR_Data.Verbose= Bool.values[1]
FAS= FilterAnaSynt()
FAS.alloc(type_sb_filter.F_MALLAT_7_9)
MR_Data.alloc(Data.nl(), Data.nc(), 4, type_transform.TO_LIFTING,FAS,sb_type_norm.NORM_L2,-1,type_undec_filter.U_B3SPLINE_2)
MR_Data.LiftingTrans=type_lift.TL_F79
MR_Data.transform(Data)


#Command Options

command="-t14 -n4 -L -v"
lstcc=[]
lstwd=[]
lstrd=[]
for k in range(1,14):
   cc=command + "-T"+str(k)
   PyMR= MRTrans(command)
   lstcc.append(cc)
   coef = PyMR.mtransform(image)
   pyfits.writeto("test_DATA.fits",Data.data[:],clobber=True)
   process = subprocess.Popen("mr_transform "+command+" test_DATA.fits test_wave_DATA.fits", shell=True)
   mr=pyfits.getdata("test_wave_DATA.fits.mr").astype("float32")
   lstwd.append(abs(mr-coef.Coefs).max())
   rec = PyMR.recons()
   lstrd.append(abs(image-rec).max())


figure()
coef = PyMR.transform(image)
imshow(coef.Coefs)
rec = PyMR.mr_recons()
imshow(abs(rec-image)/image)
plt.colorbar()

pyfits.writeto("test_DATA.fits",Data.data[:],clobber=True)
process = subprocess.Popen("mr_transform "+command+" test_DATA.fits test_wave_DATA.fits", shell=True)
mr=pyfits.getdata("test_wave_DATA.fits.mr").astype("float32")
imshow(mr-coef.Coefs)
print abs(mr-coef.Coefs).max()


import time
command="-t14 -n4 -T1 -L -v"
t0 = time.clock()
for k in range(100):
   PyMR= MRTrans(command)
   coef = PyMR.transform(image)
t1=time.clock()
print t1-t0 ;0.7
for k in range(100):
   pyfits.writeto("test_DATA.fits",Data.data[:],clobber=True)
   process = subprocess.Popen("mr_transform "+command+" test_DATA.fits test_wave_DATA.fits", shell=True)
   mr=pyfits.getdata("test_wave_DATA.fits.mr").astype("float32")
t2=time.clock()
print t2-t1, t1-t0

CppMRObj=PyMR._MRTrans__CppMRObj
imshow(PyMR.getBand(0))

cd=coef.Coefs
imshow(coef.Coefs)
CppMRObj=PyMR._MRTrans__CppMRObj
Ima=Iflt(CppMRObj.size_ima_nl(), CppMRObj.size_ima_nc())
ortho_trans_to_ima(CppMRObj,Ima)

#Test the binding
MR_Data =MultiResol()
MR_Data.Border=type_border.I_CONT
MR_Data.Verbose= Bool.values[1]
FAS= FilterAnaSynt()
FAS.alloc(type_sb_filter.F_MALLAT_7_9)
MR_Data.alloc(Data.nl(), Data.nc(), 4, type_transform.TO_MALLAT,FAS,sb_type_norm.NORM_L2,-1,type_undec_filter.U_B3SPLINE_2)
MR_Data.transform(Data)


Ima=Iflt(MR_Data.size_ima_nl(), MR_Data.size_ima_nc())
ortho_trans_to_ima(MR_Data,Ima)
wavecoeffs=Ima.data[:]
imshow(wavecoeffs)



MR_Data =MultiResol()
MR_Data.Border=type_border.I_CONT
MR_Data.Verbose= Bool.values[1]
FAS= FilterAnaSynt()
FAS.alloc(type_sb_filter.F_MALLAT_7_9)
MR_Data.alloc(Data.nl(), Data.nc(), 4, type_transform.TO_MALLAT,FAS,sb_type_norm.NORM_L2,-1,type_undec_filter.U_B3SPLINE_2)
MR_Data.transform(Data)

pyfits.writeto("test_DATA.fits",Data.data[:],clobber=True)
process = subprocess.Popen("mr_transform "+command+" test_DATA.fits test_wave_DATA.fits", shell=True)
mr=pyfits.getdata("test_wave_DATA.fits.mr").astype("float32")
imshow(mr)
imshow(wavecoeffs[0,:,:]-MR_Data.band(0).data)
print mr[0,0,0],wavecoeffs[0,0,0],MR_Data.band(0).data[0,0]
TM_TO_PYR


