#-----------------------------------------------------------------------
import numpy
import argparse
import subprocess
import pyfits
from IsapPyWrapper.sparse2d.mr_transform import *
from isap_great3_tools import *
from isap_great3_sparse2d import *
import time 

pylab


##############################
#SHORT REMARK ON ARRAY/DISPLAY CONVENTIONS IN PYTHON

#Example of numpy image, 192 lines, 96 columns
image_test=numpy.zeros((192,96), numpy.float32)+1.
for kx in range(192):
   for ky in range(96):
      image_test[kx,ky]=kx+ky*96
plt.figure(),imshow(image_test)
#This illustrates the default representation which puts origin in upper left
# but can be changed:
plt.figure(),imshow(image_test,origin='lower')


##############################
#BASIC USAGE OF THE WRAPPERS

#Now try to launch a transform
PyMR= MRTrans("") 
#Get Help
PyMR.parseTransHelp()
#Transform command
command="-t14 -n4 -v -T1 -L"
#Assign command to object PyMR
PyMR.parseTransCommand(command)
#One can also construct the transform with the constructor
PyMR= MRTrans(command) 
#Compute the coefficients (deepcopy)
coef = PyMR.transform(image)
#Reconstruct the image
t3=PyMR.recons()
plt.figure(), imshow(t3)
#Get Band0 (DEEP COPY)
kband=2
band=PyMR.getBand(kband)
plt.figure(),imshow(band.Coefs)
band.Coefs.shape
#Hard Threshold it:
n_std=numpy.std(band.Coefs[:])
band.Coefs[numpy.where(numpy.abs(band.Coefs)<3.*n_std)]=0
plt.figure(),imshow(band.Coefs),plt.colorbar()
#Note:as it is a deepcopy, the original coefs are untouched 
plt.figure(),imshow(coef.Coefs)
#set it back (checks are made to ensure it has the right MRSignature and
#the number of the band is in band0._PyMRObj__band
PyMR.setBand(band)
coefn=PyMR.getAllCoefs()
tn=PyMR.recons()
plt.figure(),imshow(np.rot90(coefn.Coefs))
plt.figure(),imshow(tn-image),title("Band {0}".format(kband))

##############################
#ADVANCED USAGE

#coef is an object, which has a signature, used in particular when the inverse
#transform is required.t
print coef.MRSignature()
#Try a different command and get the signature
PyMR.parseTransCommand("-t2 -n4 -v -u4")
coef2=PyMR.transform(image)
print coef2.MRSignature()
#Now the PyMR object and the coef object do not match anymore (two different 
#transform)
#One can manually check if coefficients where obtained with the current transform
PyMR.checkInputObject(coef)
#But Note that PyMR objects are automatically resetted to match the coefficients
#type when we set the coef for reconstruction 
PyMR.setAllCoefs(coef)
t3=PyMR.recons()
plt.figure(), imshow(t3)



#COMPARE DIRECT CALL TO C++, AND WRAPPERS
#Create gaussian noise in rectangular image
image[:]=numpy.random.randn(192,96)
Data.data[:]=image[:]
command="-t29 -v -l6 "
PyMR= MRTrans(command)
coef4 = PyMR.transform(image)
pyfits.writeto("test_DATA.fits",image,clobber=True)
process = subprocess.Popen("mr_transform "+command+
                             " test_DATA.fits test_wave_DATA.fits", shell=True)
#WE HAVE TO WAIT A LITTLE BIT TO ENSURE OUTPUT IS ACTUALLY WRITTEN 
#(AND NOT JUST PLANNED TO BE WRITTEN)
time.sleep(1)
mr=pyfits.getdata("test_wave_DATA.fits.mr").astype("float32")
plt.figure(),imshow(coef4.Coefs-mr)
numpy.max(numpy.abs(coef4.Coefs-mr)) #0.0

#HERE WE HAVE C++ OBJECTS THAT WE CAN DIRECTLY PLAY WITH
MR_Data =MultiResol()
MR_Data.Border=type_border.I_CONT
MR_Data.Verbose= Bool.values[1]
FAS= FilterAnaSynt()
FAS.alloc(type_sb_filter.F_MALLAT_7_9)
MR_Data.alloc(Data.nl(), Data.nc(), 4, type_transform.TO_LIFTING,FAS,
                         sb_type_norm.NORM_L2,-1,type_undec_filter.U_B3SPLINE_2)
MR_Data.LiftingTrans=type_lift.TL_F79
MR_Data.transform(Data)



##############################
#MORE ##############################
#MORE ADVANCED USAGE
ADVANCED USAGE

#HERE WE CHECK ON MORE FILTERS TO CHECK CONSISTENCY OF THE WRAPPERS VS C++ CODE
command="-t14 -n4 -L -v"
lstcc=[]
lstwd=[]
for k in range(1,14):
   if(k != 10):#USER DEFINED FILTER
      cc=command + " -T"+str(k)
      print(cc)
      PyMR= MRTrans(cc)
      lstcc.append(cc)
      coef = PyMR.transform(image)
      pyfits.writeto("test_DATA.fits",Data.data[:],clobber=True)
      process = subprocess.Popen("mr_transform "+cc+" test_DATA.fits test_wave_DATA.fits", shell=True)
      time.sleep(1)
      mr=pyfits.getdata("test_wave_DATA.fits.mr").astype("float32")
      lstwd.append(numpy.abs(mr-coef.Coefs).max())
print(lstwd)

#HERE WE CHECK TIME IMPROVEMENT 
import time
import cProfile
import timeit
command="-t14 -n4 -T1 -L -v"   
PyMR= MRTrans(command)
image=numpy.random.randn(90,95)
globals1=dict([("PyMR", PyMR),("image",image)])

stup='import numpy ; import subprocess; import pyfits;import IsapPyWrapper.sparse2d ; import isap_great3_tools ; import isap_great3_sparse2d; \
import IsapPyWrapper.sparse2d.mr_transform as mrt; command="-t14 -n4 -v -T1 -L "; PyMR=mrt.MRTrans(command) ; \
image=numpy.random.randn(90,95)'
timeit.timeit('Coef=PyMR.transform(image)',setup= stup, number=1000)/1000. 
#0.0021831610202789307
cProfile.runctx('Coef=PyMR.transform(image)', globals1,None)
script = """\
   pyfits.writeto("test_DATA.fits",image,clobber=True)
   process = subprocess.Popen("mr_transform "+command+" test_DATA.fits test_wave_DATA.fits", shell=True)
   time.sleep(0.01)
   mr=pyfits.getdata("test_wave_DATA.fits.mr").astype("float32")
    """
timeit.timeit(script,setup= stup, number=100)/100. 
#0.057331390380859375 OF WHICH WE SLEEP 0.01s to avoid crashes


##############################
#EXPERT USAGE

#One of the key problem for C++ interfacing of the code was conversion between
#numpy arrays and Ifloat arrays. Here we show how to (deep) copy an numpy array 
#into a Ifloat array.
image_test=numpy.random.randn(192,96)
Data=Iflt(image_test.shape[0],image_test.shape[1])
Data.data[:]=image_test[:]
#Data is actually a C++ object with many of the original IFloat members exposed
Data.nc()
Data.n_elem()
Data.maxfabs()
#The other copy around is similar
image_test_back=numpy.zeros((192,96))
image_test_back[:]=Data.data[:]
plt.figure(),imshow(image_test-image_test_back)

#The Cpp MRTrans Object is also accessible (class MultiResol)
CppMRObj=PyMR._MRTrans__CppMRObj
#And we can use it to perform orthogonal transforms directly! (We use a C++
#object in a C++ function with a C++ array, all from Python !)
Ima=Iflt(CppMRObj.size_ima_nl(), CppMRObj.size_ima_nc())
ortho_trans_to_ima(CppMRObj,Ima)

#Test the C++ binding: all the following commands operates on C++ objects
MR_Data =MultiResol()
MR_Data.Border=type_border.I_CONT
MR_Data.Verbose= Bool.values[1]
FAS= FilterAnaSynt()
FAS.alloc(type_sb_filter.F_MALLAT_7_9)
MR_Data.alloc(Data.nl(), Data.nc(), 4, type_transform.TO_MALLAT,FAS,
                         sb_type_norm.NORM_L2,-1,type_undec_filter.U_B3SPLINE_2)
MR_Data.transform(Data)
band1=MR_Data.band(1)
plt.figure(),imshow(band1.data[:])
#We can also use the functions not in the class
Ima=Iflt(MR_Data.size_ima_nl(), MR_Data.size_ima_nc())
ortho_trans_to_ima(MR_Data,Ima)
wavecoeffs=Ima.data[:]
plt.figure(),imshow(wavecoeffs)

#Note all transforms type can be visualized here
type_transform.names.keys()
#SAME THINGS for filters, e.g.
type_undec_filter.names.keys()
#EG. STARLET
MR_Data =MultiResol()
MR_Data.Border=type_border.I_CONT
MR_Data.Verbose= Bool.values[1]
MR_Data.alloc(Data.nl(), Data.nc(), 4, type_transform.TO_PAVE_BSPLINE)
MR_Data.transform(Data)
band3=MR_Data.band(3)
plt.figure(),imshow(band3.data[:])
#Equivalent C++ executable command
command=" -t2 -n4 -v "
pyfits.writeto("test_DATA.fits",Data.data[:],clobber=True)
process = subprocess.Popen("mr_transform "+command+ 
                              " test_DATA.fits test_wave_DATA.fits", shell=True)
time.sleep(0.5)
mr=pyfits.getdata("test_wave_DATA.fits.mr").astype("float32")
mr.shape()
#Compare the two approaches
plt.figure(),imshow(mr[3,:,:]-band3.data)
plt.figure(),imshow(mr[0,:,:]-MR_Data.band(0).data)
print mr[0,0,1],MR_Data.band(0).data[0,1]



