#-------------------
#Example of script to launch python binding to mr_transform
import numpy
import argparse
import subprocess
import pyfits
import IsapPyWrapper.sparse2d
import IsapPyWrapper.sparse2d.mr_transform as mrt
from isap_great3_tools import *
from isap_great3_sparse2d import *
pylab

command="-t14 -n2 -v -T1 -L -b0"
PyMR= mrt.MRTrans(command) #Need only to be initialized once

PyMR.parseTransHelp()


#FIRST CHECK CENTERING
image1=numpy.random.randn(84,84)
image=numpy.zeros((96,96))
image[7:91,7:91]= image1 #For strict equivalence, all 7 filter elements should be 0
image2=numpy.zeros((192,192))
#96x96 image centered on 48:144
image2[48:144,48:144]=image
#image2[144,144]=0
#plt.figure(),imshow(image2)

Data=Iflt(image.shape[0],image.shape[1])
Data.data[:]=image[:]

MRObj = PyMR.transform(image)#Get structure containing all coefficients
MRObj2 = PyMR.transform(image2)#Get structure containing all coefficients
#MRObj.MRSignature()
plt.figure(),imshow(MRObj2.Coefs[120:168,120:168]-MRObj.Coefs[48:96,48:96]),\
           plt.colorbar()
plt.figure(),imshow(MRObj2.Coefs[24:72,120:168]-MRObj.Coefs[0:48,48:96]),\
           plt.colorbar()
plt.figure(),imshow(MRObj2.Coefs[120:168,24:72]-MRObj.Coefs[48:96,0:48]),\
           plt.colorbar()
plt.figure(),imshow(MRObj2.Coefs[24:72,24:72]-MRObj.Coefs[0:48,0:48]),\
           plt.colorbar()

print((MRObj2.Coefs[120:125,120:125]-MRObj.Coefs[48:53,48:53]))
print((MRObj2.Coefs[164:168,164:168]-MRObj.Coefs[92:96,92:96]))

#NOW SEE WHICH COEFFICIENTS ARE AFFECTED BY DIFFERENT BORDER CONDITIONS (0 vs mirror)
image=numpy.random.randn(96,96)
image2=numpy.zeros((192,192))
#96x96 image centered on 48:144
image2[48:144,48:144]=image
MRObj = PyMR.transform(image)#Get structure containing all coefficients
MRObj2 = PyMR.transform(image2)#Get structure containing all coefficients
plt.figure(),imshow(MRObj2.Coefs[120:168,120:168]-MRObj.Coefs[48:96,48:96]),\
           plt.colorbar()
plt.figure(),imshow(MRObj2.Coefs[24:72,120:168]-MRObj.Coefs[0:48,48:96]),\
           plt.colorbar()
plt.figure(),imshow(MRObj2.Coefs[120:168,24:72]-MRObj.Coefs[48:96,0:48]),\
           plt.colorbar()
plt.figure(),imshow(MRObj2.Coefs[24:72,24:72]-MRObj.Coefs[0:48,0:48]),\
           plt.colorbar()

command="-t14 -n4 -v -T1 -L -b0"
image=numpy.random.randn(90,95)
image2=numpy.zeros((180,190))
image2[48:138,48:143]=image
PyMR= mrt.MRTrans(command) #Need only to be initialized once
MRObj = PyMR.transform(image)#Get structure containing all coefficients
MRObj2 = PyMR.transform(image2)#Get structure containing all coefficients
icoefs=PyMR.getWaveInnerCoefsLim(image.shape)
#An offset by 1 pix = non 0 diff for two borders for every plane
#icoefs=icoefs+1 
#icoefs=icoefs-1
sz_sc=numpy.asarray(image.shape)
for ksc in range( PyMR.NbrPlan-1):
   kb=3*ksc
   off1=24/(2**ksc)
   off2=off1+sz_sc//2
   sz_sc=(sz_sc+1)//2
   inner_image=MRObj.Coefs[icoefs[0,0, kb]:icoefs[1,0, kb],icoefs[0,1, kb]:icoefs[1,1, kb]]
   inner_image2=MRObj2.Coefs[icoefs[0,0, kb]+off1:icoefs[1,0, kb]+ off1,icoefs[0,1, kb]+ off2[1]:icoefs[1,1, kb]+ off2[1]]
   plt.figure(),imshow(inner_image-inner_image2),plt.colorbar()
   kb=3*ksc+1
   inner_image=MRObj.Coefs[icoefs[0,0, kb]:icoefs[1,0, kb],icoefs[0,1, kb]:icoefs[1,1, kb]]
   inner_image2=MRObj2.Coefs[icoefs[0,0, kb]+ off2[0]:icoefs[1,0, kb]+ off2[0],icoefs[0,1, kb]+ off1:icoefs[1,1, kb]+ off1]
   plt.figure(),imshow(inner_image-inner_image2),plt.colorbar()
   kb=3*ksc+2
   inner_image=MRObj.Coefs[icoefs[0,0, kb]:icoefs[1,0, kb],icoefs[0,1, kb]:icoefs[1,1, kb]]
   inner_image2=MRObj2.Coefs[icoefs[0,0, kb]+ off2[0]:icoefs[1,0, kb]+ off2[0],icoefs[0,1, kb]+ off2[1]:icoefs[1,1, kb]+ off2[1]]
   plt.figure(),imshow(inner_image-inner_image2),plt.colorbar()
kb=3*(ksc+1)
inner_image=MRObj.Coefs[icoefs[0,0, kb]:icoefs[1,0, kb],icoefs[0,1, kb]:icoefs[1,1, kb]]
inner_image2=MRObj2.Coefs[icoefs[0,0, kb]+off1:icoefs[1,0, kb]+ off1,icoefs[0,1, kb]+off1:icoefs[1,1, kb]+off1]
plt.figure(),imshow(inner_image-inner_image2),plt.colorbar()



PyMR._MRTrans__CppMRObj.size_ima_nl()
PyMR._MRTrans__CppMRObj.size_ima_nc()
#OK: a[nl,nc]

PyMR._MRTrans__CppMRObj.filter_bank_column().size_analysis()
PyMR._MRTrans__CppMRObj.filter_bank_column().analysis()
sbf_col=SubBandFilter(PyMR._MRTrans__CppMRObj.filter_bank_column())
sbf_lin=SubBandFilter(PyMR._MRTrans__CppMRObj.filter_bank_line())
s=SubBand2D(sbf_lin, sbf_col)

reload(mrt)

import timeit
import cProfile
image=numpy.random.randn(90,95)
MRObj = PyMR.transform(image)#Get structure containing all coefficients
#USE MGRID
ilist=0
a1=0
numpy.random.seed(0)
globals1=dict([("PyMR", PyMR),("MRObj",MRObj),("ilist",ilist),("a1",a1)])
cProfile.runctx('ilist=PyMR.getWaveInnerCoefsList()', globals1,None)
cProfile.runctx('a1=PyMR.getWaveInnerCoefsFromList(ilist,MRObj)', globals1,None)
cProfile.runctx('ilist=PyMR.getWaveInnerCoefsList() ; a1=PyMR.getWaveInnerCoefsFromList(ilist,MRObj)', globals1,None)
#(0.007,0.000,0.007)
stup='import numpy ; import IsapPyWrapper.sparse2d ; import isap_great3_tools ; import isap_great3_sparse2d; \
import IsapPyWrapper.sparse2d.mr_transform as mrt; command="-t14 -n4 -v -T1 -L -b0"; PyMR=mrt.MRTrans(command) ; \
image=numpy.random.randn(90,95) ;MRObj = PyMR.transform(image)'
timeit.timeit('ilist=PyMR.getWaveInnerCoefsList()',setup= stup, number=10000)/10000.
stup1=stup+' ; ilist=PyMR.getWaveInnerCoefsList()'
timeit.timeit('PyMR.getWaveInnerCoefsFromList(ilist,MRObj,method=0)',setup= stup1, number=10000)/10000.
timeit.timeit('ilist=PyMR.getWaveInnerCoefsList() ; PyMR.getWaveInnerCoefsFromList(ilist,MRObj,method=0)',setup= stup, number=10000)/10000.
#(0.0015865876913070678,9.305448532104493e-05, 0.0016791526079177857)

#USE INDICES
ilist2=0
a2=0
numpy.random.seed(0)
globals2=dict([("PyMR", PyMR),("MRObj",MRObj),("ilist2",ilist2),("a2",a2)])
cProfile.runctx('ilist2=PyMR.getWaveInnerCoefsList(method=1)', globals2,None)
cProfile.runctx('a2=PyMR.getWaveInnerCoefsFromList(ilist2,MRObj,method=1)', globals2,None)
cProfile.runctx('ilist2=PyMR.getWaveInnerCoefsList(method=1) ; a2=PyMR.getWaveInnerCoefsFromList(ilist2,MRObj,method=1)', globals3,None)
#(0.006,0.001,0.007)
timeit.timeit('ilist2=PyMR.getWaveInnerCoefsList(method=1)',setup= stup, number=10000)/10000.
stup2=stup+' ; ilist2=PyMR.getWaveInnerCoefsList(method=1)'
timeit.timeit('PyMR.getWaveInnerCoefsFromList(ilist2,MRObj,method=1)',setup= stup2, number=10000)/10000.
timeit.timeit('ilist2=PyMR.getWaveInnerCoefsList(method=1) ; PyMR.getWaveInnerCoefsFromList(ilist2,MRObj,method=1)',setup= stup, number=10000)/10000.
#(0.0015532696962356567,0.0003838843107223511,0.0019288868904113769)

#USE OGRID
ilist3=0
a3=0
numpy.random.seed(0)
globals3=dict([("PyMR", PyMR),("MRObj",MRObj),("ilist3",ilist3),("a3",a3)])
cProfile.runctx('ilist3=PyMR.getWaveInnerCoefsList(method=2)', globals3,None)
cProfile.runctx('a3=PyMR.getWaveInnerCoefsFromList(ilist3,MRObj,method=2)', globals3,None)
cProfile.runctx('ilist3=PyMR.getWaveInnerCoefsList(method=2) ; a3=PyMR.getWaveInnerCoefsFromList(ilist3,MRObj,method=2)', globals3,None)
#(0.004,0.002,0.006)
timeit.timeit('ilist3=PyMR.getWaveInnerCoefsList(method=2)',setup= stup, number=10000)/10000.
stup3=stup+' ; ilist3=PyMR.getWaveInnerCoefsList(method=2)'
timeit.timeit('PyMR.getWaveInnerCoefsFromList(ilist3,MRObj,method=2)',setup= stup3, number=10000)/10000.
timeit.timeit('ilist3=PyMR.getWaveInnerCoefsList(method=2) ; PyMR.getWaveInnerCoefsFromList(ilist3,MRObj,method=2)',setup= stup, number=10000)/10000.
#(0.0010706301927566528,0.0005711760997772217,0.0016435285091400147)

#USE MASK
ilist4=0
a4=0
numpy.random.seed(0)
globals4=dict([("PyMR", PyMR),("MRObj",MRObj),("ilist4",ilist4),("a4",a4),("numpy",numpy)])
cProfile.runctx('mask=PyMR.getWaveInnerCoefsMask() ; ilist4=numpy.nonzero(mask)', globals4,None)
cProfile.runctx('a4=MRObj.Coefs[ilist4]', globals4,None)
cProfile.runctx('mask=PyMR.getWaveInnerCoefsMask() ; ilist4=numpy.nonzero(mask) ; a4=MRObj.Coefs[ilist4]', globals4,None)
#(0.003,0.001,0.004)
timeit.timeit('mask=PyMR.getWaveInnerCoefsMask() ; ilist4=numpy.nonzero(mask)',setup= stup, number=10000)/10000.
stup4=stup+' ; mask=PyMR.getWaveInnerCoefsMask() ; ilist4=numpy.nonzero(mask)'
timeit.timeit('a4=MRObj.Coefs[ilist4]',setup= stup4, number=10000)/10000.
timeit.timeit('mask=PyMR.getWaveInnerCoefsMask() ; ilist4=numpy.nonzero(mask); a4=MRObj.Coefs[ilist4]',setup= stup, number=10000)/10000.
#(0.00066387128829956062,0.0004037280797958374, 0.0010598371028900148)

#('TIME =', 2.8998851776123047, 1.5499591827392578)
#OK EITHER MASK OR MGRID
numpy.max(numpy.abs(a1-a2)) #0.0
numpy.max(numpy.abs(a1-a3)) #0.0
numpy.max(numpy.abs(sort(a1)-sort(a4))) #0.0

#TEST CHI2 EVALUATION
#USE MASK ONLY
numpy.random.seed(0)
chi2_1=0
b= numpy.random.randn(MRObj.Coefs.shape[0],MRObj.Coefs.shape[1]) ;\
nmask= numpy.abs(numpy.random.randn(MRObj.Coefs.shape[0],MRObj.Coefs.shape[1]))*mask
cProfile.run('chi2_1=numpy.sum(nmask*(b-MRObj.Coefs)**2/0.1)')
cstup_chi2=stup + ' ;mask=PyMR.getWaveInnerCoefsMask() ; b= numpy.random.randn(MRObj.Coefs.shape[0],MRObj.Coefs.shape[1]) ;\
nmask= numpy.abs(numpy.random.randn(MRObj.Coefs.shape[0],MRObj.Coefs.shape[1]))*mask'
print(timeit.timeit('chi2_1=numpy.sum(nmask*(b-MRObj.Coefs)**2/0.1)',setup= stup_chi2, number=10000)/10000.,chi2_1)
#(0.0001558711051940918,101327.56531021901)

#USE MGRID
numpy.random.seed(0)
chi2_2=0
ilist=PyMR.getWaveInnerCoefsList()
b1d=numpy.ravel(b)[ilist]
mask1d=numpy.ravel(nmask)[ilist]
cProfile.run('chi2_2=numpy.sum(mask1d*(b1d-numpy.ravel(MRObj.Coefs)[ilist])**2/0.1)')
stup_chi2=stup + ' ;mask=PyMR.getWaveInnerCoefsMask() ; b= numpy.random.randn(MRObj.Coefs.shape[0],MRObj.Coefs.shape[1]) ;\
nmask= numpy.abs(numpy.random.randn(MRObj.Coefs.shape[0],MRObj.Coefs.shape[1]))*mask ;ilist=PyMR.getWaveInnerCoefsList() ;\
b1d=numpy.ravel(b)[ilist] ; mask1d=numpy.ravel(nmask)[ilist]'
print(timeit.timeit('chi2_2=numpy.sum(mask1d*(b1d-numpy.ravel(MRObj.Coefs)[ilist])**2/0.1)',setup= stup_chi2, number=10000)/10000.,chi2_2)
#('TIME =', 0.00021337649822235108, 101327.56531021952)

numpy.random.seed(0)
chi2_3=0
row2,col2=PyMR.getWaveInnerCoefsList(method=1)
b2d=b[row2,col2]
mask2d= nmask[row2,col2]
cProfile.run('chi2_3=numpy.sum(mask2d*(b2d-MRObj.Coefs[row2,col2])**2/0.1)')
stup_chi3=stup + ' ; import numpy;mask=PyMR.getWaveInnerCoefsMask() ; b= numpy.random.randn(MRObj.Coefs.shape[0],MRObj.Coefs.shape[1]) ;\
nmask= numpy.abs(numpy.random.randn(MRObj.Coefs.shape[0],MRObj.Coefs.shape[1]))*mask ;row2,col2=PyMR.getWaveInnerCoefsList(method=1) ;\
b2d=b[row2,col2] ; mask2d= nmask[row2,col2]'
print(timeit.timeit('chi2_3=numpy.sum(mask2d*(b2d-MRObj.Coefs[row2,col2])**2/0.1)',setup= stup_chi3, number=10000)/10000.,chi2_3)
#('TIME =',0.0005341107130050659, 101327.56531021952)

#USE MASK+indices
numpy.random.seed(0)
chi2_4=0
mask=PyMR.getWaveInnerCoefsMask() 
ilist4=numpy.nonzero(mask)
b2d=b[ilist4]
mask2d=nmask[ilist4]
start = time.time()
cProfile.run('chi2_4=numpy.sum(mask2d*(b2d-MRObj.Coefs[ilist4])**2/0.1)')
stup_chi4=stup + ' ; import numpy;mask=PyMR.getWaveInnerCoefsMask() ; b= numpy.random.randn(MRObj.Coefs.shape[0],MRObj.Coefs.shape[1]) ;\
nmask= numpy.abs(numpy.random.randn(MRObj.Coefs.shape[0],MRObj.Coefs.shape[1]))*mask ;ilist4=numpy.nonzero(mask) ;\
b2d=b[ilist4]; mask2d=nmask[ilist4]'
print(timeit.timeit('chi2_4=numpy.sum(mask2d*(b2d-MRObj.Coefs[ilist4])**2/0.1)',setup= stup_chi4, number=10000)/10000.,chi2_4)
#('TIME =', 0.0005260081052780152, 101327.56531021901)


#GET DIRECTLY THE COEFS
numpy.random.seed(0)
chi2_5=0
ilist=PyMR.getWaveInnerCoefsList()
mask1d=numpy.ravel(nmask)[ilist]
b1d=numpy.ravel(b)[ilist]
cProfile.run('chi2_5=numpy.sum(mask1d*(b1d-PyMR.getWaveInnerCoefs(MRObj))**2/0.1)')
stup_chi5=stup + ' ;mask=PyMR.getWaveInnerCoefsMask() ; b= numpy.random.randn(MRObj.Coefs.shape[0],MRObj.Coefs.shape[1]) ;\
nmask= numpy.abs(numpy.random.randn(MRObj.Coefs.shape[0],MRObj.Coefs.shape[1]))*mask ;ilist=PyMR.getWaveInnerCoefsList() ;\
b1d=numpy.ravel(b)[ilist] ; mask1d=numpy.ravel(nmask)[ilist]'
print(timeit.timeit('chi2_5=numpy.sum(mask1d*(b1d-PyMR.getWaveInnerCoefs(MRObj))**2/0.1)',setup= stup_chi5, number=10000)/10000.,chi2_5)
#('TIME =', 0.0006892024993896484, 101327.56531021952)

chi2_1-chi2_2 # -5.0931703299283981e-10
chi2_1-chi2_3 # -5.0931703299283981e-10
chi2_1-chi2_4 #0.0
chi2_1-chi2_5 #-5.0931703299283981e-10

#CCL: USE MASK ONLY DIRECTLY

#Check border effects
command="-t14 -n4 -v -T1 -L -b1" #Garantee exact reconstruction using 7/9
image=numpy.random.randn(128,128)
PyMR= mrt.MRTrans(command) #Need only to be initialized once
MRObj = PyMR.transform(image)#Get structure containing all coefficients
command2="-t14 -n4 -v -T2 -L -b3"
PyMR2=mrt.MRTrans(command2)
MRObj2 = PyMR2.transform(image)#Get structure containing all coefficients
MRObj2.coefs[:]= MRObj.coefs[:]
tt=PyMR2.recons()
plt.figure(),imshow(image-tt),plt.colorbar()
tt=PyMR.recons()
plt.figure(),imshow(image-tt),plt.colorbar()

plt.figure(),imshow(MRObj.Coefs-MRObj2.Coefs)
mask=PyMR.getWaveInnerCoefsMask() 
d=MRObj.Coefs-MRObj2.Coefs
d[where(mask==1)]=-0.5
plt.figure(),imshow(d),plt.colorbar()
where(d==0)
MRObj.Coefs*=mask
PyMR.setAllCoefs(MRObj)


