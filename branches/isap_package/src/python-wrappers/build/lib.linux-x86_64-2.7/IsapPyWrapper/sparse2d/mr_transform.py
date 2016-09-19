"""! 
   mr_transform.py - Compute Multi-Resolution transforms using ISAP 
   boost-pythons wrappers.
"""

# -- Python imports

import numpy
import argparse
from isap_great3_tools import *
from isap_great3_sparse2d import *
import copy
import cPickle




def boostEnumFindKey(typeValue,enumType):
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

def boostEnumFindObject(KeyValue,enumType):
   """! 
      Get the value associated to a boost python type in a 
      boost python enumeration
      @param KeyValue: a boost python key 
      @param enumType: a boost python enumeration
      @return value in the boost enumeration
   """
   Nobj=len(enumType.names.keys())
   return [enumType.values.values()[k] for k in range(Nobj) \
                                if  enumType.values.keys()[k]== KeyValue][0]



###############################################################################
class PyMRObj:
   """! 
      Contain either wavelet decomposition or band of wavelet coefficients along
      with the signature of the MR Transform that was used to derive them.
   """
   __signature =0
   __band=0
   __allcoeffs=0
   def __init__(self,command):
      self.__signature = { "nLines": 0, "nColumns": 0,"NbrPlan": 0,
            "Border": type_border.values[0],"Transform": 0, "SB_Filter": 0,
            "U_Filter": 0, "LiftingTrans": 0, "Norm": 0,"WaveCoeffs": 0,
            "NbrUndec":-1,"Command":command }
      self.__band=-1
      self.__allcoeffs=False

   @property
   def nLines(self):
      """! @return the nLines of input image. """
      return self.__signature["nLines"]
   @property
   def nColumns(self):
      """! @return the nColumns of input image. """
      return self.__signature["nColumns"]
   @property
   def NbrPlan(self):
      """! @return the Number of Scales. """
      return self.__signature["NbrPlan"]
   @property
   def Border(self):
      """! @return the Border Type applied in MResol. """
      return self.__signature["Border"]
   @property
   def Transform(self):
      """! @return the transform Type. """
      return self.__signature["Transform"]
   @property
   def SB_Filter(self):
      """! @return the decimated Filter Type. """
      return self.__signature["SB_Filter"]
   @property
   def U_Filter(self):
      """! @return the undecimated Filter Type. """
      return self.__signature["U_Filter"]
   @property
   def LiftingTrans(self):
      """! @return the Lifting Transform Filter Type. """
      return self.__signature["LiftingTrans"]
   @property
   def Norm(self):
      """! @return the Norm applied to the filter. """
      return self.__signature["Norm"]
   @property
   def WaveCoeffs(self):
      """! @return the Wavelet Coefficients. """
      return self.__signature["WaveCoeffs"]
   @property
   def NbrUndec(self):
      """! @return the number of undecimated scales. """
      return self.__signature["NbrUndec"]
   @property
   def Command(self):
      """! @return the command applied. """
      return self.__signature["Command"]
   @property
   def Coefs(self):
      """! @return wavelet coefficients. """
      return self.__signature["WaveCoeffs"]
      
   def CopySwitchBoostEnum(self,backward=False):
      """! @return a copy of self without boostenum \
                                          or get it back(to pickle)."""
      dc=copy.deepcopy(self)
      sign_cp=dc.__signature
      if backward:
         sign_cp["Transform"]= boostEnumFindObject(sign_cp["Transform"],
                            type_transform)
         sign_cp["Border"]= boostEnumFindObject(sign_cp["Border"],
                         type_border)
         sign_cp["SB_Filter"]= boostEnumFindObject(sign_cp["SB_Filter"],
                         type_sb_filter)
         sign_cp["U_Filter"]= boostEnumFindObject(sign_cp["U_Filter"],
                           type_undec_filter)
         sign_cp["LiftingTrans"]= boostEnumFindObject(sign_cp["LiftingTrans"],
                               type_lift)
         sign_cp["Norm"]= boostEnumFindObject(sign_cp["Norm"],sb_type_norm)
      
      else:
         sign_cp["Transform"]= boostEnumFindKey(sign_cp["Transform"],
                            type_transform)
         sign_cp["Border"]= boostEnumFindKey(sign_cp["Border"],
                         type_border)
         sign_cp["SB_Filter"]= boostEnumFindKey(sign_cp["SB_Filter"],
                         type_sb_filter)
         sign_cp["U_Filter"]= boostEnumFindKey(sign_cp["U_Filter"],
                           type_undec_filter)
         sign_cp["LiftingTrans"]= boostEnumFindKey(sign_cp["LiftingTrans"],
                               type_lift)
         sign_cp["Norm"]= boostEnumFindKey(sign_cp["Norm"],sb_type_norm)
      return dc      
   
   
#   @property
   def MRSignature(self):
      """! @return full wavelet coefficients signature. """
      sign_cp=copy.deepcopy(self.__signature)
      sign_cp["Transform"]= boostEnumFindObject(sign_cp["Transform"],
                            type_transform)
      sign_cp["Border"]= boostEnumFindObject(sign_cp["Border"],
                         type_border)
      sign_cp["SB_Filter"]= boostEnumFindObject(sign_cp["SB_Filter"],
                         type_sb_filter)
      sign_cp["U_Filter"]= boostEnumFindObject(sign_cp["U_Filter"],
                           type_undec_filter)
      sign_cp["LiftingTrans"]= boostEnumFindObject(sign_cp["LiftingTrans"],
                               type_lift)
      sign_cp["Norm"]= boostEnumFindObject(sign_cp["Norm"],sb_type_norm)
      
      return sign_cp

 #  @MRSignature.setter
   def SetMRSignature(self,SignDict):
      """! Set full wavelet coefficients signature. """
      if all(k in SignDict for k in self.__signature if not k=="WaveCoeffs"):
         if self.MRSignatureCheck(SignDict):
            for k in self.__signature :
                if not k=="WaveCoeffs":
                   self.__signature[k]=copy.deepcopy(SignDict[k])
            self.__signature["Transform"]= boostEnumFindObject(\
                                         self.Transform,type_transform)
            self.__signature["Border"]= boostEnumFindObject(\
                                         self.Border,type_border)
            self.__signature["SB_Filter"]= boostEnumFindObject(\
                                      self.SB_Filter,type_sb_filter)
            self.__signature["U_Filter"]= boostEnumFindObject(\
                                      self.U_Filter,type_undec_filter)
            self.__signature["LiftingTrans"]= boostEnumFindObject(\
                                      self.LiftingTrans,type_lift)
            self.__signature["Norm"]= boostEnumFindObject(\
                                      self.Norm,sb_type_norm)
         else:
            raise ValueError("Incorrect key value for Signature") 
      else:
            raise TypeError("Incorrect Signature Type") 
      


   def MRSignatureCheck(self,SignDict):
      """! Check if the signature corresponds to a MR object. """
      if (SignDict["nLines"] < 0):
         print "Incorrect value for key nLines: should be positive"
         return False
      elif (SignDict["nColumns"] < 0):
         print "Incorrect value for key nColumns: should be positive"
         return False
      elif (SignDict["NbrPlan"] <= 0):
         print "Incorrect value for key NbrPlan: should be strictly positive"
         return False
      elif (SignDict["Border"] not in type_border.names.values()):
         print "Incorrect value for key Border: should be in type_border"
         return False
      elif (SignDict["Transform"] not in type_transform.names.values()):
         print "Incorrect value for key Transform: should be in type_transform"
         return False
      elif (SignDict["SB_Filter"] not in type_transform.names.values()):
         print "Incorrect value for key SB_Filter: should be in type_sb_filter"
         return False
      elif (SignDict["U_Filter"] not in type_transform.names.values()):
         print "Incorrect value for key U_Filter:"\
                               "should be in type_undec_filter"
         return False
      elif (SignDict["LiftingTrans"] not in type_lift.names.values()):
         print "Incorrect value for key LiftingTrans: should be in type_lift"
         return False
      elif (SignDict["Norm"] not in  sb_type_norm.names.values()):
         print "Incorrect value for key Norm: should be in  sb_type_norm"
         return False
      elif (SignDict["NbrUndec"] < -1):
         print "Incorrect value for key NbrUndec: should be > -1"
         return False
      else:
         return True

   def MRSignatureCompare(self,SignDict):
      """! Compare the current signature with an external signature. """
      if all(k in SignDict for k in self.__signature if not k=="WaveCoeffs"):
         if all(SignDict[k] ==self.__signature[k] for k in self.__signature if \
                            (not k == "Command" ) and (not k=="WaveCoeffs")):
            return True
         else:
            return False 
      else:
         raise TypeError("Incorrect Signature Type") 
         
   @Coefs.setter
   def Coefs(self,WaveCoeffs):
      self.__signature.WaveCoeffs=copy.deepcopy(WaveCoeffs)

   def save(self, filename):
      """! Save the PyMRObj using cpickle. """
      f = open(filename,'wb')
      temp=self.CopySwitchBoostEnum()
      cPickle.dump(temp, f, 0)
      f.close() 

   @classmethod
   def load(cls,filename):
      """! Load a PyMRObj saved using cpickle. """
      f = open(filename,'rb')
      temp=cPickle.load(f)
      if isinstance(temp, PyMRObj):
         obj=temp.CopySwitchBoostEnum(True)
      f.close()
      return obj

###############################################################################
class MRTrans:
   
   """! 
      Perform Wavelet transforms on a 2D image. A Cpp multiresolution object is
      assigned and the signature is saved.
   """

   __signPyMRTrans=0
   __CppMRObj=0
   __CppFlagCompute=0
   Verbose=0
   def __init__(self,command):
      """! 
         Initialize PyMRTrans from command line arguments
         @param command mr_transform options 
     (see \sa PyMRTrans::parseTransHelp or mr_transform executable)
      """
         
      self.setTransParser()
      self.__signPyMRTrans = { "nLines": 0, "nColumns": 0,"NbrPlan": 0,
            "Border": type_border.values[0],"Transform": 0, "SB_Filter": 0,
            "U_Filter": 0, "LiftingTrans": 0, "Norm": 0,"Command": '',
            "NbrUndec":-1 }
      self.parseTransCommand(command)
      self.__CppMRObj=MultiResol()
      self.__CppFlagCompute=True
      self.Verbose = Bool.values[0]
      self.__FAS=FilterAnaSynt()

   # ~~~~~~~~~~~
   # Properties 
   # ~~~~~~~~~~~

#------------------------------------------------------------------------------
   # --- Getters
   @property
   def nLines(self):
      """! @return the number of Lines. """
      return self. __signPyMRTrans["nLines"]

   @property
   def nColumns(self):
      """! @return the number of Columns. """
      return self. __signPyMRTrans["nColumns"]

   @property
   def Transform(self):
      """! @return the transform Type. """
      return self. __signPyMRTrans["Transform"]

   @property
   def LiftingTrans(self):
      """! @return the transform Type. """
      return self. __signPyMRTrans["LiftingTrans"]

   @property
   def SBFilter(self):
      """! @return the decimated filter. """
      return self. __signPyMRTrans["SB_Filter"]

   @property
   def UFilter(self):
      """! @return the undecimated filter. """
      return self. __signPyMRTrans["U_Filter"]

   @property
   def Norm(self):
      """! @return the normalization of the filter. """
      return self. __signPyMRTrans["Norm"]

   @property
   def Command(self):
      """! @return the Helper class """
      return self. __signPyMRTrans["Command"]
   
   @property
   def Border(self):
      """! @return the Border condition. """
      return self. __signPyMRTrans["Border"]

   @property
   def NbrUndec(self):
      """! @return the number of undecimated scales """
      return self. __signPyMRTrans["NbrUndec"]

   @property
   def NbrPlan(self):
      """! @return the number of Scales """
      return self. __signPyMRTrans["NbrPlan"]

   @property
   def NbrBand(self):
      """! @return the number of Bands """
      if not self.__CppFlagCompute:
         return self.__CppMRObj.nbr_band()
      else:
         raise RuntimeError("Fist need to set the transform")

   @property
   def MRSignature(self):
      """! @return the signature of the MR transform """
      return self.__signPyMRTrans
            
   # ~~~~~~~~~~~~~~
   # Public methods 
   # ~~~~~~~~~~~~~~
#------------------------------------------------------------------------------
# --- Parser
   def setTransParser(self):
      """! 
         Set the parser for mr_transform 
      """
      #Get transform or filter size
      nUndecFilter=len(type_undec_filter.names.keys())
      nTransform=len(type_transform.names.keys())
      nLift=len(type_lift.names.keys())
      nFilter=len(type_sb_filter.names.keys())
      UFilterDef=boostEnumFindKey(type_undec_filter.U_B3SPLINE_2,
                                     type_undec_filter)
      TransformDef=boostEnumFindKey(type_transform.TO_PAVE_BSPLINE,
                    type_transform)
      LiftingDef=boostEnumFindKey(type_lift.TL_INT_HAAR,type_lift)
      SBFilterDef=boostEnumFindKey(type_sb_filter.F_MALLAT_7_9,
                                                        type_sb_filter)
      BorderDef=boostEnumFindKey(type_border.I_CONT, type_border)
      
      parser = argparse.ArgumentParser(prog='mr_transform',add_help=False,
                                   description="Process mr_transform options.")

      parser.add_argument('-U',dest='UFilter',type=int,default= UFilterDef,
                help="Type of Undecimated Filter. Default is SplineB3-Id:" 
                "H=[1,4,6,4,1]/16, Ht=H, G=Id-H*H, Gt=Id",
                choices=range(1,nUndecFilter+1))

      parser.add_argument('-u',dest='NbrUndec',default="-1",type=int,
              help="Number of Undecimated Scales", choices=range(nUndecFilter))

      parser.add_argument('-L',dest='L2Norm',action='store_true',default=False,
                          help="Use a L2 normalization. Default is L1")

      parser.add_argument('-v',dest='Verbose',action='store_true',
                default=False,help="Verbose. Default is no.")
      parser.add_argument('-t',dest='TypeTransform',type=int,help="Type of" 
            "MultiResolution Transform. Default is B-Spline Wavelet Transform" 
            "with a trous algorithm.",choices=range(1,nTransform+1),
             default= TransformDef)

      parser.add_argument('-l',dest='TypeLifting',default= LiftingDef,type=int,
                   help="Type of Lifting Transform. Default is Lifting scheme:"
                   "integer Haar WT.", choices=range(1,nLift+1))

      parser.add_argument('-T',dest='TypeFilter',type=int,help="Filter Type."
             "default is Biorthogonal 7/9 filters.",choices=range(1,nFilter+1),
             default= SBFilterDef)
                         
      parser.add_argument('-n',dest='NbrPlan',default=4,type=int,
                 help="Number of scales used in the multiresolution transform."
                "Default is 4.",choices=range(2,10))

      parser.add_argument('-c',dest='NbrIter',default=0,type=int,
        help="Iterative transformation. Iter=number of iterations.This option"
        "is valid only if the chosen multiresolution transform is pyramidal"
        "(6-11). The reconstruction is not exact and we need few iterations."
        "Generally, we take 3.",choices=range(2,20))

      parser.add_argument('-b',dest='Border',default= BorderDef,type=int,
              help="Border constraint. Generally continuous border conditions"\
              "(value 0) is taken.",choices=range(4))

      self.__transparser=parser

   def parseTransCommand(self,command):
      """! 
         Parse the command line according to \sa according to 
         PyMRTrans::setTransParser and set the MultiResol Object
         @param command: string containing the options
         @return parsed mr_transform options
      """
       #Parse arguments
      parsedOptions=self.__transparser.parse_args(command.split())
      U_Filter = type_undec_filter.values[parsedOptions.UFilter-1]
      LiftingTrans = type_lift.values[parsedOptions.TypeLifting]
      Transform = type_transform.values[parsedOptions.TypeTransform-1]
      SB_Filter = type_sb_filter.values[parsedOptions.TypeFilter]
      Border = type_border.values[parsedOptions.Border]
      if parsedOptions.Verbose:
         self.Verbose=Bool.values[1]
      if parsedOptions.L2Norm:
         Norm= sb_type_norm.NORM_L2
      else:
         Norm= sb_type_norm.NORM_L1
      #Check if it matches current MRObject signature
      test_sign = { "nLines": 0, "nColumns": 0,"NbrPlan": 0,
            "Border": type_border.values[0],"Transform": 0, "SB_Filter": 0,
            "U_Filter": 0, "LiftingTrans": 0, "Norm": 0,"Command": '',
            "NbrUndec":-1 }
      test_sign["NbrPlan"]= parsedOptions.NbrPlan
      test_sign["Transform"]= Transform
      test_sign["SB_Filter"]= SB_Filter
      test_sign["U_Filter"]= U_Filter
      test_sign["LiftingTrans"]= LiftingTrans
      test_sign["Norm"]= Norm
      test_sign["NbrUndec"]= parsedOptions.NbrUndec
      test_sign["Border"] = Border
      test_sign["nLines"] = self.__signPyMRTrans["nLines"]
      test_sign["nColumns"] = self.__signPyMRTrans["nColumns"]
      test_sign["Command"] = self.__signPyMRTrans["Command"]
      #If one of the signature does not match, need to reset MultiResol Object
      if not all([True if test_sign[ksign]==self. __signPyMRTrans[ksign] \
              else False for ksign in test_sign.keys()]):
         self.__CppFlagCompute=True
         test_sign["Command"]= command
         for ksign in test_sign.keys():
            self.__signPyMRTrans[ksign]= test_sign[ksign]
      else:
         self. __signPyMRTrans["Command"]= command

   def parseTransHelp(self):
      """! 
         Return the help of the mr_transform parser 
         \sa PyMRTrans::setTransParser
         @return mr_transform parser help
      """
      return self.__transparser.print_help()

#------------------------------------------------------------------------------
#Simple Getters: Direct interface to Cpp
   def getBand(self,kband):
      """! 
         Return a band of wavelet coefficients
         @param kband: band index
         @return wavelet coeffs in band kband
      """
      MR_Obj=self.__CppMRObj
      if not self.__CppFlagCompute:
         if(kband >= 0 and kband < self.NbrBand):
            PyObj=PyMRObj(self.Command)
            PyObj.Coefs= MR_Obj.band(kband).data[:]
            PyObj.SetMRSignature(self.__signPyMRTrans)
            PyObj._PyMRObj__band=kband
            PyObj.allcoeffs=False
            return PyObj
         else:
            raise ValueError("band should be between 0 and {1}".format(0,
                            self.NbrBand-1))
      else:
         raise RuntimeError("Fist need to set the transform")

   def getBandInfo(self,kband):
      """! 
         Return information on a band of wavelet coefficients
         @param kband: band index
         @return wavelet coeffs in band kband
      """
      if not self.__CppFlagCompute:
         if(kband >= 0 and kband < self.NbrBand):
            return self.__CppMRObj.band_to_scale(kband)
         else:
            raise ValueError("band should be between 0 and {1}".format(0,
                              self.NbrBand-1))
      else:
         raise RuntimeError("Fist need to set the transform")

   def getScaleInfo(self,kscale):
      """! 
         Return information on a band of wavelet coefficients
         @param kband: band index
         @return wavelet coeffs in band kband
      """
      if not self.__CppFlagCompute:
         if(kscale >= 0 and kscale < self.NbrPlan):
            return self.__CppMRObj.band_to_scale(kband)
         else:
            raise ValueError("band should be between 0 and {1}".format(0,
                             self.NbrBand-1))
      else:
         raise RuntimeError("Fist need to set the transform")

      
#------------------------------------------------------------------------------
#Simple Setters: Direct interface to Cpp
   # --- Setters
   def setBand(self,band):
      """! 
         Set a band of wavelet coefficients
         @param kband: band index
         @return wavelet coeffs in band kband
      """
      if not self.checkInputObject(band):
         print "Band does not match the current multiresolution transform"
         print "Update the MR transform first (self.UpdateSignatureFromObject)"
      elif band._PyMRObj__band < 0:
         raise ValueError("Band index is not valid: Is it a single band?")
      else:
         kband=band._PyMRObj__band
         if not self.__CppFlagCompute:
            if(kband >= 0 and kband < self.NbrBand):
               shapeBand=band.Coefs.shape
               if (shapeBand[0] == self.__CppMRObj.size_band_nl(kband)) and \
                        (shapeBand[1] == self.__CppMRObj.size_band_nc(kband)):
                  CppBand=Iflt(shapeBand[0], shapeBand[1])
                  CppBand.data[:]=band.Coefs[:]
                  self.__CppMRObj.insert_band(CppBand,kband)
               else:
                  raise ValueError("band should be of size {0} and {1} vs {2}"
                         " and {3}".format(self.__CppMRObj.size_band_nl(kband),
                         self.__CppMRObj.size_band_nc(kband),shapeBand[0],
                          shapeBand[1]))
            else:
               raise ValueError("band should be between 0 and {1}".format(0,
                                                             self.NbrBand-1))
         else:
            raise RuntimeError("Fist need to set the transform")

#------------------------------------------------------------------------------
#Simple I/O Routines: Direct interface - for C++ executable compatibility
   def read_fits(self,name):
      """! 
         Direct Wrapper to read MR fits obtain from the C++ executable
         @param name: name of the fitsio to read 
      """
      self.__CppMRObj.read(name)

   def write_fits(self,name):
      """! 
         Direct Wrapper to write MR fits as from the C++ executable
         @param name: name of the fitsio to write 
      """
      self.__CppMRObj.write(name)


#------------------------------------------------------------------------------
# --- Complex Getters (Rewritten in Python)

   def checkInputObject(self,InputObject):
      if not isinstance(InputObject, PyMRObj):
         raise TypeError("The input object is not a PyMRObj")
      if not InputObject.MRSignatureCompare(self.__signPyMRTrans):
         print "BEWARE: Detected a mismatch in between current MR signature"\
               "and Object signature."
         return False
      else:
         return True

   def UpdateSignatureFromObject(self,InputObject):
      if not isinstance(InputObject, PyMRObj):
         raise TypeError("The input object is not a PyMRObj")
      if not self.checkInputObject(InputObject):
         dc= InputObject.MRSignature()
         for k in self.__signPyMRTrans:
            self.__signPyMRTrans[k]=dc[k]
         self.__CppFlagCompute=True
         print "Cpp signature Updated and Cpp Object reset."
         self.create_CppMRObject()
         self.__CppFlagCompute=False

   def create_CppMRObject(self):
      """! 
         Initialize the Cpp multiresolution transform if necessary
      """
       #Create and Initialize all necessary instances if needed
      if self.__CppFlagCompute:
         if (self.Transform == type_transform.TO_MALLAT) or\
            (self.Transform == type_transform.TO_UNDECIMATED_MALLAT):
            self.__FAS.Verbose = self.Verbose
            self.__FAS.alloc(self.SBFilter)
         self.__CppMRObj.free()
         self.__CppMRObj.Verbose = self.Verbose
         self.__CppMRObj.alloc (self.nLines, self.nColumns, self. NbrPlan,
            self.Transform, self.__FAS, self.Norm, self. NbrUndec,self.UFilter)
         self.__CppMRObj.Border=self.Border
         if self.Verbose==Bool.values[1]:
            print 'Number of bands = {0}\n'.format(self.__CppMRObj.nbr_band())
         if(self.Transform == type_transform.TO_LIFTING):
            self.__CppMRObj.LiftingTrans = self.LiftingTrans
         self.__CppMRObj.Border = self.__signPyMRTrans["Border"]
         if self.Verbose==Bool.values[1]:
            print self.__CppMRObj.print_info()
         self.__CppFlagCompute=False

   def getAllCoefsDims(self):
      """! 
         Extract dimensions of the Wavelet coefficient array
         @return shape of wavelet coefficients
      """
      MR_Obj=self.__CppMRObj
      if self.__CppFlagCompute:
         raise RuntimeError("Fist need to set the transform")
      elif (MR_Obj.Set_Transform== set_transform.TRANSF_UNDECIMATED_MALLAT) or\
           (MR_Obj.Set_Transform == set_transform.TRANSF_DIADIC_MALLAT) or\
           (MR_Obj.Set_Transform == set_transform.TRANSF_PAVE) or\
           (MR_Obj.Set_Transform == set_transform. TRANSF_SEMIPYR):
         if(MR_Obj.Type_Transform  == type_transform.TC_FCT):
            naxes=numpy.zeros(1, numpy.float32)
            Nelem=1+MR_Obj.nbr_scale()
            for b in range(MR_Obj.nbr_band()):
               Nelem += 2 + MR_Obj.size_band_nl(b) * MR_Obj.size_band_nc(b)
            naxes[0] = Nelem
         else:
            naxes=numpy.zeros(3, numpy.float32)
            naxes[2] = MR_Obj.size_ima_nc()
            naxes[1] = MR_Obj.size_ima_nl()
            naxes[0] = MR_Obj.nbr_band()
      elif (MR_Obj.Set_Transform== set_transform.TRANSF_PYR):
         naxes=numpy.zeros(2, numpy.float32)
         naxes[1] = 2*MR_Obj.size_ima_nc()
         naxes[0] = 2*MR_Obj.size_ima_nl()
      elif (MR_Obj.Set_Transform== set_transform.TRANSF_MALLAT) or\
           (MR_Obj.Set_Transform== set_transform. TRANSF_FEAUVEAU):
         naxes=numpy.zeros(2, numpy.float32)
         naxes[1] = MR_Obj.size_ima_nc()
         naxes[0] = MR_Obj.size_ima_nl()
      else:
         raise TypeError("Error in gathering wavelet coeffs:" 
                         "bad Set_Transform")
      return naxes

   def getAllCoefs(self):
      """! 
         Return in the PyMRObj all wavelet coefficients stored
         @return a PyMRObj objet
      """
      #Get All coefficients
      MR_Obj=self.__CppMRObj
      if self.__CppFlagCompute:
         raise RuntimeError("Fist need to set the transform")
      naxes= self.getAllCoefsDims()
      wavecoeffs=numpy.zeros(naxes,numpy.float32)
      if (MR_Obj.Set_Transform == set_transform.TRANSF_DIADIC_MALLAT) or\
         (MR_Obj.Set_Transform == set_transform.TRANSF_PAVE):
         for kband in range(MR_Obj.nbr_band()):
            wavecoeffs[kband,:,:] = MR_Obj.band(kband).data[:]
      elif (MR_Obj.Set_Transform ==set_transform.TRANSF_UNDECIMATED_MALLAT) or\
           (MR_Obj.Set_Transform == set_transform.TRANSF_SEMIPYR):
         if(MR_Obj.Type_Transform  == type_transform.TC_FCT):
            wavecoeffs[0]=MR_Obj.nbr_scale()
            for kscale in range(MR_Obj.nbr_scale()):
               wavecoeffs[kscale+1]=MR_Obj.TabNrBandPerResol(kscale)
            
            offsetCoef=MR_Obj.nbr_scale()+1
            for kband in range(MR_Obj.nbr_band()):
               nCoefInBand=MR_Obj.nbr_coeff_in_band(kband)
               wavecoeffs[offsetCoef : offsetCoef + nCoefInBand] =numpy.ravel(
                                               MR_Obj.band(kband).data[:],'A')
               offsetCoef += nCoefInBand
         else:
            for kband in range(MR_Obj.nbr_band()):
               wavecoeffs[kband,0:MR_Obj.size_band_nl(kband),
                      0:MR_Obj.size_band_nc(kband)]=MR_Obj.band(kband).data[:]
      elif (MR_Obj.Set_Transform == set_transform.TRANSF_PYR):
         nCoefs=naxes[0]*naxes[1]
         lineStart=0
         columnStart=0
         for kband in range(MR_Obj.nbr_band()):
            wavecoeffs[lineStart:lineStart+MR_Obj.size_band_nl(kband),
                     columnStart:columnStart+MR_Obj.size_band_nc(kband)]=\
                     MR_Obj.band(kband).data[:]
            lineStart=lineStart+MR_Obj.size_band_nl(kband)
            ColumnStart=lineStart+MR_Obj.size_band_nc(kband)
      elif (MR_Obj.Set_Transform == set_transform.TRANSF_MALLAT) or\
           (MR_Obj.Set_Transform == set_transform. TRANSF_FEAUVEAU):
         Ima=Iflt(MR_Obj.size_ima_nl(), MR_Obj.size_ima_nc())
         ortho_trans_to_ima(MR_Obj,Ima)
         wavecoeffs[:]=Ima.data[:]
      PyObj=PyMRObj(self.Command)
      PyObj.Coefs=wavecoeffs
      PyObj.SetMRSignature(self.__signPyMRTrans)
      PyObj._PyMRObj__band=-1
      PyObj._PyMRObj__allcoeffs=True
      return PyObj

   def setAllCoefs(self,AllBands):
      """! 
         Set Wavelet Coefficients in a PyMRObj to the PyMRTrans
      """
      #Check if the current Cpp MR Object matches the one used to get the coefs
      #If not, assign the correct signature to the Cpp MR Object
      if not isinstance(AllBands, PyMRObj):
         raise TypeError("The input object is not a PyMRObj")
      if not (AllBands._PyMRObj__allcoeffs):
         raise ValueError("The object reported does not contain all coeffs")
      self.UpdateSignatureFromObject(AllBands)
      if self.__CppFlagCompute:
         self.create_CppMRObject()
      MR_Obj=self.__CppMRObj
      naxes= self.getAllCoefsDims()
      if any(AllBands.Coefs.shape[k] != naxes[k] for k in range(len(naxes))):
         raise ValueError("Incompatibilty in between the coefs size {0} and" 
             "expected size for transform {1}\n".format(AllBands.Coefs.shape, naxes))

      Coefs=AllBands.Coefs
      #Assign the coeffs
      if (MR_Obj.Set_Transform == set_transform.TRANSF_DIADIC_MALLAT) or\
         (MR_Obj.Set_Transform == set_transform.TRANSF_PAVE):
         band=Iflt(MR_Obj.size_ima_nl(),MR_Obj.size_ima_nc())
         for kband in range(MR_Obj.nbr_band()):
            band.data[:]= Coefs[kband,:,:]
            MR_Obj.insert_band(band,kband)
      elif (MR_Obj.Set_Transform ==set_transform.TRANSF_UNDECIMATED_MALLAT) or\
           (MR_Obj.Set_Transform == set_transform.TRANSF_SEMIPYR):
         if(MR_Obj.Type_Transform  == type_transform.TC_FCT):
            if MR_Obj.nbr_scale() != Coefs[0]:
               raise ValueError("Incompatibilty in between nbr scales {1} and" 
                  "expected {0}\n".format(MR_Obj.nbr_scale(), Coefs[0]))
            for kscale in range(MR_Obj.nbr_scale()):
               if TabNrBandPerResol(kscale)!= Coefs[1+kscale]:
                  raise ValueError("Incompatibilty in between nbr bands {0}"
                    " per scales {1} and expected {2}\n".format(
                    Coefs[1+kscale],kscale,TabNrBandPerResol(kscale)))
               
               offsetCoef=1+MR_Obj.nbr_scale()
               band=Iflt()
               for kband in range(MR_Obj.nbr_band()):
                  nCoefInBand=MR_Obj.nbr_coeff_in_band(kband)
                  band.alloc(nCoefInBand)
                  band.data[:]= Coefs[offsetCoef : offsetCoef + nCoefInBand]
                  MR_Obj.insert_band(band,kband)
                  offsetCoef += nCoefInBand
         else:
            band=Iflt()
            for kband in range(MR_Obj.nbr_band()):
               MR_Obj.band(kband).data[:]= Coefs[kband,0:MR_Obj.size_band_nl(
                                           kband),0:MR_Obj.size_band_nc(kband)]
      elif (MR_Obj.Set_Transform == set_transform.TRANSF_PYR):
         nCoefs=naxes[0]*naxes[1]
         lineStart=0
         columnStart=0
         for kband in range(MR_Obj.nbr_band()):
            MR_Obj.band(kband).data[:]= Coefs[lineStart:lineStart+\
                    MR_Obj.size_band_nl(kband), columnStart:columnStart+\
                    MR_Obj.size_band_nc(kband)]
            lineStart=lineStart+MR_Obj.size_band_nl(kband)
            ColumnStart=lineStart+MR_Obj.size_band_nc(kband)
      elif (MR_Obj.Set_Transform == set_transform.TRANSF_MALLAT) or\
           (MR_Obj.Set_Transform == set_transform. TRANSF_FEAUVEAU):
         Ima=Iflt(MR_Obj.size_ima_nl(), MR_Obj.size_ima_nc())
         Ima.data[:]= Coefs[:]
         ima_to_ortho_trans(MR_Obj,Ima)
	 
   #------------------------------------------------------------------------------
   # --- Routines to obtain inner coefficients (i.e. not affected by border effects)
   def getWaveInnerCoefsLim(self,shape_ima):
      """! 
          Return the inner coef indices limits in wavelet space 
	  NOTE: ONLY IMPLEMENTED FOR (BI)-ORTHOGONAL FILTERS WITH MALLAT TRANSFORM
          @param shape_ima shape of input image
          @return a 3D Array of type [[min:max],[line,column],[Band]] to be used for slicing
      """

      MR_Obj=self.__CppMRObj
      if self.__CppFlagCompute:
      	 self.create_CppMRObject()
      naxes= self.getAllCoefsDims()
      if (MR_Obj.Set_Transform == set_transform.TRANSF_MALLAT):
         sz_ima=numpy.array(shape_ima)
	 
         LFilt=MR_Obj.filter_bank_line()
         CFilt=MR_Obj.filter_bank_column()
	 sz_lin_bank=numpy.array([LFilt.size_analysis(),LFilt.size_synthesis()]) #[H0,G0] for line
         sz_col_bank=numpy.array([CFilt.size_analysis(),CFilt.size_synthesis()]) #[H0,G0] for column
         sbf_lin=SubBandFilter(LFilt)
         sbf_col=SubBandFilter(CFilt)
	 start_lin_filter_index=sz_lin_bank//2 #position of 0 from start in filter for lines
	 end_lin_filter_index=(sz_lin_bank-1)-start_lin_filter_index #position of 0 from end in filter for lines
	 start_col_filter_index=sz_col_bank//2 #position of 0 from start in filter for cols
	 end_col_filter_index=(sz_col_bank-1)-start_col_filter_index #position of 0 from end in filter for cols
	 start_h_pixel=numpy.array([0,0]) #start from 0 or 1 for h filters [hline,hcol]
	 start_g_pixel=numpy.array([0,0]) #start from 0 or 1 for g filters [gline,gcol]
	 #print("LIN",start_lin_filter_index,end_lin_filter_index)
	 #print("COL",start_col_filter_index,end_col_filter_index)
	 if not (sbf_lin.SubSample_H_Even):
	    start_h_pixel[0]=1
	 if not (sbf_col.SubSample_H_Even):
	    start_h_pixel[1]=1
	 if(sbf_lin.SubSample_G_Odd):
	    start_g_pixel[0]=1
	 if(sbf_col.SubSample_G_Odd):
	    start_g_pixel[1]=1
         #print("START_G=",start_g_pixel)
         #print("START_H=",start_h_pixel)
	 slice_inner_coefs=numpy.zeros((2,2,self.NbrBand)) #[[min:max],[line,column],scale]
	 lp_inner_coefs=numpy.array([[0,0],[sz_ima[0]-1,sz_ima[1]-1]])
	 for ksc in range(self.NbrPlan-1):
	    sz_ima=(sz_ima+1)//2
	    #H along lines followed by G along columns (Vertical in code)
	    slice_inner_coefs[0,0,3*ksc]=(start_lin_filter_index[0]-start_h_pixel[0]+lp_inner_coefs[0,0]+1)//2
	    slice_inner_coefs[1,0,3*ksc]=(lp_inner_coefs[1,0]-end_lin_filter_index[0]-start_h_pixel[0])//2+1
	    slice_inner_coefs[0,1,3*ksc]=(start_col_filter_index[1]-start_g_pixel[1]+lp_inner_coefs[0,1]+1)//2+sz_ima[1]
	    slice_inner_coefs[1,1,3*ksc]=(lp_inner_coefs[1,1]-end_col_filter_index[1]-start_g_pixel[1])//2+1+sz_ima[1]
	    #G along lines followed by H along columns (Vertical in the code)
	    slice_inner_coefs[0,0,3*ksc+1]=(start_lin_filter_index[1]-start_g_pixel[0]+lp_inner_coefs[0,0]+1)//2+sz_ima[0]
	    slice_inner_coefs[1,0,3*ksc+1]=(lp_inner_coefs[1,0]-end_lin_filter_index[1]-start_g_pixel[0])//2+1+sz_ima[0]
	    slice_inner_coefs[0,1,3*ksc+1]=(start_col_filter_index[0]-start_h_pixel[1]+lp_inner_coefs[0,1]+1)//2
	    slice_inner_coefs[1,1,3*ksc+1]=(lp_inner_coefs[1,1]-end_col_filter_index[0]-start_h_pixel[1])//2+1
	    #G along lines followed by G along columns (Diagonal in the code)
	    slice_inner_coefs[0,0,3*ksc+2]=(start_lin_filter_index[1]-start_g_pixel[0]+lp_inner_coefs[0,0]+1)//2+sz_ima[0]
	    slice_inner_coefs[1,0,3*ksc+2]=(lp_inner_coefs[1,0]-end_lin_filter_index[1]-start_g_pixel[0])//2+1+sz_ima[0]
	    slice_inner_coefs[0,1,3*ksc+2]=(start_col_filter_index[1]-start_g_pixel[1]+lp_inner_coefs[0,1]+1)//2+sz_ima[1]
	    slice_inner_coefs[1,1,3*ksc+2]=(lp_inner_coefs[1,1]-end_col_filter_index[1]-start_g_pixel[1])//2+1+sz_ima[1]
	    #H along lines followed by H along columns (Smooth)
	    lp_inner_coefs[0,0]=(start_lin_filter_index[0]-start_h_pixel[0]+lp_inner_coefs[0,0]+1)//2
	    lp_inner_coefs[1,0]=(lp_inner_coefs[1,0]-end_lin_filter_index[0]-start_h_pixel[0])//2
	    lp_inner_coefs[0,1]=(start_col_filter_index[0]-start_h_pixel[1]+lp_inner_coefs[0,1]+1)//2
	    lp_inner_coefs[1,1]=(lp_inner_coefs[1,1]-end_col_filter_index[0]-start_h_pixel[1])//2
	    #print(lp_inner_coefs)
	 lp_inner_coefs[1,0]+=1 #Due to Slice not reaching last number
	 lp_inner_coefs[1,1]+=1 #Due to Slice not reaching last number
         slice_inner_coefs[:,:,self.NbrBand-1]=lp_inner_coefs   
	    
      else:
         raise RuntimeError("Not Yet Implemented")
      
      return slice_inner_coefs
   
   def getWaveInnerCoefs(self,MRObj):
      """! 
          Return the inner coef of a MRObj
          @param MRObj an MRObj containing wavelet coefficients
          @return Inner wavelet coefficients 
	  Uses getWaveInnerCoefsLim()
      """
      ReInitFlag=False
      if not (self.checkInputObject(MRObj)):
         prev_signature=copy.deepcopy(self.__signPyMRTrans)
         self.UpdateSignatureFromObject(MRObj)
 	 ReInitFlag=True
      InnerLim=self.getWaveInnerCoefsLim((self.nLines,self.nColumns))
      LastBand=self.NbrBand-1
      Coefs=numpy.ndarray.flatten(MRObj.Coefs[InnerLim[0,0, LastBand]:InnerLim[1,0,LastBand],InnerLim[0,1,LastBand]:InnerLim[1,1,LastBand]])
      for kb in range(self.NbrBand-2,-1,-1):
          Coefs=numpy.concatenate((Coefs,numpy.ndarray.flatten(MRObj.Coefs[InnerLim[0,0, kb]:InnerLim[1,0, kb],InnerLim[0,1, kb]:InnerLim[1,1, kb]])))
      
      if(ReInitFlag):
         self.__signPyMRTrans=prev_signature
         self.__CppFlagCompute=True
         self.create_CppMRObject()
         self.__CppFlagCompute=False
      return Coefs

   def getWaveInnerCoefsMask(self,MRObj=None):
      """! 
          Return a mask where inner coefs are set to 1 and border coefs to 0
          @param MRObj an MRObj containing wavelet coefficients (optional)
          @return inner coefs mask
	  Uses getWaveInnerCoefsLim()
      """
      ReInitFlag=False
      if(MRObj is not None):
         if not (self.checkInputObject(MRObj)):
            prev_signature=copy.deepcopy(self.__signPyMRTrans)
            self.UpdateSignatureFromObject(MRObj)
 	    ReInitFlag=True
      naxes=self.getAllCoefsDims()
      InnerLim=self.getWaveInnerCoefsLim((self.nLines,self.nColumns))
      Mask=numpy.zeros(naxes)
      for kb in range(self.NbrBand):
          Mask[InnerLim[0,0, kb]:InnerLim[1,0, kb],InnerLim[0,1, kb]:InnerLim[1,1, kb]]=1
      
      if(ReInitFlag):
         self.__signPyMRTrans=prev_signature
         self.__CppFlagCompute=True
         self.create_CppMRObject()
         self.__CppFlagCompute=False
      return Mask
      
   def getWaveInnerCoefsList(self,MRObj=None,method=0):
      """! 
         Get list of Inner wavelet coefficients
         @param MRObj an MRObj containing wavelet coefficients (optional)
         @param method (default: 0 using mgrid, 1: numpy.indices, 2.: ogrid)
         @return 1D or 2D list of wavelet coefficient indices according to method used
      """
      ReInitFlag=False
      if(MRObj is not None):
         if not (self.checkInputObject(MRObj)):
            prev_signature=copy.deepcopy(self.__signPyMRTrans)
            self.UpdateSignatureFromObject(MRObj)
 	    ReInitFlag=True
	    
      if((method == 0) or (method >2)): #use mgrid
         InnerLim=numpy.int32(self.getWaveInnerCoefsLim((self.nLines,self.nColumns)))
         LastBand=self.NbrBand-1
         xx,yy=numpy.mgrid[InnerLim[0,0, LastBand]:InnerLim[1,0,LastBand],InnerLim[0,1,LastBand]:InnerLim[1,1,LastBand]]
         InnerList=numpy.ravel(xx*self.nColumns+yy)
         for kb in range(self.NbrBand-2,-1,-1):
            xx,yy=numpy.mgrid[InnerLim[0,0, kb]:InnerLim[1,0, kb],InnerLim[0,1, kb]:InnerLim[1,1, kb]]
	    InnerList=numpy.concatenate((InnerList,numpy.ravel(xx*self.nColumns+yy)))
      elif(method ==1): #use numpy.indices
         InnerLim=self.getWaveInnerCoefsLim((self.nLines,self.nColumns))
         LastBand=self.NbrBand-1
         row,col=numpy.indices((InnerLim[1,0,LastBand]-InnerLim[0,0, LastBand],InnerLim[1,1,LastBand]-InnerLim[0,1,LastBand]))
         row +=InnerLim[0,0, LastBand]
         col +=InnerLim[0,1,LastBand]
         row=numpy.ravel(row)
         col=numpy.ravel(col)
         for kb in range(self.NbrBand-2,-1,-1):
            rowt,colt=numpy.indices((InnerLim[1,0,kb]-InnerLim[0,0, kb],InnerLim[1,1,kb]-InnerLim[0,1,kb]))
            rowt +=InnerLim[0,0, kb]
            colt +=InnerLim[0,1,kb]
	    row=numpy.hstack((row,numpy.ravel(rowt)))
	    col=numpy.hstack((col,numpy.ravel(colt)))
         InnerList=(row,col)
      elif(method==2): #use ogrid
         InnerLim=numpy.int32(self.getWaveInnerCoefsLim((self.nLines,self.nColumns)))
         LastBand=self.NbrBand-1
         InnerList=numpy.ogrid[InnerLim[0,0, LastBand]:InnerLim[1,0,LastBand],InnerLim[0,1,LastBand]:InnerLim[1,1,LastBand]]
         for kb in range(self.NbrBand-2,-1,-1):
            InnerList=numpy.vstack((InnerList,numpy.ogrid[InnerLim[0,0, kb]:InnerLim[1,0, kb],InnerLim[0,1, kb]:InnerLim[1,1, kb]]))
      
      if(ReInitFlag):
         self.__signPyMRTrans=prev_signature
         self.__CppFlagCompute=True
         self.create_CppMRObject()
         self.__CppFlagCompute=False
      return InnerList

   def getWaveInnerCoefsFromList(self,List,MRObj,method=0):
      """! 
         Get Inner wavelet coefficients from list obtained with getWaveInnerCoefsList()
         @param List a list generated using getWaveInnerCoefsList(method=method)
         @param MRObj an MRObj containing wavelet coefficients 
         @param method (default: 0 using mgrid, 1: numpy.indices, 2.: ogrid)
         @return 1D or 2D list of wavelet coefficients according to method used
      """
     
      if((method == 0) or (method >2)): #use mgrid
         return numpy.ravel(MRObj.Coefs)[List] #1D output
      elif(method==1):
         return MRObj.Coefs[List[0],List[1]] #2D output
      elif(method==2):
         return numpy.hstack([numpy.ravel(MRObj.Coefs[List[ksc].tolist()]) for ksc in range(List.shape[0])]) #2D output
	 
	 
   #------------------------------------------------------------------------------
   # --- Forward and Backward Transforms
   def transform(self,image,nScale=None,nUndecScale=None,border=None):
      """! 
         Perform MultiResolution decomposition on an image using ISAP C++ code
         @param image a 2D image 
         @param nScale: number of decimated Scales
         @param nUndecScale: number of undecimated Scales
         @return all wavelet coefficients of the image
      """

      if nScale is None:
          nScale=self.NbrPlan
      if nUndecScale is None:
          nUndecScale =self.NbrUndec
      if border is None:
          border =self.Border

      #Get image size
      shapeImage=image.shape
      testND=len(shapeImage)
      if testND != 2:
         print '{0}{1} *** ERROR *** {2}\n'.format("MR_TRANSFORM", 
               time.strftime("%H:%M:%S", time.localtime()), 
               "Can only process 2d image")
      
      #Check if the CPP MR Object should be updated
      if (nScale != self. __signPyMRTrans["NbrPlan"]) or\
         (shapeImage[0] != self. __signPyMRTrans["nLines"]) or\
         (shapeImage[1] != self. __signPyMRTrans["nColumns"]) or\
         (nUndecScale !=self. __signPyMRTrans["NbrUndec"]) or\
         (border !=self. __signPyMRTrans["Border"]):
         self.__signPyMRTrans["NbrPlan"]= nScale
         self.__signPyMRTrans["nLines"]= shapeImage[0]
         self.__signPyMRTrans["nColumns"]= shapeImage[1]
         self.__signPyMRTrans["NbrUndec"]= nUndecScale
         self.__signPyMRTrans["Border"]= border
         self.__CppFlagCompute=True

      #Create and Initialize all necessary instances if needed
      if self.__CppFlagCompute:
         self.create_CppMRObject()

      #Allocate, Copy Data and do the transform
      Data=Iflt(shapeImage[0], shapeImage[1])
      Data.data[:]=image[:]
      self.__CppMRObj.transform(Data)
#      if self.Verbose==Bool.values[0]:
#         print self.__signPyMRTrans
      #Gather all wavelet Coefficients and the signature
      return self.getAllCoefs()

   def recons(self):
      """! 
         Perform MultiResolution reconstruction
         @return in a numpy array the synthesized image from the wavelet
         coefficients
      """
      Recons=Iflt(self.nLines,self.nColumns)
      self.__CppMRObj.recons(Recons,self.__signPyMRTrans["Border"])
      Data=numpy.zeros((self.nLines,self.nColumns), numpy.float32)
      Data[:]=Recons.data[:]
      return copy.deepcopy(Data)


