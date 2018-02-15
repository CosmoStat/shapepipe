/******************************************************************************
**                   Copyright (C) 1995 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 3.1
**
**    Author: Jean-Luc Starck
**
**    Date:  96/05/02 
**    
**    File:  MR_Filter.h
**
*******************************************************************************
**
**    DESCRIPTION  
**    ----------- 
**                 
**    PARAMETRES    
**    ----------    
** 
** 
**    RESULTS      
**    -------  
**
**
******************************************************************************/


#ifndef __MWFILTER__
#define  __MWFILTER__

const int DEF_MEM_ALPHA_CST = 1;
const int DEF_MEM_ALPHA_OPT = 2;
const int DEF_MEM_ALPHA_OPT_SCALE = 3;
const float DEF_MAX_MEM_ALPHA = 100.0;

const int DEF_MEM_FILT_MAX_ITER = 20;
const float DEF_MEM_FILT_CONGER = 0.01;

/****************************************************************************/ 

void mw_filter(MultiResol &MR_Data, MRNoiseModel & NoiseModel, 
        int TypeFilter, fltarray &TabAlpha, Bool DataSNR=False, 
	Bool DataModel=False, MultiResol *MR_Model=NULL,
	float ConvgParam = DEF_MEM_FILT_CONGER, int MaxIter=DEF_MEM_FILT_MAX_ITER, 
	Bool PositivConst=False, Bool Verbose=False, Bool UseSupport=False);
/* Apply the multiscale entropy to the multiresolution data set MR_Data
   MR_Data = in-out: multiresolution transform of an image
   NoiseModel = in: noise modeling of the image in the wavelet space
   TypeFilter = in: type of filtering
             = DEF_MEM_ALPHA_CST  ==> alpha is fixed 
	     = DEF_MEM_ALPHA_OPT  ==> alpha is globally optimized
	     = DEF_MEM_ALPHA_OPT_SCALE ==> alpha is optimized at each scale
	     
   TabAlpha = in-out: TabAlpha(b) gives the alpha regularization parameter for
                  the band b
   DataModel = in: true if we have a model
   MR_Model = in: pointer to a multiresolution model or NULL (if none)
   DataSNR  = in: true if alpha parameter are modified using the SNR
                    of the data
   PositivConst = in: apply the positivy constraint
   UseSupport= in: true if alpha parameter are modified using the multiresolution
                        support
  Warning: works well, when each wavelet can be modelized by a Gaussian		    
*/	

/****************************************************************************/ 
 
void mw_fm1(MultiResol & MR_Data, MRNoiseModel & NoiseModel, 
            fltarray &TabAlpha,  CMemWave & MemWObj, 
	    Bool UseModel,  MultiResol *MR_Model, Bool UseDataSNR, Bool UseSupport=False);
/* Apply the multiscale entropy to the multiresolution data set MR_Data
   MR_Data = in-out: multiresolution transform of an image
   NoiseModel = in: noise modeling of the image in the wavelet space
   TabAlpha = in: TabAlpha(b) gives the alpha regularization parameter for
                  the band b
   MemWObj = in: class for entropy filtering
   UseModel = in: true if we have a model
   MR_Model = in: pointer to a multiresolution model or NULL (if none)
   UseDataSNR = in: true if alpha parameter are modified using the SNR
                    of the data
   UseSupport= in: true if alpha parameter are modified using the multiresolution
                        support

*/		    	


/****************************************************************************/ 

void mw_fm2(MultiResol &MR_Data, MRNoiseModel & NoiseModel, 
        fltarray &TabAlpha, CMemWave & MemWObj, Bool DataSNR, Bool DataModel, MultiResol *MR_Model,
	float ConvgParam, int MaxIter, Bool PositivConst=False, Bool Verbose=False, Bool UseSupport=False);
/* Apply the multiscale entropy to the multiresolution data set MR_Data
   The Parameter Alpha is estimated iteratively (by dichotomy) in order to have
   a residu compatible with the noise.
   
   MR_Data = in-out: multiresolution transform of an image
   NoiseModel = in: noise modeling of the image in the wavelet space
   TabAlpha = in: TabAlpha(b) gives the alpha regularization parameter for
                  the band b
   MemWObj = in: class for entropy filtering
   UseModel = in: true if we have a model
   MR_Model = in: pointer to a multiresolution model or NULL (if none)
   UseDataSNR = in: true if alpha parameter are modified using the SNR
                    of the data
   ConvgParam = in: Convergence parameter
   MaxIter = in: maximum number of iterations
   UseSupport= in: true if alpha parameter are modified using the multiresolution
                        support
*/	

/****************************************************************************/ 


void mw_fm3(MultiResol &MR_Data, MRNoiseModel & NoiseModel, 
        fltarray &TabAlpha, CMemWave & MemWObj, Bool DataSNR, Bool DataModel, MultiResol *MR_Model,
	float ConvgParam, int MaxIter, Bool Verbose=False,
	Bool RecEachIter=False, Bool PosConstraint=False, Bool UseSupport=False);
/* Apply the multiscale entropy to the multiresolution data set MR_Data
   The Parameter Alpha is estimated iteratively at scale separately
   (by dichotomy) in order to have a residu compatible with the noise
   at all the scales
   
   MR_Data = in-out: multiresolution transform of an image
   NoiseModel = in: noise modeling of the image in the wavelet space
   TabAlpha = in: TabAlpha(b) gives the alpha regularization parameter for
                  the band b
   MemWObj = in: class for entropy filtering
   UseModel = in: true if we have a model
   MR_Model = in: pointer to a multiresolution model or NULL (if none)
   UseDataSNR = in: true if alpha parameter are modified using the SNR
                    of the data
   ConvgParam = in: Convergence parameter
   MaxIter = in: maximum number of iterations
   RecEachIter = in: if true, a reconstruction, followed by a transform
                     is performed at each iteration. This allows to resolve
		     the problem of non-orthogonality of the coefficient.
		     Set this option for orthogonal transform has no sense.
  PosConstraint: in: if set and if RecEachIter is set, the positivity contraint
                     is applied.
   UseSupport= in: true if alpha parameter are modified using the multiresolution
                        support
*/	

/****************************************************************************/ 

#endif




