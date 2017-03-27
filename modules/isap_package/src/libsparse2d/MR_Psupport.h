/******************************************************************************
**                   Copyright (C) 1995 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 3.2
**
**    Author: Jean-Luc Starck
**
**    Date:  96/06/13
**    
**    File:  MR_Psupport.h
**
*******************************************************************************
**
**    DESCRIPTION  
**    ----------- 
**                 
**    PARAMETRES    
**    ----------    
** 
**    RESULTS      
**    -------  
**
**
******************************************************************************/


#ifndef __PSUPPORT__
#define __PSUPPORT__

#include"MR_Sigma.h"

#define		CUBE_FABS(x) (fabs(x) * fabs(x) * fabs(x))

#define DEF_BORDER I_CONT
#define DEF_N_SCALE 6  /* number of scales */
#define NAME_EVENT_SUP "xx_Event"
#define NAME_WAVELET_SUP "xx_Wavelet"
#define NAME_COEF_RED "xx_Reduced"
#define NAME_SUPPORT_RED "xx_Support_Red"
#define NAME_SUPPORT "xx_Support"

#define VAL(Name,i,j) ( ((i<Nl) && (i>=0) && (j>=0) && (j<Nc))?Name(i,j):0.0 )
#define VAL0(Name,scale,i,j) ( ((i<Nl) && (i>=0) && (j>=0) && (j<Nc))?Name(scale,i,j):0.0 )

/* prototype of functions */
extern void building_imag_ascii (char *Name_Imag_In, Iint &Event_Image, 
                                 Ifloat &Image);
extern void building_imag_imag (Ifloat &Image, Iint &Event_Image);
extern void counting_event(Iint & Event_Image, MultiResol &MR_Data_Event,
                           type_border Border);
extern void reducing_coeff(const MultiResol &MR_Data_Event, 
                           MultiResol &MR_Data);

extern void event_set_support(const MultiResol &MR_Data, int CurrentScale, 
                      Iint &Event_Image, type_border Border,
                      const Ifloat &I_Abaque, MRNoiseModel & Model);

extern void mr_psupport(Iint &Event_Image, MultiResol &MR_Data,
                        Ifloat &I_Abaque, MRNoiseModel &Model, 
                        type_border Border, Bool WriteAll);

extern void mr_psupport(MultiResol &MR_Data, MRNoiseModel & Model, 
                        type_border Border);

void event_one_scale(Iint & Event_Image, int s, Iint &EventCount, type_border Border);

#endif
