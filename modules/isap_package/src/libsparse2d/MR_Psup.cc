/******************************************************************************
**                   Copyright (C) 1995 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 14:46:34
**
**    Author: Jean-Luc Starck & C. Delattre
**
**    Date:  96/06/13
**    
**    File:  MR_Psup.cc
**
**
****************************************************************************/

 
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "IM_Obj.h"
#include "IM_Math.h"
#include "IM_IO.h"
#include "MR_Obj.h"
#include "MR_Abaque.h"
#include "MR_Noise.h"
#include "MR_NoiseModel.h"
#include "MR_Psupport.h"

/************************************************************************/

void building_imag_ascii(char *Name_Imag_In,Iint &Event_Image,Ifloat &Image)
{
 int c;
 FILE *input;
 int Nl,Nc,pixel_x,pixel_y,x_cur,y_cur,x_pos,y_pos,p,q;
 float x_lu, y_lu,distance_x,distance_y,energy_x,energy_y;
	

 /* Read the size of the input image */
 input = fopen(Name_Imag_In,"r");
 if (input == NULL) {
   cout << "Error in allocation of file " <<  Name_Imag_In << " ... or file doesn't exist" << endl;
   exit(-1);
 }
 fscanf(input,"%d\t%d\n",&Nl,&Nc);
 /* create the resulting image */
 Event_Image.resize(Nl,Nc);
 Event_Image.init(0);
 Image.resize(Nl,Nc);
 Image.init(0.0);

 /* Read the input file */
 while ( (c=getc(input)) != EOF )	
 {
 	ungetc(c,input);
	if (fscanf(input,"%f\t%f\n",&x_lu,&y_lu) == 2)
	{
	if ((x_lu < 0) || (x_lu >= Nl-0.5) || (y_lu < 0) || (y_lu >= Nc-0.5))
	{
	   cerr << "Error: incorrect event coordinates: " << x_lu << " " << y_lu << endl;
	   cerr << "       x coordinates must be between 0 and " << Nl-0.5 << endl;
	   cerr << "       y coordinates must be between 0 and " << Nc-0.5 << endl;
	   exit (-1);
	}

	/* determine the coordinates of the left low corner pixel */
	/* containing (x_lu,y_lu).                                */
	pixel_x = (int) floor(x_lu);
	pixel_y = (int) floor(y_lu);
	
	/* counting the number of events */
	(ABS(x_lu - pixel_x)>0.5) ? (p=pixel_x + 1) : (p=pixel_x);
	(ABS(y_lu - pixel_y)>0.5) ? (q=pixel_y + 1) : (q=pixel_y);
	Event_Image(p,q) ++;
	
	/* examine the neighbourhood of real pixel for computing */
	/* Bspline filtering                                     */
	for (x_cur=-2;x_cur<=2;x_cur++) 
	{
	   x_pos = x_cur + pixel_x;
	   distance_x = x_pos-x_lu;
	   if (distance_x <= 2.0)
 	     energy_x = (CUBE_FABS(distance_x - 2) + CUBE_FABS(distance_x + 2)
	          - 4 * (CUBE_FABS(distance_x - 1) + CUBE_FABS(distance_x + 1))
		  + 6 * CUBE_FABS(distance_x)) / 12.;
	   else  energy_x = 0.;
	   for (y_cur=-2;y_cur<=2;y_cur++) 
	   {
	        y_pos = y_cur + pixel_y;
		distance_y = y_pos-y_lu;
 		if (distance_y <= 2.0)
		{
 	          energy_y = (CUBE_FABS(distance_y - 2) + CUBE_FABS(distance_y + 2)
	          - 4 * (CUBE_FABS(distance_y - 1) + CUBE_FABS(distance_y + 1))
		  + 6 * CUBE_FABS(distance_y)) / 12.;  	           
 		}   
		else energy_y = 0.;
		if ((x_pos >= 0) && (y_pos >= 0) && (x_pos < Nl) && (y_pos < Nc)) 
		                       Image(x_pos,y_pos) += energy_x * energy_y ;
	   }   
	}   	
      }
    }

}
   




/************************************************************************/

void building_imag_imag(Ifloat &Image, Iint &Event_Image)
{	
 int Nl,Nc,i,j,x_cur,y_cur;
 float distance_x,distance_y,energy_x,energy_y;
 float Dxr=0., Dyr=0.;
 
 // init_random();
 Nl = Image.nl();
 Nc = Image.nc();
 Event_Image.resize(Nl,Nc);
 
//io_write_ima_float ("Bef.fits", Image);

 /* the initial image already containts the number of events */
 for (i=0;i<Nl;i++)
 for (j=0;j<Nc;j++) 
      Event_Image(i,j) = (int) (Image(i,j) + 0.5);
      
 Image.init();
 
 //io_write_ima_float ("Event.fits", Event_Image);
 /* Filter with Bspline */
 for (i=0;i<Nl;i++)
 for (j=0;j<Nc;j++) 
 {
    
    
    for (int e=0; e <  Event_Image(i,j); e++)
    {
       // Dxr = get_random(-0.5,0.5);
       // Dyr = get_random(-0.5,0.5);
       for (x_cur=-2;x_cur<=2;x_cur++)
        {
	   distance_x = x_cur + Dxr;
 	   energy_x = (CUBE_FABS(distance_x - 2) + CUBE_FABS(distance_x + 2)
	          - 4 * (CUBE_FABS(distance_x - 1) + CUBE_FABS(distance_x + 1))
		  + 6 * CUBE_FABS(distance_x)) / 12.;

	   for (y_cur=-2;y_cur<=2;y_cur++) 
           {
	      distance_y = y_cur+Dyr;
	      energy_y = (CUBE_FABS(distance_y - 2) + CUBE_FABS(distance_y + 2)
	          - 4 * (CUBE_FABS(distance_y - 1) + CUBE_FABS(distance_y + 1))
		  + 6 * CUBE_FABS(distance_y)) / 12.;
 	      if ((i+x_cur >= 0)    &&  (i+x_cur < Nl)
	        &&  (j+y_cur >= 0)    &&  (j+y_cur < Nc))
	                 Image(i+x_cur,j+y_cur)  += energy_x * energy_y;
	   }  
       }
    }
 }
 
//io_write_ima_float ("After.fits", Image);
   
}


/*************************************************************************/

void event_one_scale(Iint & Event_Image, int s, Iint &EventCount, type_border Border)
{
    int Nl,Nc,i,j,k;
    int Size;
    int Total=0;
    Nl = Event_Image.nl();
    Nc = Event_Image.nc();

    Size = (int) (pow((double)2., (double)(s+2)) + 0.5);
    // printf("Scale %d: size = %d\n", s+1, 2*Size+1);
    Nl = Event_Image.nl();
    Nc = Event_Image.nc();

    for (i =0; i < Nl; i++)
    {
        Total = 0;
        for (k =-Size; k <= Size; k++)
        for (j =-Size; j <= Size; j++)
           Total += Event_Image(i+k, j, Border);

        EventCount(i,0) = Total;
        for (j =1; j < Nc; j++)
        {
            for (k = -Size; k <= Size; k++)
            {
                Total -= (int) Event_Image(i+k,j-Size-1,Border);
                Total += (int) Event_Image(i+k,j+Size, Border);
            }
            EventCount(i,j) = Total;
        }
    }
}

/*************************************************************************/

void reducing_coeff(const MultiResol &MR_Data_Event, MultiResol &MR_Data)
{
   int s,i,j,Nl,Nc;
   float sigma=0.040717,alpha;
   int N_Scale = MR_Data_Event.nbr_scale()-1;

   Nl = MR_Data_Event.size_ima_nl();
   Nc = MR_Data_Event.size_ima_nc();

   for (s=0,alpha=1.0;s<N_Scale;s++,alpha *= 4)
      for (i=0;i<Nl;i++)
      for (j=0;j<Nc;j++)
         if (MR_Data_Event(s+1,i,j) > 0.5) 
            MR_Data(s,i,j) = MR_Data(s,i,j) * alpha 
                             / (sqrt(MR_Data_Event(s+1,i,j)) * sigma);
         else MR_Data(s,i,j) = 0.;
}

/**************************************************************************/

void event_set_support(const MultiResol &MR_Data, int CurrentScale, 
                      Iint &Event_Image, type_border Border,
                      const Ifloat &I_Abaque, MRNoiseModel & Model)
{
   int i,j,s,Nl,Nc,n_event_real,n_event,power;
   float seuil_max,seuil_min, Sm;
   int ABAQUE_N_SCALE = I_Abaque.nl();
   //float sigma=0.040717,alpha;
   float sigma=0.0405078,alpha;
   Nl = MR_Data.size_ima_nl();
   Nc = MR_Data.size_ima_nc();
   Iint EventCount(Nl,Nc, "EventCount");

      for (s=0,alpha=1.0;s<CurrentScale;s++) alpha *= 4;

      s = CurrentScale;
      event_one_scale(Event_Image, s, EventCount, Border);
      
//char File[80];
//sprintf (File, "EvImOld_%d.fits", CurrentScale);
//io_write_ima_float (File, Event_Image);
//sprintf (File, "EvCntOld_%d.fits", CurrentScale);
//io_write_ima_float (File,EventCount );

      for (i=0;i<Nl;i++)
      for (j=0;j<Nc;j++) 
      {
   	  /* reading the number of events */
	  n_event_real = EventCount(i,j);
	  
	    /* compute "power" : 2^(power-1) < n_event_real <= 2^(power) */
	    power = 0;
	    n_event = 1;
	    while (n_event<n_event_real) 
	    {
	       power ++;
	       n_event *= 2;
	    }

	    /* interpolation of thresholds using abaque */
	    if (power > ABAQUE_N_SCALE) 
	    {
		power = ABAQUE_N_SCALE;
		seuil_min = I_Abaque(power,0);
		seuil_max = I_Abaque(power,1);
//seuil_min *= Model.CFewEventPoisson2d->_HistoSigma(power);		
//seuil_max *= Model.CFewEventPoisson2d->_HistoSigma(power);		
	    }
	    else if (n_event_real == n_event) 
	    {
	    
//if (power > 0) 
//cout << "Sigma computed local(Event:" << n_event_real << ",Scale:" 
//     <<  CurrentScale << ") : " <<  sqrt((float) n_event_real) * sigma / alpha
//     << "Sigma in Histo : " 
//     << Model.CFewEventPoisson2d->_HistoSigma(power) << endl;
     
		seuil_min = I_Abaque(power,0);
		seuil_max = I_Abaque(power,1);
		
//seuil_min *= Model.CFewEventPoisson2d->_HistoSigma(power);		
//seuil_max *= Model.CFewEventPoisson2d->_HistoSigma(power);			
	    }
	    else 
	    {
	     	if (power > 0)
	     	{
	     	   seuil_min = -2*(I_Abaque(power,0) - I_Abaque(power-1,0))
                          *(n_event - n_event_real)/n_event + I_Abaque(power,0);
		   seuil_max = -2*(I_Abaque(power,1) - I_Abaque(power-1,1))
                          *(n_event -n_event_real)/n_event + I_Abaque(power,1);
	     	
//float InterpSig = -2*(   Model.CFewEventPoisson2d->_HistoSigma(power) 
//                       - Model.CFewEventPoisson2d->_HistoSigma(power-1))
//                    *(n_event - n_event_real)/n_event 
//		    + Model.CFewEventPoisson2d->_HistoSigma(power);	
//seuil_min *= InterpSig;
//seuil_max *= InterpSig;		
      		
		}
	     	else
	     	{
	     	   seuil_min = I_Abaque(0,0);
	     	   seuil_max = I_Abaque(0,1);
		   
//seuil_min *= Model.CFewEventPoisson2d->_HistoSigma(0);		
//seuil_max *= Model.CFewEventPoisson2d->_HistoSigma(0);			   
	     	}
	    }
	    
	    
	    if (power > 0)
	    {
	       seuil_min *= sqrt((float) n_event_real) * sigma / alpha;
	       seuil_max *= sqrt((float) n_event_real) * sigma / alpha;
	    }
	    else
	    {
	       seuil_min *=  sigma / alpha;
	       seuil_max *=  sigma / alpha;
	    }
	    
	    

	    
	    
	     /* detection of signal : comparison with thresholds */
	  Model.support(s,i,j) = VAL_SupNull;
	  if ((MR_Data(s,i,j) <= seuil_min) || (MR_Data(s,i,j) >= seuil_max))
	  {
	     Model.support(s,i,j) = VAL_SupOK;
	     if (n_event_real < Model.MinEventNumber) 
                                 Model.support(s,i,j) = VAL_SupMinEv;
	     if ((Model.OnlyPositivDetect == True) && (MR_Data(s,i,j) < 0))
	                         Model.support(s,i,j) = VAL_SupNull;
	     if (Model.FirstDectectScale > s) 
	                         Model.support(s,i,j) = VAL_SupFirstScale;
	  }
	  
	  Sm = (ABS(seuil_min) > ABS(seuil_max)) ? ABS(seuil_min) : ABS(seuil_max);
	  Model.sigma(s,i,j) = Sm / Model.NSigma[s];
	  // Model.NSigma[s] = 1.;
	} /* end FOR (j=0 ...) */
}

/********************************************************************/

void mr_psupport(Iint &Event_Image, MultiResol &MR_Data, 
                 Ifloat &I_Abaque, MRNoiseModel & Model, 
                 type_border Border, Bool WriteAll)
{
 int s,Nc,Nl;
 Nl = Event_Image.nl();
 Nc = Event_Image.nc();
 int NScale = MR_Data.nbr_scale();

 /* writing on disk the coefficients of wavelet */
  if (WriteAll == True)  MR_Data.write(NAME_WAVELET_SUP);


 /* detecting signal : comparison with the ABAQUE of thresholds */
 cout << "Detect signal ... " << endl;

 for (s=0;s< NScale-1; s++)
 {
    event_set_support(MR_Data, s,Event_Image,  Border, I_Abaque, Model);
 }
}


/********************************************************************/

void mr_psupport(MultiResol &MR_Data, MRNoiseModel & Model, type_border Border)
{
 int i,j,s;
 int NScale = MR_Data.nbr_scale();
 Iint EventCount;
 int Nl = MR_Data.size_ima_nl();
 int Nc = MR_Data.size_ima_nc();   
	
 if (Model.TypeThreshold == T_FDR) EventCount.alloc(Nl,Nc,"ImagCount");

 for (s=0;s< NScale-1; s++)
 {
    //Model.CFewEventPoisson2dCEventPois->find_threshold((Model.TabEps)[s]);
    //event_set_support(MR_Data, s, Model.Event_Image,  Border, 
    //                  Model.CEventPois->Threshold, Model);
    if (Model.TypeThreshold != T_FDR)
    {
       if( Model.which_noise() == NOISE_EVENT_POISSON ) {
          if( Model.old_poisson() == True ){
	     Model.CFewEventPoisson2d->find_threshold((Model.TabEps)[s]);
	     event_set_support(MR_Data, s, Model.Event_Image,  Border, 
                        Model.CFewEventPoisson2d->_Threshold, Model);
	  } else {
	     Model.CFewEvent2d->find_threshold((Model.TabEps)[s]);
	     event_set_support(MR_Data, s, Model.Event_Image,  Border, 
                        Model.CFewEvent2d->mThreshold, Model);
          }
       }
    }
    else
    {
        //cout << "event_one_scale" << Model.Event_Image.nl() << " " <<  Model.Event_Image.nc() <<endl;
        //cout << "EventCount, " << EventCount.nl() << " " <<  EventCount.nc() <<endl;
	
	event_one_scale(Model.Event_Image, s, EventCount, MR_Data.Border);
        dblarray TabPValue(Nc,Nl);
        // cout << "CFewEventPoisson2d loop" << endl;
	for (i = 0; i < Nl; i++)
        for (j = 0; j < Nc; j++) 
	{
            
            double P =0;
            if( Model.which_noise() == NOISE_EVENT_POISSON ) {
	       if( Model.old_poisson() == True )
                  P = (Model.CFewEventPoisson2d)->a_trou_repartition(
	                            MR_Data(s,i,j),  EventCount(i,j) , s );
               else
	          P = (Model.CFewEvent2d)->a_trou_repartition(
	                            MR_Data(s,i,j),  EventCount(i,j) , s );
            }
             
	    if (MR_Data(s,i,j) > 0) P = 1. - P;
	    TabPValue(j,i) = P;				    
	}
        float Alpha = (Model.TabEps)[s];
	// cout << "get fdrs" << endl;
        double PDet = fdr_pvalue(TabPValue.buffer(), Nl*Nc, Alpha);
        (Model.NSigma)[s] = ABS(xerfc(0.5+(1-PDet)/2.));
	
        for (i = 0; i < Nl; i++)
        for (j = 0; j < Nc; j++) 
	{
	   
	   if (TabPValue(j,i) < PDet) 
	   {
	      (Model.support)(s,i,j) = VAL_SupOK;
	      Model.sigma(s,i,j) = ABS(MR_Data(s,i,j)) / Model.NSigma[s];
	      
	      if (EventCount(i,j) < Model.MinEventNumber)  (Model.support)(s,i,j) = VAL_SupMinEv;
	      if ((Model.OnlyPositivDetect == True) && (MR_Data(s,i,j) < 0))
	                                                    (Model.support)(s,i,j) = VAL_SupNull;
	      if (Model.FirstDectectScale > s) (Model.support)(s,i,j) = VAL_SupFirstScale;
	   }
	   else 
	   {
	      (Model.support)(s,i,j) = VAL_SupNull;
	       Model.sigma(s,i,j) =  ABS(MR_Data(s,i,j));
	   }
  	}
 	
        printf("FDR: band %d ==> Alpha = %f, PDet= %f, NSigma = %f\n", s+1, Alpha, PDet, (Model.NSigma)[s]);
    }
    // cout << "event_set_support" << endl;
    
    
   // cout << "end event_set_support" << endl;
 }
 // cout << "OUF mr_psupport" << endl;
}

/********************************************************************/


