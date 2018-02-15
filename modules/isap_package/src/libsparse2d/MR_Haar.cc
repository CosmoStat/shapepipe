/*****************************************************************
**
**    UNIT
**
**    Version: 3.1
**
**    Author: Jean-Luc Starck
**
**    Date:  96/05/02 
**    
**    File:  MR_Haar.cc
**
******************************************************************
**
**    DESCRIPTION  routines used for the haar wavelet transform algorithm 
**    -----------  
**
******************************************************************
** 
** void haar_2d_transform (Ifloat &Data, MultiResol &T_Haar)
** othogonal Haar transform
**
******************************************************************
** 
** void haar_2d_reconstruct (MultiResol & T_Haar, Ifloat &Data)
** othogonal Haar reconstruction
**
******************************************************************
**
** void haar_2d_atrou_dir2_transform (Ifloat &Data, MultiResol &T_Haar, Bool EdgeLineTransform)
** dyadic haar transform (no decimation and two bands per scale)
** if EdgeLineTransform == True then all line are line and all columns are
**                         wavelet transformed a second time
**
******************************************************************
**
** void haar_2d_atrou_dir2_recons (MultiResol & T_Haar, Ifloat &Data, Bool EdgeLineTransform)
** dyadic haar reconstruction  
** if EdgeLineTransform == True then all line are line and all columns are
**                         wavelet transformed a second time
**
*****************************************************************/ 

// static char sccsid[] = "@(#)MR_Mallat.cc 3.1 96/05/02 CEA 1994 @(#)";


#include <stdio.h>
#include <math.h>
#include <string.h>

#include "MR_Obj.h"

/****************************************************************************/

void haar_2d_reconstruct (MultiResol & T_Haar, Ifloat &Data)
{
    int Nbr_Scale = T_Haar.nbr_scale();
    int i,j,s,Nls,Ncs,Nl1,Nc1;
    int Nl = Data.nl();
    int Nc = Data.nc();
    Ifloat Buff(Nl, Nc, "Buff");
    float Sqrt2 = 2.; // sqrt(2.);
    
    s = Nbr_Scale-1;
 //   int b=T_Haar.scale_to_band(s,I_SMOOTH);
 //    cout << "stb(s) = " << b << endl;
//     cout << "b1(s) = " <<   T_Haar.size_band_nl(b) << endl;
//     cout << "b2(s) = " <<  T_Haar.size_band_nc(b) << endl;
    
    s = Nbr_Scale-1;
    Nls = size_ima_resol(s, Nl);
    Ncs = size_ima_resol(s, Nc);
    Buff.resize(Nls,Ncs);
    for (i = 0; i < Nls; i++)
    for (j = 0; j < Ncs; j++) Buff(i,j) = T_Haar(s,i,j,I_SMOOTH);
    
    for (s=Nbr_Scale-2; s >= 0; s--)
    {
        Nls = size_ima_resol(s, Nl);
        Ncs = size_ima_resol(s, Nc);    
           
        // cout << "Haar: Nls = " << Nls << "  Ncs = " << Ncs  << endl;
        Data.resize(Nls, Ncs);
        Nl1 = (Nls % 2 == 0) ? Nls : Nls - 1;  
        Nc1 = (Ncs % 2 == 0) ? Ncs : Ncs - 1;
  
        for (i = 0; i < Nl1; i++)
        for (j = 0; j < Nc1; j++)
        {
           int ii = i / 2;
           int jj = j / 2;

           if ((i % 2 == 0) && (j % 2 == 0))
           {
              
              Data(i,j)=(Buff(ii,jj) -
                      T_Haar(s,ii,jj,D_VERTICAL) -
                      T_Haar(s,ii,jj,D_HORIZONTAL) +
                      T_Haar(s,ii,jj,D_DIAGONAL)) / 2.;
           }
           else if ((i % 2 == 0) && (j % 2 == 1))
           {
              Data(i,j)=( Buff(ii,jj) -
                      T_Haar(s,ii,jj,D_VERTICAL) +
                      T_Haar(s,ii,jj, D_HORIZONTAL) -
                      T_Haar(s,ii,jj,D_DIAGONAL)) / 2.;
           }
           else if ((i % 2 == 1) && (j % 2 == 0))
           {
              Data(i,j)=(Buff(ii,jj) +
                      T_Haar(s,ii,jj,D_VERTICAL) -
                      T_Haar(s,ii,jj,D_HORIZONTAL) -
                      T_Haar(s,ii,jj,D_DIAGONAL)) / 2.;
           }
           else
           {
              Data(i,j)=(Buff(ii,jj) +
                      T_Haar(s,ii,jj,D_HORIZONTAL) +
                      T_Haar(s,ii,jj,D_VERTICAL) +
                      T_Haar(s,ii,jj,D_DIAGONAL)) / 2.;
           }
       }
       if (Nl1 != Nls)
       {
           i = Nls-1;
           for (j = 0; j < Nc1; j++) 
             if (j % 2 == 0)
                Data(i,j) = (Buff(i/2,j/2) - T_Haar(s,i/2,j/2,D_HORIZONTAL)) / Sqrt2;
             else
                Data(i,j) = (Buff(i/2,j/2) + T_Haar(s,i/2,j/2,D_HORIZONTAL)) / Sqrt2;
              
       }
       if (Nc1 != Ncs)
       {
            j = Ncs-1;
            for (i = 0; i < Nl1; i++)   
             if (i % 2 == 0) 
                 Data(i,j) = (Buff(i/2,j/2) - T_Haar(s,i/2,j/2,D_VERTICAL)) / Sqrt2;
             else Data(i,j) = (Buff(i/2,j/2) + T_Haar(s,i/2,j/2,D_VERTICAL)) / Sqrt2;
        }  
       if ((Nl1 != Nls) && (Nc1 != Ncs))
       {
            i = Nls-1;
            j = Ncs-1;
            Data(i,j) = Buff(i/2,j/2);
       }       
 
       
       
       
       
       if (s != 0) Buff = Data;
   }
}
 
/****************************************************************************/

void haar_2d_transform (Ifloat &Data, MultiResol &T_Haar)
{
    int Nl = Data.nl();
    int Nc = Data.nc();
    int Nbr_Scale = T_Haar.nbr_scale();
    int b,i,j,s,Nls=Nl,Ncs=Nc,Nlr,Ncr;
    Ifloat Buff(Nl,Nc, "haar buffer");
    Ifloat Aux( (Nl+1)/2, (Nc+1)/2, "haar buffer");
    Buff = Data;
    float Sqrt2 = 1.; // sqrt(2.);
    
    for (s=0; s < Nbr_Scale - 1; s++)
    { 
         int Nl1,Nc1;
         Nlr = size_ima_resol(s, Nl);
         Ncr = size_ima_resol(s, Nc);
         
         b = T_Haar.scale_to_band(s,D_HORIZONTAL); 
         Nls = T_Haar.size_band_nl(b);
         Ncs = T_Haar.size_band_nc(b);
         Nl1 = (Nlr % 2 == 0) ? Nls : Nls - 1;  
         Nc1 = (Ncr % 2 == 0) ? Ncs : Ncs - 1;
        // horizontal band
        for (i = 0; i < Nl1; i++)
        for (j = 0; j < Ncs; j++)
        {
           T_Haar(b,i,j) =  (-Buff(2*i,2*j)+ Buff(2*i,2*j+1)-
                 Buff(2*i+1,2*j) + Buff(2*i+1,2*j+1)) / 2.; 
         
        }
        // test non square image
        // the sigma must stay equal to 1 ==> normalization by Sqrt2
        if (Nl1 != Nls)
        {
           i = Nls-1;
           for (j = 0; j < Ncs; j++) 
            T_Haar(b,i,j) = (-Buff(2*i,2*j) + Buff(2*i,2*j+1))/ Sqrt2;  
        }        

        
        // vertical band
        b = T_Haar.scale_to_band(s,D_VERTICAL);
        Nls = T_Haar.size_band_nl(b);
        Ncs = T_Haar.size_band_nc(b);   
        Nl1 = (Nlr % 2 == 0) ? Nls : Nls - 1;  
        Nc1 = (Ncr % 2 == 0) ? Ncs : Ncs - 1;
        for (i = 0; i <  Nls; i++)
        for (j = 0; j <  Nc1; j++)
        {
           T_Haar(b,i,j) = (-Buff(2*i,2*j)- Buff(2*i,2*j+1)+
                 Buff(2*i+1,2*j) + Buff(2*i+1,2*j+1)) / 2.;
        }
        // test non square image
        // the sigma must stay equal to 1 ==> normalization by Sqrt2
       if (Nc1 != Ncs)
       {
          j = Ncs-1;
          for (i = 0; i < Nls; i++) 
              T_Haar(b,i,j) = (-Buff(2*i,2*j) + Buff(2*i+1,2*j))/ Sqrt2;
       }
           
        // diagonal band
        b = T_Haar.scale_to_band(s,D_DIAGONAL);
        for (i = 0; i <  T_Haar.size_band_nl(b); i++)
        for (j = 0; j <  T_Haar.size_band_nc(b); j++)
        {
            T_Haar(b,i,j) =  (Buff(2*i,2*j)- Buff(2*i,2*j+1)-
                 Buff(2*i+1,2*j) + Buff(2*i+1,2*j+1)) / 2.; 
         }
      
        // smooth component
        Nls =  T_Haar.size_scale_nl(s,I_SMOOTH);
        Ncs =  T_Haar.size_scale_nc(s,I_SMOOTH);
        Nl1 = (Nlr % 2 == 0) ? Nls : Nls - 1;  
        Nc1 = (Ncr % 2 == 0) ? Ncs : Ncs - 1;
        Aux.resize(Nls,Ncs);
        for (i = 0; i < Nl1; i++)
        for (j = 0; j < Nc1; j++) 
           Aux(i,j) = (Buff(2*i,2*j)+ Buff(2*i,2*j+1)+
                 Buff(2*i+1,2*j) + Buff(2*i+1,2*j+1)) / 2.;
                 
        if (Nl1 != Nls)
        {
            i = Nls-1;
            for (j = 0; j < Nc1; j++) 
                 Aux(i,j) = (Buff(2*i,2*j) + Buff(2*i,2*j+1));  
         }
         if (Nc1 != Ncs)
         {
            j = Ncs-1;
            for (i = 0; i < Nl1; i++)    
                Aux(i,j) = (Buff(2*i,2*j) + Buff(2*i+1,2*j));  
         }  
         if ((Nl1 != Nls) && (Nc1 != Ncs))
         {
            i = Nls-1;
            j = Ncs-1;
            Aux(i,j) = Buff(2*i,2*j);
         }
         
        if (s != Nbr_Scale-2) Buff = Aux;
        else
        {
           b = T_Haar.scale_to_band(s,I_SMOOTH);
           for (i = 0; i < Nls; i++)
           for (j = 0; j < Ncs; j++) T_Haar(b,i,j) = Aux(i,j);
        }
    }
}   

/****************************************************************************/
 

/*************************************************************************/

static void haar_1d_atrou_convol(int N, float *Input, float *OutH, 
                          float *OutG, int Step, double NormVal)
/* convolves the data with the filter h1 */
{
    int s,i, Index;
    int Ind =0;
    
    for (i = 0; i < N; i += 2*Step)
    for (s = i; s < MIN(i+Step,N); s ++)
    {
       Index =  (s+Step < N) ? s+Step : s;
       OutH[Ind] = (float) (NormVal*(Input[s] + Input[Index]));
       OutG[Ind] = (float) (NormVal*(Input[s] - Input[Index]));
       Ind++;
    }
}

/*************************************************************************/

static void haar_1d_atrou_rec_convol(int N, float *InH,  float *InG, 
                              float *Output, int Step, double NormVal)
{
   int s,i;
   int Ind=0;
   for (i = 0; i < N; i += 2*Step)
   for (s = i; s < MIN(i+Step,N); s ++)
   {
       float Det = (Ind < N/2) ? InG[Ind]: 0;
       float Smooth = InH[Ind];
       Output[s] = NormVal*(Det + Smooth);
       if (s+Step < N) Output[s+Step] = NormVal*(Smooth-Det);
       Ind++;
   }
}


/*************************************************************************/
 
static void haar_1d_atrou_transform(fltarray &Data, fltarray &Mallat, 
                             int Nbr_Plan, double NormVal, int WhichScale)
{
    int Np_2;
    int i,s;
    float *ImagLow, *ImagHigh, *DataResol;
    int Np = Data.nx();
    int Nps = size_resol(1, Np);
    int Step = POW2(WhichScale);
       
    DataResol    = new float[Np];
    ImagHigh     = new float[Nps];
    ImagLow      = new float[Nps];
    if (Np !=  Mallat.n_elem()) Mallat.alloc(Np);
    
    for (i=0; i< Np; i++) DataResol[i] = Data(i);
  
    // Compute the wavelet coefficients 
    for (s = 0; s < Nbr_Plan-1; s++)
    {  
         Nps = size_resol(s, Np);
        Np_2 = (Nps+1)/2;
  // cout << " step " << s << " Nps = " <<   Nps << " Np_2 = " << Np_2 << endl;
        haar_1d_atrou_convol(Nps, DataResol, ImagLow,ImagHigh, Step,  NormVal);
	for (i=0; i < Nps/2; i++)  Mallat(Np_2+i) = ImagHigh[i];
        if (s != Nbr_Plan-1)
             for (i = 0; i < Np_2; i++) DataResol[i] = ImagLow[i];

     }
     Np_2 = size_resol(Nbr_Plan-1, Np);
     for (i=0; i< Np_2; i++) Mallat(i) = ImagLow[i];
    delete [] ImagHigh;
    delete [] ImagLow;
    delete [] DataResol;
}

/****************************************************************************/

static void haar_1d_atrou_reconstruct (fltarray &Mallat, fltarray &Data, int Nbr_Plan, 
                                double NormVal, int WhichScale)
{
    int  Nps;
    register int i,s;
    float *image_h1;
    float *image_g1;
    int Np = Mallat.nx();
    int Np_2=Np;
    int Step = POW2(WhichScale);
 
    // cout << "rec mallat : " << Nbr_Plan << " np = " << Np << endl;
    
    if (Data.n_elem() != Np) Data.alloc(Np);
    /* Allocation */
    image_h1 = new float [Np];
    image_g1 = new float [Np];
  
    // last scale size  
    Np_2 = size_resol(Nbr_Plan-1, Np);
         
    // initial signal construction : image_h1  
    for (i=0; i< Np_2; i++) image_h1[i] = Mallat(i);

    for (s = Nbr_Plan-2; s >= 0; s--)
    {
        Nps = size_resol(s, Np);
        Np_2 = (Nps+1)/2;
 	float *Temp = new float [Nps];
	for (i=0; i < Nps/2; i++)  image_g1[i] = Mallat(Np_2+i);
	// cout << " step " << s+1 << " Nps = " << Nps << endl;
        haar_1d_atrou_rec_convol(Nps, image_h1, image_g1, Temp, Step,NormVal);
			      
 	/* Next iteration */
        for (i = 0; i< Nps; i++) image_h1[i] = Temp[i];
     
  	delete [] Temp;
    }
   //  cout << "out" << endl;
    
    for (i=0; i< Np; i++) Data(i) = image_h1[i];
    delete [] image_h1;
    delete [] image_g1;
}

/****************************************************************************/

void haar_2d_atrou_dir2_transform (Ifloat &Data, MultiResol &T_Haar, Bool EdgeLineTransform)
{
    int Nl = Data.nl();
    int Nc = Data.nc();
    int Nbr_Scale = T_Haar.nbr_scale();
    int i,j,s;
    Ifloat Buff(Nl,Nc, "haar buffer");
    Buff = Data;
    fltarray Line(Nc);
    fltarray TLine(Nc);
    fltarray Col(Nl);
    fltarray TCol(Nl);
    double NormVal = (T_Haar.TypeNorm == NORM_L2) ? 1. / sqrt(2.): 1. / 2.; 
    double NormVal2 = NormVal*NormVal;
    for (s=0; s < Nbr_Scale - 1; s++)
    { 
        // Line transform
        int Step = POW2(s);
        int Ns1DNL = iround((float)log((float) (Nl / 4. * 3.) / log(2.))) - s;
        int Ns1DNC = iround((float)log((float) (Nc / 4. * 3.) / log(2.))) - s;
 
        for (j = 0; j < Nc; j++)
        {
	  if ((EdgeLineTransform == False) || (Ns1DNL <= 1))
	     for (i = 0; i < Nl; i++) 
	        T_Haar(2*s,i,j) = (float) (NormVal*(Buff(i,j) - Buff(i,j+Step,T_Haar.Border)));
          else 
	  {
	     for (i = 0; i < Nl; i++)
	       Col(i) = (float) (NormVal*(Buff(i,j) - Buff(i,j+Step,T_Haar.Border)));
  	     haar_1d_atrou_transform(Col, TCol, Ns1DNL, NormVal, s);
             for (i = 0; i < Nl; i++) T_Haar(2*s,i,j) = TCol(i);    
	  }
        }
              

       for (i = 0; i < Nl; i++)
       {
	   if ((EdgeLineTransform == False) || (Ns1DNC <= 1))
	   {
	     for (j = 0; j < Nc; j++)
	       T_Haar(2*s+1,i,j) = (float) (NormVal*(Buff(i,j) - Buff(i+Step,j, T_Haar.Border)));
	   }
	   else
	   {
              for (j = 0; j < Nc; j++)
	         Line(j) = (float) (NormVal*(Buff(i,j) - Buff(i+Step,j, T_Haar.Border)));
  	      haar_1d_atrou_transform(Line, TLine, Ns1DNC, NormVal, s);
	      for (j = 0; j < Nc; j++) T_Haar(2*s+1,i,j) = TLine(j);
	   }
        }
       
       for (i = 0; i < Nl; i++)
       for (j = 0; j < Nc; j++)
         Buff(i,j) = (float) (NormVal2*(Buff(i,j)+ Buff(i,j+Step, T_Haar.Border) +
                                Buff(i+Step,j, T_Haar.Border) +
				Buff(i+Step,j+Step, T_Haar.Border)));
    }
    T_Haar.band(T_Haar.nbr_band()-1) = Buff;

}

/****************************************************************************/ 

void haar_2d_atrou_dir2_recons (MultiResol & T_Haar, Ifloat &Data, Bool EdgeLineTransform)
{
    int Nbr_Scale = T_Haar.nbr_scale();
    int i,j,s;
    int Nl = Data.nl();
    int Nc = Data.nc();
    Ifloat Buff(Nl, Nc, "Buff");
    double NormVal = (T_Haar.TypeNorm == NORM_L2) ? 1. / sqrt(2.): 1./2.; 
    double NormValRec = (T_Haar.TypeNorm == NORM_L2) ? 1. / sqrt(2.): 1.; 
    fltarray Line(Nc);
    fltarray TLine(Nc);
    fltarray Col(Nl);
    fltarray TCol(Nl);
         
    s = Nbr_Scale-1;
    Buff = T_Haar.band(T_Haar.nbr_band()-1);
    for (s=Nbr_Scale-2; s >= 0; s--)
    {
        int Step = POW2(s);
	if (EdgeLineTransform == True)
	{
           int Ns1DNL = iround((float)log((float) (Nl / 4. * 3.) / log(2.))) - s;
           int Ns1DNC = iround((float)log((float) (Nc / 4. * 3.) / log(2.))) - s;
  	 
           if (Ns1DNL > 1) 
           {
	     for (j = 0; j < Nc; j++)
 	     {
	      for (i = 0; i < Nl; i++) TCol(i) = T_Haar(2*s,i,j);
	      haar_1d_atrou_reconstruct (TCol,  Col,  Ns1DNL, NormValRec, s);  
	      for (i = 0; i < Nl; i++) T_Haar(2*s,i,j) = Col(i);
	     }
	   }

       	   if (Ns1DNC > 1)
	   {
              for (i = 0; i < Nl; i++)
	      {
	         for (j = 0; j < Nc; j++) TLine(j) = T_Haar(2*s+1,i,j);
 	         haar_1d_atrou_reconstruct (TLine,  Line,  Ns1DNC, NormValRec, s);  
	         for (j = 0; j < Nc; j++) T_Haar(2*s+1,i,j) = Line(j);
	      }
	   }
	}      
        for (i = 0; i < Nl; i++)
	for (j = 0; j < Nc; j++)
	{
  	    float u1 = Buff(i,j);
	    float u2 = T_Haar(2*s,i,j);
	    float u3 = T_Haar(2*s,i+Step,j,T_Haar.Border);
	    float u4 = T_Haar(2*s+1,i,j);
 	    Buff(i,j) =  u1 + NormVal*(u2 + u3 + 2*u4);
	    if (T_Haar.TypeNorm == NORM_L2) Buff(i,j) /= 2.;
	}
    }
    Data = Buff;
}

/****************************************************************************/


/****************************************************************************/

void haar_2d_atrou_dir3_transform (Ifloat &Data, MultiResol &T_Haar)
{
    int Nl = Data.nl();
    int Nc = Data.nc();
    int Nbr_Scale = T_Haar.nbr_scale();
    int i,j,s;
    Ifloat Buff(Nl,Nc, "haar buffer");
    Buff = Data;
    fltarray Line(Nc);
    fltarray TLine(Nc);
    fltarray Col(Nl);
    fltarray TCol(Nl);
    double NormVal = (T_Haar.TypeNorm == NORM_L2) ? 1. / sqrt(2.): 1. / 2.; 
    double NormVal2 = NormVal*NormVal;
    for (s=0; s < Nbr_Scale - 1; s++)
    { 
        // Line transform
        int Step = POW2(s);
  
        for (j = 0; j < Nc; j++)
  	for (i = 0; i < Nl; i++) 
	    T_Haar(3*s,i,j) =  (float) (NormVal*(Buff(i,j) - Buff(i,j+Step,T_Haar.Border)));
              
       for (i = 0; i < Nl; i++)
       for (j = 0; j < Nc; j++)
	    T_Haar(3*s+1,i,j) = (float) (NormVal*(Buff(i,j) - Buff(i+Step,j, T_Haar.Border)));

       for (i = 0; i < Nl; i++)
       for (j = 0; j < Nc; j++)
	    T_Haar(3*s+2,i,j) = (float) (NormVal*(Buff(i,j) - Buff(i+Step,j+Step, T_Haar.Border)));

          
       for (i = 0; i < Nl; i++)
       for (j = 0; j < Nc; j++)
         Buff(i,j) = (float) (NormVal2*(Buff(i,j)+ Buff(i,j+Step, T_Haar.Border) +
                                Buff(i+Step,j, T_Haar.Border) +
				Buff(i+Step,j+Step, T_Haar.Border)));
    }
    T_Haar.band(T_Haar.nbr_band()-1) = Buff;
}

/****************************************************************************/ 

void haar_2d_atrou_dir3_recons (MultiResol & T_Haar, Ifloat &Data)
{
    int Nbr_Scale = T_Haar.nbr_scale();
    int i,j,s;
    int Nl = Data.nl();
    int Nc = Data.nc();
    Ifloat Buff(Nl, Nc, "Buff");
    double NormVal = (T_Haar.TypeNorm == NORM_L2) ? 1. / sqrt(2.): 1./2.; 
    double NormValRec1 =  1. / (NormVal * 4.);  
    double NormValRec2 =  1. / (NormVal*NormVal * 4.);         
    s = Nbr_Scale-1;
    Buff = T_Haar.band(T_Haar.nbr_band()-1);
    for (s=Nbr_Scale-2; s >= 0; s--)
    {
        // int Step = POW2(s);
        for (i = 0; i < Nl; i++)
	for (j = 0; j < Nc; j++)
	{
  	    float u1 = Buff(i,j);
	    float u2 = T_Haar(3*s,i,j);
	    float u3 = T_Haar(3*s+1,i,j);
	    float u4 = T_Haar(3*s+2,i,j);
 	    Buff(i,j) =  NormValRec2*u1 + NormValRec1*(u2 + u3 + u4);
	}
    }
    Data = Buff;
}

/****************************************************************************/
