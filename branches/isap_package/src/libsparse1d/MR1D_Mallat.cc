/*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  19/08/98 
**    
**    File:  MR1D_Mallat.cc
**
*****************************************************************************
**
**    DESCRIPTION 1D orthogonal wavelet transform and reconstruction
**    -----------    
**   
****************************************************************************/

#include "MR1D_Obj.h"
#include "SB_Filter1D.h"

// #include "NR.h"

/****************************************************************************/


void mallat_1d_transform (fltarray &Data, fltarray &Mallat, 
                          int Nbr_Plan, SubBand1D &SBF)
{
    int i,s;
    int Np_2,Nps;
    float *ImagLow, *ImagHigh, *DataResol;
    int Np = Data.nx();
      
   // cout << "mallat : " << Nbr_Plan << " np = " << Np << endl;
    
    if (Np !=  Mallat.n_elem()) Mallat.alloc(Np);
    
    /* Allocation  */
    Nps = size_resol(1, Np);
    DataResol    = new float[Np];
    ImagHigh     = new float[Nps];
    ImagLow      = new float[Nps];
 
    for (i=0; i< Np; i++) DataResol[i] = Data(i);
  
    /* Compute the wavelet coefficients */
    for (s = 0; s < Nbr_Plan-1; s++)
    {  
  //  cout << " step " << s << endl;
        Nps = size_resol(s, Np);
        Np_2 = (Nps+1)/2;
	
        SBF.transform(Nps,DataResol,ImagLow,ImagHigh);
        for (i=0; i < Nps/2; i++)  Mallat(Np_2+i) = ImagHigh[i];
 
  	if (s != Nbr_Plan-1)
             for (i = 0; i < Np_2; i++)  
	               DataResol[i] = ImagLow[i];
     }
     Np_2 = size_resol(Nbr_Plan-1, Np);
     for (i=0; i< Np_2; i++) Mallat(i) = ImagLow[i];
 // cout << "exit" << endl;
    delete [] ImagHigh;
    delete [] ImagLow;
    delete [] DataResol;
}

/****************************************************************************/

void mallat_1d_reconstruct (fltarray &Mallat, fltarray &Data, int Nbr_Plan, 
                            SubBand1D &SBF)
{
    int  Nps;
    register int i,s;
    float *image_h1;
    float *image;
    float *image_g1;
    int Np = Mallat.nx();
    int Np_2=Np;
 
    // cout << "rec mallat : " << Nbr_Plan << " np = " << Np << endl;
    
    if (Data.n_elem() != Np) Data.alloc(Np);

    /* Allocation */
    image    = new float [Np];
    image_h1 = new float [Np];
    image_g1 = new float [Np];
  
    /* last scale size */
    Np_2 = size_resol( Nbr_Plan-1, Np);
         
     /* initial image construction : image_h1_h1. */
    for (i=0; i< Np_2; i++) image_h1[i] = Mallat(i);

    for (s = Nbr_Plan-2; s >= 0; s--)
    {
    // cout << " step " << s << endl;
        Nps = size_resol(s, Np);
        Np_2 = (Nps+1)/2;
        for (i=0; i< Nps/2; i++) image_g1[i] = Mallat(i+Np_2);
        SBF.recons (Nps, image_h1, image_g1, image);
        /* Next iteration */
        for (i = 0; i < Nps; i++) image_h1[i] = image[i];
    }
   //  cout << "out" << endl;
    
    for (i=0; i< Np; i++) Data(i) = image[i];
    delete [] image;
    delete [] image_h1;
    delete [] image_g1;
}

/*********************************************************************/
/*
void lift1_1d_transform (fltarray &Data, fltarray &Lift, int Nbr_Plan)
{
    int i,s,Index;
    int Np_2,Nps;
    float *ImagLow, *ImagHigh, *DataResol;
    int Np = Data.nx();
 
    if (Np != Lift.n_elem())  Lift.alloc(Np);
    // cout << "lift1d : " << Nbr_Plan << " np = " << Np << endl;
          
    // Allocation   
    Nps = size_resol(1, Np);
    DataResol    = new float[Np];
    ImagHigh     = new float[Nps];
    ImagLow      = new float[Nps];
 
    for (i=0; i< Np; i++) DataResol[i] = Data(i);
  
    /. Compute the wavelet coefficients  
    for (s = 0; s < Nbr_Plan-1; s++)
    {  
        Nps = size_resol(s, Np);
        Np_2 = (Nps+1)/2;

        // linear prediction
        for (i = 1; i < Nps; i += 2)
        {
            Index = test_index_mirror(i+1, Nps);
	    float Predict = 0.5*(DataResol[i-1]+DataResol[Index]);
	    ImagHigh[i/2] = DataResol[i] - Predict;
        }  
        for (i = 0; i < Nps; i += 2)
        {
            Index = test_index_mirror(i/2+1, Nps/2);
	    float Update = 0.25*(ImagHigh[i/2]+ImagHigh[Index]);
            ImagLow[i/2] = DataResol[i]  + Update;
        }
        for (i=0; i < Nps/2; i++)  Lift(i+Np_2) = ImagHigh[i];
 
  	if (s != Nbr_Plan-1)
             for (i = 0; i < Np_2; i++)  
	               DataResol[i] = ImagLow[i];
     }
     Np_2 = size_resol(Nbr_Plan-1, Np);
     for (i=0; i< Np_2; i++)  Lift(i) = ImagLow[i];

    delete [] ImagHigh;
    delete [] ImagLow;
    delete [] DataResol;
}
 
void lift1_1d_reconstruct (fltarray &Lift, fltarray &Data, int Nbr_Plan)
{
    int  Nps,Index;
    register int i,s;
    float *image_h1;
    float *image;
    float *image_g1;
    int Np = Lift.nx();
    int Np_2=Np;

    if (Data.n_elem() != Np) Data.alloc(Np);
        
    // Allocation  
    image    = new float [Np];
    image_h1 = new float [Np];
    image_g1 = new float [Np];
  
    // last scale size  
    Np_2 = size_resol( Nbr_Plan-1, Np);
         
     // initial image construction : image_h1_h1. 
    for (i=0; i< Np_2; i++)
                 image_h1[i] = Lift(i);

    for (s = Nbr_Plan-2; s >= 0; s--)
    {
        Nps = size_resol(s, Np);
        Np_2 = (Nps+1)/2;
        for (i=0; i< Nps/2; i++) image_g1[i] = Lift(i+Np_2);
        for (i = 0; i < Nps; i += 2)
        {
         Index = test_index_mirror(i/2+1, Np_2);
       	 float Update = 0.25*( image_g1[i/2]+image_g1[Index]);
         image[i] = image_h1[i/2]-Update;
        }
        for (i = 1; i < Nps; i += 2)
        {
         Index = test_index_mirror(i+1, Nps);
	 float Predict = 0.5*(image[i-1]+image[Index]);
         image[i] = image_g1[i/2]+Predict;
        }        
        // Next iteration  
        for (i = 0; i < Nps; i++) image_h1[i] = image[i];
    }
    
    for (i=0; i< Np; i++) Data(i) = image[i];

    delete [] image;
    delete [] image_h1;
    delete [] image_g1;
}
*/

/*********************************************************************/

static void wp1d_step_mallat(float *Signal, int Np, 
                      int Nstep, SubBand1D &SBF)
{
    float *ImagLow, *ImagHigh, *DataResol;
    int i,Np_2L = (Np+1)/2;
    int Np_2H = Np/2;
    
    DataResol    = new float[Np];
    ImagHigh     = new float[Np_2L];
    ImagLow      = new float[Np_2L];
    for (i=0; i< Np; i++) DataResol[i] = Signal[i];
    // cout << "tran: " <<  Np_2H << "  " << Np_2L <<  endl;

    SBF.transform(Np,DataResol,ImagLow,ImagHigh);
    for (i=0; i< Np_2L; i++) Signal[i] = ImagLow[i];
    for (i=0; i< Np_2H; i++) Signal[i+Np_2L] = ImagHigh[i];
    if (Nstep > 1)
    {
       wp1d_step_mallat(Signal+Np_2L, Np_2H, Nstep-1, SBF);
       wp1d_step_mallat(Signal, Np_2L,Nstep-1, SBF);
    }
    delete [] ImagHigh;
    delete [] ImagLow;
    delete [] DataResol;
}

/****************************************************************************/

void wp1d_mallat_transform (fltarray &Data, fltarray &Mallat,
                            int Nstep, SubBand1D &SBF)
{
    int i,Np = Data.nx();     
    if (Np != Mallat.n_elem())  Mallat.alloc(Np);

   /* Allocation  */
    float *WP = Mallat.buffer();
    for (i=0; i< Np; i++) Mallat(i) = Data(i);
    wp1d_step_mallat(WP, Np, Nstep, SBF);
}
			    
/****************************************************************************/

static void wp1d_rec_step_mallat(float *Signal, int Np, 
                      int Nstep, SubBand1D &SBF)
{
    float *ImagLow, *ImagHigh, *DataResol;
    int i,Np_2L = (Np+1)/2;
    int Np_2H = Np/2;
    
    if (Nstep > 1)
    {
       // cout << "Step: " << Nstep << endl;
       // cout << "   High resol" << endl;
       wp1d_rec_step_mallat(Signal+Np_2L, Np_2H, Nstep-1, SBF);
       // cout << "   Low resol" << endl;
       wp1d_rec_step_mallat(Signal, Np_2L,Nstep-1, SBF);
    }
    
    DataResol    = new float[Np];
    ImagHigh     = new float[Np_2L];
    ImagLow      = new float[Np_2L];
    // cout << "Rec: " <<  Np_2H << "  " << Np_2L <<  endl;
    for (i=0; i< Np_2H; i++) ImagHigh[i] = Signal[i+Np_2L];
    for (i=0; i< Np_2L; i++) ImagLow[i] = Signal[i];
    SBF.recons (Np, ImagLow, ImagHigh, DataResol);
    for (i=0; i< Np; i++) Signal[i] = DataResol[i];
    delete [] ImagHigh;
    delete [] ImagLow;
    delete [] DataResol;
}

/****************************************************************************/

void wp1d_mallat_rec (fltarray &Mallat, fltarray &Data,
                      int Nstep, SubBand1D &SBF)
{
    int i,Np = Mallat.nx();
 
    if (Np != Data.n_elem())  Data.alloc(Np);
    
   /* Allocation  */
    float *WP = Data.buffer();
    for (i=0; i< Np; i++) Data(i) = Mallat(i);
    wp1d_rec_step_mallat(WP, Np, Nstep, SBF);
} 

/****************************************************************************/
