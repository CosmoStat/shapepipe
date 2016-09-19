/******************************************************************************
**                   Copyright (C) 2004 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  30/11/04 
**    
**    File:  Pyr1D.cc
**
*******************************************************************************
**
**    DESCRIPTION   1D Pyramidal WT 
**    -----------      
**
**
******************************************************************************/ 


#include "Pyr1D.h"


/****************************************************************************/

#define C1 .00245
#define C2 -.00915
#define C3 .03414
#define C4 -0.12720
#define C5 0.60048

#define CS ((C1+C2+C3+C4+C5)*2.)
#define C1n (C1/CS)
#define C2n (C2/CS)
#define C3n (C3/CS)
#define C4n (C4/CS)
#define C5n (C5/CS)

static float Tab_Coef_Inter[10] = {C1n,C2n,C3n,C4n,C5n,C5n,C4n,C3n,C2n,C1n};

/***************************************************************************/

void PYR1D_WT::smooth_b3spline(fltarray &Data, fltarray &Res)
{
   Res.resize( Data.nx());
   for (int i = 0; i < Data.nx(); i ++) 
         Res(i) = 0.0625 *  ( Data(i - 2, Border) + Data(i + 2, Border))
                         + 0.25 *  ( Data(i - 1, Border) + Data(i + 1, Border)) + 0.375 * Data(i);
}

/***************************************************************************/

void PYR1D_WT::smooth_b3spline(fltarray &Data)
{
   smooth_b3spline(Data,Buff);
   Data = Buff;
}

/***************************************************************************/

void PYR1D_WT::extand_sample (fltarray &Data, fltarray &Data_Ext)
{
   int i,k;
   int Ns = Data_Ext.nx();
   for (i = 0; i < Ns; i+=2) Data_Ext(i) = Data(i/2);
   for (int j = 1; j < Ns; j+=2)
   {
       // i = j / 2;
       // Data_Ext(j) = 0.;
       // for (k = 0; k < 10; k++) Data_Ext(j) += Data(i - 4 + k, Border)*Tab_Coef_Inter [k];
       Data_Ext(j) = 0.5 * (Data_Ext(j-1, Border) + Data_Ext(j+1, Border));
   }
}

/***************************************************************************/

void PYR1D_WT::reduc_sample (fltarray &Data, fltarray &Data_Sub)
{
   int i;
   for (i = 0; i < Data.nx(); i+=2) Data_Sub(i/2)  = Data (i);
}


/***************************************************************************/

int PYR1D_WT::pyr_tab_pos(int NbrScale, int N)
{
    int s,Ns=N;
    int TotSize=0;
    if (TabSize.nx() != NbrScale) TabSize.reform(NbrScale);
    if (TabPos.nx() != NbrScale)  TabPos.reform(NbrScale);

    for (s = 0; s < NbrScale; s++)
    {
       TotSize += Ns;
       TabSize(s) = Ns;
       Ns = (Ns+1)/2;
    }
    TabPos(NbrScale-1) = 0;
    for (s =  NbrScale-2; s >= 0; s--) TabPos(s) = TabPos(s+1) + TabSize(s+1);

//     for (s = 0; s < NbrScale; s++)
//     {
//         cout << " Band " << s+1 << " size = " << TabSize(s) << " Pos = " << TabPos(s) << endl;
//     }
//     cout << " TotSize = " << TotSize << endl;
    return TotSize;
}

/****************************************************************************/

static void tri_1d (float *fenetre, int Window_Size)
{
    int i,j;
    float Aux;

    for (i = 1; i < Window_Size; i++)
    for (j = 0; j < Window_Size - i; j++)
        if (fenetre[j] > fenetre[j+1])
        {
           Aux = fenetre[j+1];
           fenetre[j+1] = fenetre[j];
           fenetre[j] = Aux;
        }
}
static void filt1d_mediane (const fltarray &S1, fltarray &S2, int N, int Window_Size,
                     type_border Border)
{
    int i,k,ind_fen;
    int Dim,Dim2;
    float *fenetre;

    fenetre = new float [Window_Size];
    Dim = Window_Size;
    Dim2 = Dim / 2;
    for (i = 0; i < N; i++) 
    {
        ind_fen = 0;
        for (k = i -Dim2; k <= i+Dim2; k++) fenetre[ind_fen++] = S1(k, Border);
        S2(i) = get_median(fenetre, Window_Size);
    }
    delete[] fenetre;
}

void PYR1D_WT::transform (fltarray &Data, fltarray &Trans)
{
    int i,s, N = Data.nx();
    int Depi, Ni, N_s=N;
  
    if ((Np != Data.nx()) || ( Nbr_Plan < 2))
    {
       cout << "Error: the class has not been properly allocated ... " << endl;
       exit(-1);
    }
    
    // cout << "WT  transform " << Nbr_Plan << " N = " << N << endl;
    if (Ima_b.nx() != N) Ima_b.resize (N);
    if (Data_Sub.nx() != N) Data_Sub.resize (N);
    if (Buff.nx() != N) Buff.resize (N);
    PyrPixNbr = pyr_tab_pos(Nbr_Plan, N);
    if (Trans.nx() !=  PyrPixNbr) Trans.reform(PyrPixNbr);
    Ima_h = Data;

    for (s = 0; s < Nbr_Plan - 1; s++)
    {
       Depi = pyr_pos(s);
       Ni = pyr_size(s);
       if (Ima_h.nx() != Ni) Ima_h.resize (Ni);
       if (Ima_b.nx() != Ni) Ima_b.resize (Ni);
       
       // cout << "Band " << s+1 << "  Depi = " <<  Depi << " Ni = " << Ni << endl;
       if (Ni != N_s)
       {
          cout << "Errror: Ni and Ns are not equal ... " << endl;
          exit(-1);
       }		
       switch (TypeWT)
       {
          case WT_PYR_BSPLINE:
	    smooth_b3spline(Ima_h, Ima_b);
            break;
          case WT_PYR_LINEAR :
	    for (i = 0; i < N; i ++)
                  Ima_b(i) = 0.25 *  ( Ima_h(i - 1, Border) + Ima_h(i + 1,  Border)) + 0.5 * Ima_h(i);
 	    break;
          case WT_PYR_MEDIAN:
            filt1d_mediane (Ima_h, Ima_b, Ni, WinSize, Border);
            break;
        default:
            fprintf (stderr,"Error: proc. wave_transform.\n");
            fprintf (stderr,"This transform is not computed this procedure\n");
            break;
      }
      Data_Sub.resize( pyr_size(s+1));
      reduc_sample (Ima_b, Data_Sub);
      extand_sample (Data_Sub, Ima_b);
      if (ModifiedPWT == True) smooth_b3spline(Ima_b);
      
      for (i = 0; i < Ni; i++) Trans(i+Depi) = Ima_h(i) - Ima_b(i);
      N_s = (N_s+1)/2;
      if (s != Nbr_Plan - 2) Ima_h = Data_Sub;
   }
   
   s = Nbr_Plan - 1;
   Depi = pyr_pos(s);
   Ni = pyr_size(s);
   if (Ni != N_s)
   {
       cout << "Errror: Ni and Ns are not equal .... " << endl;
       exit(-1);
   }
   for (i = 0; i < Ni; i++) Trans(i+Depi) = Data_Sub(i);
}

/****************************************************************************/

void PYR1D_WT::recons(fltarray &Trans, fltarray &Signal)
{
    if ((pyr_np() != Trans.nx()) || ( Nbr_Plan < 2))
    {
       cout << "Error: the class has not been properly allocated ... " << endl;
       exit(-1);
    }
    int Ni,i;
    int s = Nbr_Plan - 1;
    int Depi = pyr_pos(s);
    Ni = pyr_size(s);
    
    Ima_b.resize(Ni);
    for (i = 0; i < Ni; i++) Ima_b(i) = Trans(i+Depi); 
    for (s = Nbr_Plan - 2; s >= 0; s--)
    {
        int Depi = pyr_pos(s);
        int Ni = pyr_size(s);
        Ima_h.resize(Ni);
        extand_sample (Ima_b, Ima_h);
        if (ModifiedPWT == True) smooth_b3spline(Ima_h);
        for (i = 0; i < Ni; i++) Ima_h(i) += Trans(i+Depi);
	if (s != 0) Ima_b = Ima_h;
    }
    Signal = Ima_h;
}

/****************************************************************************/

void PYR1D_WT::alloc(int TWT, int NbrScale, int N)
{
    TypeWT=TWT;
    ModifiedPWT=False;
    Nbr_Plan = NbrScale;
    Np = N;
    PyrPixNbr=pyr_tab_pos(NbrScale,N);
    Border=I_MIRROR;
    WinSize=5;
}


/****************************************************************************/


