 /*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  20/08/98 
**    
**    File:  SB_Filter.cc
**
**    Modification history:
**
*******************************************************************************
**
**    DESCRIPTION  Sub-Band decomposition
**    -----------  
**                 
******************************************************************************/
 

#include "IM_Obj.h"
#include "SB_Filter.h"
 

/*************************************************************/
// 2D decimated subband decomposition
/*************************************************************/
 
void SubBand2D::transform2d (Ifloat &Data, Bool UseSameBuff, Ifloat *Horiz, 
                          Ifloat *Vert, Ifloat *Diag, Ifloat *Smooth)
{
    int Nl = Data.nl();
    int Nc = Data.nc();
    float *ImagHigh = Data.buffer();
    int Nc_2, Nl_2;
    int j, i, Index;
    float *image_h0, *image_g0, *Col_h0, *Col_g0, *Col_h0_h0;
    float *Col_h0_g0, *Col_g0_h0,*Col_g0_g0;
  
    Nl_2 = (Nl+1)/2;
    Nc_2 = (Nc+1)/2;

    // each line is convolved by h and g
    // one pixel out of two (in column) is kept 
    image_h0 = new float [Nl*Nc_2];
    image_g0 = new float [Nl*Nc_2];

    for (i = 0; i < Nl; i++) 
    {
      float *PtrHigh = ImagHigh + Nc*i;
      float *PtrLow = image_h0 + Nc_2*i;
      float *PtrDet = image_g0 + Nc_2*i;
      Ptr_SB1D_LINE->transform(Nc, PtrHigh, PtrLow, PtrDet);
   }
   // Ifloat FF;
   // FF.getbuff(image_g0, Nc_2, Nl);
   // fits_write_fltarr ("xx_tt.fits", FF);


    Col_h0    = new float [Nl];
    Col_g0    = new float [Nl];
    Col_h0_h0 = new float [Nl_2];
    Col_h0_g0 = new float [Nl_2];
    Col_g0_h0 = new float [Nl_2];
    Col_g0_g0 = new float [Nl_2];
    if (UseSameBuff ==  False) (*Smooth).resize(Nl_2, Nc_2);

    for (j = 0; j< Nc_2; j++)
    {
        for (i = 0; i< Nl; i++)
	{
            Index = Nc_2 * i + j;
            Col_h0[i] = image_h0[Index];
            if (j + Nc_2 < Nc) Col_g0[i] = image_g0[Index];
	}
	Ptr_SB1D_COL->transform(Nl, Col_h0, Col_h0_h0, Col_h0_g0);
 	if (j + Nc_2 < Nc) 
	    Ptr_SB1D_COL->transform(Nl, Col_g0, Col_g0_h0, Col_g0_g0);
        
        for (i = 0; i< Nl_2; i++)
	{
	    if (UseSameBuff == True)
	    {
 	       Data(i,j) = Col_h0_h0[i];
 	       if (i+Nl_2 < Nl) Data(i+Nl_2,j) = Col_h0_g0[i];
 	       if (j+Nc_2 < Nc) Data(i,j+Nc_2) = Col_g0_h0[i];
 	       if ((i+Nl_2 < Nl)  && (j+Nc_2 < Nc)) 
	               Data(i+Nl_2,j+Nc_2) = Col_g0_g0[i];
            }
	    else
	    {
	       (*Smooth)(i,j) = Col_h0_h0[i];
	       if (i+Nl_2 < Nl) (*Vert)(i,j) = Col_h0_g0[i];
	       if (j+Nc_2 < Nc) (*Horiz)(i,j) = Col_g0_h0[i];
	       if ((i+Nl_2 < Nl) && (j+Nc_2 < Nc)) (*Diag)(i,j) = Col_g0_g0[i];
	    }
        }
    }
    
    delete [] image_h0;
    delete [] image_g0;
    delete [] Col_h0;
    delete [] Col_g0;
    delete [] Col_h0_h0;
    delete [] Col_h0_g0;
    delete [] Col_g0_h0;
    delete [] Col_g0_g0;
}			

/****************************************************************************/

void SubBand2D::transform2d(Ifloat &Data, Ifloat &Horiz, Ifloat &Vert, 
	                  Ifloat &Diag, Ifloat &Smooth, int Step)
{
    int Nl = Data.nl();
    int Nc = Data.nc();
    float *ImagHigh = Data.buffer();
    int j, i, Index;
    float *image_h0, *image_g0, *Col_h0, *Col_g0, *Col_h0_h0;
    float *Col_h0_g0, *Col_g0_h0,*Col_g0_g0;
 
    // each line is convolved by h and g
    image_h0 = new float [Nl*Nc];
    image_g0 = new float [Nl*Nc];
    for (i = 0; i < Nl; i++) 
    {
      float *PtrHigh = ImagHigh + Nc*i;
      float *PtrLow = image_h0 + Nc*i;
      float *PtrDet = image_g0 + Nc*i;
      Ptr_SB1D_LINE->transform(Nc, PtrHigh, PtrLow, PtrDet, Step);
    }

    Col_h0    = new float [Nl];
    Col_g0    = new float [Nl];
    Col_h0_h0 = new float [Nl];
    Col_h0_g0 = new float [Nl];
    Col_g0_h0 = new float [Nl];
    Col_g0_g0 = new float [Nl];

    for (j = 0; j< Nc; j++)
    {
        for (i = 0; i< Nl; i++)
	{
            Index = Nc * i + j;
            Col_h0[i] = image_h0[Index];
            Col_g0[i] = image_g0[Index];
	}
	Ptr_SB1D_COL->transform(Nl, Col_h0, Col_h0_h0, Col_h0_g0, Step);
 	Ptr_SB1D_COL->transform(Nl, Col_g0, Col_g0_h0, Col_g0_g0, Step);
        
        for (i = 0; i< Nl; i++)
	{
 	    Smooth(i,j) = Col_h0_h0[i];
	    Vert(i,j) = Col_h0_g0[i];
	    Horiz(i,j) = Col_g0_h0[i];
	    Diag(i,j) = Col_g0_g0[i];
        }
    }
    
    delete [] image_h0;
    delete [] image_g0;
    delete [] Col_h0;
    delete [] Col_g0;
    delete [] Col_h0_h0;
    delete [] Col_h0_g0;
    delete [] Col_g0_h0;
    delete [] Col_g0_g0;
}			
/****************************************************************************/

void SubBand2D::recons2d (Ifloat &Data, Bool UseSameBuff,
	             Ifloat *Horiz, Ifloat *Vert, Ifloat *Diag, Ifloat *Smooth)
{
    int  Nc_2, Nl_2;
    register int i, j, Index;
    float *image_h1, *image_g1, *Col_h1_h1;
    float *Col_h1_g1, *Col_g1_h1, *Col_g1_g1;
    float *Col_h1, *Col_g1;
    int Nl = Data.nl();
    int Nc = Data.nc();
       
    Nl_2 = (Nl+1)/2;
    Nc_2 = (Nc+1)/2;
          
    /* Allocation  */
    image_h1 = new float [Nl*Nc_2];
    image_g1 = new float [Nl*Nc_2];

    Col_h1 = new float [Nl];
    Col_g1 = new float [Nl];
    Col_h1_h1 = new float [Nl_2];
    Col_h1_g1 = new float [Nl_2];
    Col_g1_h1 = new float [Nl_2];
    Col_g1_g1 = new float [Nl_2];


    /* Transform on the columns */
    for (j = 0; j< Nc_2; j++)
    {
        for (i = 0; i< Nl_2; i++)
	{
	    if (UseSameBuff == True)
	    {	
               Col_h1_h1[i] = Data(i,j);
               if (i+Nl_2 < Nl) Col_h1_g1[i] = Data(i+Nl_2,j);
               if (j+Nc_2 < Nc) Col_g1_h1[i] = Data(i,j+Nc_2);
               if ((i+Nl_2 < Nl)  && (j+Nc_2 < Nc))
	                  Col_g1_g1[i] =  Data(i+Nl_2,j+Nc_2);
            }
	    else
	    {
	       Col_h1_h1[i] = (*Smooth)(i,j);
	       if (i+Nl_2 < Nl) Col_h1_g1[i] = (*Vert)(i,j);
	       if (j+Nc_2 < Nc) Col_g1_h1[i] = (*Horiz)(i,j);
	       if ((i+Nl_2 < Nl)  && (j+Nc_2 < Nc)) 
	                 Col_g1_g1[i] = (*Diag)(i,j);
	    }
	}
        Ptr_SB1D_COL->recons (Nl, Col_h1_h1, Col_h1_g1, Col_h1);
        Ptr_SB1D_COL->recons (Nl, Col_g1_h1, Col_g1_g1, Col_g1);
        for (i = 0; i< Nl; i++)
	{
            Index = Nc_2 * i + j;
            image_h1[Index] = Col_h1[i];
            if (j+Nc_2 < Nc) image_g1[Index] = Col_g1[i];
	}
    }

    delete [] Col_h1;
    delete [] Col_g1;
    delete [] Col_h1_h1;
    delete [] Col_h1_g1;
    delete [] Col_g1_h1;
    delete [] Col_g1_g1;

    /* Transforms on the lines */
    float *image_Output = Data.buffer();
    for (i = 0; i< Nl; i++)  Ptr_SB1D_LINE->recons (Nc, image_h1+Nc_2*i,
                                image_g1+Nc_2*i, image_Output+Nc*i);

    delete [] image_h1;
    delete [] image_g1;
}
/************************************************************************/

void SubBand2D::recons2d(Ifloat &Horiz, Ifloat &Vert, 
	               Ifloat &Diag, Ifloat &Smooth, Ifloat &Output, int Step)
{
    register int i, j, Index;
    float *image_h1, *image_g1, *Col_h1_h1;
    float *Col_h1_g1, *Col_g1_h1, *Col_g1_g1;
    float *Col_h1, *Col_g1;
    int Nl = Output.nl();
    int Nc = Output.nc();
       
    /* Allocation  */
    image_h1 = new float [Nl*Nc];
    image_g1 = new float [Nl*Nc];

    Col_h1 = new float [Nl];
    Col_g1 = new float [Nl];
    Col_h1_h1 = new float [Nl];
    Col_h1_g1 = new float [Nl];
    Col_g1_h1 = new float [Nl];
    Col_g1_g1 = new float [Nl];


    /* Transform on the columns */
    for (j = 0; j< Nc; j++)
    {
        for (i = 0; i< Nl; i++)
	{
            Col_h1_h1[i] = Smooth(i,j);
            Col_h1_g1[i] = Vert(i,j);
            Col_g1_h1[i] = Horiz(i,j);
            Col_g1_g1[i] =  Diag(i,j);
	}
        Ptr_SB1D_COL->recons (Nl, Col_h1_h1, Col_h1_g1, Col_h1, Step);
        Ptr_SB1D_COL->recons (Nl, Col_g1_h1, Col_g1_g1, Col_g1, Step);
        for (i = 0; i< Nl; i++)
	{
            Index = Nc * i + j;
            image_h1[Index] = Col_h1[i];
            image_g1[Index] = Col_g1[i];
	}
    }

    delete [] Col_h1;
    delete [] Col_g1;
    delete [] Col_h1_h1;
    delete [] Col_h1_g1;
    delete [] Col_g1_h1;
    delete [] Col_g1_g1;

    /* Transforms on the lines */
    float *image_Output = Output.buffer();
    for (i = 0; i< Nl; i++)  Ptr_SB1D_LINE->recons (Nc, image_h1+Nc*i,
                                image_g1+Nc*i, image_Output+Nc*i, Step);

    delete [] image_h1;
    delete [] image_g1;
}
/****************************************************************************/
//                       Ortho_2D_WT
/****************************************************************************/

void Ortho_2D_WT::transform (Ifloat &Imag, Ifloat  &Imag_Out, int Nbr_Plan)
{
   // copy the image in Ima_Aux
   Imag_Out = Imag;
   transform(Imag_Out, Nbr_Plan);
}

/****************************************************************************/

void Ortho_2D_WT::transform (Ifloat &Imag, int Nbr_Plan)
{
   int i,j,k;
   int Nl = Imag.nl();
   int Nc = Imag.nc();
   int Nls = Nl;
   int Ncs = Nc;   
   SubBand1D *SB1D_LINE = get_subband_method_line();
   SubBand1D *SB1D_COL = get_subband_method_col();
   for (k = 0; k < Nbr_Plan-1; k++)
   {
      int Nl_2 = (Nls+1)/2;
      int Nc_2 = (Ncs+1)/2;
      fltarray LowLine(Nc_2);
      fltarray DetLine(Nc_2);
      for (i = 0; i < Nls; i++) 
      {
         float *PtrHigh = Imag.buffer() + Nc*i;
         SB1D_LINE->transform(Ncs, PtrHigh, LowLine.buffer(), DetLine.buffer());
         for (j = 0; j < Nc_2; j++)  PtrHigh[j] = LowLine(j);
         for (j = 0; j < Ncs/2; j++) PtrHigh[j+Nc_2] = DetLine(j);
      }
      fltarray HighCol(Nls);
      fltarray LowCol(Nl_2);
      fltarray DetCol(Nl_2);
      for (j = 0; j < Ncs; j++) 
      {
          for (i = 0; i < Nls; i++) HighCol(i) = Imag(i,j);
          SB1D_COL->transform(Nls, HighCol.buffer(), LowCol.buffer(), DetCol.buffer());
          for (i = 0; i < Nl_2; i++)  Imag(i,j) = LowCol(i);
          for (i = 0; i < Nls/2; i++) Imag(i+Nl_2,j) = DetCol(i);       
      }
      // next step on a reduced image
      Nls = (Nls+1)/2;
      Ncs = (Ncs+1)/2;
   }
}

/************************************************************************/

void Ortho_2D_WT::transform (Iint &Imag, Iint  &Imag_Out, int Nbr_Plan)
{
   // copy the image in Ima_Aux
   Imag_Out = Imag;
   transform(Imag_Out, Nbr_Plan);
}

/************************************************************************/

// void Ortho_2D_WT::transform (Iint &Imag, int Nbr_Plan)
// {
//    int i,j,k;
//    int Nl = Imag.nl();
//    int Nc = Imag.nc();
// 
//    Ifloat Ima_Aux(Nl,Nc,"aux");   
//    for (k = 0; k < Nbr_Plan-1; k++)
//    {
//        Ima_Aux.resize(Nl, Nc);
//        for (i = 0; i < Nl; i++)
//        for (j = 0; j < Nc; j++) Ima_Aux(i,j) = (float) Imag(i,j);
//        
//        // transform Ima_Aux
//        transform2d(Ima_Aux);
// 
//        // copy the transform part in Imag_Out
//        for (i = 0; i < Nl; i++)
//        for (j = 0; j < Nc; j++) 
//          Imag(i,j) = (Ima_Aux(i,j) >= 0) ? 
// 	                (int) (Ima_Aux(i,j) + 0.5) : (int) (Ima_Aux(i,j) - 0.5) ;
// 
//        // next step on a reduced image
//        Nl = (Nl+1)/2;
//        Nc = (Nc+1)/2;
//    }
// }

/************************************************************************/

static inline int to_int(float Val)
{
   return (Val >= 0) ?  (int) (Val + 0.5) : (int) (Val - 0.5);
}
 
void Ortho_2D_WT::transform (Iint &Imag, int Nbr_Plan)
{
   int i,j,k;
   int Nl = Imag.nl();
   int Nc = Imag.nc();
   int Nls = Nl;
   int Ncs = Nc;
   SubBand1D *SB1D = get_subband_method();
   for (k = 0; k < Nbr_Plan-1; k++)
   {
      int Nl_2 = (Nls+1)/2;
      int Nc_2 = (Ncs+1)/2;
      fltarray HighLine(Ncs);
      fltarray LowLine(Nc_2);
      fltarray DetLine(Nc_2);
      for (i = 0; i < Nls; i++) 
      {
         for (j = 0; j < Ncs; j++) HighLine(j) = Imag(i,j);
         SB1D->transform(Ncs,  HighLine.buffer(), LowLine.buffer(), DetLine.buffer());
         for (j = 0; j < Nc_2; j++)  Imag(i,j)  = to_int( LowLine(j));
         for (j = 0; j < Ncs/2; j++) Imag(i,j+Nc_2) = to_int(  DetLine(j));
      }
      fltarray HighCol(Nls);
      fltarray LowCol(Nl_2);
      fltarray DetCol(Nl_2);
      for (j = 0; j < Ncs; j++) 
      {
          for (i = 0; i < Nls; i++) HighCol(i) = Imag(i,j);
          SB1D->transform(Nls, HighCol.buffer(), LowCol.buffer(), DetCol.buffer());
          for (i = 0; i < Nl_2; i++)  Imag(i,j) =  to_int(  LowCol(i));
          for (i = 0; i < Nls/2; i++) Imag(i+Nl_2,j) = to_int( DetCol(i));       
      }
      // next step on a reduced image
      Nls = (Nls+1)/2;
      Ncs = (Ncs+1)/2;
   } 
}

/****************************************************************************/

void Ortho_2D_WT::recons (Ifloat &Transf_in, Ifloat &Imag_out, int Nbr_Plan) 
{
   Imag_out = Transf_in;
   recons(Imag_out, Nbr_Plan);
}
 
/****************************************************************************/

void Ortho_2D_WT::recons (Ifloat &Imag, int Nbr_Plan) 
{
   int i,j,s;
   int Nl = Imag.nl();
   int Nc = Imag.nc();
   int  Nls=Nl;
   int  Ncs=Nc;
   SubBand1D *SB1D_LINE = get_subband_method_line();
   SubBand1D *SB1D_COL = get_subband_method_col();
  
   // copy the image in Ima_Aux
   for (s = Nbr_Plan-2; s >= 0; s--)
   {
       Nls = size_resol(s, Nl);
       Ncs = size_resol(s, Nc);
       int Nl_2 = (Nls+1)/2;
       int Nc_2 = (Ncs+1)/2;
       fltarray HighCol(Nls);
       fltarray LowCol(Nl_2);
       fltarray DetCol(Nl_2);

       for (j = 0; j < Ncs; j++) 
       {
          for (i = 0; i < Nl_2; i++)  LowCol(i) = Imag(i,j);
          for (i = 0; i < Nls/2; i++) DetCol(i) = Imag(i+Nl_2,j);
          SB1D_COL->recons (Nls, LowCol.buffer(), DetCol.buffer(), HighCol.buffer());
          for (i = 0; i < Nls; i++) Imag(i,j) = HighCol(i);
       }
       // fltarray HighLine(Ncs);
       fltarray LowLine(Nc_2);
       fltarray DetLine(Nc_2);

       for (i = 0; i < Nls; i++) 
       {
          float  *PtrHigh = Imag.buffer() + Nc*i;
          for (j = 0; j < Nc_2; j++)  LowLine(j) = Imag(i,j);
          for (j = 0; j <  Ncs/2; j++) DetLine(j) = Imag(i,j+Nc_2);
          SB1D_LINE->recons(Ncs, LowLine.buffer(), DetLine.buffer(), PtrHigh);
          // for (j = 0; j < Ncs; j++) Imag(i,j) = HighLine(j);
       }
   } 
}


/****************************************************************************/

void Ortho_2D_WT::recons (Iint & Transf_in, Iint &Imag_out, int Nbr_Plan)
{
   Imag_out = Transf_in;
   recons(Imag_out, Nbr_Plan);
}

/****************************************************************************/

void Ortho_2D_WT::recons (Iint &Imag, int Nbr_Plan)
{
   int i,j,s;
   int Nl = Imag.nl();
   int Nc = Imag.nc();
   SubBand1D *SB1D = get_subband_method();
   int  Nls=Nl;
   int  Ncs=Nc;
 
   // copy the image in Ima_Aux
   for (s = Nbr_Plan-2; s >= 0; s--)
   {
       Nls = size_resol(s, Nl);
       Ncs = size_resol(s, Nc);
       int Nl_2 = (Nls+1)/2;
       int Nc_2 = (Ncs+1)/2;
       fltarray HighCol(Nls);
       fltarray LowCol(Nl_2);
       fltarray DetCol(Nl_2);
       for (j = 0; j < Ncs; j++) 
       {
          for (i = 0; i < Nl_2; i++)  LowCol(i) = Imag(i,j);
          for (i = 0; i < Nls/2; i++) DetCol(i) = Imag(i+Nl_2,j);
          SB1D->recons (Nls, LowCol.buffer(), DetCol.buffer(), HighCol.buffer());
          for (i = 0; i < Nls; i++) Imag(i,j) = to_int(HighCol(i));
       }
       fltarray HighLine(Ncs);
       fltarray LowLine(Nc_2);
       fltarray DetLine(Nc_2);
       for (i = 0; i < Nls; i++) 
       {
          for (j = 0; j < Nc_2; j++)   LowLine(j) = Imag(i,j);
          for (j = 0; j < Ncs/2; j++) DetLine(j) = Imag(i,j+Nc_2);
          SB1D->recons(Ncs, LowLine.buffer(), DetLine.buffer(), HighLine.buffer());
          for (j = 0; j < Ncs; j++) Imag(i,j) = to_int(HighLine(j));
       }
   }  
}

/****************************************************************************/

int PAVE_2D_WT::alloc (Ifloat * & TabBand, int Nl, int Nc, int Nbr_Plan)
{
    int s,NbrBand_per_Resol = 3;
    int NbrBand = NbrBand_per_Resol*(Nbr_Plan-1)+1;
    TabBand = new Ifloat [NbrBand];
    for (s = 0; s < NbrBand; s++) TabBand[s].alloc (Nl, Nc); 
    return NbrBand;
}

/****************************************************************************/

void PAVE_2D_WT::free(Ifloat *TabBand, int Nbr_Plan)
{
    if (Nbr_Plan != 0) delete [] TabBand;
}

/****************************************************************************/

void PAVE_2D_WT::one_scale_transform (Ifloat &Imag, Ifloat *TabTrans, int Step, int Pos)
{
    transform2d (Imag, TabTrans[Pos],  TabTrans[Pos+1],  
	            TabTrans[Pos+2], TabTrans[Pos+3], Step);
}

/****************************************************************************/

void PAVE_2D_WT::transform (Ifloat &Imag, Ifloat *TabTrans, int Nbr_Plan)
{
   for (int s = 0; s < Nbr_Plan-1; s++)
   {
      int Step = POW2(s);
      if (s == 0)
           transform2d (Imag, TabTrans[3*s],  TabTrans[3*s+1],  
	            TabTrans[3*s+2], TabTrans[3*s+3], Step);
      else transform2d (TabTrans[3*s], TabTrans[3*s],  TabTrans[3*s+1],  
	            TabTrans[3*s+2], TabTrans[3*s+3], Step);
    }   
}

/****************************************************************************/

void PAVE_2D_WT::one_scale_recons(Ifloat *TabTrans, Ifloat &Imag, int Step, int Pos)
{
     recons2d(TabTrans[Pos],  TabTrans[Pos+1],  TabTrans[Pos+2], 
	            TabTrans[Pos+3], Imag, Step);
}

/****************************************************************************/

void PAVE_2D_WT::recons (Ifloat *TabTrans, Ifloat &Imag, int Nbr_Plan)
{
   for (int s = Nbr_Plan-2; s >= 0; s--)
   {
      int Step = POW2(s);
      if (s == Nbr_Plan-2)
           recons2d(TabTrans[3*s],  TabTrans[3*s+1],  TabTrans[3*s+2], 
	            TabTrans[3*s+3], Imag, Step);
      else recons2d (TabTrans[3*s],  TabTrans[3*s+1], TabTrans[3*s+2], 
                     Imag, Imag, Step);
    }   
}

/****************************************************************************/


void HALF_DECIMATED_2D_WT::free(Ifloat *TabBand, int Nbr_Plan)
{
    if (Nbr_Plan != 0) delete [] TabBand;
}

void HALF_DECIMATED_2D_WT::set_tabdec(int NumUndec, Bool * & TabDec, int Nbr_Plan)
{
    int i;
    int NU = (NumUndec < 0) ? Nbr_Plan: NumUndec;
    TabDec = new Bool [Nbr_Plan];
    for (i=0; i < MIN(NU,Nbr_Plan); i++) TabDec[i] = False;
    for (i=MIN(NU,Nbr_Plan); i < Nbr_Plan; i++) TabDec[i] = True;
}


int HALF_DECIMATED_2D_WT::alloc (Ifloat * & TabBand, 
	                          int Nl, int Nc, int Nbr_Plan,int NumUndec)
{
    Bool *TabDec;
    set_tabdec(NumUndec, TabDec, Nbr_Plan);
    int NbrBand = alloc (TabBand,Nl,Nc,Nbr_Plan,TabDec);
    delete [] TabDec;
    return NbrBand;
}

int HALF_DECIMATED_2D_WT::alloc (Ifloat * & TabBand, 
	                          int Nl, int Nc, int Nbr_Plan, Bool *TabDec)
{
    char ch[80];
    int Nl_s = Nl;
    int Nc_s = Nc;
    int NbrBand_per_Resol = 3;
    int Nlb, Ncb, s ;
    int NbrBand = NbrBand_per_Resol*(Nbr_Plan-1)+1;
    TabBand = new Ifloat [NbrBand];
    for (s = 0; s < NbrBand-1; s+=3)
    {
        Bool Dec = TabDec[s/3];
	Nlb = Nl_s;
	Ncb = Nc_s;
	
        if (Dec == True)
        {
           Nlb = (Nl_s+1)/2;
           Ncb =  Nc_s/2;
        }
	sprintf (ch, "band_%d", s+1);
        TabBand[s].alloc (Nlb, Ncb, ch);   
            
	if (Dec == True)
        {
            Nlb = Nl_s/2;
            Ncb = (Nc_s+1)/2;
	}
        sprintf (ch, "band_%d", s+2);
        TabBand[s+1].alloc (Nlb, Ncb, ch);
           
	if (Dec == True)
        {        
            Nlb = Nl_s/2;
            Ncb = Nc_s/2;
	}
        sprintf (ch, "band_%d", s+3);
        TabBand[s+2].alloc(Nlb, Ncb,ch);
	
	if (Dec == True)
        {
            Nl_s = (Nl_s+1)/2;
            Nc_s = (Nc_s+1)/2;
        }
     }
     s = NbrBand-1;
     Nlb = Nl_s;
     Ncb = Nc_s;
     sprintf (ch, "band_%d", s+1);
     TabBand[s].alloc(Nlb, Ncb, ch);
     return NbrBand;
}

void HALF_DECIMATED_2D_WT::transform(Ifloat &Imag, Ifloat *TabTrans, 
	      int Nbr_Plan, int NumUndec)
{
    Bool *TabDec;
    set_tabdec(NumUndec, TabDec, Nbr_Plan);
    transform (Imag,TabTrans, Nbr_Plan, TabDec);
    delete [] TabDec;
}

static void fisz_trans(Ifloat &Coeff, Ifloat &Smooth, Bool Dec)
{
   int i,j,Nl = Coeff.nl();
   int Nc= Coeff.nc();
   // cout << "fisz_trans" << endl; 
   if (Dec == False)
   {
     for (i=0; i < Nl;i++)
     for (j=0; j < Nc;j++)
     {
        if (Smooth(i,j) > 0) Coeff(i,j) /= sqrt(Smooth(i,j));
	else Coeff(i,j) = 0;
     }
   }
   else
   {
     for (i=0; i < Nl;i++)
     for (j=0; j < Nc;j++)
     {
        if (Smooth(i/2,j/2) > 0) Coeff(i,j) /= sqrt(Smooth(i/2,j/2));
	else Coeff(i,j) = 0;
     }
   }
}
static void fisz_rec(Ifloat &Coeff, Ifloat &Smooth, Bool Dec)
{
   int i,j,Nl = Coeff.nl();
   int Nc= Coeff.nc();
   if (Dec == False)
   {
     for (i=0; i < Nl;i++)
     for (j=0; j < Nc;j++)
     {
        if (Smooth(i,j) > 0) Coeff(i,j) *= sqrt(Smooth(i,j));
	else Coeff(i,j) = 0;
     }
   }
   else
   {
     for (i=0; i < Nl;i++)
     for (j=0; j < Nc;j++)
     {
        if (Smooth(i/2,j/2) > 0) Coeff(i,j) *= sqrt(Smooth(i/2,j/2));
	else Coeff(i,j) = 0;
     }
   }
}
void HALF_DECIMATED_2D_WT::transform(Ifloat &Imag, Ifloat *TabTrans, 
	      int Nbr_Plan, Bool *TabDec)
{
   int Step = 1;
   PAVE_2D_WT PWT(*Ptr_SB1D_LINE,*Ptr_SB1D_COL);
   SubBand2D SBT(*Ptr_SB1D_LINE,*Ptr_SB1D_COL);
   int NbrBand = 3*(Nbr_Plan-1)+1;
   int Nl = Imag.nl();
   int Nc = Imag.nc();
   Ifloat Ima_Aux(Nl,Nc,"aux");

   for (int s = 0; s < Nbr_Plan-1; s++)
   {
#if SBDEBUG
       cout << "Scale " << s+1 << " Step = " << Step << endl;
       if (s == 0) cout << "  Imag: Nl = " << Imag.nl() << " Nc = " << Imag.nc() <<  " sigma = " <<  Imag.sigma()<<  endl;
       else cout << "  Imag: Nl = " <<   Ima_Aux.nl() << " Nc = " <<   Ima_Aux.nc() <<  " sigma = " << Ima_Aux.sigma() << endl;
#endif   
      if (s == 0)
      {
         if (TabDec[s] == False) 
           PWT.transform2d (Imag, TabTrans[3*s],  TabTrans[3*s+1],  
	            TabTrans[3*s+2], Ima_Aux, Step);
	 else
  	    SBT.transform2d (Imag, False, &(TabTrans[3*s]),  &(TabTrans[3*s+1]),  
	            &(TabTrans[3*s+2]), &Ima_Aux);
      }    
      else 
      {
         Ptr_SB1D_LINE->DistPix = Step;
	 Ptr_SB1D_COL->DistPix = Step;
         if (TabDec[s] == False) 
	   PWT.transform2d (Ima_Aux, TabTrans[3*s],  TabTrans[3*s+1],  
	            TabTrans[3*s+2], Ima_Aux, Step);
         else 
	    SBT.transform2d (Ima_Aux, False, &(TabTrans[3*s]),  &(TabTrans[3*s+1]),  
	            &(TabTrans[3*s+2]), &Ima_Aux);
      }
      if (FiszTrans == True)
      {
         fisz_trans(TabTrans[3*s], Ima_Aux, TabDec[s]);
 	 fisz_trans(TabTrans[3*s+1], Ima_Aux, TabDec[s]);
 	 fisz_trans(TabTrans[3*s+2], Ima_Aux, TabDec[s]);
      }
      if (TabDec[s] == False) Step *= 2;
      
#if SBDEBUG      
      cout << "  Horiz: Nl = " << TabTrans[3*s].nl() << " Nc = " << TabTrans[3*s].nc()  <<  " sigma = " << TabTrans[3*s].sigma() << endl;
      cout << "  Vert: Nl = " << TabTrans[3*s+1].nl() << " Nc = " << TabTrans[3*s+1].nc()  <<  " sigma = " << TabTrans[3*s+1].sigma() << endl;
      cout << "  Diag: Nl = " << TabTrans[3*s+2].nl() << " Nc = " << TabTrans[3*s+2].nc()  <<  " sigma = " << TabTrans[3*s+2].sigma() << endl;
      cout << "  Smooth: Nl = " << Ima_Aux.nl() << " Nc = " << Ima_Aux.nc()  <<  " sigma = " << Ima_Aux.sigma() << endl << endl;
#endif
    } 
    TabTrans[NbrBand-1] = Ima_Aux;
}

void HALF_DECIMATED_2D_WT::recons(Ifloat *TabTrans, Ifloat &Imag, int Nbr_Plan, int NumUndec)
{
    Bool *TabDec;
    set_tabdec(NumUndec, TabDec, Nbr_Plan);
    recons (TabTrans, Imag, Nbr_Plan, TabDec);
    delete []TabDec;
}

void HALF_DECIMATED_2D_WT::recons (Ifloat *TabTrans, Ifloat &Imag, int Nbr_Plan, Bool *TabDec)
{
   int s,Step = 1;
   PAVE_2D_WT PWT(*Ptr_SB1D_LINE,*Ptr_SB1D_COL);
   SubBand2D SBT(*Ptr_SB1D_LINE,*Ptr_SB1D_COL);
   int Nl = Imag.nl();
   int Nc = Imag.nc();
   Ifloat Imag_Aux(Nl,Nc,"aux");
   int *TabNl = new int [Nbr_Plan];
   int *TabNc = new int [Nbr_Plan];
   
   for (s = 0; s < Nbr_Plan-1; s++)
                    if (TabDec[s] == False) Step *= 2;
   TabNl[0] = Nl;
   TabNc[0] = Nc;
   for (s = 1; s < Nbr_Plan; s++)
   {  
      if (TabDec[s-1] == True)
      {
         TabNl[s] = (TabNl[s-1]+1)/2;
	 TabNc[s] = (TabNc[s-1]+1)/2;
      }
      else 
      {
         TabNl[s] = TabNl[s-1];
	 TabNc[s] = TabNc[s-1];
      }
      // cout << "Resol " << s-1 << " size = " << TabNl[s-1] << "  " << TabNc[s-1] << endl;
   }
   
   for (s = Nbr_Plan-2; s >= 0; s--)
   {      
      Ifloat Buff(TabNl[s], TabNc[s], "Buff");
      // Buff.getbuff(Imag.buffer(), TabNc[s], TabNl[s]);
      // cout << "resize to " <<  TabNl[s] << " " << TabNc[s]  << endl;
      if (TabDec[s] == False) Step /= 2;

#if SBDEBUG      
      cout << "Scale " << s+1 << " Step = " << Step << endl;
      cout << "  Horiz: Nl = " << TabTrans[3*s].nl() << " Nc = " << TabTrans[3*s].nc() << " sigma = " << TabTrans[3*s].sigma() << endl;
      cout << "  Vert: Nl = " << TabTrans[3*s+1].nl() << " Nc = " << TabTrans[3*s+1].nc()  << " sigma = " << TabTrans[3*s+1].sigma() << endl;
      cout << "  Diag: Nl = " << TabTrans[3*s+2].nl() << " Nc = " << TabTrans[3*s+2].nc()  << " sigma = " << TabTrans[3*s+2].sigma() << endl;
      if (s == Nbr_Plan-2) cout << "  Smooth: Nl = " << TabTrans[3*s+3].nl() << " Nc = " << TabTrans[3*s+3].nc()  << " sigma = " << TabTrans[3*s+3].sigma() << endl;
      else cout << "  Smooth: Nl = " <<  Imag_Aux.nl() << " Nc = " <<  Imag_Aux.nc()  << " sigma = " <<  Imag_Aux.sigma() << endl;
#endif
      if (FiszRec == True)
      {
         if (s == Nbr_Plan-2)
	 {
            fisz_rec(TabTrans[3*s], TabTrans[3*s+3], TabDec[s]);
	    fisz_rec(TabTrans[3*s+1], TabTrans[3*s+3], TabDec[s]);
	    fisz_rec(TabTrans[3*s+2], TabTrans[3*s+3], TabDec[s]);
         }
	 else
	 {
            fisz_rec(TabTrans[3*s], Imag_Aux, TabDec[s]);
	    fisz_rec(TabTrans[3*s+1], Imag_Aux, TabDec[s]);
	    fisz_rec(TabTrans[3*s+2], Imag_Aux, TabDec[s]);
         }
      }
      
      if (s == Nbr_Plan-2)
      {
         Ptr_SB1D_LINE->DistPix = Step;
	 Ptr_SB1D_COL->DistPix = Step;

         if (TabDec[s] == False)
           PWT.recons2d(TabTrans[3*s],  TabTrans[3*s+1],  TabTrans[3*s+2], 
	                TabTrans[3*s+3], Buff, Step);
         else SBT.recons2d (Buff, False, &(TabTrans[3*s]),  &(TabTrans[3*s+1]),  
	            &(TabTrans[3*s+2]), &(TabTrans[3*s+3]));
      }
      else 
      {
         Ptr_SB1D_LINE->DistPix =  Step;
	 Ptr_SB1D_COL->DistPix = Step;
	 if (TabDec[s] == False)
	    PWT.recons2d (TabTrans[3*s],  TabTrans[3*s+1], TabTrans[3*s+2], 
                     Imag_Aux, Buff, Step);
         else
	    SBT.recons2d(Buff, False, &(TabTrans[3*s]),  &(TabTrans[3*s+1]),  
	            &(TabTrans[3*s+2]), &Imag_Aux);
      }
      if (s==0) Imag = Buff;
      else
      {
         Imag_Aux.resize(Buff.nl(), Buff.nc());
         Imag_Aux = Buff;
      }           
#if SBDEBUG
       cout << "  Imag: Nl = " << Buff.nl() << " Nc = " << Buff.nc() << " sigma = " <<  Buff.sigma() << endl << endl;
#endif     
   }
}

/****************************************************************************/
//                       Fisz Transform
/****************************************************************************/
static int get_nbr_scale (int N) 
{  
    int Max_SIZE_LAST_SCALE = 4;
    int Nmin = N;
    int ScaleMax;
    ScaleMax=iround(log((double)Nmin/(double)Max_SIZE_LAST_SCALE) / log(2.)+ 1.);
    return (ScaleMax);
}

void fisz2d(Ifloat &Data, Ifloat &Recons,  Bool Reverse, int NbrScale, 
                  type_sb_filter Fil, int NumUndec)
{
  // Fil = F_MALLAT_7_9;
  // Fil = F_BI2HAAR;
  int Nl = Data.nl();
  int Nc = Data.nc();
  SubBandFilter  USF(Fil, NORM_L1);
  Ifloat *TabTrans;
  HALF_DECIMATED_2D_WT WT(USF);
  int Nbr_Plan = NbrScale;
  if (Nbr_Plan < 2) Nbr_Plan = get_nbr_scale(MIN(Nl,Nc));
  
  int NbrBand = WT.alloc (TabTrans, Nl, Nc, Nbr_Plan, NumUndec);
  // cout << "FISZ TRANS: Nbr_Plan = " << Nbr_Plan << " NbrBand = "  << NbrBand << endl; 				  
  if (Reverse == False) WT.FiszTrans = True;
  else WT.FiszRec = True;
  
  WT.transform(Data,TabTrans, Nbr_Plan, NumUndec);
  if ((Recons.nl() != Nl) || (Recons.nc() != Nc)) Recons.alloc(Nl,Nc);
  else Recons.init();
  WT.recons(TabTrans, Recons, Nbr_Plan, NumUndec);
  WT.free(TabTrans, Nbr_Plan);
  NbrBand=0; // to remove the warning when compiling
}
  
/****************************************************************************/


/****************************************************************************/
//                   A TROUS ALGORITHM
/****************************************************************************/

void ATROUS_2D_WT::b3spline_filtering(Ifloat & Im_in, Ifloat &Im_out,  int Step_trou, Bool GTilde, int nb_thr)
{
    int Nl = Im_in.nl();
    int Nc = Im_in.nc();
    int i,j,Step;
    double Coeff_h0 =  3. / 8.;
    double Coeff_h1 = 1. / 4.;
    double  Coeff_h2 = 1. / 16.;
    Ifloat Buff(Nl,Nc,"Buff smooth_bspline");
    type_border Type = Bord;
    double Val;
    
    Step = (int)(pow((double)2., (double) Step_trou) + 0.5);
    
    #ifdef _OPENMP
    if(nb_thr <1) nb_thr=1;
    #pragma omp parallel for private(i,j,Val) schedule(static) shared(Nl,Nc,Coeff_h0,Coeff_h1,Coeff_h2,Step,Type,Buff) num_threads(nb_thr)
    #endif
    for (i = 0; i < Nl; i ++)
        for (j = 0; j < Nc; j ++)
        {
            Val = Coeff_h0 * (double) Im_in(i,j)
            + Coeff_h1 * (  Im_in (i, j-Step, Type) 
                          + Im_in (i, j+Step, Type)) 
            + Coeff_h2 * (  Im_in (i, j-2*Step, Type) 
                          + Im_in (i, j+2*Step, Type));
            Buff(i,j) = (float) Val;
        }
    
    #ifdef _OPENMP
    #pragma omp parallel for private(i,j) schedule(static) shared(Nl,Nc,Coeff_h0,Coeff_h1,Coeff_h2,Step,Type,Buff) num_threads(nb_thr)
    #endif
    
    for (i = 0; i < Nl; i ++)
        for (j = 0; j < Nc; j ++)
            Im_out(i,j) = Coeff_h0 * Buff(i,j)
            + Coeff_h1 * (  Buff (i-Step, j, Type) 
                          + Buff (i+Step, j, Type)) 
            + Coeff_h2 * (  Buff (i-2*Step, j, Type) 
                          + Buff (i+2*Step, j, Type));
    
    if (GTilde == True)  Im_out += Im_in;
}


/****************************************************************************/

void ATROUS_2D_WT::transform(Ifloat &Image, Ifloat *TabBand, int NbrScale, int nb_thr)
{
    int s;
    TabBand[0] = Image;
    Ifloat Data_out;
    
    if (ModifiedAWT == True) Data_out.alloc(Image.nl(), Image.nc()," Buff ATW");
    
    for (s = 0; s < NbrScale-1; s++)
    {
        b3spline_filtering(TabBand[s], TabBand[s+1], s,nb_thr);
        if (ModifiedAWT == True) 
        {
            b3spline_filtering(TabBand[s+1], Data_out, s,nb_thr);
            TabBand[s] -= Data_out;
        }
        else TabBand[s] -= TabBand[s+1];
    }
}

/****************************************************************************/

void ATROUS_2D_WT::recons(Ifloat *TabBand, Ifloat &Image, int NbrScale, Bool AddLastScale, int nb_thr)
{
    int i,j,s;
    
    if(nb_thr <1) nb_thr=1;
    
    if ((ModifiedAWT == False) && (AdjointRec == False))
    { 
        int Last_Scale_Used = (AddLastScale == True) ? NbrScale : NbrScale-1;
        Image = TabBand[0];
        for (s = 1; s < Last_Scale_Used; s++) 
    #ifdef _OPENMP
    #pragma omp parallel for private(i,j) schedule(static) shared(s) num_threads(nb_thr)
    #endif
            for (i=0; i < Image.nl(); i++) {
		     	for (j=0; j < Image.nc(); j++) Image(i,j) += (TabBand[s])(i,j);
            }
    } 
    else
    {
        Ifloat temp(Image.nl(), Image.nc(), "Rec adjoint pave");
        if (AddLastScale == True) Image = TabBand[NbrScale-1];
        else Image.init();
        
        if (ModifiedAWT == True)
        {
            //cout << " MOD REC" << endl;
            for (s=NbrScale-2; s>= 0 ; s--)
            {
                b3spline_filtering (Image, temp, s,nb_thr);
                
    #ifdef _OPENMP
    #pragma omp parallel for private(i,j) schedule(static) shared(s) num_threads(nb_thr)
    #endif
                for (i=0; i < Image.nl(); i++) {
                    for (j=0; j < Image.nc(); j++) Image(i,j) = temp(i,j) +  (TabBand[s])(i,j);
                }
            }
        }
        else 
        {
            //cout << " ADJOINT REC" << endl;
            for (s=NbrScale-2; s>= 0 ; s--)
            {
                // int Step = (int)(pow((double)2., (double) s) + 0.5);
                // cout << " Scale " << s << ": Step = " << Step <<  endl;
                b3spline_filtering (Image, temp, s,nb_thr);
                b3spline_filtering (TabBand[s], Image, s, True,nb_thr);
    #ifdef _OPENMP
    #pragma omp parallel for private(i,j) schedule(static) num_threads(nb_thr)
    #endif
                for (i=0; i < Image.nl(); i++) {
                    for (j=0; j < Image.nc(); j++) Image(i,j) += temp(i,j) ;
                }
            }
        }
    }
}

/****************************************************************************/

void ATROUS_2D_WT::free(Ifloat *TabBand, int Nbr_Plan)
{
    if (Nbr_Plan != 0) delete [] TabBand;
}

/****************************************************************************/

void ATROUS_2D_WT::alloc (Ifloat * & TabBand, int Nl, int Nc, int NbrBand)
{
    char ch[80];
    TabBand = new Ifloat [NbrBand];
    for (int s = 0; s < NbrBand; s++)
    {
        sprintf (ch, "band_%d", s+1);
        TabBand[s].alloc (Nl, Nc, ch);   
    }
}

/****************************************************************************/

float ATROUS_2D_WT::norm_band(int s)
{
    static double TN[10] = 
    {0.889434,  0.200105, 0.0857724, 0.0413447, 0.0202689, 0.00995628, 0.00513504,
        0., 0., 0.};
    if ((s < 0) || (s >= 10)) return 0;
    else return TN[s];
}

/****************************************************************************/
//                   PMT ALGORITHM
/****************************************************************************/

void PMT_2D::transform(Ifloat &Image, Ifloat *TabBand, int NbrScale, type_border Bord)
{
   int s;
   int Nl = Image.nl();
   int Nc = Image.nc();
   Ifloat Buff(Nl, Nc, "buff pmt");

   TabBand[0] = Image;
   for (s = 0; s <  NbrScale-1; s++)
   {
        Buff.resize((TabBand[s]).nl(), (TabBand[s]).nc());
        smooth_mediane (TabBand[s], Buff, Bord, 0, MedianWindowSize);
        TabBand[s] -= Buff;
        im_reduce_size_2 (Buff, TabBand[s+1]);
   }
}

/****************************************************************************/

void PMT_2D::recons(Ifloat *TabBand, Ifloat &Image, int NbrScale, type_border Bord)
{
   int s,i,j;
   int Nl = (TabBand[0]).nl();
   int Nc = (TabBand[0]).nc();
   Ifloat Buffer(Nl, Nc, "buff pmt");

   if ((Nl != Image.nl()) || (Nc != Image.nc())) Image.resize(Nl,Nc);
   Nl = (TabBand[NbrScale-1]).nl();  
   Nc = (TabBand[NbrScale-1]).nc();  
   Image.resize (Nl, Nc);
   Image = TabBand[NbrScale-1];
   for (s =  NbrScale -2; s >= 0; s--)
   {
      Nl = (TabBand[s]).nl();  
      Nc = (TabBand[s]).nc();  
      Buffer.resize (Nl, Nc);
      im_increase_size_2 (Image, Buffer, Bord);
      Image.resize (Nl, Nc);
      for (i = 0; i < Nl; i++)
      for (j = 0; j < Nc; j++)
                      Image(i,j) = Buffer(i,j) + (TabBand[s])(i,j);
   }
}

/****************************************************************************/

void PMT_2D::free(Ifloat *TabBand, int Nbr_Plan)
{
    if (Nbr_Plan != 0) delete [] TabBand;
}

/****************************************************************************/

void PMT_2D::alloc (Ifloat * & TabBand, int Nl, int Nc, int NbrBand)
{
    char ch[80];
    int Nls = Nl;
    int Ncs = Nc;
    TabBand = new Ifloat [NbrBand];
    for (int s = 0; s < NbrBand; s++)
    {
	sprintf (ch, "band_%d", s+1);
        TabBand[s].alloc (Nls, Ncs, ch);   
        Nls = (Nls+1) / 2;
        Ncs = (Ncs+1) / 2;
    }
}

/****************************************************************************/

float PMT_2D::norm_band(int s)
{
   static double TN[10] = 
    {0.970291,  0.338748, 0.178373, 0.0992649, 0.0539517, 0.0317369, 0.0188998,
      0.013650, 0., 0.}; 
   if ((s < 0) || (s >= 10)) return 0;
   else return TN[s];
}

/****************************************************************************/

/****************************************************************************/
//                   Qunincunx Orthogonal transform
/****************************************************************************/

void Quincunx::band_to_ima(Ifloat *TabTrans, Ifloat &Imag, int NbrBand)
{
   int s,i,j,ind_i,ind_j;
   int Nl = Imag.nl();
   int Nc = Imag.nc();
   for (s = 0; s < NbrBand-1; s++)
   for (i = 0; i < TabTrans[s].nl(); i++)
   for (j = 0; j < TabTrans[s].nc(); j++)
   {
      ind_quinc(s,i,j,Nl,Nc,ind_i,ind_j);
      Imag(ind_i,ind_j) =  TabTrans[s](i,j);
   }
   //  s = NbrBand-1;
   //  for (i = 0; i < TabTrans[s].nl(); i++)
   //  for (j = 0; j < TabTrans[s].nc(); j++) Imag(i,j) = TabTrans[s](i,j);
}

void Quincunx::ind_quinc(int Band, int i,int j,  
                      int Nl, int Nc, int &ind_i, int &ind_j)
{
   ind_i = i;
   ind_j = j;
   int Nls = Nl; //size_ima_resol(s+1, Nl);
   int Ncs = Nc; //size_ima_resol(s+1, Nc);
   int Nscale = Band / 2;
   for (int s=0; s <= Nscale; s++)
   {
       Nls = (Nls+1)/2;
       Ncs = (Ncs+1)/2;
   }    
   if (Band % 2 == 0) ind_j += Ncs;
   else ind_i += Nls; 
}


void Quincunx::free(Ifloat *TabBand, int Nbr_Plan)
{
    if (Nbr_Plan != 0) delete [] TabBand;
}
int Quincunx::alloc (Ifloat * & TabBand, int Nl, int Nc, int Nbr_Plan )
{
    char ch[80];
    int Nl_s = Nl;
    int Nc_s = Nc;
    int NbrBand_per_Resol = 2;
    int s;
    int NbrBand = NbrBand_per_Resol*(Nbr_Plan-1)+1;
    TabBand = new Ifloat [NbrBand];
    for (s = 0; s < NbrBand-1; s+=2)
    {
  	sprintf (ch, "band_%d", s+1);
        TabBand[s].alloc (Nl_s,  Nc_s/2, ch);   
        sprintf (ch, "band_%d", s+2);
        TabBand[s+1].alloc (Nl_s/2,  (Nc_s+1)/2, ch);
        Nl_s = (Nl_s+1)/2;
        Nc_s = (Nc_s+1)/2;
     }
     s = NbrBand-1;
     sprintf (ch, "band_%d", s+1);
     TabBand[s].alloc(Nl_s, Nc_s, ch);
     return NbrBand;
}

void Quincunx::transform(Ifloat &Data, Ifloat * & TabBand, int Nbr_Plan)
{
   int s,Nb;
   Ifloat Smooth(Data.nl()/2, Data.nl()/2, "Smooth");
   if (TabBand == NULL) Nb = alloc(TabBand, Data.nl(), Data.nc(), Nbr_Plan);
   else Nb = 2*Nbr_Plan-1;
   for (s = 0; s < Nbr_Plan-1; s++)
   {
       // cout << "Step " << s+1 << endl;
       if (s == 0) 
         transform_one_step (Data, TabBand[0],  TabBand[1], Smooth);
       else 
         transform_one_step (Smooth, TabBand[2*s],  
                             TabBand[2*s+1], Smooth);
   }
   TabBand[Nb-1] = Smooth;
   // cout << "ET" << endl;
}

void Quincunx::recons(Ifloat *TabBand, Ifloat &Imag, int Nbr_Plan)
{
    int s;   
    Ifloat ImaRec;
    for (s = Nbr_Plan-2; s >= 0; s--)
    {
       if (s == Nbr_Plan-2) 
            recons_one_step(TabBand[2*s+2], TabBand[2*s+1], 
	                          TabBand[2*s],  Imag);
       else recons_one_step(ImaRec, TabBand[2*s+1], 
	                          TabBand[2*s], Imag);
      if (s != 0) ImaRec = Imag;
    }
}

void Quincunx::transform_one_step (Ifloat & Data, Ifloat &Half, Ifloat & Resol, 
                                  Ifloat & Smooth)
{
    int Nl = Data.nl();
    int Nc = Data.nc();
    int Nl1 = Nl/2;
    int Nc1 = Nc/2;
    int Nl2 = (Nl+1)/2;
    int Nc2 = (Nc+1)/2;

    // cout << " Half" << Half.nl()  << " " <<  Half.nc() << endl;
    // cout << " Resol " << Resol .nl()  << " " <<   Resol.nc() << endl;
    // cout << " Smooth " <<  Smooth.nl()  << " " << Smooth.nc() << endl;

    if ((Half.nl() != Nl) || (Half.nc() != Nc1))      Half.resize(Nl,Nc1);
    if ((Resol.nl() != Nl1) || (Resol.nc() != Nc2))   Resol.resize(Nl1,Nc2);

    int j, i;
   
    // each line is convolved by h and g
    // one pixel out of two (in column) is kept 
    Ifloat image_h0(Nl,Nc2,"image_h0");

    for (i = 0; i < Nl; i++) 
    {
       fltarray Line (Nc);
       fltarray LineH(Nc2);
       fltarray LineG(Nc2);
       float *PtrDet = Half.buffer() + Nc2*i;
       if (i % 2 == 0)
       {
          float *PtrHigh = Data.buffer() + Nc*i;
          float *PtrLow = image_h0.buffer() + Nc2*i;
          Ptr_SB1D->transform(Nc, PtrHigh, PtrLow, PtrDet);
       }
       else
       {
          Line(0) = Data(i,Nc-1);
          for (j = 1; j < Nc; j++) Line(j) = Data(i,j-1);
          Ptr_SB1D->transform(Nc, Line.buffer(), LineH.buffer(), LineG.buffer());
          for (j = 0; j < Nc2-1; j++) image_h0(i,j) = LineH(j+1);
          for (j = 0; j < Nc/2; j++) Half(i,j) = LineG(j);
          image_h0(i,Nc2-1) = LineH(0);  
       }
    }
    // cout << " col " << endl;
    // io_write_ima_float("xx_tr1.fits", image_h0);
    if ((Smooth.nl() != Nl2) || (Smooth.nc() != Nc2)) Smooth.resize(Nl2,Nc2);

    //cout << "transform_one_step" << Nl << " " << Nc << endl;
    //cout << " Half" << Half.nl()  << " " <<  Half.nc() << endl;
    //cout << " Resol " << Resol .nl()  << " " <<   Resol.nc() << endl;
    //cout << " Smooth " <<  Smooth.nl()  << " " << Smooth.nc() << endl;


    for (j = 0; j < Nc2; j++) 
    {
       fltarray Col(Nl);
       fltarray ColG(Nl/2+1);
       fltarray ColH(Nl/2+1);
       for (i = 0; i < Nl; i++) Col(i) = image_h0(i,j);
       Ptr_SB1D->transform(Nl, Col.buffer(), ColH.buffer(), ColG.buffer());
       for (i = 0; i < Nl/2; i++) Resol(i,j) = ColG(i);
       for (i = 0; i < (Nl+1)/2; i++) Smooth(i,j) = ColH(i);
    }   
    // cout << " end one step " << endl;
}


void Quincunx::recons_one_step(Ifloat &Smooth, Ifloat & Resol, 
	                      Ifloat & Half, Ifloat &ImaRec)
{
   int i,j;
   int Nl = Half.nl();
   int Nc = Resol.nc()+Half.nc();
   SubBand1D *SB1D = get_subband_method();
   int Nl_2 = (Nl+1)/2;
   int Nc_2 = (Nc+1)/2;
   if ((ImaRec.nl() != Nl) || ( ImaRec.nc() != Nc)) ImaRec.resize(Nl,Nc);

   for (j = 0; j < Nc_2; j++) 
   {
      fltarray HighCol(Nl);
      fltarray LowCol(Nl_2);
      fltarray DetCol(Nl_2);
      for (i = 0; i < Nl_2; i++) LowCol(i) = Smooth(i,j);
      for (i = 0; i < Nl/2; i++) DetCol(i) = Resol(i,j);
      SB1D->recons (Nl, LowCol.buffer(), DetCol.buffer(), HighCol.buffer());
      for (i = 0; i < Nl; i++) ImaRec(i,j) = HighCol(i);
   }

   // io_write_ima_float("xx_rec1.fits", ImaRec);
   for (i = 0; i < Nl; i++) 
   {
      fltarray LowLine(Nc_2+1);
      fltarray DetLine(Nc_2+1);
      fltarray HighLine(Nc+1);
      if (i % 2 == 0)
      {
         float *PtrHigh = ImaRec.buffer() + Nc*i;
         for (j = 0; j < Nc_2; j++) LowLine(j) = ImaRec(i,j);
         for (j = 0; j < Nc/2; j++) DetLine(j) =  Half(i,j);
         SB1D->recons(Nc, LowLine.buffer(), DetLine.buffer(), PtrHigh);
      }
      else
      {
         for (j = 1; j <= Nc_2; j++) LowLine(j) = ImaRec(i,j-1);
         LowLine(0) =  ImaRec(i,Nc_2-1);
         for (j = 0; j < Nc/2; j++) DetLine(j) = Half(i,j);
         SB1D->recons(Nc, LowLine.buffer(), DetLine.buffer(), HighLine.buffer());    
         for (j = 0; j < Nc-1; j++) ImaRec(i,j) = HighLine (j+1);
         ImaRec(i,Nc-1) = HighLine (0);
      }
   }   
}

/****************************************************************************/
