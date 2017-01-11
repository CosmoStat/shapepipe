/*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  98/05/15 
**    
**    File:  MR_Obj.cc
**
*******************************************************************************
**
**    DESCRIPTION  MR_Obj to an image 
**    -----------  
**                 
*******************************************************************************/#include "MR_Obj.h"

#include "MR_Obj.h"

static void ind_orthog_quinc(int Scale, int i,int j, details which_detail, 
                       int Nl, int Nc, int &ind_i, int &ind_j);

/************************************************************************/

void ind_orthog_transf(int Scale, int i,int j, details which_detail, 
                       int Nl, int Nc, int &ind_i, int &ind_j)
{
   ind_i = i;
   ind_j = j;
   int s,Nls,Ncs;
   
   Nls = Nl;
   Ncs = Nc;
   for (s=0; s <= Scale; s++)
   {
      Nls = (Nls+1) / 2;
      Ncs = (Ncs+1) / 2;
   }
    
    switch (which_detail)
    {
        case D_HORIZONTAL:
               ind_j += Ncs;
               break;
        case D_DIAGONAL:  
               ind_i += Nls;
               ind_j += Ncs;
             break;
        case D_VERTICAL:  
               ind_i += Nls;
               break;
        case I_SMOOTH: break;
        case D_NULL:
        default: cerr << "Error: unknown detail" << endl;
                 exit(0);
                 break;
   }
}

/************************************************************************/

void pos_orthog_transf(int Num, int & s, int & i,int & j, details & which_detail, int Nl, int Nc, int NbrScale)
{
   int Val=Num, Nls=Nl, Ncs=Nc;
   int Nls_L = (Nls+1)/2;
   int Ncs_L = (Ncs+1)/2;
   int Nls_H =  Nls_L;
   int Ncs_H =  Ncs /2;
   int Nls_V =  Nls /2;
   int Ncs_V =  Ncs_L;
   int Nls_D =  Nls /2;
   int Ncs_D =  Ncs /2;
   
   s = 0;
   int SizeB3 = Nls_H*Ncs_H + Nls_V*Ncs_V+ Nls_D*Ncs_D;
   while (Val >= SizeB3)
   {
      s++;
      Val -= SizeB3;
      Nls = (Nls+1)/2;
      Ncs = (Ncs+1)/2;
      Nls_L = (Nls+1)/2;
      Ncs_L = (Ncs+1)/2;
      Nls_H =  Nls_L;
      Ncs_H =  Ncs /2;
      Nls_V =  Nls /2;
      Ncs_V =  Ncs_L;
      Nls_D =  Nls /2;
      Ncs_D =  Ncs /2;
      SizeB3 = Nls_H*Ncs_H+ Nls_V*Ncs_V+ Nls_D*Ncs_D;      
   }
    
   if (s >= NbrScale-1)
   {
       which_detail = I_SMOOTH;
       s = NbrScale-2;
       i = Val / Ncs;
       j = Val - i*Ncs;
   }
   else
   {
       int Size_H = Nls_H*Ncs_H;
       int Size_V = Nls_V*Ncs_V;
       
       if (Val < Size_H)
       {
          which_detail = D_HORIZONTAL;
          i = Val / Ncs_H;
          j = Val - i*Ncs_H;
       }
       else if (Val <  Size_H+Size_V)
       {
          which_detail = D_VERTICAL;
          Val -=  Size_H;
          i = Val / Ncs_V;
          j = Val - i*Ncs_V;
       }
       else 
       {
          which_detail = D_DIAGONAL;
          Val -=  Size_H+Size_V;
          i = Val / Ncs_D;
          j = Val - i*Ncs_D;
       }
   }
}

/************************************************************************/

/*static void pos_pyr_transf(int Num, int & s, int & i,int & j, int Nl, int Nc,
                           Bool SemiPyr)
{   
   int Val=Num, Nls=Nl, Ncs=Nc;

   s = 0;
   while (Val - Nls*Ncs >= 0)
   {
       s++;
       Val -= Nls*Ncs;
       if ((SemiPyr != True) || (s != 1))
       {
          Nls = (Nls+1)/2;
          Ncs = (Ncs+1)/2;
       }
   }
   i = Val / Ncs;
   j = Val - i*Ncs;
}*/

/************************************************************************/

static void ind_orthog_quinc(int Scale, int i,int j, details which_detail, 
                       int Nl, int Nc, int &ind_i, int &ind_j)
{
   ind_i = i;
   ind_j = j;
   int Nls = Nl; //size_ima_resol(s+1, Nl);
   int Ncs = Nc; //size_ima_resol(s+1, Nc);
    for (int s=0; s <= Scale; s++)
    {
        Nls = (Nls+1)/2;
        Ncs = (Ncs+1)/2;
    }      
    switch (which_detail)
    {
       case D_HALF_RESOL:ind_j +=  Ncs;break;
       case D_RESOL:ind_i += Nls;break;
       case I_SMOOTH:break;
       default: cerr << "Error: unknown detail" << endl;
                    exit(0); 
                    break;
    }
}

 
/************************************************************************/

void ortho_trans_to_ima(MultiResol &MR_Data, Ifloat &Ima, int Resol)
{
   int NbrScale = MR_Data.nbr_scale();
   int Nbrband = MR_Data.nbr_band();
   int i=0,j=0,ind_i=0,ind_j=0,s = NbrScale-1;
   int b = Nbrband - 1;
   details which_detail;
   int Nli = MR_Data.size_ima_nl();
   int Nci = MR_Data.size_ima_nc();
   int Nl = size_ima_resol(Resol, Nli);
   int Nc = size_ima_resol(Resol, Nci);
   Ima.resize(Nl, Nc);
   
   if ((Resol < 0) || (Resol >= NbrScale))
   {
      cerr << "Error: bad resolution number in ortho_trans_to_ima ... " << endl;
      exit(-1);
   }
   if ((MR_Data.Set_Transform != TRANSF_MALLAT) &&
      (MR_Data.Set_Transform !=  TRANSF_FEAUVEAU))
   {
      cerr << "Error: bad transformation type in ortho_trans_to_ima ... " << endl;
      exit(-1);
   }   
   
   MR_Data.band_to_scale(b, s, which_detail);
 		    
   while ((s >= Resol) || ((s == Resol-1) && (which_detail == I_SMOOTH)))
   {
      for (i=0; i < MR_Data.size_band_nl(b); i++)
      for (j=0; j < MR_Data.size_band_nc(b); j++)
      {
         if (MR_Data.Set_Transform == TRANSF_MALLAT)
              ind_orthog_transf(s, i, j, which_detail, Nli, Nci, ind_i, ind_j);
         else ind_orthog_quinc (s, i, j, which_detail, Nli, Nci, ind_i, ind_j);        
         Ima(ind_i,ind_j) = MR_Data(b,i,j);           
      }    
      b --;
      MR_Data.band_to_scale(b, s, which_detail);
      if (b  < 0) s = -1;
   }
}

/************************************************************************/

void ima_to_ortho_trans(MultiResol &MR_Data, Ifloat &Ima, int Resol)
{
   int NbrScale = MR_Data.nbr_scale();
   int Nbrband = MR_Data.nbr_band();
   int s = NbrScale-1;
   int i=0,j=0,ind_i=0,ind_j=0,b = Nbrband - 1;
   details which_detail;
   int Nli = MR_Data.size_ima_nl();
   int Nci = MR_Data.size_ima_nc();
   int Nl = size_ima_resol(Resol, Nli);
   int Nc = size_ima_resol(Resol, Nci);
   if ((Nl != Ima.nl()) || (Nc != Ima.nc()))
   {
      cerr << "Error: bad image size in  ima_to_ortho_trans ... " << endl;
      exit(-1);
   }
    
   if ((Resol < 0) || (Resol >= NbrScale))
   {
      cerr << "Error: bad resolution number in  ima_to_ortho_trans ... " << endl;
      exit(-1);
   }
   if ((MR_Data.Set_Transform != TRANSF_MALLAT) &&
      (MR_Data.Set_Transform  !=  TRANSF_FEAUVEAU))
   {
      cerr << "Error: bad transformation type in  ima_to_ortho_trans ... " << endl;
      exit(-1);
   }   
   
   MR_Data.band_to_scale(b, s, which_detail);
   while ((s >= Resol) || ((s == Resol-1) && (which_detail == I_SMOOTH)))
   {
      for (i=0; i < MR_Data.size_band_nl(b); i++)
      for (j=0; j < MR_Data.size_band_nc(b); j++)
      {
         if (MR_Data.Set_Transform == TRANSF_MALLAT)
              ind_orthog_transf(s, i, j, which_detail, Nli, Nci, ind_i, ind_j);
         else ind_orthog_quinc (s, i, j, which_detail, Nli, Nci, ind_i, ind_j);
         MR_Data(b,i,j) = Ima(ind_i,ind_j);    
      }          
      b --;
      MR_Data.band_to_scale(b, s, which_detail);
      if (b  < 0) s = -1;
    }
}
