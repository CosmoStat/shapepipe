/*******************************************************************************
**                   Copyright (C) 2003 by CEA
*******************************************************************************
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  03/04/03 
**    
**    File:  LineCol.cc
**
**    Modification history:
**
*******************************************************************************
**
**    DESCRIPTION  Line Column Multiscale  decomposition
**    -----------  
**                 
******************************************************************************/
 

#include "IM_Obj.h"
#include "LineCol.h"
 #include "IM_IO.h"

/****************************************************************************/
//                   LINE COLUMN Orthogonal transform
/****************************************************************************/

void LineCol::transform(Ifloat &Data, Ifloat &Result, int NbrScaleLine, int NbrScaleCol)
{
   Result = Data;
   transform(Result, NbrScaleLine, NbrScaleCol);
}

/****************************************************************************/

void LineCol::recons(Ifloat &Trans, Ifloat &Result, int NbrScaleLine, int NbrScaleCol)
{
   Result = Trans;
   recons(Result, NbrScaleLine, NbrScaleCol);
}

/****************************************************************************/

void LineCol::transform(Ifloat &Data, int NbrScaleLine, int NbrScaleCol)
{
   int s, Dep=0;
   int Nl = Data.nl();
   int Nc = Data.nc();
   
   for (s = 0; s < NbrScaleLine-1; s++)
   {
      transform_one_step_line (Data,  Data.nl(), Nc, Dep);
      if (Mirror==False)  Nc = (Nc+1) / 2;
      else
      {
  	  Dep += (Nc+1) / 2;
	  Nc /= 2;
      }
   }
   
   int NScale = (NbrScaleCol < 0) ? NbrScaleLine : NbrScaleCol;
   Dep = 0;
   for (s = 0; s < NScale-1; s++)
   {
      transform_one_step_col (Data, Nl, Data.nc(), Dep);
      if (Mirror==False)  Nl = (Nl+1) / 2;
      else
      {
  	  Dep += (Nl+1) / 2;
	  Nl /= 2;
      }
   }
}

/****************************************************************************/

void LineCol::undec_transform(Ifloat &Data, Ifloat ** & TabTrans, int NbrScaleLine, int NbrScaleCol)
{
   int i,j,s, Step=1;
   int Nl = Data.nl();
   int Nc = Data.nc();
  
   TabTrans = new Ifloat * [NbrScaleLine];
   for (i=0;i < NbrScaleLine; i++) TabTrans[i] = new Ifloat [NbrScaleCol];
   for (i=0;i < NbrScaleLine; i++) 
   for (j=0;j < NbrScaleCol; j++) (TabTrans[i][j]).alloc(Nl, Nc);
   
   Ifloat Buff;
   if (NbrScaleLine == 1) TabTrans[0][0] = Data;
   else Buff = Data;
   for (s = 0; s < NbrScaleLine-1; s++)
   { 
      Step = POW2(s);
      undec_transform_one_step_line (Buff, TabTrans[s+1][0], TabTrans[s][0], Step);
      if (s != NbrScaleLine-2) Buff = TabTrans[s+1][0];
    }
   // cout << "COL " << NScale << endl;
   int NScale = (NbrScaleCol < 0) ? NbrScaleLine : NbrScaleCol;
   for (i = 0; i < NbrScaleLine; i++)
   {
      Buff=TabTrans[i][0];
      for (s = 0; s < NScale-1; s++)
      {
         Step = POW2(s);
         undec_transform_one_step_col (Buff, TabTrans[i][s+1], TabTrans[i][s], Step);
         // INFO(TabTrans[i][s+1], "smooth");
         if (s != NScale-2) Buff = TabTrans[i][s+1];
      }
   }
}

/****************************************************************************/

void LineCol::undec_transform(Ifloat &Data, Ifloat *& TabTrans, int NbrScaleLine, int NbrScaleCol)
{
   int i,s, Step=1;
//   int Nl = Data.nl();
//   int Nc = Data.nc();
//   int IndBand=0;
  //  cout << "   NbrScaleCol " << NbrScaleCol << endl;
   Ifloat Buff;
   if (NbrScaleLine == 1) TabTrans[0] = Data;
   else Buff = Data;
 //  cout << "  NbrScaleLine " << NbrScaleLine << endl;
   for (s = 0; s < NbrScaleLine-1; s++)
   { 
 //  cout << " undec_transform " << s+1 << endl;
      Step = POW2(s);
      undec_transform_one_step_line (Buff, TabTrans[(s+1)*NbrScaleCol], TabTrans[s*NbrScaleCol], Step);
      if (s != NbrScaleLine-2) Buff = TabTrans[(s+1)*NbrScaleCol];
   }
   // cout << "COL " << NScale << endl;
   int NScale = (NbrScaleCol < 0) ? NbrScaleLine : NbrScaleCol;
   for (i = 0; i < NbrScaleLine; i++)
   {
//   cout << " undec_transform i " << i+1 << endl;
      Buff=TabTrans[i*NbrScaleCol];
      for (s = 0; s < NScale-1; s++)
      {
         Step = POW2(s);
         undec_transform_one_step_col (Buff, TabTrans[i*NbrScaleCol+s+1], TabTrans[i*NbrScaleCol+s], Step);
         // INFO(TabTrans[i][s+1], "smooth");
         if (s != NScale-2) Buff = TabTrans[i*NbrScaleCol+s+1];
      }
   }
}

/****************************************************************************/

void LineCol::undec_recons(Ifloat *& TabTrans, Ifloat &Data, int NbrScaleLine, int NbrScaleCol)
{
   int i,s, Step=1;
//   int Nl = Data.nl();
//   int Nc = Data.nc();
   
   Ifloat Buff;
   
   int NScale = (NbrScaleCol < 0) ? NbrScaleLine : NbrScaleCol;
   if (NScale > 1)
   {
      for (i = 0; i < NbrScaleLine; i++)
      {
         Buff=TabTrans[i*NbrScaleCol+NScale-1];
         for (s = NScale-2; s >= 0; s--)
         {
           Step = POW2(s);
           // INFO(Buff, "rec smooth");
           undec_recons_one_step_col(Buff, TabTrans[i*NbrScaleCol+s], Data, Step);
           if (s != 0) Buff = Data;
         }
         TabTrans[i*NbrScaleCol] = Data;
       }
   }
   if (NbrScaleLine > 1)  
   {
      Buff = TabTrans[(NbrScaleLine-1)*NbrScaleCol];
      for (s = NbrScaleLine-2; s >= 0; s--)
      {
        Step = POW2(s);
        undec_recons_one_step_line (Buff, TabTrans[s*NbrScaleCol], Data, Step);
        if (s != 0) Buff = Data;
     }
  }
}

/****************************************************************************/

 
void LineCol::undec_recons(Ifloat ** & TabTrans, Ifloat &Data, int NbrScaleLine, int NbrScaleCol)
{
   int i,s, Step=1;
//   int Nl = Data.nl();
//   int Nc = Data.nc();
   
   Ifloat Buff;
   
   int NScale = (NbrScaleCol < 0) ? NbrScaleLine : NbrScaleCol;
   if (NScale > 1)
   {
      for (i = 0; i < NbrScaleLine; i++)
      {
         Buff=TabTrans[i][NScale-1];
         for (s = NScale-2; s >= 0; s--)
         {
           Step = POW2(s);
           // INFO(Buff, "rec smooth");
           undec_recons_one_step_col(Buff, TabTrans[i][s], Data, Step);
           if (s != 0) Buff = Data;
         }
         TabTrans[i][0] = Data;
       }
   }
   if (NbrScaleLine > 1)  
   {
      Buff = TabTrans[NbrScaleLine-1][0];
      for (s = NbrScaleLine-2; s >= 0; s--)
      {
        Step = POW2(s);
        undec_recons_one_step_line (Buff, TabTrans[s][0], Data, Step);
        if (s != 0) Buff = Data;
     }
  }
}

/****************************************************************************/

void LineCol::recons(Ifloat & Data, int NbrScaleLine, int NbrScaleCol)
{
    int s,s1, Dep=0;  
    int Nl = Data.nl();
    int Nc = Data.nc(); 
    int Nls,Ncs;
    int NScale = (NbrScaleCol < 0) ? NbrScaleLine : NbrScaleCol;
    for (s = NScale-2; s >= 0; s--)
    {
        if (Mirror==False)  Nls = size_resol(s, Nl);
	else 
	{
	   Nls = Nl;
	   Dep = 0;
 	   for (s1 = 0; s1 < s; s1++) 
	   {
	     Dep += (Nls+1) / 2;
	     Nls /= 2;
	   }
	}
        recons_one_step_col (Data, Nls, Nc, Dep);
    }
    for (s = NbrScaleLine-2; s >= 0; s--)
    {
        if (Mirror==False)  Ncs = size_resol(s, Nc);
 	else 
	{
	   Ncs = Nc;
	   Dep = 0;
 	   for (s1 = 0; s1 < s; s1++) 
	   {
	     Dep += (Ncs+1) / 2;
	     Ncs /= 2;
	   }
	}
        recons_one_step_line (Data, Nl, Ncs, Dep);
    }
}


/****************************************************************************/

void LineCol::transform_one_step_line (Ifloat & Data, int Nl, int Nc, int Dep)
{
    int Nc2 = (Nc+1)/2;
    int j,i;
    // cout << "Line step " << Nl << " " << Nc << endl;
    for (i = 0; i < Nl; i++) 
    {
       fltarray LineG(Nc2);
       fltarray LineH(Nc2);
       float *PtrHigh = Data.buffer() + Data.nc()*i + Dep;
       Ptr_SB1D->transform(Nc, PtrHigh,  LineH.buffer(), LineG.buffer());
       for (j = 0; j < Nc2; j++) PtrHigh[j] = LineH(j);
       for (j = 0; j < Nc/2; j++)  PtrHigh[j+Nc2] = LineG(j);
   }
}

/****************************************************************************/

void LineCol::undec_transform_one_step_line (Ifloat & Data, Ifloat & Low, Ifloat & High, int Step)
{
    int Nl=Data.nl();
    int Nc=Data.nc();
    fltarray LineG(Nc);
    fltarray LineH(Nc);
    int j,i;
    for (i = 0; i < Nl; i++) 
    {
       float *PtrHigh = Data.buffer() + Data.nc()*i;
       Ptr_SB1D->transform(Nc, PtrHigh,  LineH.buffer(), LineG.buffer(), Step);
       for (j = 0; j < Nc; j++) Low(i,j) = LineH(j);
       for (j = 0; j < Nc; j++) High(i,j) = LineG(j);
   }
}

/****************************************************************************/

void LineCol::recons_one_step_line(Ifloat &Data, int Nl, int Nc, int Dep)
{
   int i,j;
   SubBand1D *SB1D = get_subband_method();
   int Nc2 = (Nc+1)/2;
   // cout << "Rec Line step " << Nl << " " << Nc << endl;
   for (i = 0; i < Nl; i++) 
   {
      fltarray LineH(Nc2);
      fltarray LineG(Nc2);
      float *PtrHigh =  Data.buffer() + Data.nc()*i + Dep;
      for (j = 0; j < Nc2; j++)  LineH(j) = Data(i,j+ Dep);
      for (j = 0; j < Nc/2; j++) LineG(j) = Data(i,j+Nc2+ Dep);
      SB1D->recons(Nc,  LineH.buffer(), LineG.buffer(), PtrHigh);
   }  
}

/****************************************************************************/

void LineCol::undec_recons_one_step_line (Ifloat & Low, Ifloat & High, Ifloat & Data,  int Step)
{
    int Nl=Data.nl();
    int Nc=Data.nc();
    fltarray Rec(Nc);
    int j,i;
    for (i = 0; i < Nl; i++) 
    {
       float *PtrHigh = High.buffer() + Nc*i;
       float *PtrLow = Low.buffer() + Nc*i;
       Ptr_SB1D->recons(Nc, PtrLow, PtrHigh, Rec.buffer(), Step);
       for (j = 0; j < Nc; j++) Data(i,j) = Rec(j);
    }
}

/****************************************************************************/

void LineCol::transform_one_step_col (Ifloat & Data, int Nl, int Nc, int Dep)
{
    int Nl2 = (Nl+1)/2;
    int j,i;
    // cout << "Col step " << Nl << " " << Nc << endl;
    for (j = 0; j < Nc; j++) 
    {
       fltarray Col(Nl);
       fltarray ColG(Nl2);
       fltarray ColH(Nl2);
       for (i = 0; i < Nl; i++) Col(i) =  Data(i+Dep,j);
       Ptr_SB1D->transform(Nl, Col.buffer(), ColH.buffer(), ColG.buffer());
       for (i = 0; i < Nl2; i++)   Data(i+Dep,j) = ColH(i);
       for (i = 0; i < Nl/2; i++)  Data(i+Nl2+Dep,j) = ColG(i);
    }   
}

/****************************************************************************/

void  LineCol::undec_transform_one_step_col (Ifloat & Data, Ifloat & Low, Ifloat & High,  int Step)
{
    int Nl=Data.nl();
    int Nc=Data.nc();
    fltarray Col(Nl);
    fltarray ColG(Nl);
    fltarray ColH(Nl);
    int j,i;
    for (j = 0; j < Nc; j++) 
    {
       for (i = 0; i < Nl; i++) Col(i) =  Data(i,j);
       Ptr_SB1D->transform(Nl, Col.buffer(), ColH.buffer(), ColG.buffer(), Step);
       for (i = 0; i < Nl; i++) Low(i,j) = ColH(i);
       for (i = 0; i < Nl; i++) High(i,j) = ColG(i);
    }   
}

/****************************************************************************/

void LineCol::recons_one_step_col(Ifloat &Data, int Nl, int Nc, int Dep)
{
   int i,j;
   SubBand1D *SB1D = get_subband_method();
   int Nl2 = (Nl+1)/2;
    // cout << "Rec Col step " << Nl << " " << Nc << endl;

   for (j = 0; j < Nc; j++) 
   {
      fltarray  Col(Nl);
      fltarray  ColH(Nl2);
      fltarray  ColG(Nl2);
      for (i = 0; i < Nl2; i++)  ColH(i) =  Data(i+Dep,j);
      for (i = 0; i < Nl/2; i++) ColG(i) =  Data(i+Nl2+Dep,j);
      SB1D->recons (Nl,  ColH.buffer(),  ColG.buffer(), Col.buffer());
      for (i = 0; i < Nl; i++)  Data(i+Dep,j) = Col(i);
   }
}

/*****************************************************************************/

void LineCol::undec_recons_one_step_col(Ifloat & Low, Ifloat & High, Ifloat & Data,  int Step)
{
    int Nl=Data.nl();
    int Nc=Data.nc();
    fltarray  ColH(Nl);
    fltarray  ColG(Nl);
    fltarray Rec(Nl);
    int j,i;
    for (j = 0; j < Nc; j++) 
    {
       for (i = 0; i < Nl; i++) ColH(i) =  Low(i,j);
       for (i = 0; i < Nl; i++) ColG(i) =  High(i,j);
       Ptr_SB1D->recons(Nl, ColH.buffer(), ColG.buffer(), Rec.buffer(), Step);
       for (i = 0; i < Nl; i++) Data(i,j) = Rec(i);
    }
}

/*****************************************************************************/
/*****************************************************************************/
void LineCol::alloc(SubBand1D &SB1D, Bool UseMirror)
{
   Ptr_SB1D = &SB1D; 
   Mirror=UseMirror;
}

void LineCol::directional_transform(Ifloat &Data, Ifloat &DataTrans, float AngleRot, 
     int NbrScaleLine, int NbrScaleCol, Bool RadianUnitAngle)
{
    Ifloat Buff;
    int Nl = Data.nl();
    int Nc = Data.nc();
    int Nlt = 2*Nl;
    int Nct = 2*Nc;
    if ((DataTrans.nl() != Nlt) || (DataTrans.nc() != Nct))
          DataTrans.resize(Nlt, Nct);
    Buff.alloc(Nlt, Nct,"buffer");
    im_extend (Data, Buff);
    Rot.RadianUnitAngle = RadianUnitAngle;
    Rot.im_move_rotate(Buff, DataTrans, AngleRot);
    transform(DataTrans, NbrScaleLine, NbrScaleCol);
}

/****************************************************************************/

void LineCol::directional_recons(Ifloat &DataTrans, Ifloat &Data, float AngleRot, 
     int NbrScaleLine, int NbrScaleCol, Bool RadianUnitAngle)
{
    Ifloat Buff;
    int Nlt = DataTrans.nl();
    int Nct = DataTrans.nc();
    int Nl = Nlt / 2;
    int Nc = Nct / 2;
    if ((Data.nl() != Nl) || (DataTrans.nc() != Nc)) Data.resize(Nl, Nc);
    Buff.alloc(Nlt, Nct,"buffer");
    recons(DataTrans, NbrScaleLine, NbrScaleCol);
    Rot.RadianUnitAngle = RadianUnitAngle;
    Rot.im_move_rotate(DataTrans, Buff, -AngleRot);
    im_extract(Buff, Data);
}

/****************************************************************************/
/****************************************************************************/

void DirectionalLineCol::transform_col(Ifloat & Data, int Nld, int Ncd, int NScale, int Dep)
{
    int Nl = Nld;
    int Nc = Ncd;
    int Nl2 = 2*Nl;
    int s,j,i;
    fltarray Col(Nl);
    fltarray ColG(Nl2);
    fltarray ColH(Nl2);
    // cout << "Col step " << Nl << " " << Nc << endl;
    Dep = 0;
    for (s = 0; s < NScale-1; s++)
    {
       int Nl2 = (Nl+1)/2;
       for (j = 0; j < Nc; j++) 
       {
           for (i = 0; i < Nl; i++) Col(i) =  Data(i,j+Dep);
           Ptr_SB1D->transform(Nl, Col.buffer(), ColH.buffer(), ColG.buffer());
           for (i = 0; i < Nl2; i++)   Data(i,j+Dep) = ColH(i);
           for (i = 0; i < Nl/2; i++)  Data(i+Nl2,j+Dep) = ColG(i);
       }   
       Nl = (Nl+1) / 2;
    }
}

/****************************************************************************/

void DirectionalLineCol::transform(Ifloat & Data, Ifloat & Trans, int NbrScaleLine, int NbrScaleCol)
{
    int Nl = Data.nl();
    int Nc = Data.nc(); 
    int Step = 1;
    int Nct = (NbrUndecimatedScaleCol < 1) ? Nc: Nc*(NbrUndecimatedScaleCol+1);
    int Nlt = (NbrUndecimatedScaleLine < 1) ? Nl: Nl*(NbrUndecimatedScaleLine+1);
    int s,j,i;
    int Dep=0;
    
    if (Ptr_SB1D == NULL)
    {
       cout << "ERROR: class DirectionalLineCol not correctly allocated ... " << endl;
       exit(-1);
    }
    if ((Trans.nl() != Nlt) || (Trans.nc() != Nct)) Trans.resize(Nlt,Nct);

    cout << "Trans size " << Nlt << " " << Nct << endl;
    fltarray LineHigh(Nc);
    fltarray LineG(Nc);
    fltarray LineH(Nc);
    float *PtrHigh=NULL;
    
    // 1D Wavelet Transform along lines.
    for (s = 0; s < NbrScaleCol-1; s++)
    {
       int Nc2 = (Nc+1)/2;
       Dep = (NbrUndecimatedScaleCol-s)*Nc;
       cout << "Scale " << s+1 << " Nl = " << Nl << " Nc = " << Nc << " Step = " << Step << endl;
       if (is_decimated_col_scale(s) == True) cout << "DECIMATED SCALE " << Nc2 << " " << Nc/2 << endl;
       else cout << "UNDECIMATED SCALE, Dep = " << Dep << " Step = " << Step <<  endl;
       
        Ptr_SB1D->DistPix = Step;
       for (i = 0; i < Nl; i++) 
       {
          if (s == 0) PtrHigh = Data.buffer() + Data.nc()*i;
	  else PtrHigh = Trans.buffer() + Trans.nc()*i;
	  
          if (is_decimated_col_scale(s) == True) 
	  {
	     // cout << "DECIMATED SCALE " << endl;
             Ptr_SB1D->transform(Nc, PtrHigh,  LineH.buffer(), LineG.buffer());
	     if (s == 0) PtrHigh = Trans.buffer() + Trans.nc()*i;
	     for (j = 0; j < Nc2; j++) PtrHigh[j] = LineH(j);
             for (j = 0; j < Nc/2; j++) PtrHigh[j+Nc2] = LineG(j);
	  }
	  else
	  {
	     // cout << "UNDECIMATED SCALE " << endl;
             Ptr_SB1D->transform(Nc, PtrHigh,  LineH.buffer(), LineG.buffer(), Step);
	     if (s == 0) PtrHigh = Trans.buffer() + Trans.nc()*i;
	     for (j = 0; j < Nc; j++) PtrHigh[j] = LineH(j);
             for (j = 0; j < Nc; j++)  PtrHigh[Dep+j] = LineG(j);
          }
       }      
       if (is_decimated_col_scale(s) == True) Nc = (Nc+1) / 2;
       else Step *= 2;
   }
   
   // 1D Wavelet Transform along colomns.
//    Nc = Data.nc();
//    for (s = 0; s < NbrScaleCol-1; s++)
//    {
//       int Nc2 = (Nc+1)/2;
//       int NScale =
//       int Dep = (is_decimated_col_scale(s) == True) ? (NbrUndecimatedScaleCol-s)*Nc: Nc2;
//       if (is_decimated_col_scale(s) == True) Nc = (Nc+1) / 2;
//       transform_col(Trans, Nl, Nc, NScale, Dep);
// 
//    }
}


/****************************************************************************/

void DirectionalLineCol::recons(Ifloat & Trans, Ifloat & Data, int NbrScaleLine, int NbrScaleCol)
{
    int Nl = Data.nl();
    int Nc = Data.nc(); 
    int Step = 1;
//    int Nct = (NbrUndecimatedScaleCol < 1) ? Nc: Nc*NbrUndecimatedScaleCol;
    int s,j,i;
    int Dep=0;
    
    if (Ptr_SB1D == NULL)
    {
       cout << "ERROR: class DirectionalLineCol not correctly allocated ... " << endl;
       exit(-1);
    }
    
    // if ((Trans.nl() != Nl) || (Trans.nc() != Nc)) Trans.resize(Nl,Nct);

    // cout << "Line step " << Nl << " " << Nc << endl;
    cout << "DATA " << " Nl = " << Nl << " Nc = " << Nc << endl;

    fltarray LineHigh(Nc);
    fltarray LineG(Nc);
    fltarray LineH(Nc);
    float *PtrHigh=NULL;
    int *TabNl = new int [NbrScaleLine];
    int *TabNc = new int [NbrScaleCol];
   
    TabNl[0] = Nl;
    TabNc[0] = Nc;
    for (s = 1; s < NbrScaleLine; s++)
    {  
       if (is_decimated_line_scale(s-1) == True) TabNl[s] = (TabNl[s-1]+1)/2;
       else TabNl[s] = TabNl[s-1];
    }
    for (s = 1; s < NbrScaleCol; s++)
    {  
       if (is_decimated_col_scale(s-1) == True) TabNc[s] = (TabNc[s-1]+1)/2;
       else TabNc[s] = TabNc[s-1];
    }
    
    for (s = 0; s < NbrUndecimatedScaleCol; s++)
                    if (is_decimated_col_scale(s) == False) Step *= 2;   


    for (s = NbrScaleCol-2; s >= 0; s--)
    {
       char SS[100];
       int Nc2 = (TabNc[s]+1)/2;
       Dep = (NbrUndecimatedScaleCol-s)*Nc;
       if (is_decimated_col_scale(s) == False) Step /= 2;  
       cout << "Scale " << s+1 << " Nl = " << Nl << " Nc = " << TabNc[s] << " Step = " << Step << endl;
       if (is_decimated_col_scale(s) == True) cout << "DECIMATED SCALE " << Nc2 << " " << TabNc[s]/2 << endl;
       else cout << "UNDECIMATED SCALE, Dep =" << Dep << " Step = " << Step << endl;
       Ptr_SB1D->DistPix = Step;
       for (i = 0; i < Nl; i++) 
       {
          if (is_decimated_col_scale(s) == True)
          {
 	    for (j = 0; j < Nc2; j++)  LineH(j) = Trans(i,j);
            for (j = 0; j < TabNc[s]/2; j++) LineG(j) = Trans(i,j+Nc2);
	    PtrHigh = Trans.buffer() + Trans.nc()*i;
            Ptr_SB1D->recons (TabNc[s], LineH.buffer(), LineG.buffer(), PtrHigh);
	  }
          else  
          {
 	     for (j = 0; j < Nc; j++)  LineH(j) = Trans(i,j);
             for (j = 0; j < Nc; j++) LineG(j) = Trans(i,j+Dep);
	     if (s != 0) PtrHigh = Trans.buffer() + Trans.nc()*i;
	     else PtrHigh = Data.buffer() + Data.nc()*i;
             Ptr_SB1D->recons (Nc, LineH.buffer(), LineG.buffer(), PtrHigh, Step);
  	  }
       }
       sprintf(SS,"xx_rec_%d.fits", s+1);
       io_write_ima_float(SS, Trans);
   }
}
