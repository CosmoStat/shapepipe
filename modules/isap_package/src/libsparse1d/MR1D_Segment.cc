



#include "Array.h"
#include "IM_Noise.h"
#include "MR1D_Obj.h"
#include "MR1D_Filter.h"
#include "NR.h"
#include "MR1D_Segment.h"
#include "IM_IO.h"
#include "IM1D_IO.h"

#define DEBUG_SEG 0

/***********************************************************************/

void MR1D_Seg::clean()
{
   if (Nbr_Plan > 0)
   {
      if (TabObjIndStruct != NULL)
      {
          for (int s=1; s <= NbrObj; s++)
             if (TabObjIndStruct[s] != NULL) delete[] TabObjIndStruct[s];
          delete[] TabObjIndStruct;
      }
      NbrStruct = 0;
      NbrObj = 0;
      if (TabStruct != NULL) delete[] TabStruct;
      if (TabMax != NULL) delete[] TabMax;
      if (TabScaleMax != NULL) delete[] TabScaleMax;
      if (TabPosMax != NULL) delete[] TabPosMax;
      if (TabNbrStructPerObj != NULL) delete[] TabNbrStructPerObj;
   }
   TabStruct = NULL;
   TabMax = NULL;
   TabScaleMax = NULL;
   TabPosMax = NULL;
   TabNbrStructPerObj = NULL;
   TabObjIndStruct = NULL;
}


/***********************************************************************/

void MR1D_Seg::init(MR_1D & MR_Data)
{
   if (Nbr_Plan > 0) this->free();
   FirstScale = 0;
   PosDetect = False;
   NegDetect = False;
   Verbose = False;
   Np = MR_Data.size_ima_np();
   Nbr_Plan = MR_Data.nbr_scale();
   LastScale = Nbr_Plan-1;
   NbrStruct = 0;
   NbrObj = 0;
   TabStruct = NULL;
   TabMax = NULL;
   TabScaleMax = NULL;
   TabPosMax = NULL;
   TabNbrStructPerObj = NULL;
   TabObjIndStruct = NULL;
   MaxRecIter = MAX_REC_ITER_OBJ1D;
   Segment.alloc(MR_Data.size_ima_np(), MR_Data.Type_Transform, 
                "Segmentation", MR_Data.nbr_scale(),
                 MR_Data.Scale_0, MR_Data.Nu_0,
                 MR_Data.Nbr_Voie);
}

/***********************************************************************/

MR1D_Seg::MR1D_Seg()
{
   FirstScale = 0;
   PosDetect = False;
   NegDetect = False;
   Verbose = False;
   Np = 0;
   Nbr_Plan = 0;
   LastScale = 0;
   NbrStruct = 0;
   NbrObj = 0;
   TabStruct = NULL;
   TabMax = NULL;
   TabScaleMax = NULL;
   TabPosMax = NULL;
   MaxRecIter = MAX_REC_ITER_OBJ1D;
   TabNbrStructPerObj = NULL;
   TabObjIndStruct = NULL;
}

/***********************************************************************/

MR1D_Seg::MR1D_Seg(MR_1D & MR_Data)
{
   this->init(MR_Data);
}

/***********************************************************************/

void MR1D_Seg::free()
{
   if (Nbr_Plan > 0)
   {
      if (TabObjIndStruct != NULL)
      {
          for (int s=1; s <= NbrObj; s++)
             if (TabObjIndStruct[s] != NULL) delete[] TabObjIndStruct[s];
          delete[] TabObjIndStruct;
      }
      
      FirstScale = 0;
      LastScale = 0;
      PosDetect = False;
      NegDetect = False;
      Verbose = False;
      Np = 0;
      Nbr_Plan = 0;
      NbrStruct = 0;
      NbrObj = 0;
      MaxRecIter = MAX_REC_ITER_OBJ1D;
      if (TabStruct != NULL) delete[] TabStruct;
      if (TabMax != NULL) delete[] TabMax;
      if (TabScaleMax != NULL) delete[] TabScaleMax;
      if (TabPosMax != NULL) delete[] TabPosMax;
      if (TabNbrStructPerObj != NULL) delete[] TabNbrStructPerObj;
      Segment.free();
   }
}

/***********************************************************************/

MR1D_Seg::~MR1D_Seg()
{
   this->free();
}

/***********************************************************************/

void MR1D_Seg::segment(MR_1D & MR_Data, MR1DNoiseModel & NoiseModel)
{
  int i,j,s,p;
  int NbrScale = Nbr_Plan;
  int Nstruc = 0;
  Bool Stoploop;

  if (Np != MR_Data.size_ima_np())
  {
   cerr << "Error: MR1D_Seg::segment  incompatible Np in MR_Data ... " << endl;
  }
  if (Nbr_Plan != MR_Data.nbr_scale())
  {
    cerr << "Error: MR1D_Seg::segment  incompatible Nbr_Plan in MR_Data ... " << endl;
  }
  if (TabMax != NULL) delete [] TabMax;
  if (TabScaleMax != NULL) delete [] TabScaleMax;
  if (TabPosMax != NULL) delete [] TabPosMax;
  if (TabStruct != NULL) delete [] TabStruct;

  for (s=0; s < NbrScale; s++)
  for (i=0; i < Np; i++)  Segment(s,i) = 0.;

  // Segmentation scale by scale
  // we take into account change of sign ...
  float sign;

#if DEBUG_SEG
cout << " Start loop " << endl;
#endif

  for (s=0; s < NbrScale-1; s++)
  {
     TabNbrStruct[s] = 0;
     p = 0;
     while (p < Np)
     {
        Stoploop = False;
        do 
        {
           if (NoiseModel(s,p) == False) p++;
           else Stoploop = True;
        } while ((p < Np) && (Stoploop == False));
        if (Stoploop == True)
        {
            Nstruc ++;
            TabNbrStruct[s]++;
            Stoploop = False;
	    sign = (MR_Data(s,p) >= 0) ? 1 : -1;
            do 
            {
#if DEBUG_SEG
                if ((NoiseModel(s,p) == True) && (sign*MR_Data(s,p)<=0))
                {
                    cout << "SIGN" << endl;
                   cout << "s = " << s << "  p = " << p << endl;
                   cout << "MR_Data =  " << MR_Data(s,p) << "  sign = " << sign << endl;
                   cout << "NoiseModel(s,p) =  " << NoiseModel.sigma(s,p)  << endl;
                }
#endif

                if ((NoiseModel(s,p) == True) && (sign*MR_Data(s,p)>0))
                {
                   Segment(s,p) = -1;
		   p++;
                } 
                else Stoploop = True;
                //if (NoiseModel(s,p) == True) p++;
             } while ((p < Np) && (Stoploop == False));
        }
     }
  }
#if DEBUG_SEG
 cout << " step1 number of detected structures = " <<  Nstruc << endl;
 for (s=0; s < NbrScale-1; s++) cout << " Scale " << s+1 << " Nb = " << TabNbrStruct[s] << endl;
#endif

   // 1 dim object table allocation
   NbrStruct = Nstruc;
   if (NbrStruct != 0) 
   {

   TabStruct = new obj1d_struct [Nstruc+1];
   int ind_obj = 0;
   float Val;


  // fill the object table from the segmented data
  for (s=0; s < NbrScale-1; s++)
  {
     p = 0;
     while (p < Np)
     {
        Stoploop = False;
        do 
        {
           if (NoiseModel(s,p) == False) p++;
           else Stoploop = True;
        } while ((p < Np) && (Stoploop == False));
        if (Stoploop == True)
        {
           Stoploop = False;
           ind_obj ++;
           TabStruct[ind_obj].Sizex = 1;
           TabStruct[ind_obj].Scale = s;
           TabStruct[ind_obj].NumObj = ind_obj;
           TabStruct[ind_obj].Depx = p;
           Val = ABS(MR_Data(s,p));
           TabStruct[ind_obj].ValMax = MR_Data(s,p);
           TabStruct[ind_obj].Maxx = p;
           TabStruct[ind_obj].Noise = NoiseModel.sigma(s,p);
           Segment(s,p) = ind_obj;
	   sign = (MR_Data(s,p) >= 0)  ? 1 : -1;
           p++;
           if (p < Np) 
             do 
             {
                if ((NoiseModel(s,p) == True) && (sign*MR_Data(s,p)>0))
                {
                   TabStruct[ind_obj].Sizex++;
                   Segment(s,p) = ind_obj;
                   if (ABS(MR_Data(s,p)) > Val) 
                   {
                      Val = ABS(MR_Data(s,p));
                      TabStruct[ind_obj].ValMax = MR_Data(s,p);
                      TabStruct[ind_obj].Maxx = p;
                      TabStruct[ind_obj].Noise = NoiseModel.sigma(s,p);
                   }
                   p++;
                }
                else Stoploop = True;
             } while ((p < Np) && (Stoploop == False));        
        }
     }
  }
  // interscale relation taking into account subobject relations
  float *TabColor = new float [Nstruc+1];
  fltarray TabRel(Nstruc+1,Nstruc+1);
  TabRel.init();

  int posmax,indrel,colrel;
  for (i = 1; i <= Nstruc; i++) TabColor[i] = i;
  for (i = 1; i <= Nstruc; i++)
  {
     posmax = TabStruct[i].Maxx;
     s = TabStruct[i].Scale;
     if (s < NbrScale-2)
     {
         if (Segment(s+1,posmax) > FLOAT_EPSILON)
         {
             for (j = 1; j < i; j++)
               if (TabColor[j] == TabColor[i]) 
                               TabColor[j] = Segment(s+1, posmax);
             TabColor[i] = Segment(s+1, posmax);

             colrel = (int)(Segment(s+1, posmax)+ 0.5);
             TabRel(colrel,0) += 1.;
             indrel = (int)(TabRel(colrel,0) + 0.5);
             TabRel(colrel,indrel) = i;
             if (Verbose == True)
                cout << "interscale relation... " << i << " to " << colrel << endl;
         }
     }
  }
#if DEBUG_SEG
  cout << "print interscale relation table " << endl;
  for (i = 1; i <= Nstruc; i++)
  {
     colrel = (int)(TabRel(i, 0)+ 0.5);
     cout << "Rel struct " << i << "n= " << colrel << "  :" << endl;
     for (s=1; s <= colrel; s++) cout << " S" << TabRel(i, s);
  }
#endif


  // count the number of objects
  unsigned long *TabIndex = new unsigned long [Nstruc+1];
  float *TabOrderColor = new float [Nstruc+1];
  float *TabNewColor = new float [Nstruc+1];
  for (i = 0; i <= Nstruc; i++) TabNewColor[i] = i;
  for (i = 0; i <= Nstruc; i++) TabOrderColor[i] = TabColor[i];
  unsigned long No = Nstruc;
  indexx(No, TabColor, TabIndex);
  sort(No, TabOrderColor);
#if DEBUG_SEG
cout << "delete the hole in TabColor  ..." << endl;
#endif
  i=1;
  int ncol = 1;
  int Col_ind;
  int Col_obj;
  do
  {
      // delete the "hole" in TabColor by creation of a new TabColor of
      // name TabNewColor
      Col_ind = TabIndex[i];
      Col_obj = (int) TabColor[ Col_ind ];
      for (j = 1; j <= Nstruc; j++) 
                      if (TabColor[j] == Col_obj) TabNewColor[j] = ncol;

      // go to the next color
      do
      {
         i++;
         if (i > Nstruc) Stoploop = True;
         else
         {
             if (TabOrderColor[i] != TabOrderColor[i-1])
             {
                Stoploop = True;
                ncol++;
             }
             else Stoploop = False;
         }
      } while ((i <= Nstruc) && (Stoploop == False));
  } while (i <= Nstruc);
#if DEBUG_SEG
cout << "define if a structure is a starting point of a sub-object ..." << endl;
#endif

  // define if a structure is a starting point of a sub-object
  float Max_Before;
  for (i = 1; i <= Nstruc; i++) TabStruct[i].SubObj = True;
  for (i = 1; i <= Nstruc; i++)
  {
     int nbrel;

     posmax = TabStruct[i].Maxx;
     s = TabStruct[i].Scale;
     Val = TabStruct[i].ValMax;
       
     nbrel = (int) (TabRel(i,0) + 0.5);
     for (s = 1; s <= nbrel; s++)
     {
        indrel = (int) (TabRel(i,s) + 0.5);
        Max_Before = TabStruct[indrel].ValMax;
        if (ABS(Val) < ABS(Max_Before)) TabStruct[i].SubObj = False;
        else TabStruct[indrel].SubObj = False;
     }
  }
#if DEBUG_SEG
cout << "INFO obj ..." << endl;
#endif

  // search for each object its maximum value, its position, and its scale
  TabMax = new float [ncol+1];
  TabScaleMax = new int [ncol+1];
  TabPosMax = new int [ncol+1];
  NbrObj = ncol;
  for (i = 1; i <= ncol; i++) TabScaleMax[i]=-1;
  for (i = 1; i <= ncol; i++) TabMax[i]=0.;

  for (i = 1; i <= Nstruc; i++)
  {
      Col_obj = (int) TabNewColor[i];
      TabStruct[i].NumConnectObj = Col_obj;
/*
      if (TabStruct[i].SubObj == True)
      {
         if (TabStruct[i].Scale > TabScaleMax[Col_obj])
         {
             TabMax[Col_obj] = TabStruct[i].ValMax;
             TabPosMax[Col_obj] = TabStruct[i].Maxx;
             TabScaleMax[Col_obj] = TabStruct[i].Scale;
         }
      }
*/
      if (ABS(TabStruct[i].ValMax) > ABS(TabMax[Col_obj]))
      {
          TabMax[Col_obj] = TabStruct[i].ValMax;
          TabPosMax[Col_obj] = TabStruct[i].Maxx;
          TabScaleMax[Col_obj] = TabStruct[i].Scale;
      }
  }
#if DEBUG_SEG
cout << "fill the stucture array ..." << endl;
#endif

  // fill the stucture array 
  for (i = 1; i <= Nstruc; i++)
  {
     Col_obj = TabStruct[i].NumConnectObj;
     TabStruct[i].NumScale = TabScaleMax[Col_obj];
  }
#if DEBUG_SEG
cout << "last loop ..." << endl;
#endif
  for (i = 1; i <= Nstruc; i++)
  {
     int k, Dep, SizeS, ColS;
     s = TabStruct[i].Scale;
     Dep = TabStruct[i].Depx;
     SizeS = TabStruct[i].Sizex;
     ColS = TabStruct[i].NumConnectObj;
     for (k = Dep; k < Dep + SizeS; k++) Segment(s,k) = ColS;
  }
#if DEBUG_SEG
cout << "end last loop ..." << endl;
#endif
 delete[] TabColor;
 delete[] TabIndex; 
 delete[] TabOrderColor;
 delete[] TabNewColor;

#if DEBUG_SEG
cout << "set_obj_info ..." << endl;
#endif

   set_obj_info();

#if DEBUG_SEG
cout << "END SEGMENTATION..." << endl;
#endif

}
}

/***********************************************************************/

void MR1D_Seg::write_data_seg(char *FileName)
{
  extern type_1d_format IO_1D_Format;
  if (IO_1D_Format == F1D_FITS)
       fits_write_fltarr(FileName, Segment.image());
  else io_write2d_ascii(FileName, Segment.image()); 
  //   fits_write_fltarr(FileName, Segment.image());
}

/***********************************************************************/

void MR1D_Seg::print_info()
{
  int i,s;

  cout << "Number of objects = " << NbrObj << endl;
  cout << "number of structures = " <<  NbrStruct << endl;
  for (s=0; s < Nbr_Plan-1; s++) cout << "  Scale " << s+1 << " Nb = " << TabNbrStruct[s] << endl;
  cout << endl;

  for (i = 1; i <= NbrStruct; i ++) 
  {
      cout << "Structure number " << TabStruct[i].NumObj;
      cout << "  Object number " << TabStruct[i].NumConnectObj << endl;
      if (TabStruct[i].SubObj) cout << "  Sub-Object: True " << endl;
      else cout << "  Sub-Object: False " << endl;
      cout << "  Scale " << TabStruct[i].Scale;
      cout << "  NumScale " << TabStruct[i].NumScale;
      cout << "  Depx " << TabStruct[i].Depx;
      cout << "  Maxx " << TabStruct[i].Maxx<< endl;
      cout << "  ValMax " << TabStruct[i].ValMax;
      cout << "  Sizex " << TabStruct[i].Sizex; 
      cout << "  Noise " << TabStruct[i].Noise<< endl<< endl;
  }
}

/***********************************************************************/

void MR1D_Seg::set_obj_info()
{
  int i,s,IndS;

  TabNbrStructPerObj = new int [NbrObj+1];
  TabObjIndStruct = new int * [NbrObj+1];

  for (s = 0; s <= NbrObj; s ++) TabNbrStructPerObj[s] = 0;
  for (s = 1; s <= NbrObj; s ++) 
  for (i = 1; i <= NbrStruct; i ++) 
        if (TabStruct[i].NumConnectObj == s) TabNbrStructPerObj[s] ++;

  for (s = 1; s <= NbrObj; s ++)
     TabObjIndStruct[s] = new int [TabNbrStructPerObj[s]+1];

  for (s = 1; s <= NbrObj; s ++) 
  {
     IndS = 1;
     for (i = 1; i <= NbrStruct; i ++) 
     {
        if (TabStruct[i].NumConnectObj == s) TabObjIndStruct[s][IndS++]=i;
     }
  }
}

/***********************************************************************/

int MR1D_Seg::nbr_struct_in_obj(int NumObj)
{
   if ((NumObj < 1) || (NumObj > NbrObj))
   {
      cerr << "Bad object number ... " << endl;
   }
   return TabNbrStructPerObj[NumObj];
}

/***********************************************************************/

int MR1D_Seg::scale_of_obj(int NumObj)
{
   if ((NumObj < 1) || (NumObj > NbrObj))
   {
      cerr << "Bad object number ... " << endl;
   }
   return TabScaleMax[NumObj];
}

/***********************************************************************/

float MR1D_Seg::max_of_obj(int NumObj)
{
   if ((NumObj < 1) || (NumObj > NbrObj))
   {
      cerr << "Bad object number ... " << endl;
   }
   return TabMax[NumObj];
}

/***********************************************************************/

int MR1D_Seg::posmax_of_obj(int NumObj)
{
   if ((NumObj < 1) || (NumObj > NbrObj))
   {
      cerr << "Bad object number ... " << endl;
   }
   return TabPosMax[NumObj];
}

/***********************************************************************/

obj1d_struct & MR1D_Seg::struct_in_obj(int NumObj, int NumStruct)
{
   if ((NumObj < 1) || (NumObj > NbrObj))
   {
      cerr << "Bad object number ... " << endl;
      exit(-1);
   }
   if ((NumStruct < 1) || (NumStruct > NbrStruct))
   {
      cerr << "Bad structure number ... " << endl;
      exit(-1);
   }
   int Ind = TabObjIndStruct[NumObj][NumStruct];
   return TabStruct[Ind];
}

/***********************************************************************/

int MR1D_Seg::firstpos_of_obj(int NumObj)
{
   if ((NumObj < 1) || (NumObj > NbrObj))
   {
      cerr << "Bad object number ... " << endl;
      exit(-1);
   }
   int Ind,s;
   int Start = Np;

   for (s=1; s <= TabNbrStructPerObj[NumObj]; s++)
   {
      Ind = TabObjIndStruct[NumObj][s];
      if (TabStruct[Ind].Depx < Start) Start = TabStruct[Ind].Depx;
   }
   return Start;
}
/***********************************************************************/

int MR1D_Seg::lastpos_of_obj(int NumObj)
{
   if ((NumObj < 1) || (NumObj > NbrObj))
   {
      cerr << "Bad object number ... " << endl;
      exit(-1);
   }
   int s,Ind;
   int DimS, Last = 0;

   for (s=1; s <= TabNbrStructPerObj[NumObj]; s++)
   {
      Ind = TabObjIndStruct[NumObj][s];
      DimS = TabStruct[Ind].Depx+TabStruct[Ind].Sizex-1;
      if (DimS > Last) Last = DimS;
   }
   return Last;
}


/***********************************************************************/

void MR1D_Seg::print_obj_info()
{  int i,s;

  cout << "Number of objects = " << NbrObj << endl;
  cout << "number of structures = " <<  NbrStruct << endl;
  for (s=0; s < Nbr_Plan-1; s++) cout << "  Scale " << s+1 << " Nb = " << TabNbrStruct[s] << endl;
  cout << endl;

  for (s = 1; s <= NbrObj; s ++) 
  {
     cout << "Object number " << s << endl;
     cout << "  Nbr structure: " << TabNbrStructPerObj[s] << endl;
     cout << "  From " << firstpos_of_obj(s) << " to " << lastpos_of_obj(s);
     cout << "  Scale = " << scale_of_obj(s);
     cout << "  ValMax = " << max_of_obj(s) << endl;

     for (i = 1; i <= TabNbrStructPerObj[s]; i ++) 
     {
           cout << "Structure number " << struct_in_obj(s,i).NumObj;
           if (struct_in_obj(s,i).SubObj) cout << "  Sub-Object: True " << endl;
           else cout << "  Sub-Object: False " << endl;
           cout << "  Scale " << struct_in_obj(s,i).Scale;
           cout << "  NumScale " << struct_in_obj(s,i).NumScale;
           cout << "  Depx " << struct_in_obj(s,i).Depx;
           cout << "  Maxx " << struct_in_obj(s,i).Maxx<< endl;
           cout << "  ValMax " << struct_in_obj(s,i).ValMax;
           cout << "  Sizex " << struct_in_obj(s,i).Sizex;
           cout << "  Noise " << struct_in_obj(s,i).Noise<< endl<< endl;
     }
     cout << endl<< endl;
  }
}

/***********************************************************************/
/*
void MR1D_Seg::print_obj_info()
{
  int i,s;

  cout << "Number of objects = " << NbrObj << endl;
  cout << "number of structures = " <<  NbrStruct << endl;
  for (s=0; s < Nbr_Plan-1; s++) cout << "  Scale " << s+1 << " Nb = " << TabNbrStruct[s] << endl;
  cout << endl;

  for (s = 1; s <= NbrObj; s ++) 
  {
     cout << "Object number " << s << endl;
     for (i = 1; i <= NbrStruct; i ++) 
     {
        if (TabStruct[i].NumConnectObj == s)
        {
           cout << "Structure number " << TabStruct[i].NumObj;
           cout << "  Object number " << TabStruct[i].NumConnectObj << endl;
           if (TabStruct[i].SubObj) cout << "  Sub-Object: True " << endl;
           else cout << "  Sub-Object: False " << endl;
           cout << "  Scale " << TabStruct[i].Scale;
           cout << "  NumScale " << TabStruct[i].NumScale;
           cout << "  Depx " << TabStruct[i].Depx;
           cout << "  Maxx " << TabStruct[i].Maxx<< endl;
           cout << "  ValMax " << TabStruct[i].ValMax;
           cout << "  Sizex " << TabStruct[i].Sizex;
           cout << "  Noise " << TabStruct[i].Noise<< endl<< endl;
        }
     }
     cout << endl<< endl;
  }
}
*/
/***********************************************************************/

void MR1D_Seg::recons_adjoint_iter(MR_1D & MR_Data, int NumObj, fltarray & SignalRec, Bool Pos, Bool Neg)
{
   int i,Iter,s;

   MR_1D MR_Temp(Np,MR_Data.Type_Transform,"MR_Temp",Nbr_Plan,MR_Data.Interp,
                 MR_Data.Scale_0, MR_Data.Nu_0, MR_Data.Nbr_Voie);
   fltarray Resi(Np);
   MR_Temp.Border = MR_Data.Border;
   SignalRec.reform(Np);
   SignalRec.init();
    
   for (Iter = 0; Iter < MaxRecIter; Iter++)
   {
      // cout << "Obj " << NumObj << "  Iter " << Iter+1 << " sigma  = " << sigma(SignalRec) << " Max = " << max(SignalRec) << " Min = " << min(SignalRec) << endl;

      // calculate the residual wavelet coefficients 
      MR_Temp.transform(SignalRec);
      for (s=0; s < Nbr_Plan-1; s++)
      for (i=0; i < MR_Temp.size_scale_np(s); i++)
      {
         if (Segment(s,i) != NumObj) MR_Temp(s,i) = 0.;
         else MR_Temp(s,i) = MR_Data(s,i) - MR_Temp(s,i);
      }
      //for (i=0; i < MR_Temp.size_scale_np(Nbr_Plan-1); i++)
      //                                         MR_Temp(Nbr_Plan-1,i) = 0;
      // MR_Temp.recons(Resi);

      MR_Temp.rec_adjoint(Resi, False);
      SignalRec += Resi;
      if (Pos == True)
        for (i=0; i < Np; i++) if (SignalRec(i) < 0) SignalRec(i) = 0.;
      if (Neg == True)
        for (i=0; i < Np; i++) if (SignalRec(i) > 0) SignalRec(i) = 0.;
   }
}


/***********************************************************************/

void MR1D_Seg::recons_direct(MR_1D & MR_Data, int NumObj, fltarray & SignalRec, Bool Pos, Bool Neg, Bool UseAdjointRec)
{
   int i,s;

   MR_1D MR_Temp(Np,MR_Data.Type_Transform,"MR_Temp",Nbr_Plan,MR_Data.Interp,
                 MR_Data.Scale_0, MR_Data.Nu_0, MR_Data.Nbr_Voie);
   MR_Temp.Border = MR_Data.Border;
   SignalRec.alloc(Np);

   for (s=0; s < Nbr_Plan-1; s++)
   for (i=0; i < MR_Temp.size_scale_np(s); i++)
   {
      if (Segment(s,i) != NumObj) MR_Temp(s,i) = 0.;
      else MR_Temp(s,i) = MR_Data(s,i);
   }
   if (UseAdjointRec == False) MR_Temp.recons(SignalRec);
   else MR_Temp.rec_adjoint(SignalRec, False);

   if (Pos == True)
        for (i=0; i < Np; i++) if (SignalRec(i) < 0) SignalRec(i) = 0.;
   if (Neg == True)
        for (i=0; i < Np; i++) if (SignalRec(i) > 0) SignalRec(i) = 0.;
}

/***********************************************************************/

void MR1D_Seg::recons(MR_1D & MR_Data, int NumObj, fltarray & SignalRec, Bool Pos, Bool Neg)
{

   if ((NumObj < 1) || (NumObj > NbrObj))
   {
      cerr << "Error: bad object number in MR1D_Seg::recons ... " << endl;
      exit(-1);
   }

   switch (MR_Data.Type_Transform)
   {
      case TO1_PAVE_LINEAR: 
      case TO1_PAVE_B1SPLINE: 
      case TO1_PAVE_B3SPLINE: 
      case TO1_PAVE_MORLET:
      case TO1_PAVE_MEX:
      case TO1_PAVE_FRENCH:
      case TO1_PYR_LINEAR:
      case TO1_PYR_B3SPLINE: 
          recons_adjoint_iter(MR_Data, NumObj, SignalRec, Pos, Neg);
          break;
      case TM1_PAVE_MEDIAN:
      case TM1_PYR_MEDIAN:
          recons_direct(MR_Data, NumObj, SignalRec, Pos, Neg);
          break;
      default:
          cerr << "Unknown transform ..." << endl;
          break;
   }   
}

/***********************************************************************/

void MR1D_Seg::recons(MR_1D & MR_Data, fltarray & Tab_SignalRec)
{
    int ind, obj, i, NRec=0;
    fltarray SignalRec;
    int *TabIsRec = new int [NbrObj+1];

    for (obj=1; obj <= NbrObj; obj++)
    {
       if ( (TabScaleMax[obj] >= FirstScale) &&  (TabScaleMax[obj] <= LastScale) &&
            ((PosDetect == False) || (TabMax[obj] > 0)) &&
            ((NegDetect == False) || (TabMax[obj] < 0)))
       {
          TabIsRec[obj] = NRec;
          NRec++;
       }
       else TabIsRec[obj] = -1;
    }

    Tab_SignalRec.alloc(Np, NRec);

    for (obj=1; obj <= NbrObj; obj++)
    {
       if (TabIsRec[obj] >= 0)
       {
           if (Verbose == True)
           { 
            cout << "Reconstruct Obj " << obj << " Scale = " << TabScaleMax[obj] << " Max = " << TabMax[obj] << endl;
            cout << "                   From " << firstpos_of_obj(obj) << " to " << lastpos_of_obj(obj) << endl;
           }
           ind = TabIsRec[obj];
           if (TabMax[obj] >= 0)
                recons(MR_Data, obj, SignalRec, True, False);
           else recons(MR_Data, obj, SignalRec, False, True);
           // recons(MR_Data, obj, SignalRec, PosDetect, NegDetect);
           for (i=0;i < Np; i++) Tab_SignalRec(i,ind) = SignalRec(i);
       }
    }
    delete [] TabIsRec;
}
/***********************************************************************/

void MR1D_Seg::recons(MR_1D & MR_Data, fltarray & Tab_SignalRec, fltarray & Tab_Add)
{
    int i,obj;

    this->recons(MR_Data, Tab_SignalRec);
    int No = Tab_SignalRec.axis(2);
    Tab_Add.alloc(Np);
    for (obj=0; obj < No; obj++) 
    for (i=0; i < Np; i++) Tab_Add(i) += Tab_SignalRec(i,obj);
}


/***********************************************************************/

void MR1D_Seg::recons(MR_1D & MR_Data, fltarray & Tab_SignalRec, fltarray & Tab_Add, fltarray & Data, MR1DNoiseModel &NoiseModel)
{
    int i,obj;
    float Val;

    this->recons(MR_Data, Tab_SignalRec);
    int No = Tab_SignalRec.axis(2);
    Tab_Add.alloc(Np);
    for (obj=0; obj < No; obj++) 
    for (i=0; i < Np; i++) 
    {
       switch (NoiseModel.which_noise())
       {
          case NOISE_GAUSSIAN: 
          case NOISE_NON_UNI_ADD:
          case NOISE_UNDEFINED:
          case NOISE_UNI_UNDEFINED:
               // just do nothing
           break;
          case NOISE_POISSON:
          case NOISE_GAUSS_POISSON:
           Val = NoiseModel.val_transform(Data(i));
           Tab_SignalRec(i,obj) = Data(i) -
                NoiseModel.val_invtransform(Val - Tab_SignalRec(i,obj));
           break;
        case NOISE_MULTI: 
        case NOISE_NON_UNI_MULT:
        default:
           break;
       }
       Tab_Add(i) += Tab_SignalRec(i,obj);
    }
}

/***********************************************************************/
