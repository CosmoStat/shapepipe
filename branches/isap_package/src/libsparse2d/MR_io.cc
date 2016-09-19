/*******************************************************************************
**
**    UNIT
**
**    Version: 3.3
**
**    Author: Jean-Luc Starck
**
**    Date:  96/06/13 
**    
**    Modification history:
**               7-JAN-1997 R Gastaud add fits
**
**    File:  MR_io.cc
**
*******************************************************************************
**
**    DESCRIPTION  Input Output procedures
**    -----------  
**                 
**
**
******************************************************************************
**
** void MultiResol::read (char *Name)
**
** read a multiresolution transform from a file
**
******************************************************************************
**
** void MultiResol::write (char *Name)
**
** write a multiresolution transform in a file
**
******************************************************************************/

#include "IM_Obj.h"
#include "SB_Filter.h"
#include "MR_Obj.h"
#include "IM_IO.h"

static void PrintError( int status);
static void mr_io_name (char *Name, char *filename);

#define DEBUG_IO 0

/****************************************************************************/

static void mr_io_name (char *File_Name_In, char *File_Name_Out)
{
    int L;

    strcpy (File_Name_Out, File_Name_In);

    L = strlen (File_Name_In);
    if ((L < 3) || (File_Name_In[L-1] != 'r')
                || (File_Name_In[L-2] != 'm')
                || (File_Name_In[L-3] != '.'))
    {
        strcat (File_Name_Out, ".mr");
    }
}

/****************************************************************************/

int MultiResol::mr_io_fill_header( fitsfile *fptr, int Nl, int Nc, int Nbr_Plan, 
         type_transform Type_Transform,set_transform Set_Transform,
         type_format FormatInputImag, char *Name_MR, type_border Border,
         int MedianWindowSize, float Fc, int Nbr_Iter, 
	 Bool ExactPyrRec, type_lift LiftingTrans, 
	 type_sb_filter SBFilter, sb_type_norm TypeNorm, int NbrUndec, type_undec_filter U_Filter)
{
  int status = 0; // this means OK for cfitsio !!!
    /*****************************************************/
     /* write optional keyword to the header */
    /*****************************************************/


if ( ffpkyj(fptr, (char*)"Nl", (long)Nl,(char*)"Number of lines of the input image",&status))
     PrintError( status );  
if ( ffpkyj(fptr,(char*)"Nc",(long)Nc,(char*)"Number of columns of the input image",&status))
     PrintError( status );  
 if ( ffpkyj(fptr, (char*)"Nbr_Plan", (long)Nbr_Plan, (char*)"Number of scales", &status))
     PrintError( status );  
 if (ffpkyj(fptr, (char*)"Set_Tran", (long)Set_Transform,
                      (char*)StringSetTransform(Set_Transform), &status))
      PrintError( status );  
 if ( ffpkyj(fptr, (char*)"Type_Tra", (long)Type_Transform, 
                         (char*)StringTransform(Type_Transform), &status))
     PrintError( status );  
     
 if (Type_Transform  == TC_FCT)
 {
     if ( ffpkyj(fptr,(char*)"NbrDir",(long) FCT_Nbr_Angle,(char*)"NbrDir",&status))
     PrintError( status ); 
     if ( ffpkyj(fptr,(char*)"Extend",(long) ((FCT_ExtendWT == True) ? 1: 0),(char*)"Extend WT",&status))
     PrintError( status );
     Bool Iso= FCT_IsotropWT;
     if ( ffpkyj(fptr,(char*)"Isotrop",(long) ((Iso == True) ? 1: 0),(char*)"Isotropic WT",&status))
     PrintError( status );
     if ( ffpkyj(fptr,(char*)"NbrScale",(long) nbr_scale(),(char*)"Number of scales",&status))
     PrintError( status );    
     if ( ffpkyj(fptr,(char*)"Real",(long) ((FCT_RealCur == True) ? 1: 0),(char*)"Real Curv",&status))
     PrintError( status );       
 }
 
 if ((Type_Transform  == TO_MALLAT) || (Type_Transform  == TO_UNDECIMATED_MALLAT) || (Type_Transform == TO_LC))
 {
  if ( ffpkyj(fptr, (char*)"SBFilter", (long) SBFilter, (char*)"Type of filters", &status))
        PrintError( status );  
     if ( ffpkyj(fptr, (char*)"NORM", (long)  TypeNorm, (char*)"normalization", &status))
        PrintError( status );  
  if ( ffpkyj(fptr, (char*)"UNDEC", (long) NbrUndec , (char*)"nbr of undec. scales", &status))
        PrintError( status ); 
  }	
  
  if ((Type_Transform  == TO_DIV_1) || (Type_Transform  == TO_DIV_2))
  {
    if ( ffpkyj(fptr, (char*)"NORM", (long)  TypeNorm, (char*)"normalization", &status))
        PrintError( status );  
    if ( ffpkyj(fptr, (char*)"UNDEC", (long) NbrUndec , (char*)"nbr of undec. scales", &status))
        PrintError( status ); 
  }	

  if (Type_Transform  == TO_LIFTING)
  if ( ffpkyj(fptr, (char*)"LiftingT", (long)  LiftingTrans, 
                           (char*)StringLSTransform(LiftingTrans), &status))
        PrintError( status );
     
 if ( ffpkyj(fptr, (char*)"FormatIn",(long)FormatInputImag,(char*)"I dunno", &status))
     PrintError( status );  

 if ( ffpkys(fptr, (char*)"Name_MR",Name_MR ,(char*)"MR OBJ", &status))
     PrintError( status );
 if ( ffpkyj(fptr, (char*)"Border",(long)Border,(char*)"BORDER Type", &status))
     PrintError( status );  

 if ( ffpkyj(fptr, (char*)"MedianWi",(long)MedianWindowSize,
        (char*)"size of the running window for median", &status))
     PrintError( status );  

 if ( ffpkye(fptr, (char*)"Fc", Fc, 6, (char*)"cut-off frequency", &status))
     PrintError( status );  
 if ( ffpkyj(fptr, (char*)"Nbr_Iter",(long)Nbr_Iter, (char*)"Number of iterations", &status))
     PrintError( status );  
 if ( ffpkyl(fptr, (char*)"ExactPyr", ExactPyrRec, 
                     (char*)"for Pyramidal reconstruction", &status))
     PrintError( status );

 if (SBFilter == F_USER)
 {
   if (UserFilterFileName != NULL)
   {
     if ( ffpkys(fptr, (char*)"FilBank", UserFilterFileName ,(char*)"Filter", &status))
           PrintError( status );
   }
   else if ( ffpkys(fptr, (char*)"FilBank", (char*)"-" ,(char*)"Filter", &status))
          PrintError( status );
 }
 if (Type_Transform  == TO_UNDECIMATED_NON_ORTHO)
 {
    if ( ffpkyj(fptr, (char*)"UFilBank",(long) U_Filter, (char*)"Undec Filter", &status))
     PrintError( status );  
 }
 return(status);
} 


/****************************************************************************/

void MultiResol::write (char *Name)
/* new version with fits */
{
 char filename[256];
 Ifloat Ima;
 fitsfile *fptr;    
 int status;
 int i,j,s;
 float *Ptr;
 int simple;
 int bitpix;
 long naxis=0;
 long naxes[3];
 long nelements;
 long group = 1;  /* group to write in the fits file, 1= first group */
 long firstpixel = 1;    /* first pixel to write (begin with 1) */
 Ifloat Aux;
 long fpixels[3];
 long lpixels[3];
 int Nelem=0;
 
/* we keep mr as extension even if its fits ! */
 mr_io_name (Name, filename);

#if DEBUG_IO  
    cout << "Write on " << filename << endl;
#endif

 FILE *FEXIST = fopen(filename, "rb");
 if (FEXIST)
 {
    fclose(FEXIST);
    remove(filename);               /* Delete old file if it already exists */
 }

 status = 0;         /* initialize status before calling fitsio routines */

    /* open the file */
 if ( ffinit(&fptr, filename, &status) )     /* create the new FITS file */
     PrintError( status );           /* call PrintError if error occurs */
                                                                              
/* write  the header */
 simple   = True;
 bitpix   =  -32;   /* 32-bit real pixel values      */
 long pcount   =   0;  /* no group parameters */
 long gcount   =   1;  /* only a single image/group */
 int  extend   =   False;
 switch (Set_Transform)
     {
     case TRANSF_UNDECIMATED_MALLAT:
     case TRANSF_DIADIC_MALLAT:
     case TRANSF_PAVE:
     case TRANSF_SEMIPYR:
         if (Type_Transform  == TC_FCT)
         {
            naxis = 1;   
	    Nelem=1+nbr_scale();
            for (int b=0; b < NbrBand; b++)  Nelem += 2 + size_band_nl(b)*size_band_nc(b);
            naxes[0] = Nelem;
	 }
	 else
	 {
            naxis = 3;
            naxes[0] = Nc;
            naxes[1] = Nl;
            naxes[2] = NbrBand;
	 }
         break;
     case TRANSF_PYR:
         naxis = 2;
         naxes[0] = 2*Nc;
         naxes[1] = 2*Nl;
         break;
     case TRANSF_MALLAT:
     case TRANSF_FEAUVEAU:
         naxis = 2;
         naxes[0] = Nc;
         naxes[1] = Nl ;
         break;
     default:
         fprintf (stderr, "Error in mr_io_write: bad Set_Transform");
         exit(-1);
         break; 
     }




 // write first header part (parameters)
 if ( ffphpr(fptr,simple,bitpix,naxis,naxes,pcount,gcount,extend,&status) )
     PrintError( status );          /* call PrintError if error occurs */

  // cout << " version 1.26" << Nl << endl;

  /* write optional keyword to the header */

    
 // write the header of the multiresolution file
 status = mr_io_fill_header( fptr, Nl, Nc, Nbr_Plan, 
         Type_Transform, Set_Transform, FormatInputImag, Name_MR, 
         Border,MedianWindowSize, Fc, Nbr_Iter, ExactPyrRec, 
	 LiftingTrans, SBFilter, TypeNorm, NbrUndecimatedScale, U_Filter);

 // save the data
 long L1 = 1;
 switch (Set_Transform)
 {
      case TRANSF_DIADIC_MALLAT:
      case TRANSF_PAVE:
         fpixels[0] = L1;
         fpixels[1] = L1;
         fpixels[2] = L1;
         lpixels[0] = Nc;
         lpixels[1] = Nl;
         lpixels[2] = L1;
         for (s = 0; s < NbrBand; s++)
             {
              Ptr = band(s).buffer();
              if ( ffpsse(fptr, 1, naxis, naxes, fpixels, lpixels, Ptr, &status))
                 PrintError( status );  
             lpixels[2] ++; 
             fpixels[2] ++;
             }
         break; 
     case TRANSF_UNDECIMATED_MALLAT:   
     case TRANSF_SEMIPYR:
         if (Type_Transform  == TC_FCT)
         {
            fltarray Data(Nelem);
            int ind=0;
            Data(ind++) = nbr_scale();
            for  (int s=0; s < nbr_scale(); s++) Data(ind++) = TabNbrBandPerResol(s);

            // Ifloat Band;
            for (int b=0; b < nbr_band(); b++)
            {
               Data(ind++) = size_band_nl(b);
               Data(ind++) = size_band_nc(b);
               for (int i=0;  i < size_band_nl(b); i++)
               for (int j=0;  j < size_band_nc(b); j++) Data(ind++) = (*this)(b,i,j);
       
            }
            naxes[0] = ind;
            fpixels[0]=1; 
	    nelements = ind;
            if ( ffppre(fptr, group, firstpixel, nelements, Data.buffer(), &status) )
              PrintError( status );  
 	 }
	 else
	 {
         fpixels[0] = 1;
         fpixels[1] = 1;
         fpixels[2] = 1;
         lpixels[0] = Nc;
         lpixels[1] = Nl;
         lpixels[2] = 1;
	 Aux.alloc(Nl,Nc,"aux");
         for (s = 0; s < NbrBand; s++)
         {
            // cout << "Band " << s+1 << " " << size_band_nl(s) << " " << size_band_nc(s) << endl;
	    Aux.init();
	    for (i =0; i < size_band_nl(s); i++)
	    for (j =0; j < size_band_nc(s); j++) Aux(i,j) = (*this)(s,i,j);
            Ptr = Aux.buffer();
            if ( ffpsse(fptr, 1, naxis, naxes, fpixels, lpixels, Ptr, &status))
                 PrintError( status );  
            lpixels[2] ++; 
            fpixels[2] ++;
         }
	 }
         break; 
     case TRANSF_PYR:
         nelements = naxes[0] * naxes[1];
         
         // mise a zero 
         /* define the null value (must do this before writing any data) */
         // halas, it does not work with float and double
         if (ffpkyj(fptr, "BLANK", 0, "value to use for undefined pixels",
                     &status))
             PrintError( status );

         //   ffppru fills the whole file with nan, which IDL does not like
        // if(ffppru(fptr, group, firstpixel ,nelements , &status))
        //     PrintError(status);
         Ptr = new float [naxes[1]];
         for ( s=0; s < naxes[1]; s++) Ptr[s]=0;
         fpixels[0]=1;
         for ( s=0; s < naxes[0]; s++)
            {
             if(ffppre(fptr,1, fpixels[0], naxes[1], Ptr, &status))   
                PrintError(status);
             fpixels[0] += naxes[1];        
             }
          delete [] Ptr;

         fpixels[0] = 1;
         fpixels[1] = 1;
         lpixels[0] = 1;
         lpixels[1] = 1;
         for (s = 0; s < Nbr_Plan; s++)
             {
             lpixels[0] = fpixels[0] + scale(s).nc() -1; // last element exclus
             lpixels[1] = fpixels[1] + scale(s).nl() -1;
             Ptr = band(s).buffer();
             if ( ffpsse(fptr, 1, naxis, naxes, fpixels, lpixels, Ptr, &status))
                 PrintError( status );  
             fpixels[0] = lpixels[0]+1;
             fpixels[1] = lpixels[1]+1;
             }
         break;
     case TRANSF_MALLAT:
     case TRANSF_FEAUVEAU:
            // this works because everything is put in one image
         nelements = naxes[0] * naxes[1];
         Ima.alloc(size_ima_nl(),size_ima_nc(), "IO write");
         ortho_trans_to_ima((*this), Ima);
         if ( ffppre(fptr, group, firstpixel, nelements, Ima.buffer(), 
               &status) )
              PrintError( status );  
         break;
     default:
         fprintf (stderr, "Error in mr_io_write: bad Type_Transform");
         break; 
     }

 /* close the FITS file */
 if ( ffclos(fptr, &status) )  PrintError( status );  
// cout << " end of write fits " << endl;

}

/****************************************************************************/

void MultiResol::write_band (char *Name, int SelectBand)
{
    io_write_ima_float (Name, band(SelectBand));
}

/****************************************************************************/

void MultiResol::read_band (char *Name, int SelectBand)
{
    Ifloat DataScale;     
    io_read_ima_float (Name, DataScale);
    band(SelectBand) = DataScale;
}

/****************************************************************************/

static void diadic_to_ima(MultiResol & MR_Data, Ifloat & Data, int SelectScale)
{
   int i,j,Nl = MR_Data.size_ima_nl();
   int Nc = MR_Data.size_ima_nc();
   if (SelectScale <  MR_Data.nbr_scale() - 1)
   {
      Data.resize(Nl,Nc*2);
      for (i = 0; i < Nl; i++)
      for (j = 0; j < Nc; j++)
      {
         Data(i,j) = MR_Data(2*SelectScale,i,j);
	 Data(i,j+Nc) = MR_Data(2*SelectScale+1,i,j);
      }
   }
   else Data = MR_Data.band( 2* (MR_Data.nbr_scale()-1));
}

/****************************************************************************/

static void undec_3dir_to_ima(MultiResol & MR_Data, Ifloat & Data, int SelectScale)
{
   int i,j;
   int NumBand = 3*SelectScale;
   int Nl = MR_Data.size_band_nl(NumBand);
   int Nc = MR_Data.size_band_nc(NumBand);
   if (SelectScale <  MR_Data.nbr_scale() - 1)
   {
      Data.resize(Nl*2,Nc*2);
      for (i = 0; i < Nl; i++)
      for (j = 0; j < Nc; j++)
      {
         Data(i,j) = MR_Data(NumBand,i,j);
	 Data(i,j+Nc) = MR_Data(NumBand+1,i,j);
	 Data(i+Nl,j+Nc) = MR_Data(NumBand+2,i,j);
      }
   }
   else Data = MR_Data.band(MR_Data.nbr_band()-1);
}

/****************************************************************************/

void MultiResol::write (char *Name, int SelectScale)
{
   Ifloat DataScale;
   
   if ((Set_Transform == TRANSF_MALLAT) || 
             (Set_Transform == TRANSF_FEAUVEAU))
   {
       ortho_trans_to_ima((*this), DataScale, SelectScale);
      io_write_ima_float (Name, DataScale);
   }
   else if (Set_Transform == TRANSF_DIADIC_MALLAT)
   {
       diadic_to_ima ((*this), DataScale, SelectScale);
       io_write_ima_float (Name, DataScale);
   } 
   else if (Set_Transform == TRANSF_UNDECIMATED_MALLAT)
   {
      undec_3dir_to_ima ((*this), DataScale, SelectScale);
      io_write_ima_float (Name, DataScale);
   }
   else io_write_ima_float (Name, band(SelectScale));
}

/****************************************************************************/

static void ima_to_diadic(MultiResol & MR_Data, Ifloat & Data, int SelectScale)
{
   int i,j,Nl = MR_Data.size_ima_nl();
   int Nc = MR_Data.size_ima_nc();
   if (SelectScale < MR_Data.nbr_scale()  - 1)
   {
      for (i = 0; i < Nl; i++)
      for (j = 0; j < Nc; j++)
      {
         MR_Data(2*SelectScale,i,j) = Data(i,j);
	 MR_Data(2*SelectScale+1,i,j) = Data(i,j+Nc);
      }
   }
   else MR_Data.band(2* (MR_Data.nbr_scale()-1)) =  Data;
}

/****************************************************************************/

static void ima_to_undec_dir3(MultiResol & MR_Data, Ifloat & Data, int SelectScale)
{
   int i,j;
   int NumBand = 3*SelectScale;
   int Nl = MR_Data.size_band_nl(NumBand);
   int Nc = MR_Data.size_band_nc(NumBand);
   if ((Nl*2 != Data.nl()) || (Nc*2 != Data.nc()))
   {
       cout << "Error: cannot insert this image at scale " << SelectScale+1 << endl;
       exit(-1);
   }
   if (SelectScale < MR_Data.nbr_scale()  - 1)
   {
      for (i = 0; i < Nl; i++)
      for (j = 0; j < Nc; j++)
      {
         MR_Data(NumBand,i,j) = Data(i,j);
	 MR_Data(NumBand+1,i,j) = Data(i,j+Nc);
	 MR_Data(NumBand+2,i,j) = Data(i+Nl,j+Nc);
      }
   }
   else  MR_Data.band(MR_Data.nbr_band()-1) = Data;
}

/****************************************************************************/

void MultiResol::read (char *Name, int SelectScale)
{
   Ifloat DataScale;

   if ((Set_Transform == TRANSF_MALLAT) || 
             (Set_Transform == TRANSF_FEAUVEAU))
   {
      io_read_ima_float (Name, DataScale);
      ima_to_ortho_trans((*this), DataScale, SelectScale);
   }
   else if (Set_Transform == TRANSF_DIADIC_MALLAT)
   {
      io_read_ima_float (Name, DataScale);
      ima_to_diadic((*this), DataScale, SelectScale);
   } 
   else if (Set_Transform == TRANSF_UNDECIMATED_MALLAT)
   {
      io_read_ima_float (Name, DataScale);
      ima_to_undec_dir3((*this), DataScale, SelectScale);
   } 
   else 
   {
      io_read_ima_float (Name, DataScale);
      band(SelectScale) = DataScale;
   }
}

/****************************************************************************/

void MultiResol::read (char *Name)
{
    // for fits
    char filename[256];
    fitsfile *fptr;           /* pointer to the FITS file */
    int status=0, hdutype ;
    long hdunum;
    char comment[FLEN_COMMENT];
    int s,i,j,naxis;
    long naxes[3];
    long mon_long;
    int anynul = 0;
    long nulval = 0;
    long inc[3];
    void PrintError( int status);
    long nelements = 0 ; // naxes[0] * naxes[1] in the image
    Ifloat Ima;
    long fpixels[3];
    long int lpixels[3];
    
     // for multiresol
    float *Ptr;
    int my_logical; // sais pas...

     mr_io_name (Name, filename);

    inc[0]=1;  inc[1]=1; inc[2]=1;
 
#if DEBUG_IO
    cout << "Read in " << filename << endl;
#endif
   
    /* open the file */
    status = 0;         /* initialize status before calling fitsio routines */
    if ( ffopen(&fptr, filename, (int) READONLY, &status) ) 
         PrintError( status );
 // cout << " open the file " << endl;
                                   
    // ******* read the HEADER  *******

    hdunum = 1;  /*read  table */
    if ( ffmahd(fptr, hdunum, &hdutype, &status) ) /* move to the HDU */
           PrintError( status );

    int simple, bitpix, extend;
    long pcount, gcount;
    if ( ffghpr(fptr, 3, &simple, &bitpix, &naxis, naxes, &pcount,
            &gcount, &extend, &status)) /* move to the HDU */
           PrintError( status );

     nelements = naxes[0] * naxes[1];
    // cout << " begin to read " << endl;
     /* read Number of lines, columns, plans */
    if (ffgkyj(fptr,(char*)"Nl", &mon_long, comment, &status)) PrintError( status );
    int Nli = (int)mon_long; /* there is no function for reading int */
    if (ffgkyj(fptr,(char*)"Nc", &mon_long, comment, &status)) PrintError( status );
    int Nci = (int)mon_long; 
    if (ffgkyj(fptr,(char*)"Nbr_Plan", &mon_long, comment, &status)) PrintError( status );
    int NbrPlan = (int)mon_long; 
    
    type_transform  TT=DEFAULT_TRANSFORM;
    if (ffgkyj(fptr,(char*)"Type_Tra", &mon_long, comment, &status)) 
                    PrintError( status );
    else TT = (type_transform) mon_long;
    
    if (TT == TC_FCT)
    {
      if (ffgkyj(fptr,(char*)"NbrDir", &mon_long, comment, &status)) PrintError( status );
      int Ndir  = (int) mon_long;
      FCT_Nbr_Angle = Ndir;
      
//       if (ffgkyj(fptr,"NbrScale", &mon_long, comment, &status)) PrintError( status );
//       int NbrScale2D = (int) mon_long;
//       if ( ffgkyj(fptr,"Extend",&mon_long, comment, &status))
//       PrintError( status );    
//       Bool Ext  = (mon_long == 0) ? False : True;
//       if ( ffgkyj(fptr,"Isotrop",&mon_long, comment, &status))
//       PrintError( status );
//       Bool Iso  = (mon_long == 0) ? False : True;
//       if ( ffgkyj(fptr,"Real",&mon_long, comment, &status))
//       PrintError( status );    
//       Bool Real  = (mon_long == 0) ? False : True;
    }
    
    // cout <<"HH alloc" << endl;
    // alloc(Nli,Nci, NbrPlan, (type_transform) mon_long, "MR read");
    // cout <<"HH alloc 1" << endl;

    FilterAnaSynt *PtrFAS = NULL;
    int NU=-1;
    if ((TT == TO_MALLAT) || (TT == TO_UNDECIMATED_MALLAT) || (TT == TO_LC))
    {
        if (ffgkyj(fptr,(char*)"SBFilter", &mon_long, comment, &status)) SBFilter = DEF_SB_FILTER;
        else SBFilter = (type_sb_filter) mon_long;
	
	if (ffgkyj(fptr,(char*)"NORM", &mon_long, comment, &status)) TypeNorm = DEF_SB_NORM;
        else  TypeNorm = (sb_type_norm) mon_long;

       if (ffgkyj(fptr,(char*)"UNDEC", &mon_long, comment, &status))  NbrUndecimatedScale = -1;
       else NU = (int) mon_long;
 
      if (SBFilter == F_USER)
      {
          char FBName[256];
          if (ffgkys(fptr,(char*)"FilBank", FBName, comment, &status))
                                                      PrintError( status );
          if (FBName[0] == '-') UserFilterFileName = NULL; 
          else UserFilterFileName = strdup(FBName);
      }
    
       PtrFAS = new FilterAnaSynt;
       PtrFAS->Verbose = Verbose;
       PtrFAS->alloc(SBFilter);
       FilterBankAlloc=True;
    }
    if ((TT  == TO_DIV_1) || (TT  == TO_DIV_2))
    {
 	
	if (ffgkyj(fptr,(char*)"NORM", &mon_long, comment, &status)) TypeNorm = DEF_SB_NORM;
        else  TypeNorm = (sb_type_norm) mon_long;

       if (ffgkyj(fptr,(char*)"UNDEC", &mon_long, comment, &status))  NbrUndecimatedScale = -1;
       else NU = (int) mon_long;
 
       //PtrFAS = new FilterAnaSynt;
       //PtrFAS->Verbose = Verbose;
       //PtrFAS->alloc(SBFilter);
       //FilterBankAlloc=True;
    }
    
    
    type_undec_filter UF=DEF_UNDER_FILTER;
    if (TT  == TO_UNDECIMATED_NON_ORTHO)
    {
       if (ffgkyj(fptr,(char*)"UFilBank", &mon_long, comment, &status)) UF = DEF_UNDER_FILTER;
       else UF = (type_undec_filter) mon_long;
    }
    // cout <<"HH alloc" << endl;
    alloc (Nli, Nci, NbrPlan, TT, PtrFAS, TypeNorm, NU, UF);
    // cout <<"HH1" << endl;

    if (Type_Transform  == TO_LIFTING)
    {
         if (ffgkyj(fptr,(char*)"LiftingT", &mon_long, comment, &status)) 
              LiftingTrans = DEF_LIFT;
 	 else LiftingTrans = (type_lift) mon_long;
    }
 
    if (ffgkyj(fptr,(char*)"FormatIn", &mon_long, comment, &status))
         PrintError( status );
    FormatInputImag = (type_format)mon_long;

    if (ffgkys(fptr,(char*)"Name_MR", Name_MR, comment, &status))
         PrintError( status );

    if (ffgkyj(fptr,(char*)"Border", &mon_long, comment, &status))
         PrintError( status );
    Border = (type_border)mon_long; 

    if (ffgkyj(fptr,(char*)"MedianWi", &mon_long, comment, &status))
         PrintError( status );
    MedianWindowSize = (int)mon_long;
   
    if (ffgkye(fptr,(char*)"Fc", &Fc, comment, &status))
         PrintError( status ); 

    if (ffgkyj(fptr,(char*)"Nbr_Iter", &mon_long, comment, &status))
         PrintError( status );
    Nbr_Iter = (int)mon_long; 

    if (ffgkyl(fptr,(char*)"ExactPyr", &my_logical, comment, &status))
         PrintError( status ); 
    ExactPyrRec = (Bool) my_logical;
 

#if DEBUG_IO 
    cout << "Read in " << filename << endl;
    cout << "Nl = " << size_ima_nl () << endl;
    cout << "Nc = " << size_ima_nc () << endl;
    cout << "Nbr_Plan = " << nbr_scale () << endl;
    cout << "Type_Transform = " << Type_Transform << " " <<
            StringTransform(Type_Transform) << endl;
    cout << "Set_Transform = " << Set_Transform << " " <<
            StringSetTransform(Set_Transform) <<  endl;
#endif

    // ******* read the images *******
    fpixels[0] =  (long ) 1;
    fpixels[1] =  1;
    fpixels[2] =  1;
    lpixels[0] =  Nc;
    lpixels[1] =  Nl;
    lpixels[2] =  1;

    switch (Set_Transform)
    {
       case TRANSF_DIADIC_MALLAT:
       case TRANSF_PAVE:
          for (s = 0; s < NbrBand; s++)
             {
             Ptr = band(s).buffer();
             if ( ffgpve(fptr, 1, fpixels[0],  nelements, nulval, 
                         Ptr, &anynul, &status))
                 PrintError( status );
             fpixels[0] += nelements ;
             }
           break;
       case TRANSF_UNDECIMATED_MALLAT:
       case  TRANSF_SEMIPYR:
         if (Type_Transform == TC_FCT)
	 {
	      nelements = naxes[0];
	      fltarray Tab(nelements);
              Ptr = Tab.buffer();
              if ( ffgpve(fptr, 1, 1, nelements, nulval, Ptr, &anynul, &status))
                    PrintError( status );
              int ind = nbr_scale()+1;
	      for (int b=0; b < nbr_band(); b++)
              {
                int Ny = (int) Tab(ind++);
                int Nx = (int) Tab(ind++);
                for (int i=0;  i < size_band_nl(b); i++)
                for (int j=0;  j < size_band_nc(b); j++)  (*this)(b,i,j) = Tab(ind++);
		int s1 = FCT_Band2Scale(b);
		int b1 = FCT_Band2Dir(b);
		// cout << " Band " << b << "  Scake " << s1 << " Subb " << b1 << " Nl " << Ny << " Nc " << Nx << endl;
		((FCT->TabCF_Band)[s1][b1]).alloc(Ny ,Nx);
              }
	 }
	 else
	 {
	  Ima.alloc(Nl,Nc, "IO ima read");
          for (s = 0; s < NbrBand; s++)
             {
             Ptr = Ima.buffer();
             if ( ffgpve(fptr, 1, fpixels[0],  nelements, nulval, 
                         Ptr, &anynul, &status))
                 PrintError( status );
             fpixels[0] += nelements ;
	     for (i =0; i < size_band_nl(s); i++)
	     for (j =0; j < size_band_nc(s); j++) (*this)(s,i,j) = Ima(i,j);
             }
	   }
           break;
	case TRANSF_PYR:
          lpixels[0] = 1;
          lpixels[1] = 1;
          for (s = 0; s < Nbr_Plan; s++)
          {
             lpixels[0] = fpixels[0] + scale(s).nc() -1; 
             lpixels[1] = fpixels[1] + scale(s).nl() -1;
             Ptr = band(s).buffer();
            // cout << "fpixels=" << fpixels[0] <<" " << fpixels[1] << endl;
            // cout << "lpixels=" << lpixels[0] << " " << lpixels[1] << endl;
            // colnum = 0 !!! this because ffgsve calls ffgcle

             if ( ffgsve(fptr, 0 , naxis, naxes, fpixels, lpixels, inc, nulval,
                  Ptr, &anynul, &status))
                 PrintError( status );
             fpixels[0] = lpixels[0]+1;
             fpixels[1] = lpixels[1]+1;
          }
          break; 
        case TRANSF_MALLAT:
        case TRANSF_FEAUVEAU:
             // this works because everything is put in one image
             Ima.alloc(Nl,Nc, "IO ima read");
             Ptr = Ima.buffer();
             if ( ffgpve(fptr, 1, 1, nelements, nulval,Ptr, &anynul, &status))
             PrintError( status );
             ima_to_ortho_trans((*this), Ima);
             break;
        default:
          fprintf (stderr, "Error in mr_io_read: bad Set_Transform");
          break; 
    }
    strcpy(Name_MR,Name);

/* close the FITS file */
 if ( ffclos(fptr, &status) ) PrintError( status );
// cout << " end of read fits file " << endl;

#if DEBUG_IO
    cout << "Read out " << filename << endl;
#endif
}

/****************************************************************************/

/*--------------------------------------------------------------------------*/
static void PrintError( int status)
{
    /*****************************************************/
    /* Print out cfitsio error messages and exit program */
    /*****************************************************/

    char status_str[FLEN_STATUS], errmsg[FLEN_ERRMSG];
  
    if (status)
      fprintf(stderr, "\n*** Error occurred during program execution ***\n");

    ffgerr(status, status_str);        /* get the error status description */
    fprintf(stderr, "\nstatus = %d: %s\n", status, status_str);

    if ( ffgmsg(errmsg) )  /* get first message; null if stack is empty */
    {
         fprintf(stderr, "\nError message stack:\n");
         fprintf(stderr, " %s\n", errmsg);

         while ( ffgmsg(errmsg) )  /* get remaining messages */
             fprintf(stderr, " %s\n", errmsg);
    }

    exit( status );       /* terminate the program, returning error status */
}

/*--------------------------------------------------------------------------*/
