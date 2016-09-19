/******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  22/9/98 
**    
**    File:  IM1D_IO.cc
**
*******************************************************************************
**
**    DESCRIPTION   input output routines for 1D signals
**    ----------- 
**                 
**    PARAMETRES    
**    ----------    
** 
**    RESULTS      
**    -------  
**
**
**   
******************************************************************************/



// #include "IM_Obj.h"
// #include "IM_IO.h"
#include"GlobalInc.h"
#include "IM1D_IO.h"
 
type_1d_format IO_1D_Format =  F1D_UNKNOWN;


class CExelDate {
   int Nc;
   char **TabDate;
   void reset()
   { 
      if (TabDate !=  NULL)
      {
         for (int i = 0; i < Nc; i++) 
            if (TabDate[i]  !=  NULL) delete [] TabDate[i];
         delete [] TabDate;
      }
   }
   public:
   int np() { return Nc;}
   CExelDate() {TabDate = NULL;Nc=0;}
   void alloc (int N)
   {
      reset();
      TabDate = new char * [N];
      for (int i = 0; i < N; i++) TabDate[i] = new char [20];
      Nc = N;
   }
   char *date(int i) 
   { 
       if ((i < 0) || (i >= Nc))
       {
          cout << "Error: bad index in CExelDate ... " << endl;
	  exit(-1);
       }
       return TabDate[i];
   }
   ~CExelDate() {reset();}

};
CExelDate CDate;

#ifdef USELM
extern DemoLic MRDEMO;
#endif

/**************************************************************************/

type_1d_format  io_which_1d_format(const char *Format) 
{
    char Name_Up[MAXCHAR];
    char Name[MAXCHAR];
    char * Exel = ".csv";
    char * ascii = ".dat";
    char * Fits = ".fit";
    unsigned int i;

    strcpy(Name_Up, Format);
    type_1d_format Val_IO_1D_Format = F1D_UNKNOWN;
    strcpy(Name, Name_Up);
    for (i=0; i < strlen(Name); i++) Name[i] = tolower(Name[i]);
    if (strstr(Name,  Exel) != NULL)  Val_IO_1D_Format =  F_1DEXEL;
    else if (strstr(Name,  ascii) != NULL)  Val_IO_1D_Format =  F_1DASCII;
    else if (strstr(Name, Fits) != NULL)  Val_IO_1D_Format =  F1D_FITS;
 
    return  Val_IO_1D_Format;
}

/**************************************************************************/

type_1d_format  io_detect_1dformat(const char *Format) 
{
    IO_1D_Format =  F1D_UNKNOWN;
    IO_1D_Format  = io_which_1d_format(Format);
    if (IO_1D_Format  ==  F1D_UNKNOWN) IO_1D_Format  =  DEFAULT_1D_FORMAT;
    return  IO_1D_Format;
}

/**************************************************************************/

static char *asciiname(char * NameStep)
{
    char Sreturn[MAXCHAR];
    char Name[MAXCHAR];
    char *ValRet;

    strcpy(Name, NameStep);
    if (strstr(Name, ".dat") != NULL)  strcpy(Sreturn, Name);
    else sprintf(Sreturn, "%s.%s", Name, "dat");
    ValRet = strdup(Sreturn);
    return (ValRet);
}

/**************************************************************************/

static char *exelname(char * NameStep)
{
    char Sreturn[MAXCHAR];
    char Name[MAXCHAR];
    char *ValRet;

    strcpy(Name, NameStep);
    if (strstr(Name, ".csv") != NULL)  strcpy(Sreturn, Name);
    else sprintf(Sreturn, "%s.%s", Name, "csv");
    ValRet = strdup(Sreturn);
    return (ValRet);
}

/**************************************************************************/

void io_1d_set_format(const char *Format) 
{
     IO_1D_Format = io_which_1d_format(Format);
}


/****************************************************************************/

void io_read1d_ascii(char *Name_Dat_In,  fltarray & Dat, int NbrSkip, int Exell)
{
   int c;
   FILE *input=NULL;
   int Np=0;	
   float Val;
   int i=0;

   // cout << "sk" << endl;
    
   /* Read the size of the input image */
   if (input != stdin) 
   {
      if (Exell == 0) input = fopen( asciiname(Name_Dat_In),"r");
      else input = fopen(exelname(Name_Dat_In),"r");
      if (input == NULL) 
      {
        cout << "Error in allocation of file " <<  Name_Dat_In  << " ... or file doesn't exist" << endl;
        exit(-1);
      }
   }
   //fscanf(input,"%d\n",&Np);
   // Dat.alloc(Np);

    for (i=0; i < NbrSkip; i++) while ( (c=getc(input)) != '\n' );
    i = 0;
    // cout << "sk1" << endl;
    while ((c=getc(input)) != EOF )
    { 
       ungetc(c,input);
       if (Exell == 0)
       {
          if (fscanf(input,"%f",&Val) == 1)  Np++ ;
       }
       else 
       {
          while ((char) (c=getc(input)) != ',' ) if (c == EOF) break;
 	  if (fscanf(input,"%f", &Val) == 1) Np++ ;
        }
    }
    if (input != stdin) fclose(input);
    
   /* Read the size of the input image */
   if (input != stdin) 
   {
      input = fopen( Name_Dat_In,"r");
      if (input == NULL) 
      {
        cout << "Error in allocation of file " <<  Name_Dat_In  << " ... or file doesn't exist" << endl;
        exit(-1);
      }
   }
   // cout << "Data Number = " << Np << endl; 
   Dat.alloc(Np);
   for (i=0; i < NbrSkip; i++)  while ((char) (c=getc(input)) != '\n' );
   if (Exell == 1)  CDate.alloc(Np);
   for (i=0; i < Np; i++)
   {
      if (Exell == 0)  if (fscanf(input,"%f\n",&Val) == 1) Dat(i) = Val;
      if (Exell == 1) 
      {
          int ind = 0;
	  char *Ptr = CDate.date(i);
          while ((char) (c=getc(input)) != ',' ) 
 	  {
  	     Ptr[ind++] = (char) c;
          }
	  Ptr[ind] = '\0';
          if (fscanf(input,"%f\n",&Val))
          {
           Dat(i) = Val;
           // cout << i+1 << " " << CDate.date(i) << " V = " << Dat(i) << endl;
          }
      }
   }
   if (input != stdin) fclose(input);
}

/****************************************************************************/

void io_read2d_ascii(char *Name_Dat_In,  fltarray &Dat)
{
   FILE *input=NULL;
   float Val;
   int i,j;
   int Nl, Nc;
   
   input = fopen(asciiname(Name_Dat_In),"r");
   if (input == NULL) 
   {
     cout << "Error: cannot open  file " <<  Name_Dat_In << endl;
     exit(-1);
   }
   if (fscanf(input,"%d %d",&Nc, &Nl) != 2)
   {
      cout << "Error: cannot read the number of columns and the number of lines in  " << Name_Dat_In << endl;
      exit(-1);
   }
   Dat.alloc(Nc,Nl);
   for (j=0; j < Nc; j++)
   for (i=0; i < Nl; i++)
   {
      if (fscanf(input,"%f",&Val) == 1) Dat(j,i) = Val;
      else
      {
        cout << "Error: cannot read  in " <<  Name_Dat_In  << endl;
        exit(-1);
      }
   }
   fclose(input);
}

/****************************************************************************/

void io_write1d_ascii(char *Name_Dat_Out,  fltarray & Dat, Bool Exell)
{
   io_write_ascii(Name_Dat_Out, Dat, Exell);
}


/****************************************************************************/

void io_write2d_ascii(char *Name_Dat_Out,  fltarray &Dat)
{
   io_write_ascii(Name_Dat_Out, Dat);

//    FILE *PtrFile=NULL;
//    int i,j;
//   
//    PtrFile = fopen(asciiname(Name_Dat_Out),"w");
//    if (PtrFile == NULL) 
//    {
//      cout << "Error: cannot open  file " <<  Name_Dat_Out << endl;
//      exit(-1);
//    }
//    if (fprintf(PtrFile,"%d %d\n", Dat.nx(), Dat.ny()) < 0)
//    {
//       cout << "Error: cannot write the number of columns and the number of lines ... " << endl;
//       exit(-1);
//    }
//    
//    for (j=0; j < Dat.nx(); j++)
//    {
//       for (i=0; i < Dat.ny(); i++)
//       {
//          if (fprintf(PtrFile,"%f ",Dat(j,i)) < 0)  
//          {
//             cout << "Error: cannot write in " <<  Name_Dat_Out  << endl;
//             exit(-1);
// 	 }
//       }
//       fprintf(PtrFile,"\n");
//    }
//    fclose(PtrFile);
}

/****************************************************************************/

void io_write_ascii(char *Name_Dat_Out,  fltarray & Dat, Bool Exell)
{
   FILE *FileOut=NULL;
   int Np=0;	
   int i=0;
   time_t *ptp, tp;
   struct tm var;
   ptp = &tp;
   char Date[MAXCHAR];
   char FullDate[MAXCHAR];

   tp = time(ptp);
   if (tp == -1)
   {
         printf("ERROR in calling time function \n");
         exit(-1);
  }
 
   var = *gmtime(&tp);
   // strftime(Date,26,"date=%d/%b/%Y %H:%M:%S ",&var);
   strftime(Date,26,"%Y",&var);
   // sprintf(FullDate, "\"%1d/%1d/%s\",", var.tm_mday, var.tm_mon, Date);
   sprintf(FullDate, "%1d/%1d/%s,", var.tm_mday, var.tm_mon, Date);
   // cout << "Data = " << FullDate << endl;
 
   /* Read the size of the FileOut image */
   
   if (FileOut != stdout) 
   {
      if (Exell == 0) FileOut = fopen( asciiname(Name_Dat_Out),"w");
      else FileOut = fopen(exelname(Name_Dat_Out),"w");
      
      if (FileOut == NULL) 
      {
        cout << "Error: cannot open file " <<  Name_Dat_Out  << endl;
        exit(-1);
      }
   }
   if (Dat.naxis() == 1)
   {
      Np = Dat.nx();
      // cout << "Np = " << Np << endl;
      for (i=0; i < Np; i++)
      {
         if (Exell == 0) fprintf(FileOut,"%f\n", Dat(i));
         else 
	 {
	    if (i < CDate.np()) 
	      fprintf(FileOut,"%s,%f\n", CDate.date(i), Dat(i));
	    else fprintf(FileOut,"%s%f\n", FullDate, Dat(i));
	 }
      }
      if (FileOut != stdout) fclose(FileOut);
   }
   else if (Dat.naxis() == 2)
   {
      int j,Nl = Dat.ny();
      for (i=0; i < Dat.nx(); i++)
      {
         if (Exell == 1 ) 
	 {
 	    if (i < CDate.np()) 
	      fprintf(FileOut,"%s,", CDate.date(i));
	    else fprintf(FileOut,"%s", FullDate);
	 }
         for (j=0; j < Nl; j++) fprintf(FileOut,"%f  ", Dat(i,j));
         fprintf(FileOut,"\n");
      }     
   }
}


/****************************************************************************/

void io_1d_read_data(char *File_Name, fltarray & Dat, fitsstruct *Header)
{
    int NbrSkip = 0;
    io_1d_read_data(File_Name, Dat, NbrSkip, Header);
}

/****************************************************************************/

void io_1d_read_data(char *File_Name, fltarray & Dat, int NbrSkip)
{
    fitsstruct FitsHeader;    
    io_1d_read_data(File_Name, Dat, NbrSkip, &FitsHeader);
}

/****************************************************************************/

void io_1d_read_data(char *File_Name, fltarray & Dat, int NbrSkip, fitsstruct *Header)
{
    if (IO_1D_Format  == F1D_UNKNOWN)
                     IO_1D_Format  =  io_detect_1dformat(File_Name);
    // cout << "RR " << endl;
    switch (IO_1D_Format)
    {
        case  F_1DEXEL: 
	  // cout << "EXEL" << endl;
	  io_read1d_ascii(File_Name, Dat, NbrSkip, 1);
  	  break;
	case  F_1DASCII: 
	 io_read1d_ascii(File_Name, Dat, NbrSkip, 0);
  	  break;
 	case  F1D_FITS: 
          fits_read_fltarr(File_Name, Dat, Header);
	  break;
	default: cerr << "Error: unknown 1D format ... " << endl;
	         break;
   }
#ifdef USELM
    MRDEMO.test(Dat.nx());
#endif
}

/****************************************************************************/

void io_1d_write_data(char *File_Name, fltarray & Dat, fitsstruct *Header)
{
    if (IO_1D_Format  == F1D_UNKNOWN)
                     IO_1D_Format  =  io_detect_1dformat(File_Name);
    //  cout << "Format = " << String1DFormat(IO_1D_Format) << endl;
    
    switch (IO_1D_Format)
    {
        case  F_1DEXEL: 
	  io_write1d_ascii(File_Name, Dat, True);
  	  break;
	case  F_1DASCII: 
	 io_write1d_ascii(File_Name, Dat, False);
  	  break;
 	case  F1D_FITS: 
         fits_write_fltarr(File_Name, Dat, Header);
	 break;
	default: cerr << "Error: unknown 1D format ... " << endl;
	         break;
   }
}

/****************************************************************************/
