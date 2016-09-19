/*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  20/03/00 
**    
**    File: Filter.cc
**
**    Modification history:
**
*******************************************************************************
**
**    DESCRIPTION  Filter design
**    -----------  
**                 
******************************************************************************/
 
#include "Filter.h"

// Antoni 7/9 filters
/*
static float Filter7_9_h0[9] =
{
	 0.02674875741,
   	-0.0168641184 ,
	-0.0782232665 ,
	 0.26686411844,
	 0.60294901823,
	 0.26686411844,
	-0.0782232665 ,
	-0.0168641184 ,
	 0.02674875741};

static float Filter7_9_g0[7] =
{
	 0.04563588155 ,
	-0.02877176311 ,
        -0.295635881557,
	 0.557543526229,
        -0.295635881557,
   	-0.02877176311 ,
	 0.04563588155 };

static float Filter7_9_h1[7] = 
{
	-0.04563588155 , 
	-0.02877176311 , 
         0.295635881557, 
	 0.557543526229,  
         0.295635881557, 
   	-0.02877176311 ,
	-0.04563588155};

static float Filter7_9_g1[9] =
{
	 0.02674875741,
   	 0.0168641184 ,
	-0.0782232665 ,
	-0.26686411844,
	 0.60294901823,
	-0.26686411844,
	-0.0782232665 ,
	 0.0168641184 ,
	 0.02674875741};
*/
static float  AntoniniSynthesis [] = { -6.453888262893856e-02,
			      -4.068941760955867e-02,
			       4.180922732222124e-01,
			       7.884856164056651e-01,
			       4.180922732222124e-01,
			      -4.068941760955867e-02,
			      -6.453888262893856e-02 };
static float   AntoniniAnalysis[] =  {  3.782845550699535e-02,
			     -2.384946501937986e-02,
			     -1.106244044184226e-01,
			      3.774028556126536e-01,
			      8.526986790094022e-01,
			      3.774028556126537e-01,
			     -1.106244044184226e-01,
			     -2.384946501937986e-02,
			      3.782845550699535e-02 };
			      
const float  SQRT2 = sqrt(2.);
static float HaarCoeffs [3] = {1./SQRT2, 1./SQRT2, 0};

static float Daub4Coeffs [7] = {0,0,0,
                                 0.4829629131445341,  
                                 0.8365163037378077,
		                 0.2241438680420134, 
			        -0.1294095225512603};
				
// J. E. Odegard and C. S. Burrus, 
//  "Smooth biorthogonal wavelets for applications in image compression," 
// in Proceedings of DSP Workshop, Loen, Norway, September
// 1996 (http://www-dsp.rice.edu/publications). 
				
static float OdegardAnalysis[9] = {
   5.2865768532960523e-02,
  -3.3418473279346828e-02,
  -9.3069263703582719e-02,
   3.8697186387262039e-01,
   7.8751377152779212e-01,
   3.8697186387262039e-01,
  -9.3069263703582719e-02,
  -3.3418473279346828e-02,
   5.2865768532960523e-02
};

static float OdegardSynthesis[7] = {
  -8.6748316131711606e-02,
  -5.4836926902779436e-02,
   4.4030170672498536e-01,
   8.1678063499210640e-01,
   4.4030170672498536e-01,
  -5.4836926902779436e-02,
  -8.6748316131711606e-02
};

static float HaarAnalysis[3] = {
    0,
    0.707106781187,
    0.707106781187,
};

static float Haar2Synthesis[7] = {
   0,
   -0.088388347648,
   0.088388347648,
   0.707106781187 ,
   0.707106781187,
   0.088388347648,
   -0.088388347648
};
static float Haar4Synthesis[11] = {
   0,
   0.016572815184,
  -0.016572815184,
  -0.121533978016,
  0.121533978016,
  0.707106781187,
  0.707106781187,
  0.121533978016,
  -0.121533978016,
  -0.016572815184,
  0.016572815184
};

static float F5_3Analysis[7] = {
   0,
   -0.125,
   0.25,
  0.75,
  0.25,
  -0.125,
  0
};
static float F5_3Synthesis[5] = {
  0,
  0.25,
  0.5,
  0.25,
  0
};


// float h[2][5]={ {0.,-1./4.,3./4.,3./4.,-1./4.}};
 // float h1[2][5]={1./2.,1.,1./2.,0.,0.}     {{0.,-1./4.,3./4.,3./4.,-1./4.}};
 
static float DF4_4Analysis[5] = {
   0,
   -0.25,
   0.75,
  0.75,
  -0.25
};
static float F4_4Synthesis[5] = {
  0.125,
  0.375,
  0.375,
  0.125,
  0
};

static float F4_4Analysis[5] = {
   -0.25,
   0.75,
  0.75,
  -0.25,
  0,
};
static float EF4_4Synthesis[5] = {
  0,
  0.125,
  0.375,
  0.375,
  0.125
};


static float Lemarie1[23] = {
 -1.2268637e-004,
 -2.2429586e-004,
  5.1163599e-004,
  9.2337185e-004,
 -2.2019463e-003,
 -3.8832637e-003,
  9.9906062e-003,
  1.6974818e-002,
 -5.1945376e-002,
 -6.9101073e-002,
  3.9729673e-001,
  8.1764658e-001,
  3.9729673e-001,
 -6.9101073e-002,
 -5.1945376e-002,
  1.6974818e-002,
  9.9906062e-003,
 -3.8832637e-003,
 -2.2019463e-003,
  9.2337185e-004,
  5.1163599e-004,
 -2.2429586e-004,
 -1.2268637e-004
};

static float Lemarie3[41] = {
  1.4609806e-004,
 -2.3230422e-004,
 -2.8541356e-004,
  4.6209256e-004,
  5.5995183e-004,
 -9.2718607e-004,
 -1.1037477e-003,
  1.8821190e-003,
  2.1867121e-003,
 -3.8824237e-003,
 -4.3538374e-003,
  8.2014715e-003,
  8.6852877e-003,
 -1.7982279e-002,
 -1.7176319e-002,
  4.2068300e-002,
  3.2080847e-002,
 -1.1003691e-001,
 -5.0201719e-002,
  4.3392285e-001,
  7.6612988e-001,
  4.3392285e-001,
 -5.0201719e-002,
 -1.1003691e-001,
  3.2080847e-002,
  4.2068300e-002,
 -1.7176319e-002,
 -1.7982279e-002,
  8.6852877e-003,
  8.2014715e-003,
 -4.3538374e-003,
 -3.8824237e-003,
  2.1867121e-003,
  1.8821190e-003,
 -1.1037477e-003,
 -9.2718607e-004,
  5.5995183e-004,
  4.6209256e-004,
 -2.8541356e-004,
 -2.3230422e-004,
  1.4609806e-004};

static float Lemarie5[59] = {
  1.4299532e-004,
  1.5656611e-004,
 -2.2509746e-004,
 -2.4421337e-004,
  3.5556002e-004,
  3.8161407e-004,
 -5.6393016e-004,
 -5.9748378e-004,
  8.9882146e-004,
  9.3739129e-004,
 -1.4412528e-003,
 -1.4737089e-003,
  2.3286290e-003,
  2.3211761e-003,
 -3.7992267e-003,
 -3.6602095e-003,
  6.2791340e-003,
  5.7683203e-003,
 -1.0562022e-002,
 -9.0493510e-003,
  1.8208558e-002,
  1.4009689e-002,
 -3.2519969e-002,
 -2.1006296e-002,
  6.1312356e-002,
  2.9474179e-002,
 -1.2926869e-001,
 -3.7019995e-002,
  4.4246341e-001,
  7.4723338e-001,
  4.4246341e-001,
 -3.7019995e-002,
 -1.2926869e-001,
  2.9474179e-002,
  6.1312356e-002,
 -2.1006296e-002,
 -3.2519969e-002,
  1.4009689e-002,
  1.8208558e-002,
 -9.0493510e-003,
 -1.0562022e-002,
  5.7683203e-003,
  6.2791340e-003,
 -3.6602095e-003,
 -3.7992267e-003,
  2.3211761e-003,
  2.3286290e-003,
 -1.4737089e-003,
 -1.4412528e-003,
  9.3739129e-004,
  8.9882146e-004,
 -5.9748378e-004,
 -5.6393016e-004,
  3.8161407e-004,
  3.5556002e-004,
 -2.4421337e-004,
 -2.2509746e-004,
  1.5656611e-004,
  1.4299532e-004
};

char *UserFilterFileName=NULL;

/***********************************************************************/

type_sb_filter get_filter_bank(char *UserArg)
{
   int c1;
   char *ch = new char[256];
   type_sb_filter FilRet = DEF_SB_FILTER;

   int N = sscanf(UserArg,"%d,%s",&c1,ch);
   // cout << "N = " << N << endl;
   if (N < 1)
   {
      fprintf(OUTMAN, "bad type of filter: %s\n", UserArg);
      exit(-1);
   }
   if (N > 0) 
   {
       FilRet = (type_sb_filter) c1;
       if ((c1 < 1) || (c1 > NBR_SB_FILTER))
       {
	   fprintf(OUTMAN, "bad type of filter: %s\n", UserArg);
	   exit(-1);
       }
       // cout << "Filter bank = " << StringSBFilter(FilRet) << endl;
       if (N > 1) 
       {
          UserFilterFileName = ch;
          // cout << "User file name = " << UserFilterFileName  << endl;
       }
   }
   return FilRet;
}

/***********************************************************************/

const char * StringSBFilter (type_sb_filter  type)
{
    switch (type)
    {
        case F_BI2HAAR:
	      return ("Biorthogonal 2/6 Haar filters"); 
        case F_BI4HAAR:
	      return ("Biorthogonal 2/10 Haar filters");
        case  F_MALLAT_7_9:
	      return ("Biorthogonal 7/9 filters"); 
        case  F_DAUBE_4: 
              return ("Daubechies filter 4"); 
        case  F_5_3: 
              return ("5/3 filter"); 
	 case  F_3_5:
              return ("3/5 filter"); 
	case  F_MALLAT_9_7: 
              return ("Biorthogonal 9/7 filter");
	case  F_HAAR: 
              return ("Haar filter"); 
        case  F_ODEGARD_7_9:
	      return ("Odegard 9/7 filters"); 
        case  F_LEMARIE_1:
	      return ("Battle-Lemarie filters (2 vanishing moments)"); 
	case  F_LEMARIE_3:
	      return ("Battle-Lemarie filters (4 vanishing moments)"); 
	case  F_LEMARIE_5:
	      return ("Battle-Lemarie filters (6 vanishing moments)"); 
	case  F_4_4:
	      return ("4/4 Linar spline filters"); 
	case  F_USER:
	      return ("User's filters");

        default: 
	      return ("Undefined sub-band filters");
    }
}
/***********************************************************************/

void sb_usage(type_sb_filter Filter)
{
    fprintf(OUTMAN, "         [-T type_of_filters]\n");
    for (int i = 1; i <= NBR_SB_FILTER; i++)
    fprintf(OUTMAN, "              %d: %s \n",i,
                                           StringSBFilter((type_sb_filter  )i));
    fprintf(OUTMAN, "             default is %s\n\n", StringSBFilter ((type_sb_filter) Filter));
     
    fprintf(OUTMAN, "         [-L]\n");
    fprintf(OUTMAN, "              Use a L2 normalization. Default is L1.\n");
}


/***********************************************************************/
 

char *filtername(char * NameStep)
{
    char Sreturn[256];
    char Name[256];
    char *ValRet;

    strcpy(Name, NameStep);
    if (strstr(Name, ".wvf") != NULL)  strcpy(Sreturn, Name);
    else sprintf(Sreturn, "%s.%s", Name, "wvf");
    ValRet = strdup(Sreturn);
    return (ValRet);
}

/***********************************************************************/

FilterAnaSynt::FilterAnaSynt(type_sb_filter T_Filter) 
{
    reset_param();
    alloc(T_Filter);
}

/***********************************************************************/

FilterAnaSynt::FilterAnaSynt(char *FileName) 
{
    reset_param();
    FilterFileName=FileName; 
    alloc(F_USER);
}

/***********************************************************************/

FilterAnaSynt::~FilterAnaSynt() 
{
    reset_param();
}

/***********************************************************************/

void FilterAnaSynt::reset_param()
{
   Analysis = Synthesis = NULL;
   Size_Ana  = Size_Synt = 0;
   Start_Ana = Start_Synt =0;   
   TypeFilter = SB_UNKNOWN;
   Verbose = False;
   FilterFileName = NULL;
}

/***********************************************************************/

void  FilterAnaSynt::alloc(type_sb_filter T_Filter)
{
    TypeFilter = T_Filter;
    TypeNorm =  NORM_L2;
    switch(T_Filter)
    {
        case F_LEMARIE_1 :
 	   Analysis = Lemarie1;
	   Synthesis = Lemarie1;
           Size_Ana = 23;
           Size_Synt = Size_Ana;
 	   break;   
        case F_LEMARIE_3 :
 	   Analysis = Lemarie3;
	   Synthesis = Lemarie3;
           Size_Ana = 41;
           Size_Synt = Size_Ana;
 	   break;  
	case F_LEMARIE_5 :
 	   Analysis = Lemarie5;
	   Synthesis = Lemarie5;
           Size_Ana = 59;
           Size_Synt = Size_Ana;
 	   break;  
	case F_ODEGARD_7_9 :
 	   Analysis = OdegardAnalysis;
	   Synthesis = OdegardSynthesis;
           Size_Ana = 9;
           Size_Synt = 7; 
 	   break;        
        case F_MALLAT_9_7:
	   Analysis = AntoniniAnalysis;
	   Synthesis = AntoniniSynthesis;
           Size_Ana = 9;
           Size_Synt = 7; 
 	   break;
	 case F_MALLAT_7_9:
	   Analysis = AntoniniSynthesis;  
	   Synthesis = AntoniniAnalysis;
           Size_Ana = 7;
           Size_Synt = 9;
 	   break;
	 case F_HAAR:
  	    Size_Synt = 3;
	    Size_Ana = 3;
 	    Analysis = HaarCoeffs;
	    Synthesis = HaarCoeffs;
            break;	
	 case F_3_5:
	    Size_Synt = 5;
	    Size_Ana = 3; 
 	    Synthesis = F5_3Analysis;
	    Analysis = F5_3Synthesis;
	    TypeNorm =  NORM_L1;
            break;
	 case F_4_4:
  	    Size_Synt = 5;
	    Size_Ana = 5; 
 	    Analysis = F4_4Analysis;
	    Synthesis = F4_4Synthesis;
	    TypeNorm =  NORM_L1;
            break;     
	 case F_5_3:
  	    Size_Synt = 5;
	    Size_Ana = 7; 
 	    Analysis = F5_3Analysis;
	    Synthesis = F5_3Synthesis;
	    TypeNorm =  NORM_L1;
            break;       
  	 case F_BI2HAAR:
  	    Size_Synt = 7;
	    Size_Ana = 3;
 	    Analysis = HaarAnalysis;
	    Synthesis = Haar2Synthesis;
            break;	
   	 case F_BI4HAAR:
  	    Size_Synt = 11;
	    Size_Ana = 3;
 	    Analysis = HaarAnalysis;
	    Synthesis = Haar4Synthesis;
            break;
	 case F_DAUBE_4:
 	    Size_Synt = 7;
	    Size_Ana = 7;
 	    Analysis = Daub4Coeffs;
 	    Synthesis = Daub4Coeffs;
            break;
	 case F_USER:
            if (FilterFileName == NULL)
            {
               if (UserFilterFileName != NULL) 
                     FilterFileName = UserFilterFileName;
               else
               {
                  FilterFileName = DEF_USER_FILTER_FILE_NAME;
                  FILE *FileDes = fopen(FilterFileName, "r");
                  if (FileDes != NULL) fclose(FileDes);
                  else
                  {
                      FilterFileName = (char *) getenv(USER_FILER_FILE_NAME);
                      if (FilterFileName == NULL)
                      { 
                         cout << "Error: the filter bank is not defined ... " << endl;
                         exit(-1);
                      }
                  }
               }
            }
            read_from_file(FilterFileName);
            break;
         default:
           cerr << "Error: unknown filter ... " << endl;
           exit(-1);
         break;
    }   
}

/***********************************************************************/

void  FilterAnaSynt::read_from_file(char *FileName)
{
   //  [Range low] [Range high]
   //  [Analysis LP filter coefficients]
   //  .
   //  .
   //  .
   //  [Range low] [Range high]
   //  [Synthesis LP filter coefficients]
   FILE *input=NULL;
   int i,ind,AnaLow, AnaHigh, SyntLow, SyntHigh;
   float Val;
   char *FName = filtername(FileName);
   
   input = fopen(FName,"r");
   if (input == NULL) 
   {
        cout << "Error: cannot open file " <<   FileName  << " ... or file doesn't exist" << endl;
        exit(-1);
   }
   if (fscanf(input,"%d %d",&AnaLow, &AnaHigh) != 2) 
   {
      cout << "Error: bad filter file format ... " << endl;
      exit(-1);
   }
   int nLine = AnaHigh - AnaLow + 1;

   if (Verbose == True)
   {
      cout << "Read filters from file " << FileName << endl;
      cout << "  Analysis: Range low = " << AnaLow << "  Range high = " << AnaHigh << " Size = " << nLine <<  endl;
   }
   int M = MAX(AnaHigh, ABS(AnaLow));
   Size_Ana = 2*M+1;
   Start_Ana = - Size_Ana/2;
   Analysis = new float [Size_Ana];
   double Sum=0.;
   double Sum2=0.;
   for (i=0; i < Size_Ana; i++) Analysis[i] = 0.;
   for (i=0; i < nLine; i++)
   {
      if (fscanf(input,"%f",&Val) != 1) 
      {
         cout << "Error: bad filter file format ... " << endl;
         exit(-1);
      }
      ind = Size_Ana /2+i+AnaLow;
      if ((ind < 0) || (ind >= Size_Ana))
      {
        cout << "Error: bad index ind = " << ind << endl;
        exit(-1);
      }
      Analysis[ind] = Val;
      Sum += Val;
      Sum2 += Val*Val;
      // cout << "Val = " << Val << endl;
   }
   if (Verbose == True) cout << "  Sum_i H[i]^2 = " << Sum2 << endl;
   
   if (ABS(Sum2-1.) > ABS(Sum-1.)) TypeNorm = NORM_L1;
   
   if (fscanf(input,"%d %d",& SyntLow, &SyntHigh) != 2) 
   {
      cout << "Error: bad filter file format ... " << endl;
      exit(-1);
   }
   nLine = SyntHigh - SyntLow + 1;

   if (Verbose == True)
   {
       cout << "  Synthesis: Range low = " <<  SyntLow << "  Range high = " <<  SyntHigh << "  Size = " <<   nLine << endl;
   }
   M = MAX(SyntHigh, ABS(SyntLow));
   Size_Synt = 2*M+1;
   Start_Synt = - Size_Synt/2;
   Synthesis = new float [Size_Synt];
   for (i=0; i <  Size_Synt; i++) Synthesis[i] = 0.;
   for (i=0; i < nLine; i++)
   {
      if (fscanf(input,"%f",&Val) != 1) 
      {
         cout << "Error: bad filter file format ... " << endl;
         exit(-1);
      }
      ind =  Size_Synt/2+i+SyntLow;
      if ((ind < 0) || (ind >=  Size_Synt))
      {
        cout << "Error: bad index ind = " << ind << endl;
        exit(-1);
      }
       Synthesis[ind] = Val;
   }
   if (Verbose == True) cout << "  Sum_i H1[i]^2 = " << Sum2 << endl;
   fclose(input);
}

/***********************************************************************/
