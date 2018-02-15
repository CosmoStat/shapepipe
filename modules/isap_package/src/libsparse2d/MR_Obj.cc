/*******************************************************************************
**
**    UNIT
**
**    Version: 3.2
**
**    Author: Jean-Luc Starck
**
**    Date:  96/06/13 
**    
**    File:  MR_Obj.cc
**
*******************************************************************************
**
**    DESCRIPTION  multiresolution transform
**    -----------  
**                 
*******************************************************************************
**
**  float & MultiResol::operator() (int s,int i,int j, 
**                                  details which_detail) const
**
** Return  a coefficient of the transform
**               
*******************************************************************************
**
** Ifloat & MultiResol::scale (int s) const 
** 
** return a scale of the transform
**
*******************************************************************************
**
** int MultiResol::size_scale_nl (int s, details which_detail) const {
** 
** return the number of lines of the a scale
**
*******************************************************************************
**
** int MultiResol::size_scale_nc (int s) const 
**
** return the number of columns of the a scale
**
*******************************************************************************
**
** Ifloat MultiResol::extract_scale(int s, details which_detail)
** 
** extract a scale and return an image
**
*******************************************************************************
**
** void MultiResol::insert_scale(Ifloat & Image, int s, details which_detail)
**
** insert an image in the transform
**
*******************************************************************************
**
** MultiResol::MultiResol (int Nbr_Line, int Nbr_Col, int Nbr_Scale,
**                        type_transform Transform, char *Name)
**
** Constructor: allocate the memory
** initialize the varaibles
**
** Nbr_Line = number of line of the image
** Nbr_Col  = number of columns of the image
** Nbr_Scale = number of scales of the transform
** Transform = type of transform
** Name = name of the object (usefull for debugging)
** 
*******************************************************************************
**
** MultiResol::~MultiResol()
** 
** Destructor:
**
** deallocates the memory
**
*******************************************************************************
**
** void MultiResol::alloc (int Nbr_Line, int Nbr_Col, int Nbr_Scale,
**                        type_transform Transform, char *Name)
**
** allocate the memory
** initialize the varaibles
**
** Nbr_Line = number of line of the image
** Nbr_Col  = number of columns of the image
** Nbr_Scale = number of scales of the transform
** Transform = type of transform
** Name = name of the object (usefull for debugging)
** 
*******************************************************************************
**
** void MultiResol::free ()
** 
** deallocates the memory
**
******************************************************************************/


#include "MR_Obj.h"
#define EPS 1e-07
// #include "IM_Edge.h"

#define DEBUG_OBJ 0


/************************************************************************/

const char * StringSetTransform (set_transform type)
{
    switch (type)
    {
        case TRANSF_PAVE: 
              return ("cube transform");break;
        case TRANSF_PYR: 
              return ("pyramidal transform");break;
        case TRANSF_SEMIPYR: 
              return ("half-pyramidal transform");break;
        case TRANSF_MALLAT: 
              return ("MALLAT transform: MALLAT");break;
        case TRANSF_DIADIC_MALLAT: 
              return ("Diadic  transform: MALLAT");break;
        case TRANSF_FEAUVEAU: 
              return ("FEAUVEAU transform: FEAUVEAU");break;
        case TRANSF_UNDECIMATED_MALLAT: 
              return ("Undecimated mallat transform");break;
	case S_UNDEFINED: 
              return ("UNDEFINED");break;
        default:
              return ("Undefined transform");
              break;
    }
}

const char * StringTransform (type_transform type)
{
    switch (type)
    {
        case  TO_PYR_MEYER:
	 return ("Meyer's wavelets (compact support in Fourier space)");break;
	case  TO_PYR_MEYER_ISOTROP:
	 return ("Isotropic and compact support wavelet in Fourier space");break;
        case  TO_UNDECIMATED_NON_ORTHO:
	      return ("non orthogonal undecimated transform (three bands per scale)");break;
        case TO_UNDECIMATED_MALLAT:
	      return ("undecimated (bi-) orthogonal transform (three bands per scale)");break;
        case TO_PAVE_HAAR: 
              return ("undecimated Haar transform: a trous algorithm (one band per scale)");break;
        case TO_DIADIC_HAAR: 
              return ("undecimated diadic Haar transform (two bands per scale)");break;
        case TO_PAVE_LINEAR: 
              return ("linear wavelet transform: a trous algorithm");break;
        case TO_PAVE_BSPLINE: 
              return ("bspline wavelet transform: a trous algorithm");break;
        case TO_PAVE_FFT: 
              return ("wavelet transform in Fourier space");break;
        case TM_PAVE_MEDIAN: 
              return ("morphological median transform");break;
        case TM_PAVE_MINMAX:
              return ("morphological minmax transform");break;
        case TO_PYR_LINEAR: 
              return ("pyramidal linear wavelet transform");break;
        case TO_PYR_BSPLINE: 
              return ("pyramidal bspline wavelet transform");break;
        case TO_PYR_FFT_DIFF_RESOL: 
              return ("pyramidal wavelet transform in Fourier space: algo 1 (diff. between two resolutions)");break;
        case TO_PYR_FFT_DIFF_SQUARE: 
              return ("pyramidal wavelet transform in Fourier space: algo 2 (diff. between the square of two resolutions)");break;
        case TM_PYR_MEDIAN: 
              return ("pyramidal median transform (PMT)");break;
        case TM_PYR_MINMAX:
              return ("morphological pyramidal minmax transform");break;
        case TM_PYR_LAPLACIAN: 
              return ("pyramidal laplacian");break;
        case TM_PYR_SCALING_FUNCTION: 
              return ("decomposition on scaling function");break;
       case  TO_SEMI_PYR: 
              return ("half-pyramidal transform");break;
       case  TM_TO_PYR: 
              return ("mixed WT and PMT method (WT-PMT)");break; 
       case TM_TO_SEMI_PYR    : 
              return ("mixed Half-pyramidal WT and Median method (WT-HPMT)");break;   
        case TO_DIADIC_MALLAT: 
              return ("undecimated diadic wavelet transform (two bands per scale)");break;
        //case TO_MAX_DIADIC_MALLAT: 
        //      return ("Maxima diadic wavelet transform");break;
        case TO_MALLAT: 
              return ("Mallat's wavelet transform (7/9 filters)");break;
        case TM_MIN_MAX: 
              return ("G transform (morphological min-max algorithm)");break;
        case TO_FEAUVEAU: 
              return ("Feauveau's wavelet transform");break;
        case TO_HAAR: 
              return ("Haar's wavelet transform");break;
        case TO_PAVE_FEAUVEAU: 
              return ("Feauveau's wavelet transform without undersampling");
              break;
	case TO_LIFTING:
	      return ("Wavelet transform via lifting scheme");
              break;
	case TO_LC:
	      return ("Line Column Wavelet Transform (1D+1D)");
              break;
        case TC_FCT:
	      return ("Fast Curvelet Transform");
              break;
	case TO_DIV_1:
	      return ("5/3 on line and 4/4 on column");
              break;
	case TO_DIV_2:
	      return ("4/4 on line and 5/3 on column");
              break;
	case T_UNDEFINED: 
              return ("Undefined transform");
              break;
        default:
              return ("Undefined transform");
              break;
    }
} 



/************************************************************************/

int number_band_per_resol(type_transform Transform)
{
    switch (SetTransform(Transform))
    {
         case TRANSF_PAVE: 
	 case TRANSF_SEMIPYR:
         case TRANSF_PYR: 
              return 1; break;
        case TRANSF_MALLAT: 
	case TRANSF_UNDECIMATED_MALLAT:
              return 3; break;
        case TRANSF_FEAUVEAU:
	case TRANSF_DIADIC_MALLAT:
              return 2; break;
        default: return -1; break;
    }
}


/************************************************************************/

// Normalisation table
//     extern double TabNormPaveLinear[MAX_SCALE];
//     extern double TabNormPaveB3Spline[MAX_SCALE];
//     extern double TabNormPaveFFT[MAX_SCALE];
//     extern double TabNormPaveMed5[MAX_SCALE];
//     extern double TabNormPaveMed3[MAX_SCALE];
//     extern double TabNormPyrLinear[MAX_SCALE];
//     extern double TabNormPyrFFT_Diff[MAX_SCALE];
//     extern double TabNormPyrFFT_Square[MAX_SCALE];
//     extern double TabNormPyrLaplacian[MAX_SCALE];
//     extern double TabNormMallat[MAX_SCALE];
//     extern double TabNormFeauveau[MAX_SCALE];
//     extern double TabNormPyrB3Spline[MAX_SCALE];
//     extern double TabNormPyrMedian5[MAX_SCALE];
//     extern double TabNormPyrMedian3[MAX_SCALE];
//     extern double TabNorm1 [MAX_SCALE];
//     extern double TabNormMinMax [MAX_SCALE];
//     extern double TabNormDiadicMallat [MAX_SCALE];
//     extern double TabNormSemiPyr [MAX_SCALE];  

const double VALISQRT2 = 1. / sqrt(2.);

double TabNormPaveLinear[MAX_SCALE] = 
    {0.800368,  0.272765, 0.119628, 0.0575768, 0.0292096, 0.0147667, 0.00780443,
     0., 0., 0.};  
     
double TabNormPaveB3Spline[MAX_SCALE] = 
    {0.889434,  0.200105, 0.0857724, 0.0413447, 0.0202689, 0.00995628, 0.00513504,
     0., 0., 0.}; 

double TabNormPaveFFT[MAX_SCALE] = 
    {0.968976,  0.112575, 0.049161, 0.023409, 0.010840, 0.005771, 0.002001,
     0., 0., 0.};            
     	     
double TabNormPaveMed5[MAX_SCALE] = 
    {0.989353,  0.208126, 0.0930072, 0.047539, 0.0240318, 0.0137446, 0.00771696,
     0.00393747, 0.0047435, 0.};    
     
double TabNormPaveMed3[MAX_SCALE] = 
    {0.970459,  0.337561, 0.176449, 0.0979767, 0.0536651, 0.0292211, 0.0170532,
     0.010303, 0.00732031, 0.};         
  
double TabNormPyrLinear[MAX_SCALE] = 
    {0.799733,  0.27278, 0.119974, 0.0581741, 0.028896, 0.0144647, 0.0075949,
     0., 0., 0.};

// double TabNormPyrUnserSpline[MAX_SCALE] = 
//    {0.867366, 0.951371, 0.629301, 0.3356, 0.173236, 0.0881835, 0.0443031,
//     0.022, 0.011, 0.}; 
     
double TabNormPyrFFT_Diff[MAX_SCALE] = 
    {0.967509,  0.113718, 0.0481781, 0.0232149, 0.011826, 0.00577454, 0.00296613,
     0., 0., 0.}; 

double TabNormPyrFFT_Square[MAX_SCALE] = 
    {0.984522,  0.135770, 0.057616, 0.028759, 0.015023, 0.007133, 0.001372,
     0., 0., 0.}; 
 
double TabNormPyrLaplacian[MAX_SCALE] = 
    {0.916072,  0.253218, 0.0712446, 0.021007, 0.00899726, 0.004002, 0.002001,
     0., 0., 0.}; 

double TabNormMallat[MAX_SCALE] = 
          {0.5,  1./4., 1./8., 1./16., 1./32., 1./64., 1./128., 
	   1./256., 1./512., 1./1024.}; 

double TabNormFeauveau[MAX_SCALE] = 
    {0.70,  0.501087, 0.353383, 0.251154, 0.178927, 0.126607, 0.087654,
     0.0608978, 0.0441607, 0.0302819}; 
     
double TabNormPyrB3Spline[MAX_SCALE] = 
    {0.88937,  0.200924, 0.08544, 0.0410743, 0.0208055, 0.0108577, 0.00528016,
     0., 0., 0.};     
     
double TabNormPyrMedian5[MAX_SCALE] = 
    {0.970291,  0.338748, 0.178373, 0.0992649, 0.0539517, 0.0317369, 0.0188998,
      0.013650, 0., 0.}; 
//This following table has been compute using the distrinution histogram
// of an image noise image, and an eps corresponding to Nsigma=3
//double TabNormPyrMedian5[MAX_SCALE] = 
//    {0.938473,  0.211919, 0.095457, 0.0493821, 0.0265218, 0.0141343, 0.0118281,
//      0.0118281, 0., 0.}; 

double TabNormPyrMedian3[MAX_SCALE] = 
    {0.999234,  0.438327, 0.2347, 0.132756, 0.0757592, 0.0466552, 0.0326198,
     0.0401325, 0., 0}; 

double TabNorm1[MAX_SCALE] = {1., 1., 1., 1., 1., 1., 1., 1., 1., 1.}; 

double TabNormMinMax[MAX_SCALE] = {1., 0.4, 0.28, 0.25, 0.25, 0., 0., 0., 0., 0.}; 
    
double TabNormDiadicMallat[MAX_SCALE] = {2.82494, 0.741742, 0.319126, 
                                         0.156064, 0.0811592, 0.0464773, 
					  0.0349472,  0.012, 0., 0.}; 

double TabNormGradDiadicMallat[MAX_SCALE] = {1.95818,0.498473,0.213932,
                                             0.107204,0.0590752,0.0386296,
                                             0.030227,0.0282632,0.012381,0.};

double TabNormSemiPyr[MAX_SCALE] = {0.889505, 0.201211 , 0.0876197,  0.043199, 
                                  0.0250364,  0.0123624, 0.00731319 , 
				  0.00646761, 0., 0.}; 

double TabNormPaveHaar[MAX_SCALE] = {0.866227, 0.433013, 0.216506,  0.108253, 
                                  0.0541266,  0.0270633, 0.0156250, 
				  0.00781250, 0.00390625, 0.00195312};

double TabNormDiadicHaar[MAX_SCALE] = {VALISQRT2,VALISQRT2/2,VALISQRT2/4,
                                       VALISQRT2/8,VALISQRT2/16,VALISQRT2/32,
				       VALISQRT2/64,VALISQRT2/128,VALISQRT2/256,
				        VALISQRT2/512};
				       
double TabNormUndecimatedMallat[MAX_SCALE] = {VALISQRT2,VALISQRT2/2,VALISQRT2/4,
                                       VALISQRT2/8,VALISQRT2/16,VALISQRT2/32,
				       VALISQRT2/64,VALISQRT2/128,VALISQRT2/256,
				        VALISQRT2/512};


static float CST_UNOT = 0.0437185;
double TabNormUndecNonOrth_B3SPLINE_HOR[MAX_SCALE] = {0.378322,0.100304 ,0.0437185,0.0211748,0.0105051,0.00524224, 
                                0.00261960, CST_UNOT / 32., CST_UNOT / 64., CST_UNOT / 128.};
     
double TabNormUndecNonOrth_B2SPLINE_HOR[MAX_SCALE] = {0.443045,0.166417,0.0755669,0.0367922,0.0182710,
                                                   0.00911981,0.00455783,0.00227830,0.00260405,
						   0.00260405};

static float CST_UNOTD = 0.0316653;
double TabNormUndecNonOrth_B3SPLINE_DIAG[MAX_SCALE] = {0.523438 ,0.0814819,0.0316653, 
   0.0149385, 0.00736354, 0.00366867, 0.00183254, 
   CST_UNOTD / 32., CST_UNOTD / 64., CST_UNOTD / 128.};
   
double TabNormUndecNonOrth_B2SPLINE_DIAG[MAX_SCALE] = {0.523438,0.161133,0.0679931,0.0324249,
                                                  0.0160160,0.00798335,0.00398832,0.00199306,
						  0.00199306,0.00199306};

double TabNormUndecNonOrth_B3SPLINE_HOR_2[MAX_SCALE] = {0.421425,0.121383,0.0541457,0.0263682,
                                                     0.0130990,0.00653883,0.00653883/2., 0.00653883/4.,
						     0.00653883/8., 0.00653883/16.};
double TabNormUndecNonOrth_B3SPLINE_DIAG_2[MAX_SCALE] = {0.649505,0.119327,0.0485716,0.0231647,
                                                     0.0114490,0.00570776,0.00570776/2., 0.00570776/4,
						     0.00570776/8., 0.00570776/16.};
double TabNormUndecNonOrth_B3SPLINE_ISOTROP_2[MAX_SCALE] = {1.09995, 0.272904, 0.119071, 0.0577063,
                                                     0.0286334, 0.0142893, 0.0142893/2., 0.0142893/4.,
						     0.0142893/16., 0.0142893/32};

/*******************************************************************/


set_transform SetTransform (type_transform Transform)
{
    set_transform Val_Return=TRANSF_PAVE;

    // if (test_type_transform(Transform, "SetTransform"))
    {
       switch (Transform)
       {
         case TO_LC:
         case TO_UNDECIMATED_NON_ORTHO:
	 case TO_UNDECIMATED_MALLAT:
	 case TC_FCT:
	 case TO_DIV_1:
	 case TO_DIV_2:
	     Val_Return = TRANSF_UNDECIMATED_MALLAT;
	     break;
         case TO_DIADIC_MALLAT:
	 case TO_DIADIC_HAAR:
		 //case TO_MAX_DIADIC_MALLAT:
	    Val_Return = TRANSF_DIADIC_MALLAT ; 
	    break;
         case TO_SEMI_PYR:
	 case TM_TO_SEMI_PYR:
	    Val_Return = TRANSF_SEMIPYR; 
	    break;
        case TO_PAVE_LINEAR:
        case TO_PAVE_BSPLINE:
        case TO_PAVE_FFT:
        case TM_PAVE_MEDIAN:
        case TM_PAVE_MINMAX:
        case TO_PAVE_FEAUVEAU:
	case TO_PAVE_HAAR:
               Val_Return = TRANSF_PAVE;
               break;
	case TO_PYR_MEYER:      
	case TO_PYR_MEYER_ISOTROP:
        case TO_PYR_LINEAR:
        case TO_PYR_BSPLINE:
        case TM_PYR_MEDIAN:
        case TM_PYR_LAPLACIAN:
        case TO_PYR_FFT_DIFF_RESOL:
        case TO_PYR_FFT_DIFF_SQUARE:
        case TM_PYR_MINMAX:
        case TM_PYR_SCALING_FUNCTION:
	case TM_TO_PYR:
               Val_Return = TRANSF_PYR;
               break;
        case TO_MALLAT:
	case TO_LIFTING:
        case TM_MIN_MAX:
        case TO_HAAR:
               Val_Return = TRANSF_MALLAT;
               break;
        case TO_FEAUVEAU:
               Val_Return = TRANSF_FEAUVEAU;
               break;
        default:
               cerr << "Error: bad parameter Transform" << endl;
	       exit(-1);
               break;
       }
    }
    return (Val_Return);
}
/************************************************************************/

double mr_scale_norm(int s, type_transform Transform, int MedianWindowSize, details which_detail,
                     type_undec_filter U_Filter)
{
    double Val=0.;
    if (s >= MAX_SCALE) return 1;
    else
    {
    switch (Transform)
    {
        case TO_PYR_MEYER:
	       if (s == 0) Val = 0.925;
	      else Val = 0.65;
	      break;
	case TO_PYR_MEYER_ISOTROP:
	      if (s == 0) Val = 0.94;
	      else Val = 0.58;
	      break;
        case TO_PAVE_LINEAR: Val = TabNormPaveLinear[s]; break;
        case TO_PAVE_BSPLINE: Val = TabNormPaveB3Spline[s]; break;
        case TO_PAVE_FFT: Val = TabNormPaveFFT[s]; break;
	case TO_PAVE_HAAR: Val = TabNormPaveHaar[s]; break;
        case TM_PAVE_MEDIAN:
           if (MedianWindowSize == 5) Val = TabNormPaveMed5[s];
	   else if (MedianWindowSize == 3) Val = TabNormPaveMed3[s];
           else 
           {
           cerr << "Error Bad window size for the median transform ..." << endl;
           exit(-1);
           }
           break;
        case TO_PYR_LINEAR: Val = TabNormPyrLinear[s]; break;
	case TM_TO_PYR:
        case TO_PYR_BSPLINE: Val = TabNormPyrB3Spline [s]; break;    
        case TO_PYR_FFT_DIFF_RESOL:Val = TabNormPyrFFT_Diff[s]; break;
        case TO_PYR_FFT_DIFF_SQUARE: Val = TabNormPyrFFT_Square[s]; break;
        case TM_PYR_LAPLACIAN: Val =  TabNormPyrLaplacian[s]; break;
        case TO_MALLAT: Val =  TabNormMallat [s]; break;
        case TO_FEAUVEAU:
        case TO_PAVE_FEAUVEAU: Val = TabNormFeauveau[s]; break;
	case TO_LIFTING:
        case TM_MIN_MAX:
	case TO_HAAR: Val =   TabNorm1 [s]; break;
        case TM_PYR_MEDIAN:  
	   if (MedianWindowSize == 5) Val = TabNormPyrMedian5[s];
	   else if (MedianWindowSize == 3) Val = TabNormPyrMedian3[s];
           else 
           {
           cerr << "Error Bad window size for the median transform ..." << endl;
           exit(-1);
           }
           break;
        case TM_PAVE_MINMAX:
        case TM_PYR_MINMAX:
        case TM_PYR_SCALING_FUNCTION: Val = TabNormMinMax[s]; break;
 	case  TO_DIADIC_MALLAT: Val = TabNormDiadicMallat [s]; break;
	case  TO_DIADIC_HAAR: Val = TabNormDiadicHaar[s]; break;
	case  TO_UNDECIMATED_MALLAT: Val = TabNormMallat[s]; break;
	case  TO_DIV_1: Val = 1.; break;
	case  TO_DIV_2: Val = 1.; break;
	case  TO_LC: Val = 1.; break;
	case  TC_FCT: Val = 1.; break;
	case  TO_UNDECIMATED_NON_ORTHO:
	       switch  (U_Filter)
	       {
	         case  U_B2SPLINE:
 	             if (which_detail == D_DIAGONAL) Val = TabNormUndecNonOrth_B2SPLINE_DIAG[s];
	             else Val = TabNormUndecNonOrth_B2SPLINE_HOR[s]; 
	             break;
	         case  U_B3SPLINE:
   	             if (which_detail == D_DIAGONAL) Val = TabNormUndecNonOrth_B3SPLINE_DIAG[s];
	             else Val = TabNormUndecNonOrth_B3SPLINE_HOR[s];
		     break; 
	          case U_B3SPLINE_2:
		     if (which_detail == D_DIAGONAL) Val = TabNormUndecNonOrth_B3SPLINE_DIAG_2[s];
	             else Val = TabNormUndecNonOrth_B3SPLINE_HOR_2[s]; 
		      break; 
	           case U_HAAR_B3S:
		       Val = 1./ pow((double) 2., (double)(s+1));
 		      break; 
		   case U_HAAR_B3S_POS:
		       if (which_detail == D_DIAGONAL) Val = 1.06 / (pow((double) 2., (double)(s+1)) * sqrt(2.));
	               else  Val = 1.22 / (pow((double) 2., (double)(s+1)) * sqrt(2.));
 		      break; 
	       } 
	       break;
 	case  TO_SEMI_PYR: Val = TabNormSemiPyr[s]; break;
	case TM_TO_SEMI_PYR: Val = TabNormSemiPyr[s]; break;
         default:
               fprintf (stderr, "scale_norm: noise_compute => Not implemented\n");
               exit(-1);
	       break;  
    }}
    return Val;
}
/************************************************************************/

double MultiResol::scale_norm(int s, details which_detail)
{
   if ((Type_Transform == TO_LC) || (Type_Transform == TC_FCT))  return 1.;
   else if ((Type_Transform == TO_DIV_1) || (Type_Transform == TO_DIV_2))  return 1.;
   else  if ((Type_Transform == TO_MALLAT) && (TypeNorm == NORM_L2)) return 1.;
   else if ((Type_Transform == TO_DIADIC_HAAR) && (TypeNorm == NORM_L2)) return 1.;
   else if ((Type_Transform == TO_UNDECIMATED_MALLAT) && (TypeNorm == NORM_L2)) return 1.;
   else return mr_scale_norm(s, Type_Transform, MedianWindowSize, which_detail, U_Filter);
} 

/************************************************************************/

double MultiResol::band_norm(int b)
{
   details Detail;
   int s;
   band_to_scale(b,s, Detail);
   if ((Type_Transform == TO_LC) || (Type_Transform == TC_FCT)) return 1.;
   else if ((Type_Transform == TO_DIV_1) || (Type_Transform == TO_DIV_2))  return 1.;
   else if ((Type_Transform == TO_MALLAT) && (TypeNorm == NORM_L2)) return 1.;
   else if ((Type_Transform == TO_DIADIC_HAAR) && (TypeNorm == NORM_L2)) return 1.;
   else if ((Type_Transform == TO_UNDECIMATED_MALLAT) && (TypeNorm == NORM_L2)) return 1.;
   // else if ((Type_Transform == TO_UNDECIMATED_NON_ORTHO) && (TypeNorm == NORM_L2)) return 1.;
   else return mr_scale_norm(s,  Type_Transform, MedianWindowSize, Detail, U_Filter);
}

/************************************************************************/

double MultiResol::gauss_detect_level(int b, float Sigma, float NSigma)
{ 
   double Val;
   details Detail;
   int s;
   band_to_scale(b,s, Detail);   
   if (s == 0) Val = (NSigma+1) * scale_norm(s,Detail) * Sigma;
   else Val = NSigma *  scale_norm(s) * Sigma;
   
   if ( (Type_Transform == TM_PAVE_MINMAX) 
          || (Type_Transform == TM_PYR_MINMAX)
          || (Type_Transform == TM_PYR_SCALING_FUNCTION))  Val *= 2.;
   return Val;
}

/************************************************************************/

void MultiResol::threshold(int b, float Level, Bool UseAbsVal)
{
   int i,j;
   
   for (i = 0; i < size_band_nl(b); i++)
   for (j = 0; j < size_band_nc(b); j++) 
   {
      if (UseAbsVal == True)
      {
         if (ABS((*this)(b,i,j)) < Level)  (*this)(b,i,j) = 0.;
      }
      else if ((*this)(b,i,j) < Level)  (*this)(b,i,j) = 0.;
   }
}

/************************************************************************/

void MultiResol::threshold(float Sigma, float NSgima, Bool UseAbsVal)
{
   int b;
   float Level;
   
   for (b=0; b < NbrBand-1; b ++)
   {
       Level = (float) gauss_detect_level(b, Sigma, NSgima);
       threshold(b, Level, UseAbsVal);
   }
}

/************************************************************************/

void MultiResol::pos_coeff(int NumCoef, int &s, int &i, int &j, details & which_detail)
{
   int b=0;
   int PosCoef = NumCoef;
   
   while (PosCoef >= nbr_coeff_in_band(b)) 
   {
      PosCoef -= nbr_coeff_in_band(b);
      b++;
   }   
   int Ncb = size_band_nc(b);
   i = PosCoef / Ncb;
   j = PosCoef - i * Ncb;
   band_to_scale(b, s, which_detail);
}

/************************************************************************/

float & MultiResol::operator() (int NumCoef) 
{
   int i,j;
   int b=0;
   int PosCoef = NumCoef;
   
   while (PosCoef >= nbr_coeff_in_band(b)) 
   {
      PosCoef -= nbr_coeff_in_band(b);
      b++;
   }   
   int Ncb = size_band_nc(b);
   i = PosCoef / Ncb;
   j = PosCoef - i * Ncb;
   return TabBand[b](i,j);
}

/************************************************************************/

float & MultiResol::operator() (int b, int NumCoef)  
{
   int Ncb = size_band_nc(b);
   int i = NumCoef / Ncb;
   int j = NumCoef - i*Ncb;
   return (*this)(b,i,j);
}

/************************************************************************/

float MultiResol::operator() (int b, int NumCoef, type_border Bord)
{
   int Ncb = size_band_nc(b);
   int i = NumCoef / Ncb;
   int j = NumCoef -  i*Ncb;
   return (*this)(b,i,j,Bord);
}

/************************************************************************/

float MultiResol::operator() (int b, int i, int j, type_border Bord) const
{
   float Val=-1;
#ifdef TESTIND
   test_band(b);
#endif
   if ((i < 0) || (i >= size_band_nl(b)) || (j < 0) || (j >= size_band_nc(b))) 
   {
      Val = TabBand[b](i,j,Bord);
   }
   else Val = TabBand[b](i,j);
   
   return Val;   
}

/************************************************************************/

float MultiResol::operator() (int s, int i, int j, details which_detail, type_border Bord)  
{
   int b = scale_to_band(s, which_detail); 
   return (*this)(b,i,j,Bord);
}

/************************************************************************/

float & MultiResol::operator() (int b, int i, int j) const
{
#ifdef TESTIND
   test_band(b);
#endif
   return TabBand[b](i,j);
}

/************************************************************************/

float & MultiResol::operator() (int s,int i,int j, details which_detail)  
{
   int b = scale_to_band(s, which_detail); 
   return TabBand[b](i,j);
}

/****************************************************************************/

Ifloat & MultiResol::band (int b) const 
{
   test_band(b);
   return (TabBand[b]);
}

/****************************************************************************/

void band2scale(const int NumBand, type_transform Type_Transform, int NbrBand,
                   int & s, details &which_detail)
{
  switch (SetTransform(Type_Transform))
   {       
        case TRANSF_SEMIPYR:
        case TRANSF_PAVE:
        case TRANSF_PYR:
               s = NumBand;
               which_detail = D_NULL;
               break;
        case TRANSF_FEAUVEAU:
               s = NumBand / 2;
               if (NumBand-s*2 == 0) which_detail = D_HALF_RESOL;
               else if (NumBand-s*2 == 1) which_detail = D_RESOL;
               if (NumBand >= NbrBand-1)
               {
                  which_detail = I_SMOOTH;
                  s =  (NbrBand-2)/2;
               }
               break;
	case TRANSF_DIADIC_MALLAT:
	       if (NumBand >= NbrBand-1)
               {
                  which_detail = I_SMOOTH;
                  s =  (NbrBand-1)/2 + 1;
               } 
	       else
	       {
	          s = NumBand  / 2;
		  if (NumBand-s*2 == 0) which_detail = D_HORIZONTAL;
		  else which_detail =  D_VERTICAL;
 	       }
 	      break;
        case TRANSF_MALLAT:
	case TRANSF_UNDECIMATED_MALLAT:
               s = NumBand / 3;
               if (NumBand-s*3 == 0) which_detail = D_HORIZONTAL;
               else if (NumBand-s*3 == 1) which_detail = D_VERTICAL;
               else if (NumBand-s*3 == 2) which_detail = D_DIAGONAL;
               if (NumBand >=  NbrBand-1)
               {
                   which_detail = I_SMOOTH;
                   s =  (NbrBand-2)/3;
               }
               break;
        default:
               fprintf (stderr,"Error: unknown transform\n");
               exit (-1);
               break;
   }  
}


/****************************************************************************/

int scale2band(int s, type_transform Type_Transform, 
                  int NbrBand, details which_detail)   
{   
   int NumBand = 0;
   switch (SetTransform(Type_Transform))
   {       
        case TRANSF_PAVE:
        case TRANSF_PYR:
	case TRANSF_SEMIPYR:
               NumBand=s;
               break;
        case TRANSF_FEAUVEAU:
               if (which_detail == I_SMOOTH) NumBand = NbrBand-1;
               else
               {
                  if (which_detail ==  D_HALF_RESOL) NumBand = 2 * s;
                  else if (which_detail == D_RESOL) NumBand = 2 * s + 1;
                  else
                  {
                      cerr << "Error: band type details is not correct ... " << endl;
                      exit(-1);
                  }
               }
               break;
	case TRANSF_DIADIC_MALLAT:
	       if ((which_detail == I_SMOOTH) 
	           || (s >= (NbrBand-1)/2+1 ))  NumBand = NbrBand-1;
               else
               {
                  if (which_detail == D_HORIZONTAL) NumBand  = 2*s;
                  else if (which_detail == D_VERTICAL) NumBand = 2*s+1;
                  else 
                  {
                      cerr << "Error: band type details is not correct ... " << endl;
                      exit(-1);
                  }
               }
	       break;
        case TRANSF_MALLAT:
	case TRANSF_UNDECIMATED_MALLAT:
               if (which_detail == I_SMOOTH) NumBand = NbrBand-1;
               else
               {
                  if (which_detail == D_HORIZONTAL) NumBand  = 3*s;
                  else if (which_detail == D_VERTICAL) NumBand = 3*s+1;
                  else if (which_detail == D_DIAGONAL ) NumBand = 3*s+2;
                  else 
                  {
                      cerr << "Error: band type details is not correct ... " << endl;
                      exit(-1);
                  }
               }
               break;
        default:
               fprintf (stderr,"Error: unknown transform\n");
               exit (-1);
               break;
   }  
   return NumBand;
}

 
/****************************************************************************/

int MultiResol::size_scale_nl (int ScaleNumber, details which_detail)   
{
   int Val_Return = Nl;
    
   if ((which_detail != I_SMOOTH) || (ScaleNumber >= Nbr_Plan-1))
   {
      int b = scale_to_band(ScaleNumber,which_detail);
      Val_Return = size_band_nl(b);
   }
   else  for (int i=0; i <= ScaleNumber; i++)  Val_Return  = (Val_Return+1)/2;
     
   return Val_Return;
}

/****************************************************************************/

int MultiResol::size_scale_nc (int ScaleNumber, details which_detail)   
{
   int Val_Return = Nc;
   if ((which_detail != I_SMOOTH) || (ScaleNumber >= Nbr_Plan-1))
   {
      int b = scale_to_band(ScaleNumber,which_detail);
      Val_Return = size_band_nc(b);
   }
   else  for (int i=0; i <= ScaleNumber; i++)  Val_Return  = (Val_Return+1)/2;
   return Val_Return;
}

/****************************************************************************/

Ifloat MultiResol::extract_scale(int s, details which_detail)
{
    int b = scale_to_band(s, which_detail); 
    int Nls = TabBand[b].nl();
    int Ncs = TabBand[b].nc();
    Ifloat *Image_Return = NULL;

    test_band(b);
    Image_Return = new Ifloat(Nls,Ncs, "extract");;
    *Image_Return = TabBand[b];
    return (*Image_Return);
}

/****************************************************************************/

void MultiResol::insert_scale(Ifloat & Image, int s, details which_detail)
{
    int b = scale_to_band(s, which_detail);  
    insert_band(Image, b);
}


/****************************************************************************/

Ifloat MultiResol::extract_band(int b)
{
    int Nls = TabBand[b].nl();
    int Ncs = TabBand[b].nc();
    Ifloat *Image_Return = NULL;

    test_band(b);
    Image_Return = new Ifloat(Nls,Ncs, "extract");;
    *Image_Return = TabBand[b];
    return (*Image_Return);
}

/****************************************************************************/

void MultiResol::insert_band(Ifloat & Image, int b)
{
      test_band(b);
     if ( (TabBand[b].nl() != Image.nl()) ||  (TabBand[b].nc() != Image.nc()))
     {
         cerr << "Error: unable to insert band ..." << endl;
         cerr << "       band and image haven't the same size" <<endl;
         exit(0);
     }
     TabBand[b] = Image;
}

/****************************************************************************/

MultiResol::MultiResol (int Nbr_Line, int Nbr_Col, int Nbr_Scale,
                        type_transform Transform, char *Name)
{
    init_param();
    (*this).alloc (Nbr_Line, Nbr_Col, Nbr_Scale, Transform, Name);
}

 
/****************************************************************************/

MultiResol::~MultiResol()
{
   (*this).free();
}

/****************************************************************************/

void MultiResol::init_param()
{
   extern type_format Format_Imag;

   Nl=0;
   Nc=0;
   Nbr_Plan=0;   
   Nbr_MrCoef = 0;
   NbrBand_per_Resol=0;
   NbrBand=0;
   Set_Transform=S_UNDEFINED;
   Type_Transform=T_UNDEFINED;
   Name_MR[0]='\0';
   ExactPyrRec = False;
   MedianWindowSize = DEFAULT_MEDIAN_WINDOW_SIZE;
   Border = DEFAULT_BORDER;
   FormatInputImag = Format_Imag;
   Nbr_Iter = DEFAULT_ITER_SCALING;
   Fc = DEFAULT_CUT_OFF_FREQ;
   TabBand = NULL;
   TabBandNl = NULL;
   TabBandNc = NULL;
   TMTO_NSigBSpline = DEFAULT_N_SIGMA;
   TMTO_NSigMedian = 5.;
   SigmaNoise = 0.;
   LiftingTrans = DEF_LIFT;
   SBFilter = DEF_SB_FILTER;
   TypeNorm = DEF_SB_NORM;
   NormHaar=4.;
   EdgeLineTransform=False;
   Verbose = False;
   FilterBank = NULL;
   FilterBank_Line = NULL;
   FilterBank_Column = NULL;
   FilterBankAlloc = False;
   NbrUndecimatedScale= -1;
   ModifiedATWT=False;
   UndecFilterBank=NULL;
   U_Filter = DEF_UNDER_FILTER;  
   MeyerWT = NULL;
   LC_NbrScaleCol = -1;
   LC_NbrScaleLine = -1;
   LC = NULL;
   FCT = NULL;
   FCT_Nbr_Angle = DEF_FCT_NDIR;       
   FCT_ExtendWT=False;
   FCT_IsotropWT=False;
   FCT_RealCur=True;
}

/****************************************************************************/

MultiResol::MultiResol ()
{
   init_param();
}

/****************************************************************************/

void MultiResol::alloc (int Nbr_Line, int Nbr_Col, int Nbr_Scale,
                        type_transform Transform, FilterAnaSynt *FAS, 
                        sb_type_norm Norm, int NbrUndec, type_undec_filter Undec_Filter)
{
  //cout << "IN " <<  NbrUndec << endl;
   char *Name = strdup("multi trans");
   FilterBank = FAS;
   if (FAS != NULL) SBFilter = FAS->type_filter();
   NbrUndecimatedScale = NbrUndec;
   TypeNorm = Norm;
   U_Filter = Undec_Filter;
   alloc (Nbr_Line, Nbr_Col, Nbr_Scale, Transform, Name);
}

/****************************************************************************/

void MultiResol::alloc (int Nbr_Line, int Nbr_Col, int Nbr_Scale,
                        type_transform Transform, char *Name)
{
    int s, Nl_s=Nbr_Line, Nc_s=Nbr_Col;
    char ch[80];
// cout << "ALLOC " << endl;
    
    Type_Transform = Transform;
    free();
    Set_Transform = SetTransform (Transform);
    Nl = Nbr_Line;
    Nc = Nbr_Col;
    if (Nbr_Scale <= 1) Nbr_Scale = get_nbr_scale( MIN(Nl,Nc) );
    Nbr_Plan = Nbr_Scale;
    strcpy(Name_MR, Name);
    NbrBand = Nbr_Plan;
    Nbr_MrCoef=0;
    NbrBand_per_Resol= number_band_per_resol(Transform);
    TabNbrBandPerResol.alloc(NbrBand);
    TabNbrBandPerResol.init(NbrBand_per_Resol);
    FilterBankAlloc=False;   
    switch (Transform)
    {
       case TO_DIV_1: 
           // cout << "TO_DIV_1" << endl;
	   if ((FilterBank_Line == NULL) || (FilterBank_Column == NULL))
	   {
  	      TypeNorm = DEF_SB_NORM;
              FilterBank_Line = new FilterAnaSynt;
              FilterBank_Line->alloc(F_5_3);
	      FilterBank_Column = new FilterAnaSynt;
              FilterBank_Column->alloc(F_4_4);
	      FilterBankAlloc=True;
	   }
	   // cout << "  ok TO_DIV_1" << endl;
           break;	
      case TO_DIV_2:
       	   // cout << "TO_DIV_2" << endl;
	   if ((FilterBank_Line == NULL) || (FilterBank_Column == NULL))
	   {
	      TypeNorm = DEF_SB_NORM;
              FilterBank_Line = new FilterAnaSynt;
              FilterBank_Line->alloc(F_4_4);
	      FilterBank_Column = new FilterAnaSynt;
              FilterBank_Column->alloc(F_5_3);
	      FilterBankAlloc=True;
	   }
	   // cout << "  ok TO_DIV_2" << endl;
           break;	
      case TO_MALLAT:
      case TO_UNDECIMATED_MALLAT:

      case TO_LC:

          if (FilterBank == NULL)
          {
	    TypeNorm = DEF_SB_NORM;
            SBFilter = DEF_SB_FILTER;
           // cout << "FB1 " << SBFilter << endl;
           // cout << StringSBFilter(SBFilter) << endl;
           FilterBank = new FilterAnaSynt;
	   	
           FilterBank->alloc(SBFilter);
	   
	   FilterBank_Line = FilterBank;
	   FilterBank_Column = FilterBank;
           // cout << "FB2 " << endl;
           FilterBankAlloc=True;
         }
         else
         {
           FilterBank_Line = FilterBank;
	   FilterBank_Column = FilterBank;
	  
          }
#if DEBUG_OBJ
          cout << "ALLOC 1 " << endl;

         cout << "NbrBand_per_Resol = " << NbrBand_per_Resol << endl;
         cout << "NbrBand_per_Resol = " << NbrBand_per_Resol << endl;
         cout << "MultiResol: create " << Name << endl;
         cout << "Type_Transform =  " << StringTransform(Type_Transform) << endl;
         if ((Transform == TO_MALLAT) || (Transform == TO_UNDECIMATED_MALLAT))
              cout << StringSBFilter(SBFilter) << endl;
         cout << "Nbr_Scale = " << Nbr_Scale << ", Nl = " << Nbr_Line << endl;
#endif
         if (Transform == TO_LC)
         {
           if (LC_NbrScaleCol  <= 1) LC_NbrScaleCol = get_nbr_scale(Nc);
           if (LC_NbrScaleLine <= 1) LC_NbrScaleLine = get_nbr_scale(Nl);
           // cout << LC_NbrScaleCol << " " << LC_NbrScaleLine << endl;
           Nbr_Plan = LC_NbrScaleCol;
           NbrBand = LC_NbrScaleCol * LC_NbrScaleLine;
           NbrBand_per_Resol= LC_NbrScaleLine;
           TabNbrBandPerResol.resize(NbrBand);
           TabNbrBandPerResol.init(NbrBand_per_Resol);
           FilterAnaSynt *FAS = filter_bank();
           if (FAS == NULL)
           {
             cout << "Error: filter bank is not defined ... " << endl;
             exit(-1);
           }
           SBF_LC = new SubBandFilter(*FAS, NORM_L2);
	   SBF_LC->Border = I_PERIOD;
	   // cout <<  " Nbr_Plan = "  << " " <<    LC_NbrScaleCol << endl;
           LC = new  LineCol;
	   (LC->alloc)(*SBF_LC, False);
	   //  cout <<  ":ALLOC OK "  << " " <<    NbrBand_per_Resol << endl;
	   if (Verbose == True) cout << "NbrScale in Y = " << LC_NbrScaleCol << ", NbrScale in X = "  << LC_NbrScaleLine << endl; 
          }  
     	  break;
      case  TO_UNDECIMATED_NON_ORTHO:
         // TypeNorm = NORM_L2;
        // UndecFilterBank = new UndecSubBandFilter(U_Filter, TypeNorm);
	  UndecFilterBank = new UndecSubBandFilter(U_Filter);
         break;
       case TO_PYR_MEYER_ISOTROP:
       case TO_PYR_MEYER:
           {
              Bool IsotropWT = (Transform == TO_PYR_MEYER_ISOTROP) ? True: False;
              MeyerWT = new MEYER_WT;
              MeyerWT->init(Nbr_Plan,Nl,Nc,False,IsotropWT);
           }
	   break;
      case TC_FCT:
          {
             FCT = new FCUR;
             FCT->alloc_from_fine(Nbr_Scale, Nbr_Line, Nbr_Col,  FCT_Nbr_Angle, FCT_ExtendWT, FCT_IsotropWT, FCT_RealCur);
            NbrBand = FCT->nbr_tot_band();
             TabNbrBandPerResol.resize(Nbr_Scale);
            // Ifloat X(Nbr_Line, Nbr_Col);
            // FCT->cur_trans(X); // IT IS THE ONLY WAY I FOUND TO GET THE BANDS SIZES !!!
            for (s = 0; s < Nbr_Scale; s++) TabNbrBandPerResol(s) = FCT->nbr_band(s);
            FCT_Band2Scale.alloc(NbrBand);       
            FCT_Band2Dir.alloc(NbrBand);   
            FCT_Scale2Band.alloc(FCT->nbr_scale());  
           // for (s = 0; s < Nbr_Scale; s++) cout << "ALLOC FCT " << Nbr_Scale << " " << TabNbrBandPerResol(s) << endl;    
    	   // cout << "ALLOC FCT END 1" <<  endl;   
          }
      default:
         break;
   }
   

    if (SetTransform(Type_Transform) == TRANSF_DIADIC_MALLAT)
    {
	 // we allocate one buffer more than the real need
	 // real number of band = (Nbr_Plan-1)*2 + 1    
	 // but we allocate more band more for implementation purpose    
	 NbrBand = Nbr_Plan*2;
    }
	else if ((SetTransform(Type_Transform) == TRANSF_UNDECIMATED_MALLAT) && (Type_Transform != TO_LC) && (Type_Transform != TC_FCT))
         NbrBand = (Nbr_Plan-1)*3 + 1;

    switch (Set_Transform)
    {
	case TRANSF_UNDECIMATED_MALLAT:
             {
	      // cout << "IN " << endl;
                 if ((NbrUndecimatedScale < 0) || (NbrUndecimatedScale > Nbr_Plan))
                      NbrUndecimatedScale = Nbr_Plan;
		 if  ((Type_Transform != TO_UNDECIMATED_NON_ORTHO) && (Type_Transform != TO_LC)&& (Type_Transform != TC_FCT))
		 {
  	     //  cout << "WT1D " << endl;
                     SubBandFilter WT1D(*FilterBank_Line, TypeNorm);
                       HALF_DECIMATED_2D_WT HDWT(WT1D);
                       NbrBand = HDWT.alloc(TabBand,Nl,Nc,Nbr_Plan,NbrUndecimatedScale);
 	     //  cout << "WT1D " << endl;
		 }
		 else
		 {
		    TabBand = new Ifloat [NbrBand];
		    if (Transform != TC_FCT)
		    {
		       for (s = 0; s < NbrBand; s++)
                       {
                        sprintf (ch, "band_%d", s+1);
                        TabBand[s].alloc (Nl, Nc, ch);
                        }
		    }
		    else 
		    {
		       int IndBand = 0;
		       // intarray TN1,TN2;
                       // int Nb = fct_real_get_band_size(Nbr_Scale, Nbr_Line, Nbr_Col, FCT_Nbr_Angle, TN1, TN2);

		       for (s = 0; s < FCT->nbr_scale(); s++)
		       {
		       // cout << "B" << endl;
		          FCT_Scale2Band(s) = IndBand;
		          for (int d = 0; d < FCT->nbr_band(s);  d++)
		          {
 		           sprintf (ch, "band_%d", IndBand+1);
			   int Ny = (FCT->size_band_nl)(s,d);
			   int Nx = (FCT->size_band_nc)(s,d);
                          // cout << "Band " << IndBand << " " << "(" << s+1 << ","  <<  d+1 << "), Size = " <<  Ny << " " << Nx <<  " " << TN1(IndBand) << " " << TN2(IndBand) << endl;

		           TabBand[IndBand].alloc(Ny, Nx, ch);
			   FCT_Band2Scale(IndBand) = s;
			   FCT_Band2Dir(IndBand) = d;
			   IndBand++;      
  		          }
		       }
		    }
		 }
		 
                 TabBandNl = new int [NbrBand];
                 TabBandNc = new int [NbrBand];
                 Nbr_MrCoef = 0;
                 for (s = 0; s < NbrBand; s++)
                 {
                   TabBandNl[s] = TabBand[s].nl();
                   TabBandNc[s] = TabBand[s].nc();
                   if (s != NbrBand-1) Nbr_MrCoef += TabBandNl[s]*TabBandNc[s];
                //  cout << "Band " << s+1 << " " << size_band_nl(s) << " " << size_band_nc(s) << endl;
                 }
             }
             break;
        case TRANSF_DIADIC_MALLAT:
	case TRANSF_PAVE:                      
               NbrUndecimatedScale = Nbr_Plan;
               TabBandNl = new int [NbrBand];
               TabBandNc = new int [NbrBand];
               TabBand = new Ifloat [NbrBand];
               Nbr_MrCoef = (NbrBand-1)*Nl*Nc;
               for (s = 0; s < NbrBand; s++)
               {
                   sprintf (ch, "band_%d", s+1);
                   TabBand[s].alloc (Nl, Nc, ch);
                   TabBandNl[s] = Nl;
                   TabBandNc[s] = Nc;
               }
	       // correct the number of band
	       if (SetTransform(Type_Transform) == TRANSF_DIADIC_MALLAT)
	                                                        NbrBand --;
               break;
	case TRANSF_SEMIPYR:
        case TRANSF_PYR:
               if (Set_Transform == TRANSF_SEMIPYR)
                                        NbrUndecimatedScale = 1;
               else NbrUndecimatedScale = 0;

               TabBand = new Ifloat [NbrBand];
               TabBandNl = new int [NbrBand];
               TabBandNc = new int [NbrBand];
               for (s = 0; s < NbrBand; s++)
               {
                   sprintf (ch, "scale_%d", s+1);
                   TabBand[s].alloc (Nl_s, Nc_s, ch);
                   if (s != NbrBand-1) Nbr_MrCoef += Nl_s*Nc_s;
                   TabBandNl[s] = Nl_s;
                   TabBandNc[s] = Nc_s;
		   if ((s != 0) || 
		        (SetTransform(Type_Transform) != TRANSF_SEMIPYR))
	           {
                      Nl_s = (Nl_s+1)/2;
                      Nc_s = (Nc_s+1)/2;
		   }
               }
               break;
        case TRANSF_MALLAT:
               NbrUndecimatedScale = 0;
               NbrBand =  NbrBand_per_Resol*(Nbr_Plan-1)+1;
             // cout << "Nbr_Plan = " << Nbr_Plan << " NbrBand = " << NbrBand << endl;
               TabBandNl = new int [NbrBand];
               TabBandNc = new int [NbrBand];
               TabBand = new Ifloat [NbrBand];
               for (s = 0; s < NbrBand-1; s+=3)
               {
                   TabBandNl[s] = (Nl_s+1)/2;
                   TabBandNc[s] = Nc_s/2;
                   sprintf (ch, "scale_%d", s+1);
                   TabBand[s].alloc ( TabBandNl[s], TabBandNc[s], ch);   
                   Nbr_MrCoef += TabBandNl[s]*TabBandNc[s];
                   
                   TabBandNl[s+1] = Nl_s/2;
                   TabBandNc[s+1] = (Nc_s+1)/2;
                   sprintf (ch, "scale_%d", s+2);
                   TabBand[s+1].alloc ( TabBandNl[s+1], TabBandNc[s+1], ch);
                   Nbr_MrCoef += TabBandNl[s+1]*TabBandNc[s+1];
                   
                   TabBandNl[s+2] = Nl_s/2;
                   TabBandNc[s+2] = Nc_s/2;
                   sprintf (ch, "scale_%d", s+3);
                   TabBand[s+2].alloc(TabBandNl[s+2],TabBandNc[s+2],ch);
                   Nbr_MrCoef += TabBandNl[s+2]*TabBandNc[s+2];
                   Nl_s = (Nl_s+1)/2;
                   Nc_s = (Nc_s+1)/2;
               }
               s = NbrBand-1;
               TabBandNl[s] = Nl_s;
               TabBandNc[s] = Nc_s;
               sprintf (ch, "scale_%d", s+1);
               TabBand[s].alloc(TabBandNl[s],TabBandNc[s],ch);
               //for (s = 0; s < NbrBand; s++)
               //    cout << " Band " << s <<  ": " << TabBandNl[s] << "x" << TabBandNc[s] << endl;
                break;
        case TRANSF_FEAUVEAU:               
               NbrUndecimatedScale = 0;
               NbrBand = NbrBand_per_Resol*(Nbr_Plan-1)+1;
               TabBandNl = new int [NbrBand];
               TabBandNc = new int [NbrBand];
	       TabBand = new Ifloat [NbrBand];
               for (s = 0; s < NbrBand-1; s+=2)
               {
                  int Nlw = Nl_s;
                  int Ncw = Nc_s/2;
                  TabBandNl[s] = Nlw;
                  TabBandNc[s] = Ncw;
                  Nbr_MrCoef += TabBandNl[s]*TabBandNc[s];
                  sprintf (ch, "scale_%d", s+1);
                  TabBand[s].alloc ( TabBandNl[s], TabBandNc[s], ch);   
                  Nlw /= 2;
                  Ncw = (Nc_s+1)/2;
                  TabBandNl[s+1] = Nlw;
                  TabBandNc[s+1] = Ncw;
                  Nbr_MrCoef += TabBandNl[s+1]*TabBandNc[s+1];
                  sprintf (ch, "scale_%d", s+2);
                  TabBand[s+1].alloc ( TabBandNl[s+1], TabBandNc[s+1], ch); 
                  Nl_s = (Nl_s+1)/2;
                  Nc_s = (Nc_s+1)/2;
               }
               s= NbrBand-1;
               TabBandNl[s] = Nl_s;
               TabBandNc[s] = Nc_s;
               sprintf (ch, "scale_%d", s+1);
               TabBand[s].alloc ( TabBandNl[s], TabBandNc[s], ch);                 
               break;
        default:
               fprintf (stderr, "Not implemented\n");
               exit (0);
               break;
    }
}
 
/****************************************************************************/
int MultiResol::computeNbBand(int Nbr_Line, int Nbr_Col,int Nbr_Scale,type_transform Transform,int NbrUndec)
{
  NbrUndecimatedScale = NbrUndec;
    char ch[80];
    int s, Nl_s=Nbr_Line, Nc_s=Nbr_Col;
    Nl = Nbr_Line;
    Nc = Nbr_Col;
    if (Nbr_Scale <= 1) Nbr_Scale = get_nbr_scale( MIN(Nl,Nc) );
    Type_Transform = Transform;
    Set_Transform = SetTransform (Transform);
    Nbr_Plan = Nbr_Scale;
    NbrBand = Nbr_Plan;
    Nbr_MrCoef=0;
    NbrBand_per_Resol= number_band_per_resol(Transform);
    
    
    switch (Transform)
    {
          case TO_MALLAT:
      case TO_UNDECIMATED_MALLAT:

      case TO_LC:

          if (FilterBank == NULL)
          {
	    TypeNorm = DEF_SB_NORM;
            SBFilter = DEF_SB_FILTER;
           FilterBank = new FilterAnaSynt;	   	
           FilterBank->alloc(SBFilter);	   
	   FilterBank_Line = FilterBank;
	   FilterBank_Column = FilterBank;
           FilterBankAlloc=True;
         }
         else
         {
           FilterBank_Line = FilterBank;
	   FilterBank_Column = FilterBank;
	  
          }
         if (Transform == TO_LC)
         {
           if (LC_NbrScaleCol  <= 1) LC_NbrScaleCol = get_nbr_scale(Nc);
           if (LC_NbrScaleLine <= 1) LC_NbrScaleLine = get_nbr_scale(Nl);
           Nbr_Plan = LC_NbrScaleCol;
           NbrBand = LC_NbrScaleCol * LC_NbrScaleLine;
          }  
	    
      case TC_FCT:
          {
             FCT = new FCUR;
             FCT->alloc_from_fine(Nbr_Scale, Nbr_Line, Nbr_Col,  FCT_Nbr_Angle, FCT_ExtendWT, FCT_IsotropWT, FCT_RealCur);
            NbrBand = FCT->nbr_tot_band();
          }
      default:
         break;
   }
   

    if (SetTransform(Type_Transform) == TRANSF_DIADIC_MALLAT)
    {
	 // we allocate one buffer more than the real need
	 // real number of band = (Nbr_Plan-1)*2 + 1    
	 // but we allocate more band more for implementation purpose    
	 NbrBand = Nbr_Plan*2;
    }
	else if ((SetTransform(Type_Transform) == TRANSF_UNDECIMATED_MALLAT) && (Type_Transform != TO_LC) && (Type_Transform != TC_FCT))
         NbrBand = (Nbr_Plan-1)*3 + 1;


    switch (Set_Transform)
    {
	case TRANSF_UNDECIMATED_MALLAT:
             {
	       
	      // cout << "IN " << endl;
                 if ((NbrUndecimatedScale < 0) || (NbrUndecimatedScale > Nbr_Plan))
                      NbrUndecimatedScale = Nbr_Plan;
		 if  ((Type_Transform != TO_UNDECIMATED_NON_ORTHO) && (Type_Transform != TO_LC)&& (Type_Transform != TC_FCT))
		 {
  	     //  cout << "WT1D " << endl;
		  
                     SubBandFilter WT1D(*FilterBank_Line, TypeNorm);
		      
                       HALF_DECIMATED_2D_WT HDWT(WT1D);
                       NbrBand = HDWT.alloc(TabBand,Nl,Nc,Nbr_Plan,NbrUndecimatedScale);
		       
 	     //  cout << "WT1D " << endl;
		 }
		
		 
                 
             }
             break;
        case TRANSF_DIADIC_MALLAT:
	case TRANSF_PAVE:                      
	       // correct the number of band
	       if (SetTransform(Type_Transform) == TRANSF_DIADIC_MALLAT)
	                                                        NbrBand --;
               break;

        case TRANSF_MALLAT:
              
               NbrBand =  NbrBand_per_Resol*(Nbr_Plan-1)+1;
                break;
        case TRANSF_FEAUVEAU:               
               NbrBand = NbrBand_per_Resol*(Nbr_Plan-1)+1;        
               break;
        default:
               fprintf (stderr, "Not implemented\n");
               exit (0);
               break;
    }
    return NbrBand;
}
void MultiResol::print_info()
{
    cout << "Transform = " << StringTransform(Type_Transform) << endl;
    if ((Type_Transform  == TO_MALLAT) || (Type_Transform  == TO_UNDECIMATED_MALLAT))
    {
       cout << StringSBFilter(SBFilter) << endl;
       if (TypeNorm  == NORM_L2) cout << "L2 normalization" << endl;
    }
    if ( Type_Transform == TO_LIFTING)
       cout << StringLSTransform(LiftingTrans) << endl;
    cout << "Number of scales = " << Nbr_Plan << endl;
    cout << "Input image size: " << Nl << " " << Nc << endl;
    cout << "Number of bands per resolution: " << NbrBand_per_Resol << endl;
    cout << "Number of undecimated scales = " << NbrUndecimatedScale << endl;

}

/****************************************************************************/

void MultiResol::free ()
{
#if DEBUG_OBJ
    if (Set_Transform!=S_UNDEFINED)
                cout << "~MultiResol: delete " << Name_MR << endl;
    else cout << "free of undefined object : " << Name_MR << endl;
#endif
    
    if ((NbrBand != 0) && (TabBand != NULL)){
	NbrBand = 0;
	delete [] TabBand;
	TabBand == NULL;
    }
    
    if (TabBandNl != NULL) {
       delete[] TabBandNl;
       TabBandNl=NULL;
    }
    if (TabBandNc != NULL){
        delete[] TabBandNc;
	TabBandNc=NULL;
     }
    if ((FilterBankAlloc == True) && (FilterBank != NULL)) 
    {
        delete FilterBank;
        FilterBank = NULL;
    }
    
    if((FilterBankAlloc == True) && (Type_Transform == TO_DIV_1 || Type_Transform == TO_DIV_2)){
      
       if(FilterBank_Line != NULL){
	   delete FilterBank_Line;
	   FilterBank_Line = NULL;
       }
    
       if(FilterBank_Column != NULL){
	   delete FilterBank_Column;
	   FilterBank_Column = NULL;
       }
    }
   
   
    if (UndecFilterBank != NULL) 
    {
       delete UndecFilterBank;
       UndecFilterBank = NULL;
    }    
    
    if ( ((Type_Transform == TO_PYR_MEYER_ISOTROP) || (Type_Transform == TO_PYR_MEYER))
         && (MeyerWT != NULL))  {
	   delete MeyerWT;
	   MeyerWT =NULL; 
    }
    if (Type_Transform == TC_FCT) {
        delete FCT;
	FCT=NULL;
     }
    
    if( Type_Transform == TO_LC){
      if(LC != NULL){
	delete LC;
	LC = NULL;
      }
      
      if(SBF_LC != NULL){
	delete SBF_LC;
	SBF_LC = NULL;
      }
    }

//    init_param();
#if DEBUG_OBJ
cout << " END MR " <<  endl;
#endif
}

/****************************************************************************/

static void norm_abs_imag(Ifloat &Imag, float NormVal, Bool Abs)
{
   int i,j;
   float Val,Max;
   
   Max = Imag(0,0);
   for (i = 0; i < Imag.nl(); i++)
   for (j = 0; j < Imag.nc(); j++)
   {
       Val = ABS(Imag(i,j));
       if (Max < Val) Max = Val;
   }
   for (i = 0; i < Imag.nl(); i++)
   for (j = 0; j < Imag.nc(); j++)
   {
       if (Abs == True)
               Imag(i,j) = ABS(Imag(i,j)) / Max * NormVal;
       else Imag(i,j) = Imag(i,j) / Max * NormVal;
   }   
}

/****************************************************************************/

void MultiResol::norm(float NormVal, Bool Abs)
{
  for (int s = 0; s < NbrBand; s++) norm_abs_imag(band(s),  NormVal,  Abs);
}

/****************************************************************************/

void MultiResol::compute_mod_phase ( Ifloat *&prpo_Modulus,
                                     Ifloat *&prpo_Phase) {

	int s,i,j;

	switch (Type_Transform) {
	case TO_DIADIC_MALLAT :
	case TO_DIADIC_HAAR:
	//case TO_MAX_DIADIC_MALLAT :

        // warning : all scales must have same dimension
        //           dir x is at scale 2*i and y at 2*i+1

        //Compute the wavelet coefficients
        for (s=0; s<nbr_scale()-1; s++) { // each scale
            for (i=0; i<size_ima_nl(); i++) {          // each line
                for (j=0; j<size_ima_nl(); j++) {      // each point
	                prpo_Modulus[s](i,j) = 
                        sqrt(  (TabBand[2*s])(i,j) * (TabBand[2*s])(i,j) + 
                               (TabBand[2*s+1])(i,j) * (TabBand[2*s+1])(i,j));
	                ARG ( TabBand[2*s](i,j), TabBand[2*s+1](i,j),
                          prpo_Phase[s](i,j));
                }
            }
        }
        break;


	default:
        fprintf (stderr, "Not implemented\n");
        exit (0);
        break;
    }
}









// interpolate between two consecutive maxima modulus at row pi_Row
void Ortho_Proj_Operator (int pi_Scale,            // current scale 
			  int pi_Line,             // current line 
                          int pi_Row ,             // current rowe
                          int pi_Begining,         // ind of begining
	                  int pi_End,              // ind of end
		          MultiResol& pro_Mr2dData,       // Data IN
			  MultiResol& pro_Mr2dRecData) {  // Rec Data OUT


  int pi_NumScale = pi_Scale/2 +1;   // pi_NumScale goes from 1 to Max, 
                                     // and pi_Scale from 0 to max-1
  pi_End = pi_End - pi_Begining;
  int ai_Dec = pi_Begining;
  pi_Begining = 0;

  // coef of interpolate function
  double af_Step1 = pow((double)2., (double)(-pi_NumScale))*pi_Begining;
  double af_Step2 = pow((double)2., (double)(-pi_NumScale))*pi_End;
  double ad_Delta = exp(af_Step1) * exp(-af_Step2) 
                  - exp(af_Step2) * exp(-af_Step1);

  int ai_FirstIndLine  = (pi_Line != -1 ? pi_Line : pi_Begining+ai_Dec);
  int ai_SecondIndLine = (pi_Line != -1 ? pi_Line : pi_End+ai_Dec);
  int ai_FirstIndRow   = (pi_Row  != -1 ? pi_Row : pi_Begining+ai_Dec);
  int ai_SecondIndRow  = (pi_Row  != -1 ? pi_Row : pi_End+ai_Dec);

  double ad_DeltAlpha = 
                (pro_Mr2dData (pi_Scale, ai_FirstIndLine, ai_FirstIndRow) - 
                 pro_Mr2dRecData (pi_Scale, ai_FirstIndLine,ai_FirstIndRow))
                 * exp(-af_Step2) -  exp(-af_Step1) *
                (pro_Mr2dData (pi_Scale, ai_SecondIndLine, ai_SecondIndRow) - 
                 pro_Mr2dRecData(pi_Scale, ai_SecondIndLine, ai_SecondIndRow));


  double ad_DeltBeta  = 
                (pro_Mr2dData (pi_Scale, ai_SecondIndLine, ai_SecondIndRow) - 
                 pro_Mr2dRecData (pi_Scale, ai_SecondIndLine, ai_SecondIndRow))
                 * exp(af_Step1) - exp(af_Step2) *
                (pro_Mr2dData (pi_Scale, ai_FirstIndLine, ai_FirstIndRow) - 
                 pro_Mr2dRecData (pi_Scale, ai_FirstIndLine, ai_FirstIndRow));

  double ad_Alpha = ad_DeltAlpha / ad_Delta;
  double ad_Beta  = ad_DeltBeta  / ad_Delta;

  // interpolate between pi_Begining and pi_End 
  int ai_IndLine=0, ai_IndRow=0;
  for (int i=pi_Begining+1; i<pi_End; i++) {
      double ad_inter = exp (pow((double)2., (double)(-pi_NumScale))*i)
                        * ad_Alpha + ad_Beta *
                        exp (- pow((double)2., (double)(-pi_NumScale))*i);
      ai_IndLine = (pi_Line != -1 ? pi_Line : i+ai_Dec);
      ai_IndRow  = (pi_Row  != -1 ? pi_Row : i+ai_Dec);
      pro_Mr2dRecData (pi_Scale, ai_IndLine, ai_IndRow) += ad_inter;
  }
}






void interpolate  (intarray& pro_NumberMaxInLine,   // Number of max in line
                   intarray** pppo_MaxModInLine,    // loc of max mod on line
                   intarray& pro_NumberMaxInRow,    // Number of max in line
                   intarray** pppo_MaxModInRow,     // loc of max mod on line
                   MultiResol& pro_Mr2dData,        // Data IN
                   MultiResol& pro_Mr2dRecData) {   // Rec Data OUT



  int ai_Nl = pro_Mr2dData.size_ima_nl();   // number of line of image
  int ai_Nr = pro_Mr2dData.size_ima_nc();   // number of row of image
  int si_NbrPlan = pro_Mr2dData.nbr_scale();

  // x plane
  for (int s=0; s<si_NbrPlan-1; s+=2) {

    // get a line y=cste
    for (int i=0; i<ai_Nl; i++) {

      // interpolate between max of line i at scale s
      for (int l=0; l<pro_NumberMaxInLine(s,i)-1; l++){

	Ortho_Proj_Operator (2*s, i, -1, 
                             (**(pppo_MaxModInLine+s*ai_Nl+i))(l), 
                             (**(pppo_MaxModInLine+s*ai_Nl+i))(l+1),
                             pro_Mr2dData, pro_Mr2dRecData);
      }
    }

    // get a row
    for (int j=0; j<ai_Nr; j++) {

      // interpolate between max of line i at scale s
      for (int l=0; l<pro_NumberMaxInRow(s,j)-1; l++){

	Ortho_Proj_Operator (2*s+1, -1, j, 
                             (**(pppo_MaxModInRow+s*ai_Nr+j))(l), 
			     (**(pppo_MaxModInRow+s*ai_Nr+j))(l+1),
			     pro_Mr2dData, pro_Mr2dRecData);
      }
    }
  }
}


void init_last_scale (int pi_NbrPlan,                // number of plan
                      int pi_Nl,                     // number of line
                      int pi_Nr,                     // number of row
                      MultiResol& po_Mr2dData,       // Data IN
                      MultiResol& po_Mr2dRecData){   // Rec Data OUT
 
   for (int i=0; i<pi_Nl; i++)
     for (int j=0; j<pi_Nr; j++) 
      po_Mr2dRecData (2*(pi_NbrPlan-1),i,j)=po_Mr2dData(2*(pi_NbrPlan-1),i,j);
   //po_Mr2dRecData.band(2*(pi_NbrPlan-1)) = po_Mr2dData.band(2*(pi_NbrPlan-1));
}

void init_max (int pi_NbrPlan,                // number of plan
               int pi_Nl,                     // number of line
               int pi_Nr,                     // number of row
	       intarray& pro_NumberMaxInLine, // Number of max in line
               intarray** pppo_MaxModInLine,  // loc of max mod on line
               intarray& pro_NumberMaxInRow,  // Number of max in line
               intarray** pppo_MaxModInRow,   // loc of max mod on lin
               MultiResol& pro_Mr2dData,      // Data IN
               MultiResol& pro_Mr2dRecData) { // max modulus of Data IN

  
  for (int s=0; s<pi_NbrPlan-1; s++) {

    for (int i=0; i<pi_Nl; i++) 
      for (int l=0; l<pro_NumberMaxInLine(s,i); l++)
	if (pro_NumberMaxInLine(s,i) > 0)
	  pro_Mr2dRecData (2*s, i, (**(pppo_MaxModInLine+s*pi_Nl+i))(l))
	    = pro_Mr2dData(2*s, i, (**(pppo_MaxModInLine+s*pi_Nl+i))(l));
  
	
    for (int j=0; j<pi_Nr; j++)
      for (int l=0; l<pro_NumberMaxInRow(s,j); l++)
	if (pro_NumberMaxInRow(s,j) > 0)
	  pro_Mr2dRecData (2*s+1, (**(pppo_MaxModInRow+s*pi_Nr+j))(l), j)
	    = pro_Mr2dData(2*s+1, (**(pppo_MaxModInRow+s*pi_Nr+j))(l), j);
  }
}



void search_number_max (Ifloat*& prpo_ModData, 
                        intarray& pro_NumberMaxInLine,
                        intarray& pro_NumberMaxInRow,
                        int pi_NBrBand, Bool Verbose) {
  int s;
  for (s=0; s<pi_NBrBand-1; s++) 
    for (int i=0; i<prpo_ModData[s].nl(); i++)
      for (int j=0; j<prpo_ModData[s].nc(); j++) 
	if (prpo_ModData[s](i,j) < -EPS || prpo_ModData[s](i,j) > EPS) {
	  pro_NumberMaxInLine(s,i)++;
	  pro_NumberMaxInRow(s,j)++;
	}

  int ai_MaxLine=0, ai_MaxRow=0;
  for (s=0; s<pi_NBrBand-1; s++) {
    for (int i=0; i<prpo_ModData[s].nl(); i++)
      ai_MaxLine += pro_NumberMaxInLine(s,i);
    for (int j=0; j<prpo_ModData[s].nc(); j++)
      ai_MaxRow += pro_NumberMaxInRow(s,j);
    if (Verbose == True)
    {
       cout << "Number max at scale " << s << " on line : " << ai_MaxLine << endl;
       cout << "Number max at scale " << s << " on row  : " << ai_MaxRow << endl;
    }
  }
}


void search_max (Ifloat*& prpo_ModData, 
                 intarray& pro_NumberMaxInLine,
                 intarray& pro_NumberMaxInRow,
                 int pi_NBrBand) {
 
  //cout << "search_max" << endl;	 
  //cout << NOT_USED << endl;
  int s;
  for (s=0; s<pi_NBrBand-1; s++) 
    for (int i=0; i<prpo_ModData[s].nl(); i++)
      for (int j=0; j<prpo_ModData[s].nc(); j++) 
	    if (prpo_ModData[s](i,j) > - NOT_USED ) {
	      pro_NumberMaxInLine(s,i)++;
	      pro_NumberMaxInRow(s,j)++;
	}

  int ai_MaxLine=0, ai_MaxRow=0;
  for (s=0; s<pi_NBrBand-1; s++) {
    for (int i=0; i<prpo_ModData[s].nl(); i++)
      ai_MaxLine += pro_NumberMaxInLine(s,i);
    for (int j=0; j<prpo_ModData[s].nc(); j++)
      ai_MaxRow += pro_NumberMaxInRow(s,j);
    cout << "Number max at scale " << s << " on line : " << ai_MaxLine << endl;
    cout << "Number max at scale " << s << " on row  : " << ai_MaxRow << endl;
  }
}

// void compute_max_modulus (Ifloat*& prpo_ModData, 
//                           Ifloat*& prpo_PhaData, 
//                           Bool     pe_Verbose,
// 			  int      pi_NbrPlan) {
//   char atc_FileName[80];
//   for (int s=0; s<pi_NbrPlan-1;s++) {
//     EDGE ao_EdgeDetect (ED_CANNY);
//     ao_EdgeDetect.edge_from_1deriv (prpo_ModData[s], prpo_PhaData[s]);
//     if (pe_Verbose) {
//       sprintf (atc_FileName, "mr2d_mod_%d", s);      
//       io_write_ima_float (atc_FileName, prpo_ModData[s]);      
//       sprintf (atc_FileName, "mr2d_pha_%d", s);
//       io_write_ima_float (atc_FileName, prpo_PhaData[s]);
//     }
//   }
// }


void init_coord_max (Ifloat*&    prpo_ModData, 
                     int         pi_Nl, 
		     int         pi_Nr,
		     intarray&   pro_NumberMaxInLine, 
		     intarray**  pppo_MaxModInLine,
                     intarray&   pro_NumberMaxInRow, 
		     intarray**  pppo_MaxModInRow,
		     int         pi_NbrPlan) {
		     
  for (int s=0; s<pi_NbrPlan-1; s++) {
    for (int i=0; i<pi_Nl; i++) {
      if (pro_NumberMaxInLine(s,i) > 0) {
	(*(pppo_MaxModInLine+s*pi_Nl+i)) = new intarray (pro_NumberMaxInLine(s,i));
	int l=0;
	for (int j=0; j<pi_Nr; j++) 
	  if (prpo_ModData[s](i,j) < -EPS || prpo_ModData[s](i,j) > EPS)
	    (**(pppo_MaxModInLine+s*pi_Nl+i))(l++) = j;
      }
    }
    for (int j=0; j<pi_Nr; j++) {
      if (pro_NumberMaxInRow(s,j) > 0) {
	(*(pppo_MaxModInRow+s*pi_Nr+j)) = new intarray (pro_NumberMaxInRow(s,j));
	int l=0;
	for (int i=0; i<pi_Nl; i++) 
	  if (prpo_ModData[s](i,j) < -EPS || prpo_ModData[s](i,j) > EPS)
	    (**(pppo_MaxModInRow+s*pi_Nr+j))(l++) = i;
      }
    }
  } 
}

void init_multiresol_out (int         pi_Nl, 
                          int         pi_Nr, 
			  int         pi_NbrPlan, 
			  intarray&   pro_NumberMaxInLine, 
			  intarray**  pppo_MaxModInLine,
                          intarray&   pro_NumberMaxInRow,  
			  intarray**  pppo_MaxModInRow, 
			  MultiResol& pro_Mr2dData, 
			  MultiResol& pro_Mr2dRecData) {
			  
			  
  // init ao_Mr2dRecData
  for (int b=0; b<pro_Mr2dData.nbr_band(); b++) 
    for (int i=0; i<pi_Nl; i++) 
      for (int j=0; j<pi_Nr; j++)
        pro_Mr2dRecData (b,i,j) = - NOT_USED;


  // init max
  for (int s=0; s<pi_NbrPlan-1; s++) {

    for (int i=0; i<pi_Nl; i++) 
      for (int l=0; l<pro_NumberMaxInLine(s,i); l++)
	    if (pro_NumberMaxInLine(s,i) > 0)
	      pro_Mr2dRecData (2*s, i, (**(pppo_MaxModInLine+s*pi_Nl+i))(l))
	        = pro_Mr2dData(2*s, i, (**(pppo_MaxModInLine+s*pi_Nl+i))(l));

    for (int j=0; j<pi_Nr; j++)
      for (int l=0; l<pro_NumberMaxInRow(s,j); l++)
	    if (pro_NumberMaxInRow(s,j) > 0)
	      pro_Mr2dRecData (2*s+1, (**(pppo_MaxModInRow+s*pi_Nr+j))(l), j)
	        = pro_Mr2dData(2*s+1, (**(pppo_MaxModInRow+s*pi_Nr+j))(l), j);
  }

  // init last scale
  for (int i=0; i<pi_Nl; i++)
    for (int j=0; j<pi_Nr; j++) 
      pro_Mr2dRecData (2*(pi_NbrPlan-1),i,j)=pro_Mr2dData(2*(pi_NbrPlan-1),i,j);
			  
}			  
			  
			  
			  
			  
			  
			  
			  
			  
			  
			  
			  
			  
			  
			  
			  
			  
			  
			  
			  
			  
			  
			  	     
		     
