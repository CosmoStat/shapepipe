/******************************************************************************
**                   Copyright (C) 2000 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  10/02/00 
**    
**    File:  IM_Rot.h
**
**    Modification history :
**
******************************************************************************/

#ifndef _ROT_H_
#define _ROT_H_

#include"IM_Obj.h"
 
inline float REM(float Arg, float Denominator)
{
   int Val = (int)(Arg / Denominator);
   float Vali = (float) Val * Denominator;
   return (Arg - Vali);
}

void translate1d(fltarray & x, float Dx);
// 1d signal translation from right to left

void translate2d(Ifloat & Ima, float Dx, float Dy);
// 2d image translation from left to right and bottom to top



class BSPLINE_DEC {
   double InitialCausalCoefficient(double c[],long DataLength, double z, double	Tolerance);
				 // c[] = coefficients  
				 // DataLength = number of coefficients  
				 // z = actual pole  
				 // Tolerance = admissible relative error  
 
 
   double InitialAntiCausalCoefficient(double c[], long DataLength, double z)
				 // c[] = coefficients  
				 // DataLength = number of coefficients 
				 // z = actual pole  			 
   {  /* this initialization corresponds to mirror boundaries */
	return((z / (z * z - 1.0)) * (z * c[DataLength - 2L] + c[DataLength - 1L]));
   } 
   
   void ConvertToInterpolationCoefficients(double c[], long DataLength, 
                                         double z[], long NbPoles, double Tolerance);
				 // c[] = coefficients  
				 // DataLength = number of coefficients 					
				 // z[] = poles  
				 // NbPoles = number of poles  
				 // Tolerance = admissible relative error  



   void  GetColumn(float *Image, long Width, long x, double Line[], long Height);
		// Image = input image array 
		// Width = width of the image 
		// x=  x coordinate of the selected line  
		// Line[] = output linear array  
		// Height= length of the line  
   void GetRow (float *Image, long y, double Line[], long Width);
		// y=  y coordinate of the selected line */
		// Line[] =  output linear array  
		// long	Width = length of the line 
   void PutColumn (float *Image, long	Width, long  x,	 double	Line[],	long Height);
   void PutRow(float	*Image,	long y,	double	Line[],	long Width);
 
				 				 
 public:
  int SplineDegree;
  BSPLINE_DEC() {SplineDegree=3;RotationMasking=False;}
  int SamplesToCoefficients(float *Image, int Width, int Height);

  void SamplesToCoefficients(fltarray & Data)
  {
     int Err = SamplesToCoefficients(Data.buffer(), Data.nx(),  Data.ny());
     if (Err) 
     {
       printf("Error: Change of spline basis failed ...\n");
       exit(-1);
     }
  }
  void SamplesToCoefficients(Ifloat & Data)
  {
     int Err = SamplesToCoefficients(Data.buffer(), Data.nc(), Data.nl());
     if (Err) 
     {
       printf("Error: Change of spline basis failed ...\n");
       exit(-1);
     }	
  }
  
  double InterpolatedValue(fltarray &Data, double x,  double y)
  {  
      return InterpolatedValue(Data.buffer(),  Data.nx(), Data.ny(), x, y);
  }
  double InterpolatedValue(Ifloat &Data, double x,  double y)
  {  
      return InterpolatedValue(Data.buffer(),  Data.nc(), Data.nl(), x, y);
  }
  double InterpolatedValue(float *Bcoeff, int Width, int Height, double x, double y);
    		//  Bcoeff = input B-spline array of coefficients */
 		//  x = x coordinate where to interpolate  
		//  y = y coordinate where to interpolate  
		// SplineDegree = degree of the spline model  

  Bool RotationMasking;
  void shift_rotation_from_spline_coeff(float *Data, float *Result, 
                           int Nx, int Ny, double AngleParam,
                           double xShift, double yShift,  
			   double  xOrigin, double  yOrigin);

  void shift_rotation(Ifloat &Data, Ifloat & Result,  
                   double AngleParam, Bool DataSplineDecomposed = False,
		   double xShift=0, double yShift=0, 
                   int ParamxOrigin=-1, int ParamyOrigin=-1)
  {
    int Nx = Data.nc();
    int Ny = Data.nl();
    double xOrigin=(ParamxOrigin==-1) ? Nx/2: ParamxOrigin;
    double yOrigin=(ParamyOrigin==-1) ? Ny/2: ParamyOrigin;
    if ((Result.nl() != Ny) || (Result.nc() != Nx)) Result.resize(Ny,Nx);
    if (DataSplineDecomposed == False) SamplesToCoefficients(Data);
    
    shift_rotation_from_spline_coeff(Data.buffer(),  Result.buffer(),Nx, 
                           Ny, AngleParam, xShift, yShift,  
			   xOrigin, yOrigin);
  }
  void shift_rotation(fltarray &Data, fltarray & Result,  
                   double AngleParam, Bool DataSplineDecomposed = False,
		   double xShift=0, double yShift=0, 
                   int ParamxOrigin=-1, int ParamyOrigin=-1)
  {
    int Nx = Data.ny();
    int Ny = Data.nx();
    double xOrigin=(ParamxOrigin==-1) ? Nx/2: ParamxOrigin;
    double yOrigin=(ParamyOrigin==-1) ? Ny/2: ParamyOrigin;
    if ((Result.ny() != Ny) || (Result.nx() != Nx)) Result.resize(Nx,Ny);
    if (DataSplineDecomposed == False) SamplesToCoefficients(Data);
    
    shift_rotation_from_spline_coeff(Data.buffer(),  Result.buffer(),Nx, 
                           Ny, AngleParam, xShift, yShift,  
			   xOrigin, yOrigin);
  }  
};
void im_zoom(Ifloat &Ima1, Ifloat &Ima2, double ZoomX, double ZoomY);
// Ima2 = zoom(Ima1) with a zoom factor of ZoomX and ZoomY in X and Y direction
void im_zoom(Ifloat &Ima1, Ifloat &Ima2, double Zoom);
// Ima2 = zoom(Ima1) with a zoom factor of Zoom in both directions
void im_zoom(Ifloat &Ima1, Ifloat &Ima2);
// Ima2 = zoom(Ima1) with a zoom factor defined by the ratio between Ima2 and Ima1 size
// Ima2 must be allocated.

class Rotation {
   void im_shear(Ifloat &Data, char Dir, float beta); 
   public:
   BSPLINE_DEC SplineDec;
   Bool Verbose;
   Bool RadianUnitAngle; // if False, angle are given in degree

   Rotation() {Verbose=False;RadianUnitAngle=False;} 
   void im_rot90(Ifloat &DataIn, Ifloat & DataRot, int NRot);
               // rotate the image by a angle multiple of 90 degrees.
               // NRot = 1 ==> Angle = 90
               // NRot = 2 ==> Angle = 180
               // NRot = 3 ==> Angle = 270

   void im_rotate(Ifloat &Data, Ifloat &  DataRot, float Angle);
              // Rotate an image by the Three Pass Algorithm of Unser 1995
              // this rotation is reversible, because it does not 
              // introduce any averaging during the rotation process.
              // As it uses the FFT, image size must be power of 2.

  void im_move_rotate_bilinear(Ifloat &Data1, Ifloat & Data2, float AngleRot, 
                              float Dx=0, float Dy=0);
             // standard rotation method using a bilinear interpolation

  void im_move_rotate_spline(Ifloat &Data1, Ifloat & Data2, float AngleRot, 
                              float Dx=0, float Dy=0, int SlineDegree=3)
  {
     double Angle = (RadianUnitAngle == False) ? -AngleRot: -AngleRot / PI * 180.;			    
     SplineDec.SplineDegree = SlineDegree;
     SplineDec.shift_rotation(Data1, Data2, Angle, False, (double) Dx, (double) Dy);
  }
  void im_move_rotate(Ifloat &Data1, Ifloat & Data2, float AngleRot, 
                              float Dx=0, float Dy=0)
  {
      im_move_rotate_spline(Data1, Data2, AngleRot, (double) Dx, (double) Dy);
  }
  ~Rotation (){}
};


#endif
