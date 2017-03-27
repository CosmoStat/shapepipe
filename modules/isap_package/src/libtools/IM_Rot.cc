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
**    File:  IM_Rot.cc
**
**     
**
******************************************************************************/

#include "IM_Rot.h"
#include "IM_Obj.h"
#include "FFTN_1D.h"
#include "FFTN_2D.h"
using namespace std;

/**********************************************************************/
/*********************************************************************/

double BSPLINE_DEC::InitialCausalCoefficient(
					double	c[],		/* coefficients */
					long	DataLength,	/* number of coefficients */
					double	z,			/* actual pole */
					double	Tolerance	/* admissible relative error */
				)
{ /* begin InitialCausalCoefficient */
	double	Sum, zn, z2n, iz;
	long	n, Horizon;
	/* this initialization corresponds to mirror boundaries */
	Horizon = DataLength;
	if (Tolerance > 0.0) {
		Horizon = (long)ceil(log(Tolerance) / log(fabs(z)));
	}
	if (Horizon < DataLength) {
		/* accelerated loop */
		zn = z;
		Sum = c[0];
		for (n = 1L; n < Horizon; n++) {
			Sum += zn * c[n];
			zn *= z;
		}
		return(Sum);
	}
	else {
		/* full loop */
		zn = z;
		iz = 1.0 / z;
		z2n = pow(z, (double)(DataLength - 1L));
		Sum = c[0] + z2n * c[DataLength - 1L];
		z2n *= z2n * iz;
		for (n = 1L; n <= DataLength - 2L; n++) {
			Sum += (zn + z2n) * c[n];
			zn *= z;
			z2n *= iz;
		}
		return(Sum / (1.0 - zn * zn));
	}
} /* end InitialCausalCoefficient */



/*********************************************************************/

void BSPLINE_DEC::ConvertToInterpolationCoefficients
				(
					double	c[],		/* input samples --> output coefficients */
					long	DataLength,	/* number of samples or coefficients */
					double	z[],		/* poles */
					long	NbPoles,	/* number of poles */
					double	Tolerance	/* admissible relative error */
				)
{ /* begin ConvertToInterpolationCoefficients */
	double	Lambda = 1.0;
	long	n, k;
	/* special case required by mirror boundaries */
	if (DataLength == 1L) {
		return;
	}
	/* compute the overall gain */
	for (k = 0L; k < NbPoles; k++) {
		Lambda = Lambda * (1.0 - z[k]) * (1.0 - 1.0 / z[k]);
	}
	/* apply the gain */
	for (n = 0L; n < DataLength; n++) {
		c[n] *= Lambda;
	}
	/* loop over all poles */
	for (k = 0L; k < NbPoles; k++) {
		/* causal initialization */
		c[0] = InitialCausalCoefficient(c, DataLength, z[k], Tolerance);
		/* causal recursion */
		for (n = 1L; n < DataLength; n++) {
			c[n] += z[k] * c[n - 1L];
		}
		/* anticausal initialization */
		c[DataLength - 1L] = InitialAntiCausalCoefficient(c, DataLength, z[k]);
		/* anticausal recursion */
		for (n = DataLength - 2L; 0 <= n; n--) {
			c[n] = z[k] * (c[n + 1L] - c[n]);
		}
	}
} /* end ConvertToInterpolationCoefficients */


/*********************************************************************/

void BSPLINE_DEC::GetColumn(float	*Image,		/* input image array */
			long	Width,		/* width of the image */
			long	x,			/* x coordinate of the selected line */
			double	Line[],		/* output linear array */
			long	Height		/* length of the line */
		       )
{ 
	long	y;
	Image = Image + (ptrdiff_t)x;
	for (y = 0L; y < Height; y++) {
		Line[y] = (double)*Image;
		Image += (ptrdiff_t)Width;
	}
}  
/*********************************************************************/

void BSPLINE_DEC::GetRow (float	*Image,		/* input image array */
		    long	y,			/* y coordinate of the selected line */
		    double	Line[],		/* output linear array */
		    long	Width		/* length of the line */
		   )
{  
	long	x;
	Image = Image + (ptrdiff_t)(y * Width);
	for (x = 0L; x < Width; x++) {
		Line[x] = (double)*Image++;
	}
}  
/*********************************************************************/


void BSPLINE_DEC::PutColumn (float	*Image,		/* output image array */
					long	Width,		/* width of the image */
					long	x,			/* x coordinate of the selected line */
					double	Line[],		/* input linear array */
					long	Height		/* length of the line and height of the image */
				)
{ 
	long	y;
	Image = Image + (ptrdiff_t)x;
	for (y = 0L; y < Height; y++) {
		*Image = (float)Line[y];
		Image += (ptrdiff_t)Width;
	}
}

/*********************************************************************/

void BSPLINE_DEC::PutRow(
					float	*Image,		/* output image array */
					long	y,			/* y coordinate of the selected line */
					double	Line[],		/* input linear array */
					long	Width		/* length of the line and width of the image */
				)
{ 
	long	x;
	Image = Image + (ptrdiff_t)(y * Width);
	for (x = 0L; x < Width; x++) {
		*Image++ = (float)Line[x];
	}
}  

/*********************************************************************/

int BSPLINE_DEC::SamplesToCoefficients(float *Image, int Width, int Height)

{ /* begin SamplesToCoefficients */
   double	*Line;
   double	Pole[2];
   long	NbPoles;
   long	x, y;
   double DBL_Epsilon = 0.;
	/* recover the poles from a lookup table */
   switch (SplineDegree) {
		case 2L:
			NbPoles = 1L;
			Pole[0] = sqrt(8.0) - 3.0;
			break;
		case 3L:
			NbPoles = 1L;
			Pole[0] = sqrt(3.0) - 2.0;
			break;
		case 4L:
			NbPoles = 2L;
			Pole[0] = sqrt(664.0 - sqrt(438976.0)) + sqrt(304.0) - 19.0;
			Pole[1] = sqrt(664.0 + sqrt(438976.0)) - sqrt(304.0) - 19.0;
			break;
		case 5L:
			NbPoles = 2L;
			Pole[0] = sqrt(135.0 / 2.0 - sqrt(17745.0 / 4.0)) + sqrt(105.0 / 4.0)
				- 13.0 / 2.0;
			Pole[1] = sqrt(135.0 / 2.0 + sqrt(17745.0 / 4.0)) - sqrt(105.0 / 4.0)
				- 13.0 / 2.0;
			break;
		default:
			printf("Invalid spline degree\n");
			return(1);
	}
	/* convert the image samples into interpolation coefficients */
	/* in-place separable process, along x */
	Line = (double *)malloc((size_t)(Width * (long)sizeof(double)));
	if (Line == (double *)NULL) {
		printf("Row allocation failed\n");
		return(1);
	}
	for (y = 0L; y < Height; y++) {
		GetRow(Image, y, Line, Width);
		ConvertToInterpolationCoefficients(Line, Width, Pole, NbPoles, DBL_Epsilon);
		PutRow(Image, y, Line, Width);
	}
	free(Line);
	/* in-place separable process, along y */
	Line = (double *)malloc((size_t)(Height * (long)sizeof(double)));
	if (Line == (double *)NULL) {
		printf("Column allocation failed\n");
		return(1);
	}
	for (x = 0L; x < Width; x++) {
		GetColumn(Image, Width, x, Line, Height);
		ConvertToInterpolationCoefficients(Line, Height, Pole, NbPoles, DBL_Epsilon);
		PutColumn(Image, Width, x, Line, Height);
	}
	free(Line);
	return(0);
} /* end SamplesToCoefficients */

/*********************************************************************/

double BSPLINE_DEC::InterpolatedValue(float *Bcoeff,  /* input B-spline array of coefficients */
 			  int Width, 
			  int Height,
			  double	x,			/* x coordinate where to interpolate */
			  double	y			/* y coordinate where to interpolate */)
{ /* begin InterpolatedValue */
 	float	*p;
	double	xWeight[6], yWeight[6];
	double	interpolated;
	double	w, w2, w4, t, t0, t1;
	long	xIndex[6], yIndex[6];
	long	Width2 = 2L * Width - 2L, Height2 = 2L * Height - 2L;
	long	i, j, k;
	/* compute the interpolation indexes */
	if (SplineDegree & 1L) {
		i = (long)floor(x) - SplineDegree / 2L;
		j = (long)floor(y) - SplineDegree / 2L;
		for (k = 0L; k <= SplineDegree; k++) {
			xIndex[k] = i++;
			yIndex[k] = j++;
		}
	}
	else {
		i = (long)floor(x + 0.5) - SplineDegree / 2L;
		j = (long)floor(y + 0.5) - SplineDegree / 2L;
		for (k = 0L; k <= SplineDegree; k++) {
			xIndex[k] = i++;
			yIndex[k] = j++;
		}
	}
	/* compute the interpolation weights */
	switch (SplineDegree) {
		case 2L:
			/* x */
			w = x - (double)xIndex[1];
			xWeight[1] = 3.0 / 4.0 - w * w;
			xWeight[2] = (1.0 / 2.0) * (w - xWeight[1] + 1.0);
			xWeight[0] = 1.0 - xWeight[1] - xWeight[2];
			/* y */
			w = y - (double)yIndex[1];
			yWeight[1] = 3.0 / 4.0 - w * w;
			yWeight[2] = (1.0 / 2.0) * (w - yWeight[1] + 1.0);
			yWeight[0] = 1.0 - yWeight[1] - yWeight[2];
			break;
		case 3L:
			/* x */
			w = x - (double)xIndex[1];
			xWeight[3] = (1.0 / 6.0) * w * w * w;
			xWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) - xWeight[3];
			xWeight[2] = w + xWeight[0] - 2.0 * xWeight[3];
			xWeight[1] = 1.0 - xWeight[0] - xWeight[2] - xWeight[3];
			/* y */
			w = y - (double)yIndex[1];
			yWeight[3] = (1.0 / 6.0) * w * w * w;
			yWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) - yWeight[3];
			yWeight[2] = w + yWeight[0] - 2.0 * yWeight[3];
			yWeight[1] = 1.0 - yWeight[0] - yWeight[2] - yWeight[3];
			break;
		case 4L:
			/* x */
			w = x - (double)xIndex[2];
			w2 = w * w;
			t = (1.0 / 6.0) * w2;
			xWeight[0] = 1.0 / 2.0 - w;
			xWeight[0] *= xWeight[0];
			xWeight[0] *= (1.0 / 24.0) * xWeight[0];
			t0 = w * (t - 11.0 / 24.0);
			t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
			xWeight[1] = t1 + t0;
			xWeight[3] = t1 - t0;
			xWeight[4] = xWeight[0] + t0 + (1.0 / 2.0) * w;
			xWeight[2] = 1.0 - xWeight[0] - xWeight[1] - xWeight[3] - xWeight[4];
			/* y */
			w = y - (double)yIndex[2];
			w2 = w * w;
			t = (1.0 / 6.0) * w2;
			yWeight[0] = 1.0 / 2.0 - w;
			yWeight[0] *= yWeight[0];
			yWeight[0] *= (1.0 / 24.0) * yWeight[0];
			t0 = w * (t - 11.0 / 24.0);
			t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
			yWeight[1] = t1 + t0;
			yWeight[3] = t1 - t0;
			yWeight[4] = yWeight[0] + t0 + (1.0 / 2.0) * w;
			yWeight[2] = 1.0 - yWeight[0] - yWeight[1] - yWeight[3] - yWeight[4];
			break;
		case 5L:
			/* x */
			w = x - (double)xIndex[2];
			w2 = w * w;
			xWeight[5] = (1.0 / 120.0) * w * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			xWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) - xWeight[5];
			t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
			t1 = (-1.0 / 12.0) * w * (t + 4.0);
			xWeight[2] = t0 + t1;
			xWeight[3] = t0 - t1;
			t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
			t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
			xWeight[1] = t0 + t1;
			xWeight[4] = t0 - t1;
			/* y */
			w = y - (double)yIndex[2];
			w2 = w * w;
			yWeight[5] = (1.0 / 120.0) * w * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			yWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) - yWeight[5];
			t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
			t1 = (-1.0 / 12.0) * w * (t + 4.0);
			yWeight[2] = t0 + t1;
			yWeight[3] = t0 - t1;
			t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
			t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
			yWeight[1] = t0 + t1;
			yWeight[4] = t0 - t1;
			break;
		default:
			printf("Invalid spline degree\n");
			return(0.0);
	}
	/* apply the mirror boundary conditions */
	for (k = 0L; k <= SplineDegree; k++) {
		xIndex[k] = (Width == 1L) ? (0L) : ((xIndex[k] < 0L) ?
			(-xIndex[k] - Width2 * ((-xIndex[k]) / Width2))
			: (xIndex[k] - Width2 * (xIndex[k] / Width2)));
		if (Width <= xIndex[k]) {
			xIndex[k] = Width2 - xIndex[k];
		}
		yIndex[k] = (Height == 1L) ? (0L) : ((yIndex[k] < 0L) ?
			(-yIndex[k] - Height2 * ((-yIndex[k]) / Height2))
			: (yIndex[k] - Height2 * (yIndex[k] / Height2)));
		if (Height <= yIndex[k]) {
			yIndex[k] = Height2 - yIndex[k];
		}
	}
	/* perform interpolation */
	interpolated = 0.0;
	for (j = 0L; j <= SplineDegree; j++) {
		p = Bcoeff + (ptrdiff_t)(yIndex[j] * Width);
		w = 0.0;
		for (i = 0L; i <= SplineDegree; i++) {
			w += xWeight[i] * p[xIndex[i]];
		}
		interpolated += yWeight[j] * w;
	}
	return(interpolated);
} /* end InterpolatedValue */

/*********************************************************************/
 
void BSPLINE_DEC::shift_rotation_from_spline_coeff(float *Data, float *Result, 
                           int Nx, int Ny, double AngleParam,
                           double xShift, double yShift,  
			   double xOrigin, double yOrigin)
{ 
    float	*p;
    double	a11, a12, a21, a22;
    double	x0, y0, x1, y1;
    double	Angle=AngleParam;
    long	Width=Nx, Height=Ny;
    long	x, y;
    int		Masking=(RotationMasking == True) ? 1: 0;
     
    /* prepare the geometry */
    Angle *= PI / 180.0;
    a11 = cos(Angle);
    a12 = -sin(Angle);
    a21 = sin(Angle);
    a22 = cos(Angle);
    x0 = a11 * (xShift + xOrigin) + a12 * (yShift + yOrigin);
    y0 = a21 * (xShift + xOrigin) + a22 * (yShift + yOrigin);
    xShift = xOrigin - x0;
    yShift = yOrigin - y0;
    
    p = Result;
    for (y = 0L; y < Height; y++) 
    {
       x0 = a12 * (double)y + xShift;
       y0 = a22 * (double)y + yShift;
       for (x = 0L; x < Width; x++) 
       {
           x1 = x0 + a11 * (double)x;
           y1 = y0 + a21 * (double)x;
           if (Masking) 
	   {
               if ((x1 <= -0.5) || (((double)Width - 0.5) <= x1)
					|| (y1 <= -0.5) || (((double)Height - 0.5) <= y1)) 
 		     *p++ = 0.0F;
 	       else  *p++ = (float)InterpolatedValue(Data, Nx, Ny, x1, y1);
 	   }
	   else  *p++ = (float) InterpolatedValue(Data, Nx, Ny, x1, y1);
       }
    }
}
		
/**********************************************************************/
//                     END BSPLINE  FUNCTION
/**********************************************************************/
		
void translate1d(fltarray & x, float Dx)
{
   // Translate
   FFTN_1D FFT;
   int i,N = x.nx();
   float t = REM(Dx, N);
   double Pos,Sre,Sim,Rre,Rim;
 
   cfarray Line(N);
   for (i=0; i< N; i++) Line(i) = complex_f(x(i), 0.);
   FFT.fftn1d (Line);
   for (i=0; i< N; i++) 
   {
      Pos = (double) i - (double) N / 2.; //  + 0.5;
      Sre =  cos(2.*PI*Pos *t/(double) N);
      Sim = -sin(2.*PI*Pos *t/(double) N);
      Rre =  Line(i).real() * Sre -  Line(i).imag() * Sim;
      Rim =  Line(i).real() * Sim +  Line(i).imag() * Sre;
      Line(i) = complex_f((float) Rre, (float) Rim);
   }
   FFT.fftn1d (Line, True);
   for (i=0; i < N; i++)  x(i) = (Line(i)).real();
}

/**********************************************************************/

void im_zoom(Ifloat &Ima1, Ifloat &Ima2, double ZoomX, double ZoomY)
{
   int Nl2 = (int) (Ima1.nl()*ZoomY+0.5);
   int Nc2 = (int) (Ima1.nc()*ZoomX+0.5);
   Ima2.resize(Nl2,Nc2);
   BSPLINE_DEC BSD;
   BSD.SamplesToCoefficients(Ima1.buffer(), Ima1.nc(), Ima1.nl());
   for (int i=0; i < Ima2.nl(); i++) 
   for (int j=0; j < Ima2.nc(); j++) 
         Ima2(i,j) = BSD.InterpolatedValue(Ima1, (double) j / ZoomX+0.5, (double) i / ZoomY+0.5);
}

/**********************************************************************/

void im_zoom(Ifloat &Ima1, Ifloat &Ima2, double Zoom)
{
   im_zoom( Ima1,  Ima2,  Zoom,  Zoom);
}

/**********************************************************************/

void im_zoom(Ifloat &Ima1, Ifloat &Ima2)
{
    double ZoomX = (float) Ima2.nc() / (float) Ima1.nc();
    double ZoomY = (float) Ima2.nl() / (float) Ima1.nl();
    BSPLINE_DEC BSD;
    BSD.SamplesToCoefficients(Ima1.buffer(), Ima1.nc(), Ima1.nl());
    for (int i=0; i < Ima2.nl(); i++) 
    for (int j=0; j < Ima2.nc(); j++) 
        Ima2(i,j) = BSD.InterpolatedValue(Ima1, (double) j / ZoomX, (double) i / ZoomY);
}
 
/**********************************************************************/

void translate2d(Ifloat & Ima, float Dx, float Dy)
{
   // Translate
   int i,j;
   int Nl = Ima.nl();
   int Nc = Ima.nc();
   float ty = REM(-Dy, Nl);
   float tx = REM(-Dx, Nc);
   double Posi, Posj,Sre,Sim,Rre,Rim;
   Icomplex_f Data_cf(Nl, Nc, "cf");
   FFTN_2D FFT;
   
   // fft2d (Ima, Data_cf,1);
   FFT.fftn2d(Ima, Data_cf);
   for (i=0; i< Nl; i++) 
   for (j=0; j< Nc; j++) 
   {
      Posi = (double) i - (double) Nl / 2.;
      Posj = (double) j - (double) Nc / 2.;
      double u = 2.*PI*(Posj *tx/(double) Nc + Posi *ty/(double) Nl);
      Sre =  cos(u);
      Sim = -sin(u);
      Rre =   Data_cf(i,j).real() * Sre -   Data_cf(i,j).imag() * Sim;
      Rim =   Data_cf(i,j).real() * Sim +   Data_cf(i,j).imag() * Sre;
       Data_cf(i,j) = complex_f((float) Rre, (float) Rim);
   }
   // fft2d (Data_cf, -1);
   FFT.ifftn2d(Data_cf);
   for (i=0; i< Nl; i++) 
   for (j=0; j< Nc; j++)  Ima(i,j) = ( Data_cf(i,j)).real();
}

/**********************************************************************/

// void translate1d_dbl(fltarray & x, float Dx)
// {
//    // Translate
//    int i,N = x.nx();
//    double t = (double) REM(Dx, N);
//    double Pos,Sre,Sim,Rre,Rim;
//  
//    complex_double *Line = new  complex_double [N];
//    for (i=0; i< N; i++) 
//    {
//       Line[i].re = (double) x(i);
//       Line[i].im = 0.;
//    }
//    ft_cd_1d (Line, 1, N);
//    // cout << " N/2=" << x(i) << " " << Line[N/2].re << endl;
// 
//    for (i=0; i< N; i++) 
//    {
//       Pos = (double) i - (double) N / 2.; //  + 0.5;
//       Sre =  cos(2.*PI*Pos *t/(double) N);
//       Sim = -sin(2.*PI*Pos *t/(double) N);
//       Rre = Line[i].re * Sre -  Line[i].im * Sim;
//       Rim = Line[i].re * Sim +  Line[i].im * Sre;
//       Line[i].re  = Rre;
//       Line[i].im  = Rim;
//    }
//    ft_cd_1d (Line, -1, N);
// // cout << " -1 N/2=" << Line[N/2].re << endl;
//    for (i=0; i < N; i++)  x(i) = Line[i].re;
//    delete [] Line;
// }

/**********************************************************************/

void Rotation::im_rot90(Ifloat &DataIn, Ifloat & DataRot, int NRot)
{
    int i,j;
    int Nl = DataIn.nl();
    int Nc = DataIn.nc();
    int Nr = NRot % 4;
    if (Nr < 0) Nr += 4;
    if ((Nr == 1) || (Nr == 3)) DataRot.resize(Nc, Nl);
    else DataRot.resize(Nl, Nc);

    for (i=0; i< Nl; i++)
    for (j=0; j< Nc; j++)
    {
        switch (Nr)
        {
           case 1: DataRot(j,i) = DataIn(Nl-i-1,j); break;
           case 2: DataRot(i,j) = DataIn(Nl-i-1, Nc-j-1); break;
           case 3: DataRot(j,i) = DataIn(i,Nc-j-1); break;
           default: DataRot(i,j) = DataIn(i,j); break;
        }
    } 
}

/**********************************************************************/

void Rotation::im_shear(Ifloat &Data, char Dir, float beta)
{
   int i,k;
   int Nl = Data.nl();
   int Nc = Data.nc();
   float x,y;

   if (Dir == 'x')
   {
      // cout << "x " << Nl << " " << Nc <<  endl;
      int n = Nc;
      fltarray R(Nl);
      for (i=0; i < n; i++)
      {
         for (k=0; k < Nl; k++) R(k) = Data(k,i);
         y = (float) i - iround((float)n/2.) + 1./2.;
         translate1d(R,beta*y);
         for (k=0; k < Nl; k++)  Data(k,i) = R(k);
      }
   }
   if (Dir == 'y')
   {
      // cout << "y" << endl;
      int n = Nl;
      fltarray R(Nc);
      for (i=0; i < n; i++)
      {
         for (k=0; k < Nc; k++) R(k) = Data(i,k);
         x = (float) i - iround((float)n/2.) + 1./2.;
         translate1d(R,beta*x);
         for (k=0; k < Nc; k++)  Data(i,k) = R(k);
      }
   }
}

/**********************************************************************/

void Rotation::im_rotate(Ifloat &Data, Ifloat &  DataRot, float Angle)
// Rotate_3 - Image Rotation by the Three Pass Algorithm of Unser 1995
//  Inputs
//    Data     image size n by n, n even
//    Angle  angle of rotation
//  Outputs
//    DataRot    rotated image
// 
// 
// imr = rot90(im , fix( theta / (pi/2) ));
// beta = rem( theta , (pi/2) );
// % Three pass algorithm for rotation
// imr = Shear(imr,  'x', -tan(beta/2));
// imr = Shear(imr, 'y', sin(beta));
// imr = Shear(imr, 'x', -tan(beta/2));
{
   float Theta = Angle;
   if (RadianUnitAngle == False) Theta = Theta / 180. * PI;
   int NRot = (int) (Theta /(PI/2.));
   // cout << "NRot = " << NRot << endl;
   // cout << "Data  = " <<  Data.nl() << " " <<  Data.nc()<< endl;
   // cout << "DataRot  = " <<  DataRot.nl() << " " <<  DataRot.nc() << endl;

   if (Verbose == True) cout << "NRot 90 degrees = " << NRot << endl;
   im_rot90(Data, DataRot, NRot);
   float Beta = REM(Theta, (PI/2));
   // cout << "Theta = " << Theta << " Beta  = " <<  Beta << endl;
   if (Verbose == True) cout << "im_shear x" << endl;
   im_shear(DataRot, 'x', -tan(Beta/2));
   if (Verbose == True) cout << "im_shear y" << endl;
   im_shear(DataRot, 'y', sin(Beta));
   if (Verbose == True) cout << "im_shear x" << endl;
   im_shear(DataRot, 'x', -tan(Beta/2));
}

/**********************************************************************/

void Rotation::im_move_rotate_bilinear(Ifloat &Data1, Ifloat & Data2, float AngleRot, 
                              float Dx, float Dy)
{
   float  v1,v2,v3,v4;
   int Nc=Data1.nc();
   int Nl=Data1.nl();
   float Angle = -AngleRot;
   float SinAngle, CosAngle, X, Y;
   float w1,w2,w3,w4;
   float Xscale  = 1.0;
   float Yscale  = 1.0;
   int i,j,Xhalf, Yhalf;
   Data2.resize(Nl,Nc);
 
   if (RadianUnitAngle == False) Angle = PI * Angle / 180.;
   SinAngle = (float) sin( (double) Angle);
   CosAngle = (float) cos( (double) Angle);
   if (Xscale != 0.0) Xscale = 1.0 / Xscale;
   if (Yscale != 0.0) Yscale = 1.0 / Yscale;

   Yhalf = Nl / 2;
   Xhalf = Nc / 2;
   for (i = 0; i < Nl; i++) 
   for (j = 0; j < Nc; j++) 
   {
       X = Xscale * ((j-Xhalf) * CosAngle 
           - (i-Yhalf) * SinAngle) + Xhalf - Dx;
       Y = Yscale * ((i-Yhalf) * CosAngle 
           + (j-Xhalf) * SinAngle) + Yhalf - Dy;
       if (X < 0) X += Nc;
       if (X >= Nc) X -= Nc;
       if (Y < 0) Y += Nl;
       if (Y >= Nl) Y -= Nl;
         
       if ((X < 0) || (X >= Nc) || (Y < 0) || (Y >= Nl)) Data2(i,j) = 0.0;
       else 
       {
          v1 = v2 = v3 = v4 = Data1((int)Y,(int)X);
          if (X+1 < Nc)  v2 = Data1((int)Y,(int)X+1);
          if (Y+1 < Nl)  v3 = Data1((int)Y+1,(int)X);
          if ((X+1 < Nc) && (Y+1 < Nl))  v4 = Data1((int)Y+1,(int)X+1);
          w4 = (X-(int)X)* (Y-(int)Y);
          w3 = (1-X+(int)X)*(Y-(int)Y);
          w2 = (X-(int)X)*(1-Y+(int)Y);
          w1 = (1-X+(int)X)*(1-Y+(int)Y);
          Data2(i,j) = w1*v1 + w2*v2 + w3*v3 + w4*v4;
       }
    }
}

/**********************************************************************/
