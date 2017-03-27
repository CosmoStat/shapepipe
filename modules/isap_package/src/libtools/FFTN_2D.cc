/******************************************************************************
**                   Copyright (C) 2000 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 3.2
**
**    Author: Jean-Luc Starck
**
**    Date:  28/08/00
**    
**    File:  FFTN_2D.h
**
*******************************************************************************
**
**    DESCRIPTION  2D FFT Class routine for any size of images
**    ----------- 
**
******************************************************************************/

#include "FFTN_2D.h"
#include "IM_Rot.h"

/**********************************************************************/

void im_shift(Ifloat &Data1, Ifloat & Data2, int Dx, int Dy)
{
   int Nc=Data1.nc();
   int Nl=Data1.nl();
   int X,Y,i,j;
   Data2.resize(Nl,Nc);
  
   for (i = 0; i < Nl; i++) 
   for (j = 0; j < Nc; j++) 
   {
       X = j + Dx;
       Y = i + Dy;
       if (X < 0) X += Nc;
       if (X >= Nc) X -= Nc;
       if (Y < 0) Y += Nl;
       if (Y >= Nl) Y -= Nl;
       if ((X < 0) || (X >= Nc) || (Y < 0) || (Y >= Nl)) Data2(i,j) = 0.0;
       else Data2(i,j) = Data1(Y,X);
   }
}

/**********************************************************************/

void im_shift(Icomplex_f &Data1, Icomplex_f & Data2, int Dx, int Dy)
{
   int Nc=Data1.nc();
   int Nl=Data1.nl();
   int X,Y,i,j;
   Data2.resize(Nl,Nc);
  
   for (i = 0; i < Nl; i++) 
   for (j = 0; j < Nc; j++) 
   {
       X = j + Dx;
       Y = i + Dy;
       if (X < 0) X += Nc;
       if (X >= Nc) X -= Nc;
       if (Y < 0) Y += Nl;
       if (Y >= Nl) Y -= Nl;
       if ((X < 0) || (X >= Nc) || (Y < 0) || (Y >= Nl)) 
            Data2(i,j) = complex_f(0.,0.);
       else Data2(i,j) = Data1(Y,X);
   }
}

/**********************************************************************/

void im_shift(Icomplex_d &Data1, Icomplex_d & Data2, int Dx, int Dy)
{
	int Nc=Data1.nc();
	int Nl=Data1.nl();
	int X,Y,i,j;
	Data2.resize(Nl,Nc);
	
	for (i = 0; i < Nl; i++) 
		for (j = 0; j < Nc; j++) 
		{
			X = j + Dx;
			Y = i + Dy;
			if (X < 0) X += Nc;
			if (X >= Nc) X -= Nc;
			if (Y < 0) Y += Nl;
			if (Y >= Nl) Y -= Nl;
			if ((X < 0) || (X >= Nc) || (Y < 0) || (Y >= Nl)) 
				Data2(i,j) = complex_d(0.,0.);
			else Data2(i,j) = Data1(Y,X);
		}
}

/**********************************************************************/

void FFTN_2D::center(Ifloat &Image)
{
   int Nc=Image.nc();
   int Nl=Image.nl();
   Ifloat Dat(Nl,Nc,"dat");
   Dat = Image;
   im_shift(Dat, Image, (Nc+1)/2, (Nl+1)/2);
}
/**********************************************************************/

void FFTN_2D::uncenter(Ifloat &Image)
{
   int Nc=Image.nc();
   int Nl=Image.nl();
   Ifloat Dat(Nl,Nc,"dat");
   Dat = Image;
   im_shift(Dat, Image, -(Nc+1)/2, -(Nl+1)/2);
}

/**********************************************************************/

void FFTN_2D::center(Icomplex_f &Image)
{
   int Nc=Image.nc();
   int Nl=Image.nl();
   Icomplex_f Dat(Nl,Nc,"dat");
   Dat = Image;
   im_shift(Dat, Image, (Nc+1)/2, (Nl+1)/2);
}

/**********************************************************************/

void FFTN_2D::center(Icomplex_d &Image)
{
	int Nc=Image.nc();
	int Nl=Image.nl();
	Icomplex_d Dat(Nl,Nc,"dat");
	Dat = Image;
	im_shift(Dat, Image, (Nc+1)/2, (Nl+1)/2);
}

/**********************************************************************/

void FFTN_2D::uncenter(Icomplex_f &Image)
{
   int Nc=Image.nc();
   int Nl=Image.nl();
   Icomplex_f Dat(Nl,Nc,"dat");
   Dat = Image;
   im_shift(Dat, Image, -(Nc+1)/2, -(Nl+1)/2);
}

/**********************************************************************/

void FFTN_2D::uncenter(Icomplex_d &Image)
{
	int Nc=Image.nc();
	int Nl=Image.nl();
	Icomplex_d Dat(Nl,Nc,"dat");
	Dat = Image;
	im_shift(Dat, Image, -(Nc+1)/2, -(Nl+1)/2);
}

/**********************************************************************/

void FFTN_2D::swap_buff(Icomplex_f &Ima, Bool Reverse)
{
    int i,j,Indi,Indj;
    int Nl = Ima.nl();
    int Nc = Ima.nc();
    complex_f *Buff = Ima.buffer();
    complex_f temp;
    if (Reverse == False)
    {
       if (Nc % 2 != 0) 
        for (i = 0; i < Nl; i++)
        {
          temp = Ima(i, Nc/2);
          for (j = Nc/2; j < Nc-1; j++) Ima(i, j) = Ima(i, j+1);
          Ima(i, Nc-1) = temp;
        }
       if (Nl % 2 != 0)
        for (j = 0; j < Nc; j++)
         {
          temp = Ima(Nl/2, j);
          for (i = Nl/2; i < Nl-1; i++) Ima(i, j) = Ima(i+1, j);
          Ima(Nl-1, j) = temp;
         }
    }

    for (i = 0; i < Nl/2; i++)
    for (j = 0; j < Nc/2; j++)
    {
       Indi = Nl/2+i;
       Indj = Nc/2+j;
       pix_swap(Buff[i*Nc+j], Buff[Indi*Nc+Indj]);
       pix_swap(Buff[Indi*Nc+j], Buff[i*Nc+Indj]);
    }  

    if (Reverse == True)
    {
       if (Nl % 2 != 0)
        for (j = 0; j < Nc; j++)
         {
          temp = Ima(Nl-1, j);
          for (i = Nl-1; i > Nl/2; i--) Ima(i, j) = Ima(i-1, j);
          Ima(Nl/2, j) = temp;
         }
       if (Nc % 2 != 0) 
        for (i = 0; i < Nl; i++)
        {
          temp = Ima(i, Nc-1);
          for (j = Nc-1; j > Nc/2; j--) Ima(i, j) = Ima(i, j-1);
          Ima(i, Nc/2) = temp;
        }

    }
}


/**********************************************************************/

void FFTN_2D::swap_buff(Icomplex_d &Ima, Bool Reverse)
{
    int i,j,Indi,Indj;
    int Nl = Ima.nl();
    int Nc = Ima.nc();
    complex_d *Buff = Ima.buffer();
    complex_d temp;
    if (Reverse == False)
    {
		if (Nc % 2 != 0) 
			for (i = 0; i < Nl; i++)
			{
				temp = Ima(i, Nc/2);
				for (j = Nc/2; j < Nc-1; j++) Ima(i, j) = Ima(i, j+1);
				Ima(i, Nc-1) = temp;
			}
		if (Nl % 2 != 0)
			for (j = 0; j < Nc; j++)
			{
				temp = Ima(Nl/2, j);
				for (i = Nl/2; i < Nl-1; i++) Ima(i, j) = Ima(i+1, j);
				Ima(Nl-1, j) = temp;
			}
    }
	
    for (i = 0; i < Nl/2; i++)
		for (j = 0; j < Nc/2; j++)
		{
			Indi = Nl/2+i;
			Indj = Nc/2+j;
			pix_swap(Buff[i*Nc+j], Buff[Indi*Nc+Indj]);
			pix_swap(Buff[Indi*Nc+j], Buff[i*Nc+Indj]);
		}  
	
    if (Reverse == True)
    {
		if (Nl % 2 != 0)
			for (j = 0; j < Nc; j++)
			{
				temp = Ima(Nl-1, j);
				for (i = Nl-1; i > Nl/2; i--) Ima(i, j) = Ima(i-1, j);
				Ima(Nl/2, j) = temp;
			}
		if (Nc % 2 != 0) 
			for (i = 0; i < Nl; i++)
			{
				temp = Ima(i, Nc-1);
				for (j = Nc-1; j > Nc/2; j--) Ima(i, j) = Ima(i, j-1);
				Ima(i, Nc/2) = temp;
			}
		
    }
}

/******************************************************************************/

void FFTN_2D::convolve(Ifloat &Ima1, Ifloat &Ima2, Ifloat &Result)
{
    int i,j,Nl = Ima1.nl();
    int Nc = Ima1.nc();
    
    if ((Ima2.nl() != Nl) || (Ima2.nc() != Nc))
    {
       cout << "Error in FFTN_2D::convolve: images have not the same size ... " << endl;
       cout << "   Ima 1: Nl = " << Nl << " Nc = " << Nc << endl;
       cout << "   Ima 2: Nl = " << Ima2.nl() << " Nc = " << Ima2.nc() << endl;
       exit(-1);
    }
    if ((Result.nl() != Nl) || (Result.nc() != Nc)) 
                                            Result.alloc(Nl,Nc,"convol");
    
   Icomplex_f Ima1_cf (Nl, Nc, "Buffer1 conv");
   Icomplex_f Ima2_cf (Nl, Nc, "Buffer2 conv");
   fftn2d(Ima1,Ima1_cf);
   fftn2d(Ima2,Ima2_cf);
   Ima1_cf *= Ima2_cf;
   fftn2d(Ima1_cf, True);
   for (i = 0; i < Nl; i++)
   for (j = 0; j < Nc; j++) Result(i,j) = Ima1_cf(i,j).real();   
}

/******************************************************************************/

void FFTN_2D::convolve(Ifloat &Ima1, Icomplex_f & Ima2_cf, Ifloat &Result)
{
    int i,j,Nl = Ima1.nl();
    int Nc = Ima1.nc();
    
    if ((Ima2_cf.nl() != Nl) || (Ima2_cf.nc() != Nc))
    {
       cout << "Error in FFTN_2D::convolve: images have not the same size ... " << endl;
       cout << "   Ima 1: Nl = " << Nl << " Nc = " << Nc << endl;
       cout << "   Ima 2: Nl = " << Ima2_cf.nl() << " Nc = " << Ima2_cf.nc() << endl;
       exit(-1);
    }
    if ((Result.nl() != Nl) || (Result.nc() != Nc)) 
                                            Result.alloc(Nl,Nc,"convol");
    
   Icomplex_f Ima1_cf (Nl, Nc, "Buffer1 conv");
   fftn2d(Ima1,Ima1_cf);
   Ima1_cf *= Ima2_cf;
   fftn2d(Ima1_cf, True);
   for (i = 0; i < Nl; i++)
   for (j = 0; j < Nc; j++) Result(i,j) = Ima1_cf(i,j).real();   
}

/******************************************************************************/

void FFTN_2D::convolve_conj(Ifloat &Ima1, Icomplex_f & Ima2_cf, Ifloat &Result)
{
    int i,j,Nl = Ima1.nl();
    int Nc = Ima1.nc();
    
    if ((Ima2_cf.nl() != Nl) || (Ima2_cf.nc() != Nc))
    {
       cout << "Error in FFTN_2D::convolve: images have not the same size ... " << endl;
       cout << "   Ima 1: Nl = " << Nl << " Nc = " << Nc << endl;
       cout << "   Ima 2: Nl = " << Ima2_cf.nl() << " Nc = " << Ima2_cf.nc() << endl;
       exit(-1);
    }
    if ((Result.nl() != Nl) || (Result.nc() != Nc)) 
                                            Result.alloc(Nl,Nc,"convol");
    
   Icomplex_f Ima1_cf (Nl, Nc, "Buffer1 conv");
   fftn2d(Ima1,Ima1_cf);
   for (i = 0; i < Nl; i++)
   for (j = 0; j < Nc; j++) Ima1_cf(i,j) *=  conj(Ima2_cf(i,j));
   fftn2d(Ima1_cf, True);
   for (i = 0; i < Nl; i++)
   for (j = 0; j < Nc; j++) Result(i,j) = Ima1_cf(i,j).real();   
}

/******************************************************************************/

void FFTN_2D::fftn2d(Ifloat &Image, Icomplex_f &Im_out, Bool Reverse, bool normalize)
{
    int Ind=0;
    int Nl = Image.nl();
    int Nc = Image.nc();

    // for square image, we use fft2d (run faster with a better precision)
    // fft2d return the FFT whith the zero frequency in the middle
    if ((is_power_of_2(Nl) == True) && (is_power_of_2(Nc) == True) && (Nl==Nc))
    {
        if ((Reverse == True) && (CenterZeroFreq == False))
                                             swap_buff(Im_out, Reverse);
        int Dir = (Reverse == False) ? 1: -1; 
        fft2d(Image, Im_out, Dir, normalize);
        if ((Reverse == False) && (CenterZeroFreq == False)) 
                                            swap_buff(Im_out, Reverse);
    }
    else
    {
       float *FFT_Dat = (float *)  &(Im_out(0,0));
       for (int i = 0; i < Nl; i++)
       for (int j = 0; j < Nc; j++)
       {
          FFT_Dat[Ind++] = Image(i,j);
          FFT_Dat[Ind++] = 0.;
       }
       if (CenterZeroFreq == True) uncenter(Im_out);
       transform2d(FFT_Dat, Nc, Nl, Reverse, normalize);
       if (CenterZeroFreq == True) center(Im_out);
   }
}
 
/******************************************************************************/
 
void FFTN_2D::fftn2d(Ifloat &Image, Icomplex_d &Im_out, Bool Reverse)
{
    int Ind=0;
    int Nl = Image.nl();
    int Nc = Image.nc();
	
    double *FFT_Dat = (double *)  &(Im_out(0,0));
    for (int i = 0; i < Nl; i++)
	for (int j = 0; j < Nc; j++)
	{
       FFT_Dat[Ind++] = Image(i,j);
       FFT_Dat[Ind++] = 0.;
	}
    if (CenterZeroFreq == True) uncenter(Im_out);
	transform2d(FFT_Dat, Nc, Nl, Reverse);
    if (CenterZeroFreq == True) center(Im_out);
}

/******************************************************************************/

void FFTN_2D::fftn2d(Icomplex_f &Im_out, Bool Reverse, bool normalize)
{
    int Nl = Im_out.nl();
    int Nc = Im_out.nc();

    if ((is_power_of_2(Nl) == True) && (is_power_of_2(Nc) == True) && (Nl==Nc))
    {
        int Dir = (Reverse == False) ? 1: -1; 
        if ((Reverse == True) && (CenterZeroFreq == False))
                                             swap_buff(Im_out, Reverse);
        fft2d(Im_out, Dir, normalize);
        if ((Reverse == False) && (CenterZeroFreq == False)) 
                                             swap_buff(Im_out, Reverse);
    }
    else
    {
       float *FFT_Dat = (float *)  &(Im_out(0,0));
        if (CenterZeroFreq == True) uncenter(Im_out);
       transform2d(FFT_Dat, Nc, Nl, Reverse, normalize);
       if (CenterZeroFreq == True) center(Im_out);
    }
}


/******************************************************************************/
 
void FFTN_2D::fftn2d(Icomplex_d &Im_out, Bool Reverse)
{
    int Nl = Im_out.nl();
    int Nc = Im_out.nc();
	
    double *FFT_Dat = (double *)  &(Im_out(0,0));
    if (CenterZeroFreq == True) uncenter(Im_out);
    transform2d(FFT_Dat, Nc, Nl, Reverse);
	if (CenterZeroFreq == True) center(Im_out);
}

/******************************************************************************/

void get_isotropic_spectrum(Ifloat & Imag, fltarray & Spectrum, float Resol, float Dens)
{
   int i,j;
   int Nl = Imag.nl();
   int Nc = Imag.nc();
   float Np = Nl*Nc;
   Ifloat Data(Nl,Nc, "data");
   Icomplex_f Ima1_cf (Nl, Nc, "Buffer1 conv");
   FFTN_2D FFT;
   
   FFT.fftn2d(Imag,Ima1_cf);
   for (i=0; i < Nl; i++)
   for (j=0; j < Nc; j++) Data(i,j) =  (float) (norm(Ima1_cf(i,j)) / (double) Np);
   int N2 = (Nl+1)/2*sqrt(2.);
   int NpSpec = (int) ((N2 - 1) / Resol);

//    int x=Nl/2;
//    int y=Nc/2;
//    Data(x,y) = 0.;
//    for (i=x-1; i <= x+1; i++)
//    for (j=y-1; j <= y+1; j++)  if ((i!=x) || (j!=y)) Data(x,y) += Data(i,j);
//    Data(x,y) /= 8.;
//    
   Spectrum.alloc(NpSpec,3);
   Spectrum(0,1) = Data(Nl/2,Nc/2);
   Spectrum(0,0) = 0;
   Spectrum(0,2) = 0;
   
   BSPLINE_DEC BD;
   BD.SamplesToCoefficients(Data);

   // Spectrum Calculation
   for (int r = 1; r < NpSpec; r++)
   {
      float Rad = r*Resol;
      float L = Rad / Nl;
      int Nu=0;
      int Np = (int) (Dens*2*PI*Rad+0.5);  // number of points on the circle of radius r
      for (i=0; i < Np; i++)
      {
         float Angle = 2. *PI / (float) Np * i;
	 double x = Nc/2 + Rad * cos(Angle);
	 double y = Nl/2 + Rad * sin(Angle);
	 if ((x >= 0) && (y >= 0) && (x < Nc) && (y < Nl))
	 {
	    float I = MAX(0, BD.InterpolatedValue(Data,x,y));
	    Spectrum(r,1) += I;
	    Nu ++;
 	 }
	 // printf("r=%5.3f, A=%5.3f, x=%5.3f, y=%5.3f, P=%5.3f\n",Rad,Angle,x,y,Spectrum(r,1));
      }
      if (Nu > 0) Spectrum(r,1) /= (float) Nu;
      Spectrum(r,0) = L;
      Spectrum(r,2) = Spectrum(r,1) * L*(L+1);
   }
}

/****************************************************************************/
