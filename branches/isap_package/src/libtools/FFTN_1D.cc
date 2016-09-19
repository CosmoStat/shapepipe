/******************************************************************************
**                   Copyright (C) 2001 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 3.2
**
**    Author: Jean-Luc Starck
**
**    Date:  20/03/01
**    
**    File:  FFTN_1D.h
**
*******************************************************************************
**
**    DESCRIPTION  1D FFT Class routine for any size of images
**    ----------- 
**
******************************************************************************/

#include "FFTN_1D.h"
#include "DefMath.h"
using namespace std;

/**********************************************************************/

void im1d_shift(fltarray &Data1, fltarray & Data2, int Dx)
{
   int N=Data1.nx();
   int X,i;
   // Data2.reform(N);
  
   for (i = 0; i < N; i++) 
   {
       X = i + Dx;
       if (X < 0) X += N;
       if (X >= N) X -= N;
       if ((X < 0) || (X >= N)) Data2(i) = 0.0;
       else Data2(i) = Data1(X);
   }
}
/**********************************************************************/

void im1d_shift(dblarray &Data1, dblarray & Data2, int Dx)
{
   int N=Data1.nx();
   int X,i;
   // Data2.reform(N);
  
   for (i = 0; i < N; i++) 
   {
       X = i + Dx;
       if (X < 0) X += N;
       if (X >= N) X -= N;
       if ((X < 0) || (X >= N)) Data2(i) = 0.0;
       else Data2(i) = Data1(X);
   }
}
/**********************************************************************/

void im1d_shift(complex_f *Data1, complex_f *Data2, int N, int Dx)
{
   int X,i;
   
   for (i = 0; i < N; i++) 
   {
       X = i + Dx;
       if (X < 0) X += N;
       if (X >= N) X -= N;
       if ((X < 0) || (X >= N))  Data2[i] = complex_f(0.,0.);
       else Data2[i] = Data1[X];
   }
}

/**********************************************************************/

void im1d_shift(complex_d *Data1, complex_d *Data2, int N, int Dx)
{
   int X,i;
   
   for (i = 0; i < N; i++) 
   {
       X = i + Dx;
       if (X < 0) X += N;
       if (X >= N) X -= N;
       if ((X < 0) || (X >= N))  Data2[i] = complex_d(0.,0.);
       else Data2[i] = Data1[X];
   }
}

/**********************************************************************/

// void FFTN_1D::center(fltarray &Dat)
// {
//    int N=Dat.nx();
//    fltarray Buf(N);
//    Buf = Dat;
//    im1d_shift(Buf, Dat, N/2);
// }
// /**********************************************************************/
// 
// void FFTN_1D::center(dblarray &Dat)
// {
//    int N=Dat.nx();
//    dblarray Buf(N);
//    Buf = Dat;
//    im1d_shift(Buf, Dat, N/2);
// }
// /**********************************************************************/
// 
// void FFTN_1D::uncenter(fltarray &Dat)
// {
//    int N=Dat.nx();
//    fltarray Buf(N);
//    Buf = Dat;
//    im1d_shift(Buf, Dat, -N/2);
// }
// /**********************************************************************/
// 
// void FFTN_1D::uncenter(dblarray &Dat)
// {
//    int N=Dat.nx();
//    dblarray Buf(N);
//    Buf = Dat;
//    im1d_shift(Buf, Dat, -N/2);
// }

/**********************************************************************/

void FFTN_1D::center(complex_f *Dat, int N)
{
   int i;
   complex_f *Buff = new complex_f[N];
   for (i=0; i < N; i++) Buff[i] = Dat[i];
   im1d_shift( Buff,  Dat, N, (N+1)/2);
   delete [] Buff;
}

/**********************************************************************/

void FFTN_1D::uncenter(complex_f *Dat, int N)
{
   int i;
   complex_f *Buff = new complex_f[N];
   for (i=0; i < N; i++) Buff[i] = Dat[i];
   im1d_shift(Buff,  Dat, N, -(N+1)/2);
   delete [] Buff;
}
/**********************************************************************/

void FFTN_1D::center(complex_d *Dat, int N)
{
   int i;
   complex_d *Buff = new complex_d[N];
   for (i=0; i < N; i++) Buff[i] = Dat[i];
   im1d_shift( Buff,  Dat, N, (N+1)/2);
   delete [] Buff;
}

/**********************************************************************/

void FFTN_1D::uncenter(complex_d *Dat, int N)
{
   int i;
   complex_d *Buff = new complex_d[N];
   for (i=0; i < N; i++) Buff[i] = Dat[i];
   im1d_shift(Buff,  Dat, N, -(N+1)/2);
   delete [] Buff;
}
/**********************************************************************/

void FFTN_1D::swap_buff(complex_f *Buff, int N, Bool Reverse)
{
    int i,Indi;
    complex_f temp;

    if (Reverse == False)
    {
       if (N % 2 != 0) 
       {
          temp = Buff[N/2];
          for (i = N/2; i < N-1; i++) Buff[i] = Buff[i+1];
          Buff[N-1] = temp;
       }
    }

    for (i = 0; i < N/2; i++)
    {
       Indi = N/2+i;
       pix_swap(Buff[i], Buff[Indi]);
    }  

    if (Reverse == True)
    {
       if (N % 2 != 0)
       {
          temp = Buff[N-1];
          for (i = N-1; i > N/2; i--) Buff[i] = Buff[i-1];
          Buff[N/2] = temp;
       }
    }
}
/******************************************************************************/

void FFTN_1D::swap_buff(complex_d *Buff, int N, Bool Reverse)
{
    int i,Indi;
    complex_d temp;

    if (Reverse == False)
    {
       if (N % 2 != 0) 
       {
          temp = Buff[N/2];
          for (i = N/2; i < N-1; i++) Buff[i] = Buff[i+1];
          Buff[N-1] = temp;
       }
    }

    for (i = 0; i < N/2; i++)
    {
       Indi = N/2+i;
       pix_swap(Buff[i], Buff[Indi]);
    }  

    if (Reverse == True)
    {
       if (N % 2 != 0)
       {
          temp = Buff[N-1];
          for (i = N-1; i > N/2; i--) Buff[i] = Buff[i-1];
          Buff[N/2] = temp;
       }
    }
}

/******************************************************************************/

void FFTN_1D::swap_buff(dblarray &Data,  Bool Reverse)
{
    int N = Data.nx();
    int i,Indi;
    double temp;
    double *Buff = Data.buffer();

    if (Reverse == False)
    {
       if (N % 2 != 0) 
       {
          temp = Buff[N/2];
          for (i = N/2; i < N-1; i++) Buff[i] = Buff[i+1];
          Buff[N-1] = temp;
       }
    }

    for (i = 0; i < N/2; i++)
    {
       Indi = N/2+i;
       pix_swap(Buff[i], Buff[Indi]);
    }  

    if (Reverse == True)
    {
       if (N % 2 != 0)
       {
          temp = Buff[N-1];
          for (i = N-1; i > N/2; i--) Buff[i] = Buff[i-1];
          Buff[N/2] = temp;
       }
    }
}

/******************************************************************************/

void FFTN_1D::fftn1d(dblarray &Data, complex_d *Data_Out, Bool Reverse, bool normalize)
{
    int Ind=0;
    int N = Data.nx();
 
    if ((Reverse == True) && (CenterZeroFreq == True))
                                            swap_buff(Data_Out, N, Reverse);
    double *FFT_Dat = (double *) Data_Out;
    for (int i = 0; i < N; i++)
    {
       FFT_Dat[Ind++] = Data(i);
       FFT_Dat[Ind++] = 0.;
    }
    if (CenterZeroFreq == True) uncenter(Data_Out, N);
    transform1d(FFT_Dat, N, Reverse, normalize);
    if (CenterZeroFreq == True) center(Data_Out, N);
    //  transform1d return the FFT whith the zero frequency at the left
}

/******************************************************************************/

void FFTN_1D::fftn1d(fltarray &Data, complex_d *Data_Out, Bool Reverse, bool normalize)
{
    int Ind=0;
    int N = Data.nx();
 
    // if ((Reverse == True) && (CenterZeroFreq == True))
    //                                        swap_buff(Data_Out, N, Reverse);
    double *FFT_Dat = (double *) Data_Out;
    for (int i = 0; i < N; i++)
    {
       FFT_Dat[Ind++] = Data(i);
       FFT_Dat[Ind++] = 0.;
    }
    if (CenterZeroFreq == True) uncenter(Data_Out, N);
    transform1d(FFT_Dat, N, Reverse, normalize);
    if (CenterZeroFreq == True) center(Data_Out, N);
    //  transform1d return the FFT whith the zero frequency at the left
}

/******************************************************************************/

void FFTN_1D::fftn1d(fltarray &Data, complex_f *Data_Out, Bool Reverse, bool normalize)
{
    int Ind=0;
    int N = Data.nx();
 
    // if ((Reverse == True) && (CenterZeroFreq == True))
    //                                        swap_buff(Data_Out, N, Reverse);
    float *FFT_Dat = (float *) Data_Out;
    for (int i = 0; i < N; i++)
    {
       FFT_Dat[Ind++] = Data(i);
       FFT_Dat[Ind++] = 0.;
    }
    if (CenterZeroFreq == True) uncenter(Data_Out, N);
    transform1d(FFT_Dat, N, Reverse, normalize);
    if (CenterZeroFreq == True) center(Data_Out, N);
    //  transform1d return the FFT whith the zero frequency at the left
}
 
/******************************************************************************/

void FFTN_1D::dct1d(fltarray &Data, fltarray &Data_Out, Bool Reverse)
{
    int i,N = Data.nx();
    dblarray TabData(N);
    for (i = 0; i < N; i++) TabData(i) = Data(i);
    dct1d(TabData,Reverse);
    for (i = 0; i < N; i++) Data_Out(i) = (float) TabData(i);
}

/******************************************************************************/

void FFTN_1D::dct1d(dblarray &Data,  Bool Reverse)
{
    int N = Data.nx();
    void ddct(int n, int isgn, double *a, int *ip, double *w);
    int Nw = 2*N*5/4-1;
    int Nip = 2*N;
    // double i,N2 = 2. / (double) N;
    int isgn = (Reverse == False) ?  1: -1;
    if (BuffW.nx() != Nw) BuffW.alloc(Nw);
    if (BuffIp.nx() != Nip) 
    {
       BuffIp.alloc(Nip);
       BuffIp(0) = BuffIp(1) = 0;
    }
    BuffIp.init();
    BuffW.init();
     
    if ((Reverse == True) && (CenterZeroFreq == True))
                                            swap_buff(Data, Reverse);
    // if (Reverse == True) Data(0) *= 0.5;
    ddct(N, isgn,  Data.buffer(),  BuffIp.buffer(),  BuffW.buffer());
    if ((Reverse == False) && (CenterZeroFreq == True)) 
                                           swap_buff(Data, Reverse);
//     if (Reverse == True)
//     {
//        Data(0) *= 0.5;
//        for (i = 0; i < N; i++) Data(i) *= N2;
//     }
}

/******************************************************************************/

void FFTN_1D::fftn1d(complex_f *Dat, int N, Bool Reverse, bool normalize)
{
   float *FFT_Dat = (float *) Dat;
   if (CenterZeroFreq == True) uncenter(Dat, N);
   transform1d(FFT_Dat, N, Reverse, normalize);
   if (CenterZeroFreq == True) center(Dat, N);
}

/******************************************************************************/

void FFTN_1D::fftn1d(complex_d *Dat, int N, Bool Reverse, bool normalize)
{
   double *FFT_Dat = (double *) Dat;
   if (CenterZeroFreq == True) uncenter(Dat, N);
   transform1d(FFT_Dat, N, Reverse, normalize);
   if (CenterZeroFreq == True) center(Dat, N);
}
 
/******************************************************************************/

void FFTN_1D::fftn1d (cfarray &Signal, Bool Reverse, bool normalize)
{
   int N = Signal.nx();
   fftn1d(Signal.buffer(), N, Reverse, normalize);
}

/******************************************************************************/

void FFTN_1D::fftn1d(fltarray &SignalIn, cfarray &SignalOut, Bool Reverse, bool normalize)
{
   int N = SignalIn.nx();
   for (int i=0;i<N;i++)
	   SignalOut(i) = complex<float>(SignalIn(i),0);
   fftn1d(SignalOut.buffer(), N, Reverse, normalize);
}

/******************************************************************************/

void FFTN_1D::fftn1d (cfarray &SignalIn, cfarray &SignalOut, Bool Reverse, bool normalize)
{
   int N = SignalIn.nx();
   for (int i=0; i < N; i++) SignalOut(i) = SignalIn(i);
   fftn1d(SignalOut.buffer(), N, Reverse, normalize);
}

/******************************************************************************/

void FFTN_1D::fftn1d (cdarray &Signal, Bool Reverse, bool normalize)
{
   int N = Signal.nx();
   fftn1d(Signal.buffer(), N, Reverse, normalize);
}

/******************************************************************************/

void FFTN_1D::fftn1d (cdarray &SignalIn, cdarray &SignalOut, Bool Reverse, bool normalize)
{
   int N = SignalIn.nx();
   for (int i=0; i < N; i++) SignalOut(i) = SignalIn(i);
   fftn1d(SignalOut.buffer(), N, Reverse, normalize);
}

/******************************************************************************/

void FFTN_1D::convolve(fltarray &D1, fltarray &D2, fltarray &Result)
{
    int i,N = D1.nx();
     
    if (D2.nx() != N)  
    {
       cout << "Error in FFTN_1D::convolve: arrays have not the same size ... " << endl;
       cout << "   Array 1: N = " <<  N << endl;
       cout << "   Array 2: N = " << D2.nx() << endl;
       exit(-1);
    }
    if (Result.nx() != N) Result.reform(N);
    
   complex_d *B1_cf = new complex_d [N];
   complex_d *B2_cf = new complex_d [N];
   fftn1d(D1,B1_cf);
   fftn1d(D2,B2_cf);
   for (i = 0; i < N; i++) B1_cf[i] *= B2_cf[i];
   fftn1d(B1_cf, N, True);
   for (i = 0; i < N; i++) Result(i) = (float) (B1_cf[i]).real();   
}

/******************************************************************************/

void STD_Window::rect(fltarray & Data, float W)
{
   int N = Data.nx();
   float N2 = (float) N/2.;
   int i;
   for (i = 0; i < N; i++)  
   {
      float u = ((float) i - N2) / N;
      if (ABS(u) < W) Data(i) = 1.;
      else Data(i) = 0;
   }  
}

/******************************************************************************/

void STD_Window::hamming(fltarray & Data, float W)
{
   int N = Data.nx();
   float N2 = (float) N/2.;
   float C = PI/W;
   int i;
   for (i = 0; i < N; i++)
   {
      float u = ((float) i - N2) / N;
      if (ABS(u) < W) Data(i) = 0.54 + 0.46 * cos(C*u);
      else Data(i) = 0;
   }
}


/*********************************************************************/

void STD_Window::hanning(fltarray & Data, float W)
{
   int N = Data.nx();
   float N2 = (float) N/2.;
   float C = PI/W;
   int i;
   for (i = 0; i < N; i++)
   {
      float u = ((float) i - N2) / N;
      if (ABS(u) < W) Data(i) = 0.5 + 0.5 * cos(C*u);
      else Data(i) = 0;
   }
}

/*********************************************************************/

void STD_Window::gaussian(fltarray & Data, float W)
{
   int N = Data.nx();
   float N2 = (float) N/2.;
   float C = -4.5/W/W;
   int i;
   for (i = 0; i < N; i++)
   {
      float u = ((float) i - N2) / N;
      if (ABS(u) < W) Data(i) = exp((double)C*u*u);
      else Data(i) = 0;
   }
}


/*********************************************************************/

void STD_Window::blackman(fltarray & Data, float W)
{
   int N = Data.nx();
   float N2 = (float) N/2.;
   float C1 = PI/W;
   float C2 = 2.*PI/W;
   int i;
   for (i = 0; i < N; i++)
   {
      float u = ((float) i - N2) / N;
      if (ABS(u) < W) Data(i) = 0.42+0.5*cos(C1*u) + 0.08*cos(C2*u);
      else Data(i) = 0;
   }
}

/*********************************************************************/

void STD_Window::make_win(fltarray & Data, type_std_win type, float W)
{
   switch (type)
   {
      case W_HAMMING: hamming(Data,W);break;
      case W_HANNING: hanning(Data,W); break;
      case W_GAUSS: gaussian(Data,W); break;
      case W_BLACKMAN: blackman(Data,W);break;
      case W_RECT: rect(Data,W);break;
   }
}
   
/*********************************************************************/

void ST_FFTN::alloc(int N, type_std_win type, float WinPar, int Size, int StepAna)
{
   STD_Window STDW;
   Step = StepAna;
   SizeWin = Size;
   TypeWin = type;
   WinParam = WinPar;
   Win.alloc(SizeWin);
   STDW.make_win(Win, TypeWin, WinParam);
   Nlt = SizeWin;
   Nct = N/Step;
   if (N % Step != 0) Nct ++;
}

/*********************************************************************/

void ST_FFTN::transform(fltarray & Data, complex_f *Trans)
{
   int N = Data.nx();
   if (SizeWin == 0) alloc(N, DEF_STD_WIN, DEF_PARAM_STD_WIN, N/2);
   int Nw = Win.nx();
   int i,k,l;
   complex_f *Buff = new complex_f [Nw];
   int Nct = nc(N);
   int Nlt = nl(N);

   if (Step > Nw/2)
   {
      cout << "Warning: the reconstruction is impossible when the step is larger " << endl;
      cout << "         than the half window size. " << endl;
   }
//    cout << "Min = " << Win.min() << endl;
//    cout << "Max = " << Win.max() << endl;
//    cout << "Min = " << Data.min() << endl;
//    cout << "Max = " << Data.max() << endl;  
//     cout << "N = " << N << endl;   
//     cout << "Nw = " << Nw << endl;   
//     cout << "Nct = " << Nct << endl;   
//     cout << "Step = " << Step << endl;   

   for (l=0; l < Nct; l++)
   {
      k = l*Step;
      for (i=0; i < Nw; i++) 
      {
         float Re;
	 int IndData = MIN(k,N-1);
	 int u = IndData - Nw / 2 + i;
 	 if ((u >= 0) && (u < N)) Re = Data(u)*Win(i);
	 else Re = 0.;
 	 Buff[i] = complex_f(Re, (float) 0.);
      }
      fftn1d(Buff, Nw);
      for (i=0; i < Nw; i++) 
      {
          int Ind = l+i*Nct;
	  if (Ind >= Nlt*Nct)
	  {
	     cout << "Error: too large index ... " << endl;
	     cout << "Ind = " << Ind << " Nlt = " << Nlt << " Nct = " << Nct << endl;
	     cout << "Step = " << Step << " k/Step = " << k/Step << endl;
	     exit(-1);
	  }
          Trans[Ind] = Buff[i];
      }
   }
   delete [] Buff;
}

/*********************************************************************/

void ST_FFTN::dct_transform(fltarray & Data, fltarray & Trans)
{
   int N = Data.nx();
   if (SizeWin == 0) alloc(N, DEF_STD_WIN, DEF_PARAM_STD_WIN, N/2);
   int Nw = Win.nx();
   int i,k,l;
   dblarray Buff(Nw);
   int Nct = nc(N);
   int Nlt = nl(N);
   
   if ((Nct != Trans.nx()) || (Nlt != Trans.ny()))
   {
      Trans.alloc(Nct,Nlt);
      cout << "Transform size: Nx = " << Nct << " Ny = " << Nlt << endl;
   }
   if (Step > Nw/2)
   {
      cout << "Warning: the reconstruction is impossible when the step is larger " << endl;
      cout << "         than the half window size. " << endl;
   }
//    cout << "Min = " << Win.min() << endl;
//    cout << "Max = " << Win.max() << endl;
//    cout << "Min = " << Data.min() << endl;
//    cout << "Max = " << Data.max() << endl;  
//     cout << "N = " << N << endl;   
//     cout << "Nw = " << Nw << endl;   
//     cout << "Nct = " << Nct << endl;   
//     cout << "Step = " << Step << endl;   

   for (l=0; l < Nct; l++)
   {
      k = l*Step;
      for (i=0; i < Nw; i++) 
      {
         float Re;
	 int IndData = MIN(k,N-1);
	 int u = IndData - Nw / 2 + i;
 	 if ((u >= 0) && (u < N)) Re = Data(u)*Win(i);
	 else Re = 0.;
 	 Buff(i) = Re;
      }
      dct1d(Buff);
      for (i=0; i < Nw; i++) 
      {
          int Ind = l+i*Nct;
	  if (Ind >= Nlt*Nct)
	  {
	     cout << "Error: too large index ... " << endl;
	     cout << "Ind = " << Ind << " Nlt = " << Nlt << " Nct = " << Nct << endl;
	     cout << "Step = " << Step << " k/Step = " << k/Step << endl;
	     exit(-1);
	  }
          Trans(l,i)  = Buff(i);
      }
   }
}

/*********************************************************************/

void ST_FFTN::recons_win(complex_f *Trans, fltarray & Data)
{
   int N = Data.nx();
   int i,k,l;
   if (SizeWin == 0) alloc(N, DEF_STD_WIN, DEF_PARAM_STD_WIN, N/2);
   int Nw = Win.nx();
   int Nct = nc(N);
   complex_f *Buff = new complex_f [Nw];
   // double Norm = Win.energy();
   fltarray TabNorm(N);
   
   Data.init();
   TabNorm.init();
   for (l=0; l < Nct; l++)
   {
      k = l*Step;
      for (i=0; i < Nw; i++)  Buff[i] = Trans[l+i*Nct];
      fftn1d(Buff, Nw, True);
      for (i=0; i < Nw; i++) 
      {
         float Re = Buff[i].real();
 	 int IndData = MIN(k,N-1);
         int u = IndData + i - Nw / 2;
	 if ((u >= 0) && (u < N)) 
	 {
	    Re *= Win(i);
	    TabNorm(u) += Win(i)*Win(i);
	    Data(u) += Re;
	 }
      }      
   }
   for (k=0; k < N; k++) if (TabNorm(k) > FLOAT_EPSILON) Data(k) /= TabNorm(k);
   // for (k=0; k < N; k++)  Data(k) /= Norm;
   
   delete [] Buff;
}

/*********************************************************************/

void ST_FFTN::dct_recons(fltarray &Trans, fltarray & Data)
{
   int N = Data.nx();
   int i,k,l;
   if (SizeWin == 0) alloc(N, DEF_STD_WIN, DEF_PARAM_STD_WIN, N/2);
   int Nw = Win.nx();
   int Nct = nc(N);
   dblarray Buff(Nw);
   // double Norm = Win.energy();
   fltarray TabNorm(N);
   
   Data.init();
   TabNorm.init();
   for (l=0; l < Nct; l++)
   {
      k = l*Step;
      for (i=0; i < Nw; i++)  Buff(i) = Trans(l,i); // [l+i*Nct];
      dct1d(Buff, True);
      for (i=0; i < Nw; i++) 
      {
         float Re = Buff(i);
 	 int IndData = MIN(k,N-1);
         int u = IndData + i - Nw / 2;
	 if ((u >= 0) && (u < N)) 
	 {
	    Re *= Win(i);
	    TabNorm(u) += Win(i)*Win(i);
	    Data(u) += Re;
	 }
      }      
   }
   for (k=0; k < N; k++) if (TabNorm(k) > FLOAT_EPSILON) Data(k) /= TabNorm(k);
   // for (k=0; k < N; k++)  Data(k) /= Norm;
}

/*********************************************************************/

void ST_FFTN::recons_direct(complex_f *Trans, fltarray & Data)
{
   int N = Data.nx();
   int i,k;
   complex_f *Buff = new complex_f [N];

   if (Step != 1)
   {
      cout << "Errror: direct reconstruction requires a step equal to 1 ... " << endl;
      exit(-1);
   }
   
   for (k=0; k < N; k++)
   {
      for (i=0; i < N; i++)  Buff[i] = Trans[k+i*N];
      fftn1d(Buff, N, True);
      Data(k) = Buff[k].real();
   }
   delete [] Buff;
}

/*********************************************************************/

void ST_FFTN::recons(complex_f *Trans, fltarray & Data, Bool UseWindow)
{
   if (UseWindow == False) recons_direct(Trans, Data);
   else
   {
      int N = Data.nx();
      if (SizeWin == 0) alloc(N,DEF_STD_WIN, DEF_PARAM_STD_WIN, N);
      recons_win(Trans, Data);
   } 
}

/*********************************************************************/

// void ST_FFTN::spectogram(complex_f *Buff_STF, int Nx, int Ny, fltarray &Spec)
// {
//    int i,j ;
//    if ((Spec.nx() != Nx/2) || (Spec.ny() != Nx))
//    for (i=0; i <= Ny/2; i++) 
//    for (j=0; j < Nx; j++)
//    {
//        float Re = Buff_STF[i*Nx+j].real();
//        float Im = Buff_STF[i*Nx+j].imag();
//        Spec(j,i) = Re*Re+Im*Im;
//    } 
// }

/*********************************************************************/

void ST_FFTN::spectogram(fltarray & Data, fltarray &Spec)
{
   Icomplex_f STFDat;
   int Nx = Data.nx();   
   int Nlt = nl(Nx);
   int Nct = nc(Nx);
   int i,j ;

   STFDat.alloc(Nlt, Nct, "STF");
   complex_f *Buff_STF = STFDat.buffer();
   
   //Short term Fourier transform
   transform(Data, Buff_STF);
   
   // cout << "transform OK" << endl;
   
   if ((Spec.nx() != Nct) || (Spec.ny() != Nlt/2))
                                Spec.alloc(Nct, Nlt/2);

   for (i=0; i < Spec.ny(); i++) 
   for (j=0; j < Spec.nx(); j++)
   {
       float Re = Buff_STF[i*Nct+j].real();
       float Im = Buff_STF[i*Nct+j].imag();
       Spec(j,i) = Re*Re+Im*Im;
   } 
   // cout << "end spec" << endl;
}

/*********************************************************************/

void ST_FFTN::wigner_wille(fltarray & Data, fltarray & TabWV)
{
   int NbrFreq = SizeWin;
   int Tau,i,Nx = Data.nx();
   int Nf2 =  NbrFreq / 2 - 1;
   int StepAna = (Step <= 0) ? 1: Step;
   int IndTime=0;
   if ((TabWV.nx() != Data.nx() / StepAna) && (TabWV.ny() != NbrFreq/2))
                                 TabWV.alloc(Data.nx() / StepAna,  NbrFreq/2);
   complex_f *Buff = new complex_f [NbrFreq];
   for (IndTime = 0; IndTime < TabWV.nx(); IndTime ++)
   {
       int t = IndTime*StepAna;
       int Taumax = MIN (t, (Nx - t - 1));
       if (Taumax > Nf2) Taumax = Nf2;
       for (i = 0; i < NbrFreq; i++) Buff[i] =  complex_f(0.,0.);
       for (Tau = -Taumax; Tau <= Taumax; Tau++)
       {
          i = (Tau >= 0) ? Tau: Tau + NbrFreq;
          Buff[i] =  complex_f(Data(t+Tau) * Data(t-Tau), 0);
       }
       fftn1d(Buff, NbrFreq);
       for (i = 0; i < NbrFreq / 2; i++) 
       {
          if (IndTime >= TabWV.nx())
	  {
	     cout << "Error: nx array bounds write in wigner_wille ... " << endl;
	     cout << "       IndTime = " << IndTime << " TabWV.nx() = "  <<  TabWV.nx() << endl;
	     exit(-1);
	  }
	  else
	  if (i >= TabWV.ny())
	  {
	     cout << "Error: ny array bounds write in wigner_wille ... " << endl;
	     cout << "       i = " << i << " TabWV.ny() = "  <<  TabWV.ny() << endl;
	     exit(-1);
	  }
          TabWV(IndTime,i) = Buff[i].real();
      }
      // IndTime++;
   }
   delete [] Buff;
}

/*********************************************************************/

void ST_FFTN::choi_williams(fltarray & Data, fltarray &TabCW,  
                            type_std_win TWinF, int WindowF_Length, 
			    double sigma, float WinParamF)
{
   int Nx = Data.nx();
   int column, row, time, WinTime2, WinFreq2;
   int taumax, tau, mumin, mumax, mu;
   double normK, normF;
   double spreadfac = 16.0/sigma;
   double R1_real,R2_real;
   int StepAna = (Step <= 0) ? 1: Step;
   int WindowT_Length = SizeWin;
   int N_freq = SizeWin;
   float *WindowT = Win.buffer();
   cfarray Buff(N_freq);
   fltarray TabChoi(Data.nx() / StepAna, N_freq);
   
   if ((TabCW.nx() != Data.nx() / StepAna) && (TabCW.ny() != N_freq/2))
                                 TabCW.alloc(Data.nx() / StepAna,  N_freq/2);
				 
   if (WindowF_Length <= 1) WindowF_Length = N_freq / 4;
   if (WindowF_Length % 2 == 0) WindowF_Length++;
   fltarray TabWindowF(WindowF_Length);
   float *WindowF=TabWindowF.buffer();
   STD_Window STDW;
   STDW.make_win(TabWindowF, TWinF, WinParamF);
   
   WinTime2 = (WindowT_Length - 1) / 2;
   WinFreq2 = (WindowF_Length - 1) / 2;
   normF=WindowF[WinFreq2];
   for(row = 0; row < WindowF_Length; row++) WindowF[row]=WindowF[row]/normF;
 
   int MinNfreq = MIN(N_freq/2,WinFreq2);
   taumax=MinNfreq;
   dblarray TabKernel(taumax, WindowT_Length); 
 
   for(tau=1;tau<=taumax;tau++)
   {
      double Sig2 = -spreadfac*tau*tau;
      for(mu=-WinTime2;mu<=+WinTime2;mu++)
  	  TabKernel(tau-1,WinTime2+mu)=exp(mu*mu/Sig2) *WindowT[WinTime2 + mu];
   }
   for (row = 0; row < N_freq ; row++) Buff(tau) = complex_f(0., 0.);
 
   for (time = 0; time < Nx; time += StepAna)
   {
      column = time / StepAna;
      taumax = MIN((time+WinTime2), (Nx-time-1+WinTime2));
      taumax = MIN(taumax,MinNfreq);
      Buff(tau) = complex_f(Data(time)*Data(time), 0.);
      for (tau = 1; tau <= taumax; tau++)
      {
         R1_real = R2_real = 0.0;
	 mumin=MIN(WinTime2, (Nx-time-1-tau));
	 mumax=MIN(WinTime2,time-tau);
	 normK=0;
	 for(mu=-mumin;mu<=mumax;mu++) normK += TabKernel(tau-1,WinTime2+mu);
	 for(mu=-mumin;mu<=mumax;mu++)
	 {
	    double D2 = Data(time+tau-mu)*Data(time-tau-mu);
            R1_real += D2 * TabKernel(tau-1,WinTime2+mu);
            R2_real += D2 * TabKernel(tau-1,WinTime2-mu);
         }
	 R1_real *= WindowF[WinFreq2+tau]/normK;
	 R2_real *= WindowF[WinFreq2-tau]/normK;
 	 Buff(tau) = complex_f(R1_real, 0);
	 Buff(N_freq-tau) = complex_f(R2_real,0.);
      }

      tau=(int) floor(N_freq/2.);
      if ((time<=Nx-tau-1) && (time>=tau) && (tau<=WinFreq2))
      {
         R1_real = R2_real = 0;
	 mumin = MIN(WinTime2, (Nx-time-1-tau));
	 mumax = MIN(WinTime2,time-tau);
	 normK=0;
	 for(mu=-mumin;mu<=mumax;mu++) normK += TabKernel(tau-1,WinTime2+mu);
 
	 for(mu=-mumin;mu<=mumax;mu++)
	 {
  	    double D2 = Data(time+tau-mu) * Data(time-tau-mu);
            R1_real += D2 * TabKernel(tau-1,WinTime2+mu);
            R2_real += D2 * TabKernel(tau-1,WinTime2-mu);
 	 }
	 R1_real *= WindowF[WinFreq2+tau]/normK;
	 R2_real *= WindowF[WinFreq2-tau]/normK;
 	 Buff(tau) = complex_f(0.5*(R1_real+ R2_real), 0.);
      }
      fftn1d(Buff.buffer(), Buff.nx());
      for (row = 0; row < N_freq; row++)
      {
 	  TabChoi(column,row) = Buff(row).real();
 	  Buff(row) = complex_f(0.,0.);
      }
   }
   
   for (row = 0; row < TabCW.ny(); row++)
   for (time = 0; time < TabCW.nx(); time ++) 
                             TabCW(time,row) = TabChoi(time,row);
}
/***************************************************************************/
