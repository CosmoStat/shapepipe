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
**    Date:  29/03/00 
**    
**    File:  WT1D_FFT.cc
**
*******************************************************************************
**
**    DESCRIPTION   routines concerting filters used by the 1D wavelet 
**    -----------   transform which need the FFT   
**
**
******************************************************************************/ 


#include "WT1D_FFT.h"


/****************************************************************************/

int WT1D_FFT::pyr_tab_pos(int NbrScale, int N)
{
    int s,Ns=N;
    int TotSize=0;
    if (TabSize.nx() != NbrScale) TabSize.reform(NbrScale);
    if (TabPos.nx() != NbrScale)  TabPos.reform(NbrScale);

    for (s = 0; s < NbrScale; s++)
    {
       TotSize += Ns;
       TabSize(s) = Ns;
       Ns = (Ns+1)/2;
    }
    TabPos(NbrScale-1) = 0;
    for (s =  NbrScale-2; s >= 0; s--) TabPos(s) = TabPos(s+1) + TabSize(s+1);

//     for (s = 0; s < NbrScale; s++)
//     {
//         cout << " Band " << s+1 << " size = " << TabSize(s) << " Pos = " << TabPos(s) << endl;
//     }
//     cout << " TotSize = " << TotSize << endl;
    return TotSize;
}

/****************************************************************************/

void WT1D_FFT::transform (fltarray &Signal, fltarray &Trans, int Nbr_Plan)
{
    int N = Signal.nx();
  
    if (Ima_b_cf.nx() != N) Ima_b_cf.resize (N);
    if (Buff.nx() != N) Buff.resize (N);
    // for (i = 0; i < N; i ++) Ima_b_cf(i) = complex_f( Signal(i),0.);
    // fft1d (Buff, Ima_b_cf, 1);
    
    FFT.fftn1d(Signal, Ima_b_cf.buffer());
    transform_cf(Ima_b_cf, Trans, Nbr_Plan);
}

/****************************************************************************/

void WT1D_FFT::transform_cf(cfarray & Data, fltarray &Trans, int Nbr_Plan)
{
    int i,s, N = Data.nx();
    int Depi, Ni, N_s=N;
    float Coef;

    // cout << "WT FFT transform " << Nbr_Plan << " N = " << N << endl;
    if (Ima_h_cf.nx() != N) Ima_h_cf.resize (N);
    if (Buff.nx() != N) Buff.resize (N);
    switch (TypeWT)
    {
        case WT_PAVE_FFT:
            if ((Trans.nx() != N) || (Trans.ny() !=  Nbr_Plan))
                                                  Trans.reform(N, Nbr_Plan);
            for (s = 0; s < Nbr_Plan - 1; s ++)
            {
                for (i = 0; i < N; i ++) Ima_h_cf(i) = Data(i);
                wave_cf_mult (Data, Ima_h_cf, s, N, True);
		FFT.fftn1d(Ima_h_cf, Buff, True);
                // fft1d (Ima_h_cf, Buff, -1);
                for (i = 0; i < N; i ++) Trans(i,s) = Buff(i).real();
            }
            // fft1d (Data, Buff, -1);
	    FFT.fftn1d(Data, Buff, True);
            for (i = 0; i < N; i ++) Trans(i,Nbr_Plan-1) = Buff(i).real();
            break;
        case WT_PYR_FFT_DIFF_RESOL:
        case WT_PYR_FFT_DIFF_SQUARE:
            PyrPixNbr = pyr_tab_pos(Nbr_Plan, N);
            if (Trans.nx() !=  PyrPixNbr) Trans.reform(PyrPixNbr);
            for (s = 0; s < Nbr_Plan - 1; s ++)
            {
                Depi = pyr_pos(s);
                Ni = pyr_size(s);
             // cout << "Band " << s+1 << "  Depi = " <<  Depi << " Ni = " << Ni << endl;
               if (Ni != N_s)
                {
                   cout << "Errror: Ni and Ns are not equal ... " << endl;
                   exit(-1);
                }
                Coef = (float)(N) / (float)(Ni);              
                Ima_h_cf.resize (N_s);
                for (i = 0; i < N_s; i ++) Ima_h_cf(i) = Data(i);
                wave_cf_mult(Data, Ima_h_cf, s, N, True);
                Buff.resize (N_s);
                // fft1d (Ima_h_cf, Buff, -1);
		FFT.fftn1d(Ima_h_cf, Buff, True);
                for (i = 0; i < Ni; i ++) Trans(i+Depi) = Buff(i).real()/Coef;
                under_sampling_2 (Data);
                N_s = (N_s+1)/2;
            }
            s = Nbr_Plan - 1;
            Depi = pyr_pos(s);
            Ni = pyr_size(s);
            Coef = (float)(N) / (float)(Ni);              
         // cout << "  Depi = " <<  Depi << " Ni = " << Ni << " coef = " << Coef << endl;
            if (Ni != N_s)
            {
                cout << "Errror: Ni and Ns are not equal .... " << endl;
                exit(-1);
            }
            Buff.resize (Ni);
            // fft1d (Data, Buff, -1);
	    FFT.fftn1d(Data, Buff, True);
            for (i = 0; i < Ni; i++) Trans(i+Depi) = Buff(i).real()/Coef;
	    Data.resize (N);
            break;
        default:
            fprintf (stderr,"Error: proc. wave_cf_transform.\n");
            fprintf (stderr,"This transform is not computed this procedure\n");
            break;
    }
}

/****************************************************************************/

void WT1D_FFT::recons(fltarray &Trans, fltarray &Signal, int Nbr_Plan)
{
    int N = Signal.nx();
    int s,i;

    if (TypeWT == WT_PAVE_FFT)
    {
        Signal.init();
        for (s = Nbr_Plan -1; s >= 0; s--)
        for (i = 0; i < N; i ++) Signal(i) += Trans(i,s);
    }
    else 
    {
        if (Ima_b_cf.nx() != N) Ima_b_cf.resize (N);
        if (Buff.nx() != N) Buff.resize (N);         
        recons_cf(Trans, Ima_b_cf, Nbr_Plan);
        // fft1d(Ima_b_cf, Buff,-1);
	FFT.fftn1d(Ima_b_cf, Buff, True);
        for (i=0; i < N; i++) Signal(i) = Buff(i).real(); 
    }
}

/****************************************************************************/

void WT1D_FFT::recons_cf(fltarray &Trans, cfarray & Data, int Nbr_Plan)
{
    int N = Data.nx();
    Bool Odd = (N % 2 == 0) ? False: True;
    int N_s, s,i,Depi;
    float Coef;

    if (Ima_h_cf.nx() != N) Ima_h_cf.resize (N);
    if (Data.nx() != N) Data.resize (N);
    if (Buff.nx() != N) Buff.resize (N);
    // cout << "NN = " << N << endl;
    switch (TypeWT)
    {
        case WT_PAVE_FFT:
           {
              fltarray Signal(N);
              Signal.init();
              for (s = Nbr_Plan -1; s >= 0; s--)
              for (i = 0; i < N; i ++) Signal(i) += Trans(i,s);
              for (i=0; i <  N_s; i++) Buff(i)=complex_f(Signal(i),0.);
              // fft1d(Buff, Data, 1);
	      FFT.fftn1d(Buff, Data);
           }
            break;
        case WT_PYR_FFT_DIFF_RESOL:
        case WT_PYR_FFT_DIFF_SQUARE:
            PyrPixNbr = pyr_tab_pos(Nbr_Plan, N);
            if (Trans.nx() !=  PyrPixNbr)
            {
               cout << "Error: PyrPixNbr and Trans.nx() Pb ... " << endl;
               cout << "  Trans.nx() = " << Trans.nx() << endl;
	       cout << "  PyrPixNbr = " << PyrPixNbr << endl;
	       cout << "  N = " << N << "  Nbr_Plan = " << Nbr_Plan << endl;
	       exit(-1);
            }
            s = Nbr_Plan - 1;
            Depi = pyr_pos(s);
            N_s = pyr_size(s);
            Coef = (float)(N) / (float)(N_s);              
            Data.resize (N_s);
            Buff.resize (N_s);  
            for (i=0; i <  N_s; i++) Buff(i) = complex_f(Trans(i+Depi)*Coef,0.);       
            // fft1d(Buff, Data, 1);
	    FFT.fftn1d(Buff, Data);
            for (s = Nbr_Plan - 2; s >= 0; s --)
            {
                Depi = pyr_pos(s);
                N_s = pyr_size(s);
		Odd = (N_s % 2 == 0) ? False: True;
               // cout << "  Band " << s+1 << " N = " << Data.nx() <<  endl;
                Ima_h_cf.resize (N_s);
                Buff.resize (N_s);  
                Coef = (float)(N) / (float)(N_s); 
               //  cout << "  Depi = " <<  Depi << " Ni = " << N_s << " coef = " << Coef << endl;

                for (i=0; i <  N_s; i++) Buff(i) = complex_f(Trans(i+Depi)*Coef,0.);
                // fft1d(Buff, Ima_h_cf, 1);
		FFT.fftn1d(Buff, Ima_h_cf);
                over_sampling_2 (Data, Odd);
                wave_cf_mult (Data, Ima_h_cf, s, N, False);
                for (i=0; i < Data.nx() ; i++) Data(i) += Ima_h_cf(i);
            }
	    // cout  << " Nend = " << Data.nx() <<  endl;
            // Data.resize (N);
            break;
        default:
            fprintf (stderr,"Error: proc. wave_cf_transform.\n");
            fprintf (stderr,"This transform is not computed this procedure\n");
            break;
    }
}

/****************************************************************************/

void WT1D_FFT::under_sampling_2 (cfarray &Imag)
{
    int N = Imag.nx();
    int N_2 = (N+1)/2;
    int i,i1,i0,Depi;
    cfarray Buff;
    Buff.resize (N_2);
    Depi = N/4;
    for (i0 = 0; i0 < N_2; i0++)
    {
       i1 = i0 + Depi;
       Buff (i0) = Imag (i1);
    }
    Imag.resize(N_2);
    for (i = 0; i <  N_2; i++) Imag(i) = Buff(i);
}

/****************************************************************************/

void WT1D_FFT::over_sampling_2 (cfarray &Imag, Bool Odd)
{
    int N = Imag.nx();
    int N_2 = (Odd == True) ? 2*N-1: 2*N;
    int i,i1,i0,Depi;
    cfarray Buff;

    Buff.resize(N_2);
    Depi = N_2/4;
    for (i0 = 0; i0 < N_2; i0++) Buff (i0) = complex_f(0.,0.);
    for (i1 = 0; i1 < N; i1++)
    {
       i0 = i1 + Depi;
       Buff (i0) = Imag (i1);
    }
    Imag.resize(N_2);
    for (i = 0; i < N_2; i++) Imag(i) = Buff(i);
}

/****************************************************************************/

void WT1D_FFT::wave_cf_mult (cfarray &Ima_b, cfarray &Ima_h, int s, int N,
                             Bool Down)
{
    int i,u,Dep,N_s = Ima_b.nx();
    Dep = (int) (pow((double)2., (double) s) + 0.5);
    for (i = 0; i < N_s; i++)
    {
        u = Dep * (i - N_s/2);
        if ((u+N/2 < 0) || (u+N/2 >= N)) Ima_b(i) = complex_f ((float)0.,(float)0.);
        else
        {
            if (Down == True)
            {
                 Ima_b (i) *= filter (WT_H, (float) u, N);
                 Ima_h (i) *= filter (WT_G, (float) u, N);
            }
            else
            {
                 Ima_b (i) *= filter (WT_H_TILDE, (float) u,N);
                 Ima_h (i) *= filter (WT_G_TILDE, (float) u,N);
            }

        }
    }
}

/****************************************************************************/

void WT1D_FFT::create_filter (fltarray &Filter, type_wt_filter Which)
{
    int i;
    float u;
    int N = Filter.nx();
    for (i = 0; i < N; i++)
    {
        u = (float) i - (float) N / 2.;
        Filter(i) = filter (Which,u,N);
    }
}

/****************************************************************************/

double WT1D_FFT::scaling_function (float u, int N)
{
    double Fc, Scaling_Function;
    Fc = Freq_Coup * (float) N; // we fixe the cut-off frequency to N / 2
            //  . we multiply by 3/2 in order that: T(0) = 1
            //  . we multiply F by 2./Fc (b3_spline(2) = 0) 
            //    in order that: T(Fc) = 0
    Scaling_Function = 3. / 2. * b3_spline ((double)(2. * u / Fc));
    return (Scaling_Function);
}

/***************************************************************************/

double WT1D_FFT::filter_h (float Freq_u, int N)
{
    double Filter_H;
    int N_2 = N / 2;
    double u2 = 2. * Freq_u;
    double Freq1,Freq2;
    if ((u2 < - N_2 ) || (u2 >= N_2)) Filter_H = 0.;
    else
    {
         Freq1 = scaling_function (Freq_u,  N);
         Freq2 = scaling_function (u2, N);  
         if (ABS(Freq1) < FLOAT_EPSILON)  Filter_H = 0.;
         else  Filter_H = Freq2 / Freq1;
    }
    return (Filter_H);
}

/***************************************************************************/   

double WT1D_FFT::filter_g (float Freq_u, int N)
{
    double Filter_G=0.;
    double Filter_H;
     Filter_H = filter_h (Freq_u, N);
     switch (TypeWT)
     {
         case WT_PYR_FFT_DIFF_RESOL:
         case WT_PAVE_FFT: 
                      Filter_G = 1 - Filter_H;
                      break;
         case WT_PYR_FFT_DIFF_SQUARE:
                      Filter_G  = sqrt (1. - Filter_H * Filter_H);
                      break;
        default:
            fprintf (stderr, "Error: bad wavelet type in filter_g\n");
            exit (-1);
            break;
     }

    return (Filter_G);
}

/***************************************************************************/

double WT1D_FFT::filter_h_tilde (float Freq_u, int N)
{
    double Filter_H_Tilde=0, Filter_H, Filter_G, Den;

    switch (TypeWT)
    {
        case WT_PYR_FFT_DIFF_RESOL: 
        case WT_PAVE_FFT:
                Filter_H = filter_h (Freq_u,N);
                Filter_G = filter_g (Freq_u,N);
                Den = Filter_H * Filter_H + Filter_G * Filter_G;
                if (Den < FLOAT_EPSILON) Filter_H_Tilde = 0.;
                else Filter_H_Tilde = Filter_H / Den;
            break;
        case WT_PYR_FFT_DIFF_SQUARE:
            Filter_H_Tilde = filter_h(Freq_u,N);
            break;
        default:
            fprintf (stderr, "Error: bad wave in filter_h_tilde\n");
            exit (-1);
            break;
    }
    return (Filter_H_Tilde);
}

/***************************************************************************/

double WT1D_FFT::filter_g_tilde (float Freq_u, int N)
{
    double  Filter_G_Tilde=0, Filter_H, Filter_G, Den;
    switch (TypeWT)
    {
        case WT_PYR_FFT_DIFF_RESOL:
        case WT_PAVE_FFT:
                Filter_H = filter_h(Freq_u,N);
                Filter_G = filter_g(Freq_u,N);
                Den = Filter_H * Filter_H + Filter_G * Filter_G;
                if (Den < FLOAT_EPSILON) Filter_G_Tilde = 0.;
                else Filter_G_Tilde = Filter_G / Den;
                break;
        case WT_PYR_FFT_DIFF_SQUARE:
                Filter_G_Tilde = filter_g(Freq_u,N);
            break;
        default:
            fprintf (stderr, "Error: bad wavelet type in filter_g_tilde\n");
            exit (-1);
            break;
    }
    return (Filter_G_Tilde);
}
   
/***************************************************************************/   

double WT1D_FFT::filter_wavelet (float Freq_u, int N)
{
    double Filter_Psi=0;
    double u2;
    double Freq1,Freq2;

    Freq1 = scaling_function (Freq_u,N);
    u2 = Freq_u / 2.;
    Freq2 = scaling_function (u2,N);
    switch (TypeWT)
    {
        case WT_PYR_FFT_DIFF_RESOL: 
        case WT_PAVE_FFT:
            Filter_Psi = Freq2 - Freq1;
            break;
        case WT_PYR_FFT_DIFF_SQUARE:
            Filter_Psi = Freq2*Freq2 - Freq1*Freq1;
            break;
        default:
            fprintf (stderr, "Error: bad wavelet type in filter_wavelet\n");
            exit (-1);
            break;
    }
    return (Filter_Psi);
}

/***************************************************************************/

float WT1D_FFT::filter (type_wt_filter Which_Filter, float u,  int N)
{
    double Val_Return=0.;
    switch (Which_Filter)
    {
        case WT_PHI : Val_Return = scaling_function (u, N); break;
        case WT_H :         Val_Return = filter_h (u, N); break;
        case WT_H_TILDE :   Val_Return =  filter_h_tilde (u, N); break;
        case WT_G :         Val_Return = filter_g (u,N); break;
        case WT_G_TILDE :   Val_Return = filter_g_tilde (u,N); break;
        case WT_PSI :          Val_Return = filter_wavelet (u,N); break;
        default:
              fprintf (stderr, "Error: bad filter in filter\n");
              exit (0);
              break;
    }
    assert ((Val_Return >= -FLOAT_EPSILON) && (Val_Return < 2));
    return ((float) Val_Return);
}

/***************************************************************************/
/***************************************************************************/

// static void print_info_cf (cfarray &Image, int s)
// {
//     int Nl = Image.nl();
//     int Nc = Image.nx();
//     char Buff_Name[100];
//     sprintf (Buff_Name, "%s (%d, %d)", "FFT Ima_Buff", Nl, Nc);
//     Ifloat Ima_Buff(Nl, Nc, Buff_Name);
//     Ima_Buff = real(Image);
//     cout << Image.Name_Obj << "(" << s << ") : " << endl; 
//     cout << "       Sigma reel = " << sigma (Ima_Buff) << endl;
//     cout << "       Min reel = " << min (Ima_Buff) << endl;
//     cout << "       Max reel = " << max (Ima_Buff) << endl;
//     cout << "       Moy reel = " << average (Ima_Buff) << endl;
// }

/****************************************************************************/
 
