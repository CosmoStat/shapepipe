/******************************************************************************
**                   Copyright (C) 1999 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author:  J.L. Starck
**
**    Date:  18/02/99
**    
**    File:  MatOper.cc
**
*******************************************************************************
**
**    DESCRIPTION  Matrix Operation 
**    ----------- 
**
**
**    PARAMETRES    
**    ----------    
** 
**
******************************************************************************/

#include<cmath>
#include "Array.h"
#include "MatrixOper.h"
#include "NR.h"

/*********************************************************************/

void MatOper::mat_print (dblarray &CorrelMat, const char *Mes)
{
   int i,j;
      
   // print intercorrelation matrix
//    cout << endl << Mes << " :" << CorrelMat.ny() << " " << CorrelMat.nx() << endl;
//    for (i=0; i<  CorrelMat.ny(); i++)
//    {
//  	 for (j=0; j < CorrelMat.nx(); j++) 
// 	 {
// 	    cout.width(8);  // see p343 c++ stroustrup
// 	    cout.fill(' ');
// 	    cout.setf(ios::right,ios::adjustfield);
// 	    //cout.setf(ios::fixed,ios::floatfield);
// 	    cout.setf(ios::scientific,ios::floatfield);
// 	    cout.precision(4);
// 	    cout << CorrelMat(j,i) << " " ;
// 	 }
// 	 cout << endl;
//    }
   cout << endl << Mes << " :" << CorrelMat.ny() << " " << CorrelMat.nx() << endl;
   if (CorrelMat.ny() != 0)
   {
   for (i=0; i<  CorrelMat.ny(); i++)
   {
 	 for (j=0; j < CorrelMat.nx(); j++) 
	 {
            printf("%5.3f ", CorrelMat(j,i));
	 }
	 cout << endl;
   }}
   else  for (j=0; j < CorrelMat.nx(); j++) 
	 {
            printf("%5.3f ", CorrelMat(j));
	 }
	 cout << endl;
}

void MatOper::mat_print (fltarray &CorrelMat, const char *Mes)
{
   int i,j;
      
   // print intercorrelation matrix
   cout << endl << Mes << " :" << CorrelMat.ny() << " " << CorrelMat.nx() << endl;
   for (i=0; i<  CorrelMat.ny(); i++)
   {
 	 for (j=0; j < CorrelMat.nx(); j++) 
	 {
	    printf("%5.3f ", CorrelMat(j,i));
	 }
	 cout << endl;
   }
}
/*********************************************************************/
void MatOper::dblarray2fltarray (dblarray &Ud, fltarray &Uf)
{
  int i,j;
  int Nx = Ud.nx();
  int Ny = Ud.ny();
  Uf.alloc(Nx,Ny);
  for (i=0; i < Nx; i++)
  for (j=0; j < Ny; j++) Uf(i,j) = (float) Ud(i,j);
}


void MatOper::fltarray2dblarray (fltarray &Uf, dblarray &Ud)
{
  int i,j;
  int Nx = Uf.nx();
  int Ny = Uf.ny();
  Ud.alloc(Nx,Ny);
  for (i=0; i < Nx; i++)
  for (j=0; j < Ny; j++) Ud(i,j) = (double) Uf(i,j);
}

/*********************************************************************/

void MatOper::transpose (dblarray &U, dblarray &Ut)
{
  int i,j;
  int Nx = U.nx();
  int Ny = U.ny();
  Ut.alloc(Ny,Nx);
  for (i=0; i < Nx; i++)
  for (j=0; j < Ny; j++) Ut(j,i) = U(i,j);
}
void MatOper::transpose (fltarray &U, fltarray &Ut)
{
  int i,j;
  int Nx = U.nx();
  int Ny = U.ny();
  Ut.alloc(Ny,Nx);
  for (i=0; i < Nx; i++)
  for (j=0; j < Ny; j++) Ut(j,i) = U(i,j);
}

/*********************************************************************/

void MatOper::mat_mult(dblarray &U, dblarray &V, dblarray &R)
//
//  MATRIX multiplication
//  R = U # V
{
  int Q = U.nx();
  int P = U.ny();
  int Q1 = V.nx();
  int P1 = V.ny();
  int i,j,k;
  double Val;
  if (Q != P1)
  {
     cerr << "Error: cannot multiply the matrix ... " << endl;
     exit(-1);
  }
  if ((R.naxis() != 2) || (R.nx() != Q1) || (R.ny() != P)) R.alloc(Q1, P);
  
  for (i = 0; i < Q1; i++)
  for (j = 0; j < P; j++)
  { 
     Val = 0.;;
     for (k = 0; k < Q; k++) 
              Val +=  (double) U(k,j) * (double) V(i,k);
     R(i,j) = Val;
  }     
}
   
void MatOper::mat_mult(fltarray &U, fltarray &V, fltarray &R)
//
//  MATRIX multiplication
//  R = U # V
{
  int Q = U.nx();
  int P = U.ny();
  int Q1 = V.nx();
  int P1 = V.ny();
  int i,j,k;
  float Val;
  if (Q != P1)
  {
     cerr << "Error: cannot multiply the matrix ... " << endl;
     exit(-1);
  }
  if ((R.naxis() != 2) || (R.nx() != Q1) || (R.ny() != P)) R.alloc(Q1, P);
  
  for (i = 0; i < Q1; i++)
  for (j = 0; j < P; j++)
  { 
     Val = 0.;;
     for (k = 0; k < Q; k++) 
              Val +=  U(k,j) *  V(i,k);
     R(i,j) = Val;
  }     
}
   
/*********************************************************************/

void MatOper::apply_mat(dblarray &B, dblarray &Data, dblarray &Result)
{
   int i,j,k;
   dblarray V,R;
   int Nx = Data.axis(1);
   int Ny = Data.axis(2);
   int Nz = Data.axis(3); 
   int Nc = B.axis(2);

   if (B.axis(1) != Data.axis(3))
   {
      cerr << "Error: matrix have bad dimensions. Ny(first matrix) = " << B.axis(2) << endl;
      cerr << "       Second dimension must be equal to " << Nc << endl;
      exit(-1);
   }
     
   V.alloc(1,Nz);
   R.alloc(1,Nc);
   if ((Result.n_elem() == 0) || (Result.naxis() != 3)) Result.alloc(Nx,Ny,Nc);
   else if (Result.axis(1) != Nx) Result.alloc(Nx,Ny,Nc);
   else if (Result.axis(2) != Ny) Result.alloc(Nx,Ny,Nc);
   else if (Result.axis(3) != Nc) Result.alloc(Nx,Ny,Nc);
   
   for (i=0; i < Nx; i++)
   for (j=0; j < Ny; j++)
   {
      for (k=0; k < Nz; k++) V(0,k) = Data(i,j,k);
      mat_mult(B,V,R);
      for (k=0; k < Nc; k++) Result(i,j,k) = R(0,k);
   }
}

void MatOper::apply_mat(fltarray &B, fltarray &Data, fltarray &Result)
{
   int i,j,k;
   fltarray V,R;
   int Nx = Data.axis(1);
   int Ny = Data.axis(2);
   int Nz = Data.axis(3); 
   int Nc = B.axis(2);

   if (B.axis(1) != Data.axis(3))
   {
      cerr << "Error: matrix have bad dimensions. Ny(first matrix) = " << B.axis(2) << endl;
      cerr << "       Second dimension must be equal to " << Nc << endl;
      exit(-1);
   }
     
   V.alloc(1,Nz);
   R.alloc(1,Nc);
   if ((Result.n_elem() == 0) || (Result.naxis() != 3)) Result.alloc(Nx,Ny,Nc);
   else if (Result.axis(1) != Nx) Result.alloc(Nx,Ny,Nc);
   else if (Result.axis(2) != Ny) Result.alloc(Nx,Ny,Nc);
   else if (Result.axis(3) != Nc) Result.alloc(Nx,Ny,Nc);
   
   for (i=0; i < Nx; i++)
   for (j=0; j < Ny; j++)
   {
      for (k=0; k < Nz; k++) V(0,k) = Data(i,j,k);
      mat_mult(B,V,R);
      for (k=0; k < Nc; k++) Result(i,j,k) = R(0,k);
   }
}
	 
/*********************************************************************/

void MatOper::inv_mat_svd(dblarray &Mat, dblarray &InvMat)
// MATRIX inversion: number of lines >= number of columns
// InvMat = Mat^-1
// Mat: IN = matrix with P = Mat.ny  and Q = Mat.nx   (P >= Q)
//  
{
  void dsvdcmp(double **a, int m, int n, double w[], double **v);
  int i,j;
  dblarray  M1,tMat;

  int P = Mat.ny(); // number of lines
  int Q = Mat.nx(); // number of columns
  dblarray U(Q,P);  // composed of eigen value of  B B^t 
  dblarray V(Q,Q);  // composed of eigen vector of B^t B
  dblarray Ut(P,Q); // composed of eigen value of  B B^t 
  dblarray Vt(Q,Q); // composed of eigen vector of B^t B
  
  dblarray Diag(Q,Q);
  Diag.init();
  double **a; // matrix P lines, Q columns
  double **v; // eigenvectors matrix 
  double *w;  // eigenvalues vector
  EV.alloc(Q);
   
  a=dmatrix((long) 1, P, (long)1, Q);   
  v=dmatrix(1, Q, 1, Q);
  w=dvector(1, Q);
  
  for(i=1; i <= P; i++)
  for(j=1; j <= Q; j++)  a[i][j]= Mat(j-1,i-1);
   
  dsvdcmp(a,P,Q,w,v);
  
  for(i=0; i < Q; i++) EV(i) =   w[i+1];
  for(i=0; i < Q; i++) Diag(i,i) = EV(i);
  // double wmax= EV.max();
  // double minT= wmax*EpsEigen;
  // for(i=0; i < Q; i++) if (EV(i) < minT) EV(i) = 0.;

  if (Verbose == True) 
  {
     double CondNumb = condition_nbr();
     cout << "Matrix Condition number = " << CondNumb << endl;
     if (1. / CondNumb < 1e-12) cout << "WARNING: Singular matrix" << endl;
  }

  if (Verbose == True)
  {
     cout << "Eigen values: ";
     for(i=0; i < Q; i++) cout << Diag(i,i) << " ";
     cout << endl;
  }
  for(i=0; i < P; i++)
  for(j=0; j < Q; j++) U(j,i) =  a[i+1][j+1];
  
  for(i=0; i < Q; i++)
  for(j=0; j < Q; j++) V(j,i) =  v[i+1][j+1];
   
  free_dmatrix(a, 1, P, 1, Q); 
  free_dmatrix(v, 1, Q, 1, Q);
  free_dvector(w, 1, Q);

  transpose(U,Ut);
  transpose(V,Vt);
  
  for(i=0; i < Q; i++)  if (Diag(i,i) != 0) Diag(i,i) = 1. / Diag(i,i);

  mat_mult(Diag, Ut, M1);
  mat_mult(V, M1, InvMat);  
  if (Verbose == True) 
  {
     mat_mult(InvMat, Mat, M1);
     mat_print(M1, "B^-1 # B");
  }
}


void MatOper::inv_mat_svd(fltarray &Mat, fltarray &InvMat)
{   
	dblarray Matd;  
	dblarray InvMatd;  
    
	
	fltarray2dblarray( Mat, Matd);
	inv_mat_svd( Matd,  InvMatd);
	dblarray2fltarray( InvMatd, InvMat);
} 


/*********************************************************************/

void MatOper::lin_eq_svd(dblarray & Mat, dblarray & MAT_B, dblarray & MAT_X)
{
  void dsvdcmp(double **a, int m, int n, double w[], double **v);
  void dsvbksb(double **u, double w[], double **v, int m, int n, 
               double b[], double x[]);

  int i,j;
  dblarray  M1,tMat;

  int P = Mat.ny(); // number of lines (m)
  int Q = Mat.nx(); // number of columns (n)
  double **a; // matrix P lines, Q columns
  double **v; // eigenvectors matrix 
  double *w,*b,*x;  // eigenvalues vector

  a=dmatrix((long) 1, P, (long)1, Q);   
  v=dmatrix(1, Q, 1, Q);
  w=dvector(1, Q);
  b=dvector(1, P);
  x=dvector(1, Q);

  for(i=1; i <= P; i++)
  for(j=1; j <= Q; j++)  a[i][j]= Mat(j-1,i-1);
  for(i=1; i <= P; i++) b[i] = MAT_B(i-1);

  dsvdcmp(a,P,Q,w,v);
  double wmax=0.;
  for(i=1; i <= Q; i++) if (wmax < w[i]) wmax = w[i];
  EV.alloc(Q);
  for(i=0; i < Q; i++) EV(i) =   w[i+1];

  double minT= wmax*EpsEigen;
  for(i=1; i <= Q; i++) if (w[i] < minT) w[i] = 0.;

  dsvbksb(a,w,v,P,Q,b,x);

  for(i=1; i <= Q; i++)  MAT_X(i-1) = x[i];

  free_dmatrix(a, 1, P, 1, Q); 
  free_dmatrix(v, 1, Q, 1, Q);
  free_dvector(w, 1, Q);
  free_dvector(b, 1, P);
  free_dvector(x, 1, Q);
}

void MatOper::lin_eq_svd(fltarray & Mat, fltarray & MAT_B, fltarray & MAT_X)
{   
	dblarray Matd;   
	dblarray MAT_Bd; 
	dblarray MAT_Xd; 

	
	fltarray2dblarray( Mat, Matd);
	fltarray2dblarray( MAT_B, MAT_Bd);
	fltarray2dblarray( MAT_X, MAT_Xd);

	lin_eq_svd(Matd,  MAT_Bd,  MAT_Xd);
	
	dblarray2fltarray( Matd, Mat);
	dblarray2fltarray( MAT_Bd, MAT_B);
	dblarray2fltarray( MAT_Xd, MAT_X);
} 



/*********************************************************************/

void MatOper::inv_mat_iter(dblarray &Mat, dblarray &InvMat, int NIter, Bool Verbose)
// MATRIX inversion
// InvMat = Mat^-1
{
  int i,j,k;
  dblarray  M1, M2;
    
  // iterative improvment
  //   Start the iterative inversion from the matrix transposition
  double Eps=0.;
  transpose(Mat,InvMat);
  for (k=0; k < InvMat.n_elem(); k++) Eps += InvMat(k)*InvMat(k);
  Eps = 1. / Eps;
  if (Verbose == True) cout << "Eps = " << Eps << endl;
  for (k=0; k < InvMat.n_elem(); k++) InvMat(k) *= Eps;

  for (k=0; k < NIter; k++)
  {
     mat_mult(InvMat, Mat, M1);
     mat_mult(M1, InvMat, M2);
     for(i=0; i < InvMat.nx(); i++)
     for(j=0; j < InvMat.ny(); j++) InvMat(i,j)  +=  InvMat(i,j) - M2(i,j);
  }
  if (Verbose == True) 
  {
      mat_mult(InvMat, Mat, M1);
      mat_print( InvMat, "Inverse matrix");
      mat_print(M1, "B^-1 # B: after iterating");
  }
}

void MatOper::inv_mat_iter(fltarray &Mat, fltarray &InvMat, int NIter, Bool Verbose)

{
	dblarray Matd;  
	dblarray InvMatd;  
    	
	fltarray2dblarray( Mat, Matd);
	inv_mat_iter(Matd, InvMatd,  NIter,  Verbose);
	dblarray2fltarray( InvMatd, InvMat);
}

/*********************************************************************/






/***************************************************************************/

void AR_PREDICTION::mse_predict(fltarray & Signal, int Np,  dblarray &ArModel, 
                        dblarray &TabErr,  int Step, int ScaleNumber)
{
   int MaxNbrAR = ArModel.nx();
   int FirstPix = (MaxNbrAR+1)*Step;  // POW2(ScaleNumber);
   int NPixUsed = Np - FirstPix;
   int a,i,t;
   double  Pred, ErrPred;
   double Nr = (double) NPixUsed;

   TabErr.alloc(MaxNbrAR+1);
//    cout << " ScaleNumber = " << ScaleNumber << " MaxNbrAR = " << MaxNbrAR << endl;
//    cout << " FirstPix = " << FirstPix << " Np = " << Np << endl;
//    cout << " Signal = " << Signal.nx() << " Np = " << Np << endl;

   ErrPred = 0.;
   for (i = FirstPix; i < Np; i++) ErrPred += Signal(i) * Signal(i);
   TabErr(0) = ErrPred / Nr;

   for (a=0; a < MaxNbrAR; a++)
   {
       ErrPred = 0.;
       for (i = FirstPix; i < Np; i++) 
       {
          Pred = 0.;
          for (t=0; t <= a; t++) 
	  {
	     // int Pos = i-(t+1)*Step;
	     int Pos = i-1-(Step*t);
	     // Pred += ARModel(i-1) * Signal(Pos-1-(Step*(i-1)));
	     if ((Pos < 0) || (Pos >= Signal.nx()))
	     {
	        cout << "Error: Pos = " << Pos << " t = " << t << " a = " << a << " Step = " << Step << endl;
		exit(-1);
	     }
 	     Pred += ArModel(t,a) * Signal(Pos);
	  }
          ErrPred += (Signal(i)-Pred) * (Signal(i)-Pred);
       }
       TabErr(a+1) = ErrPred / Nr;
   }
}

/*********************************************************************/

void AR_PREDICTION::get_best_ar_model(fltarray & Signal, int Np,  
             int &BestNbrAR, fltarray &BestARModel, int Step, int ScaleNumber)
{
   int i,a;
   int MaxNbr_of_AR = MaxNbrAR;
   double Yamma0=0.;

   int FirstPix = MaxNbr_of_AR*Step;
   int NPixUsed = Np - FirstPix;
   if (NPixUsed  < 2)
   {
       MaxNbr_of_AR = Np / (2*Step);
       FirstPix = MaxNbr_of_AR*Step;
       NPixUsed = Np - FirstPix;
   }
   if (NPixUsed < 2)
   {
      cout << "Error: not enough point for AR estimation ..." << endl;
      cout << "  MaxNbr_of_AR = " << MaxNbr_of_AR << " Step = " << Step << endl;
      cout << "  Np = " << Np << " FirstPix = " << FirstPix << endl;
      exit(-1);
   }
   dblarray ArModel(MaxNbr_of_AR,MaxNbr_of_AR);
   dblarray ErrModel(MaxNbr_of_AR);
   dblarray Yamma(MaxNbr_of_AR+1);
   dblarray TabErr(MaxNbr_of_AR+1);
   dblarray TabFunc(MaxNbr_of_AR+1);
   dblarray PenalFunc(MaxNbr_of_AR+1);

   // Initialization for AR(1)
   for (i = FirstPix; i < Np; i++)
   {
      Yamma(0) += Signal(i)*Signal(i-1);
      Yamma0 += Signal(i)*Signal(i);
   }
   if (Yamma0 > 0)  ArModel(0,0) = Yamma(0) / Yamma0;
   ErrModel(0)  = Yamma0 *(1- ArModel(0,0)*ArModel(0,0));

   // from AR(2) to AR(MaxNbr_of_AR)
   for (a=1; a < MaxNbr_of_AR; a++)
   {
      double PhiYamma=0.;

      // Yamma(a) calculation
      for (i = FirstPix; i < Np; i++) 
                 Yamma(a) += Signal(i)*Signal(i-1-(a)*Step);

      // Phi_a,a calculation
      for (i = 0; i < a; i++)  PhiYamma += ArModel(i,a-1)*Yamma(a-i-1);
      if (ABS(ErrModel(a-1)) > 0) 
                   ArModel(a,a) = (Yamma(a) - PhiYamma) / ErrModel(a-1);

       // Phi_i,a calculation
       for (i = 0; i < a; i++)
             ArModel(i,a) = ArModel(i,a-1) - ArModel(a,a)*ArModel(a-i-1,a-1);

       // Error calculation
       ErrModel(a) =  ErrModel(a-1) * (1- ArModel(a,a)*ArModel(a,a));

//        cout << " MODELE AR(" << a+1 << ")" << endl;
//        cout << "   AR = ";
//        for (i = 0; i <= a; i++) cout << ArModel(i,a) << " ";
//        cout << endl << "   Err = " << ErrModel(a) << endl;

    }

   mse_predict(Signal, Np, ArModel, TabErr, Step, ScaleNumber);
   double Min=0.;
   int NumAr=0;
   // double Sig = Signal.sigma();
   double AIC, BIC, AICC;
   for (a=0; a <= MaxNbr_of_AR; a++)
   {
      // double Nr = (double) Np / (double ) Step;
      double Nr =  (double) Np / (double ) POW2(ScaleNumber+1);
      switch (AR_Order_Detect)
      {
         case AR_ORDER_AIC:  
              AIC = 2.*a / Nr;
              PenalFunc(a) =  AIC; break;
         case AR_ORDER_AICC: 
              AICC = (float)(a+Nr)  / (float)(Nr  - a  - 2.);
              PenalFunc(a) =  AICC ; break;
         case AR_ORDER_BIC:  
              BIC = a *log(Nr) / Nr;
              PenalFunc(a) =  BIC; break;
         default: cout << "Error: unknown criterion ... " << endl;
                  exit(-1);
      }
      TabFunc(a) = log(TabErr(a)) + PenalFunc(a);
 
      // cout << a << ": Err = " << TabErr(a) << " Log = " << log(TabErr(a)) << " Pen = " << 2*(a+1) << endl;
      // cout <<  "  Sig = " << Sig << " AIC = " << AIC << " BIC  = " <<  BIC << " AICC  = " << AICC  << endl;
      // cout << a << ": H = " << AIC << " AIC = " <<  TabFunc(a) << endl;
      if (a == 0) Min = TabFunc(a);
      else if (Min > TabFunc(a))
      {
         Min = TabFunc(a);
         NumAr = a;
      }
   }
   BestNbrAR = NumAr;
   if ( BestARModel.nx() != BestNbrAR) BestARModel.alloc(BestNbrAR);
   get_ar_model(Signal, Np, BestNbrAR, 0, BestARModel, Step);
   // for (a=0; a < BestNbrAR; a++) BestARModel(a) = (float)  ArModel(a, BestNbrAR-1); 

//    for (a=0; a < MaxNbr_of_AR; a++)
//    {
//        cout << " MODELE AR(" << a << ")" << endl;
//        cout << "   AR = ";
//        if (a > 0) for (i = 0; i < a; i++) cout << ArModel(i,a-1) << " ";
//        cout << endl << "   ErrPred = " << TabErr(a) << " Penal = " <<  PenalFunc(a) <<  endl;
//        cout <<  "   NPixUsed = " << NPixUsed << " Func = " << TabFunc(a) << endl;
//    }

//    a = BestNbrAR;
//    cout << " MODELE AR(" << a << ")" << endl;
//    cout << "   AR = ";
//    cout <<  " Estimated AR order = " << NumAr << endl;
//   if (a > 0) for (i = 0; i < a; i++) cout << ArModel(i,a-1) << " ";
//   cout << endl;
}

/*********************************************************************/

void AR_PREDICTION::get_ar_model(fltarray & Signal, int Np, int NbrAR, int MaxNbrTraining, 
                  fltarray & ARModel, int Step)
{
   int i,j,IndData;
   int t0 = Np-1-PredDistance;
   // int N = Signal.nx();
   int NbrTraining;
   NbrTraining = t0 - NbrAR*Step;
   if ((MaxNbrTraining > 0) && (NbrTraining >  MaxNbrTraining)) NbrTraining = MaxNbrTraining;

   if (NbrTraining < 0)
   {
       NbrAR = t0 / Step;
       NbrTraining = t0 - NbrAR*Step;
       if (NbrTraining < 0)
       {
       cout << "Error: Number of scales is to high ... " << endl;
       exit(-1);
       }
   }
   // cout << "NbrTraining = " << NbrTraining << endl;
   // cout << " Step = " << Step  << "t0 = " << t0 << " NbrAR = " << NbrAR  << " NbrTraining = " << NbrTraining << endl;

   // Learning values: A X = B
   dblarray MAT_X(1,NbrAR);
   dblarray MAT_A(NbrAR, NbrTraining);
   dblarray InvMAT_A;
   dblarray MAT_B(1,NbrTraining);

   ARModel.alloc(NbrAR);
 
   MAT_X.init();
   for (i=0; i < NbrTraining; i++)
   {
      for (j=0; j < NbrAR; j++)
      {
         IndData = t0 - i - 1 - j *Step;
         if (IndData < 0)
         {
            cout << endl;
            cout << "ERROR: IndData = " << IndData << endl;
            exit(-1);
         }
         MAT_A (j,i) = Signal(IndData);
      }
      IndData = t0 - i + PredDistance;
      MAT_B(i) = Signal(IndData);
   } 
   MatOper CO;
   CO.lin_eq_svd(MAT_A,  MAT_B,  MAT_X);
 
   for (j=0; j < NbrAR; j++) ARModel(j) = (float) MAT_X(j);
}

/***************************************************************************/

void fit1d_pol_deg2(fltarray &Data, dblarray  & MAT_X, dblarray &Weight, int N, int NL) 
{
   int i,Np = N;
   int NLast = NL;
   if (N == 0) Np = Data.nx();
   if (NLast == 0) NLast = Np;
   dblarray MAT_A(3,NLast);
   dblarray MAT_B(1,NLast);

   for (i=0; i < NLast; i++)
   {
      MAT_A(0,i) = i*i*Weight(i);
      MAT_A(1,i) = i*Weight(i);
      MAT_A(2,i) = 1*Weight(i);
      MAT_B(i) = Data(Np-NLast+i)*Weight(i);
   }
   
   // Learning values: A X = B
   MAT_X.init();
   MatOper CO;
   CO.lin_eq_svd(MAT_A,  MAT_B,  MAT_X);
}


/***************************************************************************/
 
void fit1d_pol_deg2(fltarray &Data, dblarray  & MAT_X, int N, int NL)
// Calcule the polynome of degree 2 fitting the Nlast points
// of the data
// Data(i) = a x^2 + b x + 1    
//           and i = [Np-NLast,Np-1]
// Np: IN = Number of pixels in Data
// Data: IN = input data
// Nlast: IN = number of points in Data to be used
// DistPred: IN = Calculate the prediction at a distance DistPred+1 from
//                the last pixel (DistPred=0) for the next pixel prediction)
//                and return this value      
{
   int i,Np = N;
   int NLast = NL;
   if (N == 0) Np = Data.nx();
   if (NLast == 0) NLast = Np;
   dblarray MAT_A(3,NLast);
   dblarray MAT_B(1,NLast);

   for (i=0; i < NLast; i++)
   {
      MAT_A(0,i) = i*i;
      MAT_A(1,i) = i;
      MAT_A(2,i) = 1;
      MAT_B(i) = Data(Np-NLast+i);
   }

   // Learning values: A X = B
   MAT_X.init();
   MatOper CO;
   CO.lin_eq_svd(MAT_A,  MAT_B,  MAT_X);
}

/***************************************************************************/
 
void fit1d_pol_deg1(fltarray &Data, dblarray  & MAT_X, int N, int NL)
{
   int i,Np = N;
   int NLast = NL;
   if (N == 0) Np = Data.nx();
   if (NLast == 0) NLast = Np;
   dblarray MAT_A(2,NLast);
   dblarray MAT_B(1,NLast);
   for (i=0; i < NLast; i++)
   {
      MAT_A(0,i) = i;
      MAT_A(1,i) = 1;
      MAT_B(i) = Data(Np-NLast+i);
   }
   // Learning values: A X = B
   MAT_X.init();
   MatOper CO;
   CO.lin_eq_svd(MAT_A,  MAT_B,  MAT_X);
}

/***************************************************************************/

void tendancy_est(fltarray & Data, fltarray & Tend, int WindowSize, int N)
{
   int i,k;
   int Np = N;
   if (N == 0) Np = Data.nx();
   dblarray MAT_X(1,3);  
   if (Tend.nx() != Np) Tend.alloc(Np);
   
   for (i=0; i<Np; i++)
   {
       dblarray W;
       int Nx = WindowSize;
       int Pmin = MAX(0,i-WindowSize/2);
       int Pmax = MIN(Np-1,i+WindowSize/2);
       Nx = Pmax-Pmin+1;
       W.alloc(Nx);
       for (k=0; k<Nx; k++) 
       {
         double D = MAX(Pmax-i,i-Pmin);
         double x = (Pmin + k - i) / D;
         W(k) = sqrt(1. - x*x);
       }
       fit1d_pol_deg2(Data, MAT_X, W, Pmax+1, Nx);
       double P, x = i-Pmin;
       P = MAT_X(0)*x*x+ MAT_X(1)*x +MAT_X(2);      
       Tend(i) = (float) P;
   }
} 

/*********************************************************************/

void autocor1d(fltarray  & Data, fltarray &TabAuto, int NBShift, int N)
{
   int i,s;
   double XY,X2,Y2;
   int Np = N;
   if (N == 0) Np = Data.nx();
   if (TabAuto.nx() != NBShift) TabAuto.alloc(NBShift);
   TabAuto(0) = 1;
   for (s=1; s < NBShift; s++)
   {
      XY = X2 = Y2 = 0.;
      for (i=s; i < Np; i++)
      {
         XY += Data(i) * Data(i-s);
         X2 += Data(i) * Data(i);
         Y2 += Data(i-s)*Data(i-s);
      }
      TabAuto(s) = (float) ( XY / (sqrt(X2)*sqrt(Y2)));
   }
}

/*********************************************************************/
