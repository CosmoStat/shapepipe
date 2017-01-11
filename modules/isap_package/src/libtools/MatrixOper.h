#ifndef _MATRIXOPER_H_
#define _MATRIXOPER_H_

#define MAX_ITER_INV_MAT_ITER 1000

class MatOper {

  public:
    dblarray EV; // eigen values vector
                 // calculated by lin_eq_svd or inv_mat_svd
    double condition_nbr()
    {
       if (EV.nx () == 0) 
       {
           cout << "Error: condition number cannot be calculated ... " << endl;
	   exit(-1);
       }
       double wmin = EV.min();
       double CondNumb = (wmin > 0) ? EV.max() / EV.min() : 0;
       if (Verbose == True)
	    cout << "Matrix Condition number = " << CondNumb << endl; 
       return CondNumb;
    }
    		 
    Bool Verbose;
    double EpsEigen; // when Cond. nbr  <  EpsEigen
                     // the reconstruction must be performed
		     // used by lin_eq_svd
    MatOper () {Verbose=False;EpsEigen=1e-6;};
    ~MatOper () {};
    
	void dblarray2fltarray (dblarray &Ud, fltarray &Uf);
    void fltarray2dblarray (fltarray &Uf, dblarray &Ud);
    // convert the type of an array

	
    void mat_print (dblarray & Mat, const char *Mes);
    void mat_print (fltarray & Mat, const char *Mes);
    // Print a matrix on stdout

    void transpose (dblarray &U, dblarray &Ut);
    void transpose (fltarray &U, fltarray &Ut);
    // transpose a matrix

    void mat_mult(dblarray &U, dblarray &V, dblarray &R);
    void mat_mult(fltarray &U, fltarray &V, fltarray &R);
    // Caculate the matrix product: R  = U # V

    void apply_mat(dblarray &B, dblarray &Data, dblarray &Result);
    void apply_mat(fltarray &B, fltarray &Data, fltarray &Result);
    // Apply the matrix B to all vectors of Data
    // B: IN = 2D array
    // Data: IN = 3D array
    // Result: OUT: 3D array  Result(i,j,*) = B # Data(i,j,*)

    void inv_mat_svd(dblarray &Mat, dblarray &InvMat);
    void inv_mat_svd(fltarray &Mat, fltarray &InvMat);
    // Inverse a matrix by singular value decomposition
    // MATRIX inversion when the number of lines >= number of columns
    // Mat: IN = matrix with P = Mat.ny  and Q = Mat.nx   (P >= Q)
    // InvMat: OUT  = Mat^-1

    void  inv_mat_iter(dblarray &Mat, dblarray &InvMat, int NIter=MAX_ITER_INV_MAT_ITER, Bool Verbose=False);
    void  inv_mat_iter(fltarray &Mat, fltarray &InvMat, int NIter=MAX_ITER_INV_MAT_ITER, Bool Verbose=False);
    // iterative matrix inversion

    void lin_eq_svd(dblarray & Mat, dblarray & MAT_B, dblarray & MAT_X);
    void lin_eq_svd(fltarray & Mat, fltarray & MAT_B, fltarray & MAT_X);
    // Resolve the linear equation: A X = B
};
    
	
	
	
	
#define NBR_AR_DETECT_ORDER_METHOD 3
enum  ar_order_detect_type {AR_ORDER_AIC,AR_ORDER_AICC,AR_ORDER_BIC};
      // type of criterion which is used to determine the best AR model
inline const char * StringARDetectType (ar_order_detect_type  type)
{
    switch (type)
    {
        case  AR_ORDER_AIC:  return ("AIC"); 
        case  AR_ORDER_AICC: return ("AICC");
        case  AR_ORDER_BIC:  return ("BIC");
        default:  return ("Error: undefined method.");
    }
}	

class AR_PREDICTION
{
   void mse_predict(fltarray & Signal, int Np,  dblarray &ArModel, 
                    dblarray &TabErr,  int Step=1, int ScaleNumber=0);
   // Calculate the prediction error for a set of AR models
   // Signal: IN = input signal
   // Np: IN = consider only the Np first in the data set (t0 = Np)
   // ArModel: IN = set of AR model 
   //               ArModel(*,a) = AR coefficients of ath AR model.
   // Step: IN = Step to use between two succesive points

  public:
   Bool Verbose;       // if true, print information

   ar_order_detect_type AR_Order_Detect;
                       // indicates which criterion to use to find 
                       // the AR order. Default is BIC

   int MaxNbrAR;       // AR maximum order value. Default is 10

   int PredDistance;  // Distance from the last point to be predicted
                      // if the input data have Np points Data(0..Np-1)
		      // prediction value at Np+PredDistance
		      // by default, PredDistance=0 (next pixel)
		      
   AR_PREDICTION() {Verbose=False;MaxNbrAR=10;
                    AR_Order_Detect=AR_ORDER_BIC;PredDistance=0;}
   ~AR_PREDICTION() {}

   void get_ar_model(fltarray & Signal, int Np, int NbrAR, int MaxNbrLearning, 
                    fltarray & ARModel, int Step=1);
   // Find the AR model of order NbrAR
   // Signal: IN = input signal
   // Np: IN = consider only the Np first in the data set (t0 = Np)
   // NbrAR: IN = AR order
   // MaxNbrLearning: use only at maximum MaxNbrLearning points for the AR model
   //                 estimation: Signal(Np-MaxNbrLearning:Np-1) is used
   //                 if MaxNbrLearning == 0, no limitations on the number
   //                 of points.
   // Step: IN = Step to use between two succesive points
   // ARModel: OUT = AR model: ArModel(0:NbrAR-1)= AR coefficients
   
   void get_best_ar_model(fltarray & Signal, int Np,
                int &BestNbrAR, fltarray &BestARModel, 
                int Step=1, int ScaleNumber=0);
   // Calculate the best order model for a signal
   // Signal: IN = input signal
   // Np: IN = consider only the Np first in the data set (t0 = Np)
   // BestNbrAR: OUT = best AR order
   // BestARModel: OUT = best AR Model: 
   //              BestARModel(0..BestNbrAR-1) = AR coefficients
   // Step: IN = Step to use between two succesive points

   inline float ar_prediction(fltarray & Signal, int NbrAR, int Pos,
                              fltarray &ARModel, int Step=1)
   {
       float Pred=0.;
       for (int i=1; i<= NbrAR; i++) 
       {
          int Ind = Pos-1-(Step*(i-1));
	  if (Ind < 0)
	  {
	      cout << "Error: Ind = " << Ind << " Step = " << Step << " Pos = " << Pos << " i = " << i << endl;
	      exit(-1);
	  }
          // Pred += ARModel(i-1) * Signal(Pos-Step*i);
	  Pred += ARModel(i-1) * Signal(Ind);
          // cout << i << " AR " << ARModel(i-1) << " S " << Signal(Pos-Step*i) << " " << Pred << endl;
       }
       return Pred;
   }
};

void autocor1d(fltarray  & Data, fltarray &TabAuto, int NBShift, int Np=0);
// Calculate the autocorrelation function of a 1D signal
// Data: IN = input Data
// NBShift: IN = Number of lags used in the AF calculation
// TabAuto: OUT = Array[0..NBShift-1]
// Np: IN = only the first Np pixels are used for the AF calculation

void tendancy_est(fltarray & Data, fltarray & Tend, int WindowSize, int Np=0);
// Estimate the tendancy in a one dimensional signal by fitting for each pixel
// a polynome of order 2 in a weighted window centered at the pixel position.
// Data: IN = input Data
// Np: IN = only the first Np pixels are used for the tendancy calculation
//          by defaullt, Np = Data.nx()
// Tend: OUT = Array[0..Np-1]
// WindowSize: IN = Window size in which the fitting is done.

void fit1d_pol_deg1(fltarray &Data, dblarray  & MAT_X, int Np=0, int NLast=0);
// Linear regression on pixels Data[Np-NLast:Np-1])
// by default all pixels are used in the linear regression (Np=NLast=Data.nx())
// Data: IN = input Data
// MAT_X: OUT = Array[2]
//        FIT:   Data(i) = MAT_X[0]*i + MAT_X[1]

void fit1d_pol_deg2(fltarray &Data, dblarray  & MAT_X, int Np=0, int NLast=0);
// Polynomial fit on pixels Data[Np-NLast:Np-1])
// by default all pixels are used in the linear regression (Np=NLast=Data.nx())
// Data: IN = input Data
// MAT_X: OUT = Array[2]
//        FIT:   Data(i) = MAT_X[0]*i^2 + MAT_X[1]*i + MAT_X[2]

void fit1d_pol_deg2(fltarray &Data, dblarray  & MAT_X, dblarray &Weight, int Np=0, int NLast=0); 
// Weighted polynomial fit on pixels Data[Np-NLast:Np-1])
// by default all pixels are used in the linear regression (Np=NLast=Data.nx())
// Data: IN = input Data
// Weight: IN = Array[0:NLast-1] 
//                Weight[i] = weight of pixel Data[Np-NLast+i]

#endif
