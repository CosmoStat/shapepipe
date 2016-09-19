/******************************************************************************
**                   Copyright (C) 2007 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck - Jerome Bobin
**
**    Date:  12/10/07
**    
**    File:  mr_gmca.cc
**
*******************************************************************************
**
**    DESCRIPTION  Generalized Morphological Component Analysis // unknown unmber of sources
**    ----------- 
**                 
**    Usage: mr_gmca options cube output
**
******************************************************************************/

#include "GMCA.h"

/****************************************************************************/

void GMCA::run_gmca(fltarray &TabCannels, fltarray & TabSource)
{ 
 // Declaring some local variables

   int Nmax = Max_GMCA_Iter;
   float KThrd;
   if (Inpainting == False) KThrd = 10;  // For K-Mad filtering  -- UN PEU ELEVE ?
   if (Inpainting == True) KThrd = 25; 
   float KThrd_min = KMin;  // May be an option ?
   float DiffKThrd = (KThrd-KThrd_min)/(Nmax-1);//,lsparse1=0,L1Crit=0;  // For K-Mad filtering - defines the speed of convergence
   fltarray WhiteMat,RecSource,Data,DataRec;    
   KStart = 1;
      
   init_mixmat(TabCannels);  // May use a random initialization
   
   //if (MatNbrScale1D < 2) MatNbrScale1D = get_nbr_scale (NbrCannels);
   FilterAnaSynt *PtrFAS;
   FAS.alloc(F_MALLAT_7_9);
   PtrFAS = &FAS;
   MR_Mat.alloc(NbrCannels, TO1_MALLAT, MatNbrScale1D, PtrFAS, NORM_L2, False);
      
   if ((Verbose == True) && (L1_Matrix_Constraint == True)) 
        cout << " WT1D of the mixing matrix: nbr of scales = " << MatNbrScale1D << endl;
        
   if (L1_Matrix_Constraint == True) KStart = (int) NbrCannels/pow((float) 2, (float) MatNbrScale1D-1);
   
  // cout << "KStart" << KStart << endl;
//   {
//    mat_wt1d_trans(TabCannels);
//    if (SVConst == True)
//    {
//    	pca(TabCannels,SingVec); // Compute the Singular Vectors
//    }
//  }
//   
   NbrCoef = TabCannels.nx();
   
   
   if (PositiveMatrix == True)
   {
        if (Verbose == True) cout << "Positivity of the mixing matrix" <<endl;
   }
   
   if (Inpainting == True)
   {
        if (Verbose == True)  cout << "Use mask and inpaint" << endl;
        recons_data(TabCannels,Data);
        fits_write_fltarr ("DataIn.fits",Data);
   }
   
   
   float SigmaData = TabCannels.sigma();
  
// Main loop
   if (Verbose == True)  cout << "Number of sources : " << NbrSources << endl;
   if (Verbose == True)  cout << "Last K-Mad : " << KMin << endl;

   for (int i=0; i < Nmax; i++)
   {
      // Update the sources :

      // Normalizing the mixing matrix                              
            
      Normal_Mixmat();  
               
      Update_Sources(TabCannels,TabSource, KThrd);
            
      Sort_Sources(TabSource);  // A retirer eventuellement

      if (PositiveSource == True) 
      {
            transrecons_sources(TabSource,RecSource);
            positive_cube_constraint(RecSource);
            transform_sources(RecSource,TabSource);
      }                                         
                        
      // Update the mixing matrix :
            
       Update_MixingMat(TabCannels,TabSource,KThrd);
       
       // Compute the L2 error :
       
      L2err =  l2_error(TabCannels,TabSource);
      L2err = -20.*log10(L2err/SigmaData + 1e-100);
      
      if (PositiveMatrix == True) positive_mixmat_constraint();
      
      if (SVConst == True) PCAMixMat(); 
                                      
       // Applying the mask
       
       if (Inpainting == True)
      {
         reform_data(TabSource,DataRec);
         fits_write_fltarr("DataRec.fits",DataRec);
         apply_mask(Data,Mask,DataRec); 
         transform_data(DataRec,TabCannels);
      } 
     
     
                                      
      // Decrease the threshold
      
      KThrd -= DiffKThrd;

     if (Verbose == True) cout << "Iteration : "<< i <<" - L2-error : " <<  L2err <<endl; // See the threshold decrease
      
   }

   //fits_write_fltarr ("cc_EstMixmat.fits", MixingMat);
   //fits_write_fltarr ("cc_EstSources.fits", TabSource);
   
   if (Inpainting == True)
   {
    transrecons_sources(TabSource,RecSource);
    fits_write_fltarr ("xx_InpaintedSources.fits", RecSource);
    reform_data(TabSource,DataRec);
    fits_write_fltarr ("xx_InpaintedData.fits", DataRec);
    cout << "Inpaint Fin" << endl;
   }
   
   //if (L1_Matrix_Constraint == True)
   //{
   // mat_wt1d_recons(MixingMat);
   //}
   
   // Calc_InvMixingMat(TabSource);
   
   //L2err =  l2_error(TabCannels,TabSource);
    
}

/****************************************************************************/

void GMCA::PCAMixMat()
{
	fltarray Mat;
	fltarray tMat;
	fltarray Res;
	fltarray Q(NbrCannels,NbrSources);
	int i,k;
	for (i=0;i<NbrCannels;i++) for (k=0;k<NbrSources;k++) Q(i,k) = SingVec(i,k);
	
	MAT.transpose(Q,tMat);
    MAT.mat_mult(Q,tMat,Mat);
    MAT.mat_mult(Mat,MixingMat,Res);
	MixingMat = Res;
}

/****************************************************************************/

void GMCA::positive_mixmat_constraint()
{
 int i,k,s;
 fltarray V(NbrSources),W(NbrSources),Q(NbrCannels);
 
 // Reconstruction
 
 if (L1_Matrix_Constraint == True) mat_wt1d_recons(MixingMat);
 
 // Positivity
 
 for (i=0;i<NbrSources;i++)    
 {
 for (k=0;k<NbrCannels;k++) Q(k) = MixingMat(k,i);  // Take care of the sign of the max ?
 V(i) = Q.max();
 W(i) = Q.min();
 }
 
 for (i=0;i<NbrCannels;i++)    
 for (k=0;k<NbrSources;k++) 
 {
    if (abs(V(k)) < abs(W(k))) MixingMat(i,k) = -MixingMat(i,k);
    if (MixingMat(i,k) < 0) MixingMat(i,k) = 0;  // Take care of the sign of the max ?
 }
 RecMixingMat = MixingMat;
 
    if (L1_Matrix_Constraint == True)  // Transformation if necessary
    {
        
        
        fltarray Vect(NbrCannels);
        
        for (k=0; k < NbrSources; k++)
        {
            for (i=0;i<NbrCannels; i++) Vect(i) = MixingMat(i,k);  
            
            MR_Mat.transform(Vect);
            
            int ind=0;
            
            for (s=0; s < MR_Mat.nbr_band (); s++)
            {
                for (i=0; i < MR_Mat.size_scale_np (s); i++) 
                {
                    MixingMat(ind+i,k) = MR_Mat(s,i);
                    
                }
                
                ind += MR_Mat.size_scale_np (s);
            }
        }
    }
    
}

/****************************************************************************/

void GMCA::positive_cube_constraint(fltarray &Data)
{
    int i,j,z;
    fltarray V(Data.nz()),W(Data.nz()),Frame(Data.nx(),Data.ny());
    
    for (z=0;z < Data.nz();z++)
    {
    for (i=0;i < Data.nx();i++)
    for (j=0;j < Data.ny();j++) Frame(i,j) = Data(i,j,z);
    V(z) = Frame.max(); 
    W(z) = Frame.min();
    }
    
    
    for (z=0;z < Data.nz();z++)
    for (i=0;i < Data.nx();i++)
    for (j=0;j < Data.ny();j++)
    {
     if (abs(V(z)) < abs(W(z))) Data(i,j,z) = -Data(i,j,z);
     if (Data(i,j,z) < 0) Data(i,j,z) = 0;
    }
    
}

/****************************************************************************/

void GMCA::ortho_hypcube(fltarray &TabSource)
{
	int i,k,l;
	int NbrCoef = TabSource.nx();
	fltarray P(NbrCoef*NbrCannels,NbrSources);
	
	int Indo = 0;
	
	for (i=0;i<NbrSources;i++)
	{
		Indo=0;
		for (k=0;k<NbrCannels;k++)
		{
			for (l=0;l<NbrCoef;l++) P(Indo,i)=MixingMat(k,i)*TabSource(l,i);
			Indo++;
			
		}	
		
	}
	
	fltarray MMV,tQ,Q;
	
	ortho_sources_svd(P,Q);
	
	MAT.transpose(Q,tQ);
	MAT.mat_mult(tQ,TabSource,MMV);
	
	TabSource = MMV;
	
}

/****************************************************************************/

void GMCA::pca(fltarray &Mat,fltarray &VMat)
// Compute the Singular Vectors of the covariance matrix of Mat
// Mat = TabCannels

{
  void dsvdcmp(double **a, int m, int n, double w[], double **v);
  int i,j;
  fltarray  M1,tMat;
  fltarray CovMatSources;
  float MMean;
  
  // 1) retirer la moyenne
  
  fltarray NMat = Mat;
  fltarray MeanMat(NMat.nx());
  for(i=0; i < NMat.nx(); i++)
  { 
  		for(j=0; j < NMat.ny(); j++) MeanMat(j) =   Mat(i,j);
  	    MMean = MeanMat.mean();
  	    // cout << MMean << endl;
  	    for(j=0; j < NMat.ny(); j++) NMat(i,j) =   Mat(i,j) - MMean;
  }
  
  //cout << NMat.nx() << endl;
  
  MAT.transpose(NMat,tMat);
  MAT.mat_mult(NMat,tMat,CovMatSources);

  //cout << "Orthogonalization" << endl;
  //cout << CovMatSources.nx() << endl;
  //cout << CovMatSources.ny() << endl;

  int P = CovMatSources.ny(); // number of lines
  int Q = CovMatSources.nx(); // number of columns
  fltarray U(Q,P);  // composed of eigen value of  B B^t 
  fltarray V(Q,Q);  // composed of eigen vector of B^t B
  fltarray Ut(P,Q); // composed of eigen value of  B B^t 
  fltarray Vt(Q,Q); // composed of eigen vector of B^t B
  fltarray EV;
  
  fltarray Diag(Q,Q);
  Diag.init();
  double **a; // matrix P lines, Q columns
  double **v; // eigenvectors matrix 
  double *w;  // eigenvalues vector
  EV.alloc(Q);
   
  a=dmatrix((long) 1, P, (long)1, Q);   
  v=dmatrix(1, Q, 1, Q);
  w=dvector(1, Q);
  
  for(i=1; i <= P; i++)
  for(j=1; j <= Q; j++)  a[i][j]= CovMatSources(j-1,i-1);
   
  dsvdcmp(a,P,Q,w,v);
  
  for(i=0; i < Q; i++) EV(i) =   w[i+1]; // Eigenvalues
  for(i=0; i < Q; i++) Diag(i,i) = EV(i);
  // double wmax= EV.max();
  // double minT= wmax*EpsEigen;
  // for(i=0; i < Q; i++) if (EV(i) < minT) EV(i) = 0.;

  
  for(i=0; i < P; i++)
  for(j=0; j < Q; j++) U(j,i) =  a[i+1][j+1];
  
  for(i=0; i < Q; i++)
  for(j=0; j < Q; j++) V(j,i) =  v[i+1][j+1];
   
  free_dmatrix(a, 1, P, 1, Q); 
  //free_dmatrix(NMat); 
  free_dmatrix(v, 1, Q, 1, Q);
  free_dvector(w, 1, Q);
  
  VMat = V;
}

/****************************************************************************/

void GMCA::ortho_sources_svd(fltarray &Mat,fltarray &VMat)
// Orthogonalization by svd

{
  void dsvdcmp(double **a, int m, int n, double w[], double **v);
  int i,j;
  fltarray  M1,tMat;
  fltarray CovMatSources;

  MAT.transpose(Mat,tMat);
  MAT.mat_mult(Mat,tMat,CovMatSources);

  //cout << "Orthogonalization" << endl;
  //cout << CovMatSources.nx() << endl;
  //cout << CovMatSources.ny() << endl;

  int P = CovMatSources.ny(); // number of lines
  int Q = CovMatSources.nx(); // number of columns
  fltarray U(Q,P);  // composed of eigen value of  B B^t 
  fltarray V(Q,Q);  // composed of eigen vector of B^t B
  fltarray Ut(P,Q); // composed of eigen value of  B B^t 
  fltarray Vt(Q,Q); // composed of eigen vector of B^t B
  fltarray EV;
  
  fltarray Diag(Q,Q);
  Diag.init();
  double **a; // matrix P lines, Q columns
  double **v; // eigenvectors matrix 
  double *w;  // eigenvalues vector
  EV.alloc(Q);
   
  a=dmatrix((long) 1, P, (long)1, Q);   
  v=dmatrix(1, Q, 1, Q);
  w=dvector(1, Q);
  
  for(i=1; i <= P; i++)
  for(j=1; j <= Q; j++)  a[i][j]= CovMatSources(j-1,i-1);
   
  dsvdcmp(a,P,Q,w,v);
  
  for(i=0; i < Q; i++) EV(i) =   w[i+1]; // Eigenvalues
  for(i=0; i < Q; i++) Diag(i,i) = EV(i);
  // double wmax= EV.max();
  // double minT= wmax*EpsEigen;
  // for(i=0; i < Q; i++) if (EV(i) < minT) EV(i) = 0.;

  
  for(i=0; i < P; i++)
  for(j=0; j < Q; j++) U(j,i) =  a[i+1][j+1];
  
  for(i=0; i < Q; i++)
  for(j=0; j < Q; j++) V(j,i) =  v[i+1][j+1];
   
  free_dmatrix(a, 1, P, 1, Q); 
  free_dmatrix(v, 1, Q, 1, Q);
  free_dvector(w, 1, Q);
  
  fltarray MMat;
  fltarray tV;
  MAT.transpose(V,tV);
  MAT.mat_mult(tV,Mat,MMat);
  Mat = MMat;
  VMat = V;
  
  //MAT.transpose(Mat,tMat);
  //MAT.mat_mult(Mat,tMat,CovMatSources);
  //int k;
  
  //if (Verbose == True)
  //{
  //   cout << "CovMat: ";
  //   for(i=0; i < P; i++) for(k=0; k < Q; k++) cout << CovMatSources(i,k) << " ";
  //   cout << endl;
  //}
}

/****************************************************************************/

void GMCA::ortho_spectra_svd()
// Orthogonalization by svd

{
  void dsvdcmp(double **a, int m, int n, double w[], double **v);
  int i,j;
  fltarray  M1,tMat;
  fltarray CovMatSources;

  MAT.transpose(MixingMat,tMat);
  MAT.mat_mult(MixingMat,tMat,CovMatSources);

  //cout << "Orthogonalization" << endl;
  //cout << CovMatSources.nx() << endl;
  //cout << CovMatSources.ny() << endl;

  int P = CovMatSources.ny(); // number of lines
  int Q = CovMatSources.nx(); // number of columns
  fltarray U(Q,P);  // composed of eigen value of  B B^t 
  fltarray V(Q,Q);  // composed of eigen vector of B^t B
  fltarray Ut(P,Q); // composed of eigen value of  B B^t 
  fltarray Vt(Q,Q); // composed of eigen vector of B^t B
  fltarray EV;
  
  fltarray Diag(Q,Q);
  Diag.init();
  double **a; // matrix P lines, Q columns
  double **v; // eigenvectors matrix 
  double *w;  // eigenvalues vector
  EV.alloc(Q);
   
  a=dmatrix((long) 1, P, (long)1, Q);   
  v=dmatrix(1, Q, 1, Q);
  w=dvector(1, Q);
  
  for(i=1; i <= P; i++)
  for(j=1; j <= Q; j++)  a[i][j]= CovMatSources(j-1,i-1);
   
  dsvdcmp(a,P,Q,w,v);
  
  for(i=0; i < Q; i++) EV(i) =   w[i+1]; // Eigenvalues
  for(i=0; i < Q; i++) Diag(i,i) = EV(i);
  // double wmax= EV.max();
  // double minT= wmax*EpsEigen;
  // for(i=0; i < Q; i++) if (EV(i) < minT) EV(i) = 0.;

  
  for(i=0; i < P; i++)
  for(j=0; j < Q; j++) U(j,i) =  a[i+1][j+1];
  
  for(i=0; i < Q; i++)
  for(j=0; j < Q; j++) V(j,i) =  v[i+1][j+1];
   
  free_dmatrix(a, 1, P, 1, Q); 
  free_dmatrix(v, 1, Q, 1, Q);
  free_dvector(w, 1, Q);
  
  fltarray MMat;
  fltarray tV;
  MAT.transpose(V,tV);
  MAT.mat_mult(tV,MixingMat,MMat);
  MixingMat = MMat;
  
  //MAT.transpose(Mat,tMat);
  //MAT.mat_mult(Mat,tMat,CovMatSources);
  //int k;
  
  //if (Verbose == True)
  //{
  //   cout << "CovMat: ";
  //   for(i=0; i < P; i++) for(k=0; k < Q; k++) cout << CovMatSources(i,k) << " ";
  //   cout << endl;
  //}
}


/****************************************************************************/

void GMCA::force_disjspec()
{
    int i,k,m;
    fltarray V(NbrSources),Q(NbrCannels),P(NbrSources);
        
   // int MaxSpecIndex = (int) NbrCannels*(1-1/((float) MatNbrScale1D));
    
   // cout << MaxSpecIndex << endl;
   
   for (i=0;i<NbrSources;i++) 
    {
        for (m=0;m<NbrCannels;m++) Q(m) = ABS(MixingMat(m,i));  //  Mieux eviter pour les basses echelles
    
        float My = Q.max();
    
        P(i) = 1;
        
        if (My > ErrRate) P(i) = My;
        
    }   
   
        
    for (i=0;i<NbrCannels;i++) 
    {
        for (k=0;k<NbrSources;k++) V(k) = ABS(MixingMat(i,k))/P(k);  //  Mieux eviter pour les basses echelles
        
        float Mx = V.max();
    
        for (k=0;k<NbrSources;k++)
        {
            float STh = ABS(MixingMat(i,k))/P(k);
            
            if (STh <  Mx) MixingMat(i,k) = 0;
        }
    }   
}

/****************************************************************************/

void GMCA::force_disjsources(fltarray &TabSource)
{
    int i,k;
    int NbrCoef = TabSource.nx();
    fltarray V(NbrSources);

        
    for (i=0;i<NbrCoef;i++) 
    {
        for (k=0;k<NbrSources;k++) V(k) = ABS(TabSource(i,k));  //  Mieux eviter pour les basses echelles
    
        float Mx = V.max();
    
        for (k=0;k<NbrSources;k++) if (ABS(TabSource(i,k)) < Mx) TabSource(i,k) = 0;
    }   
}

/****************************************************************************/
void GMCA::subs_mean(fltarray &Data)
{
	int i,k;
	int NbC = Data.ny();
	int NbS = Data.nx();
	
	fltarray V(NbS);
	
	for (i=0;i<NbC;i++)
	{
		for (k=0;k<NbS;k++) V(k) = Data(k,i);
		
		float Me = V.mean();
				
		for (k=0;k<NbS;k++) Data(k,i)-=Me;	
		
	}
	
}

/****************************************************************************/

void GMCA::print_mixmatnorm()
{
    int i,k;
    fltarray V(NbrCannels);
        
    for (i=0;i<NbrSources;i++) 
    {
        for (k=0;k<NbrCannels;k++) V(k) = MixingMat(k,i);  //  Mieux eviter pour les basses echelles
    
        cout << "Std of the columns : " << V.sigma() << endl;
    }   
}

/****************************************************************************/

void GMCA::Calc_InvMixingMat(fltarray &TabSource)
{
	fltarray LeftA;
	fltarray TabMMTemp,IndexA(NbrSources);
    int i,k,NbNzS = NbrSources,NLin_Mat = MixingMat.axis(1),NCol_Mat = MixingMat.axis(2);

    // Evalve the indices of the non-zero columns of the mixing matrix :
    Where_MixingMat(MixingMat,IndexA,NbNzS);

     TabMMTemp.alloc(NLin_Mat,NbNzS);
     
     int Indi = 0;
     
     for (i=0; i < NCol_Mat;i++)
     {
        if (IndexA(i) == 1)
        { 
                for (k=0; k < NLin_Mat;k++) TabMMTemp(k,Indi) = MixingMat(k,i);
                Indi++;
        }
     }

     // Compute the pseudo-inverse:

    leftpi_mat(TabMMTemp, LeftA);
	TabMMTemp.alloc(MixingMat.axis(2),MixingMat.axis(1));
			
	for (i=0; i < NbNzS;i++)
    for (k=0; k < MixingMat.axis(1);k++) TabMMTemp(i,k) = LeftA(i,k);
    
    if (NbNzS < MixingMat.axis(2))
    for (i=NbNzS; i < MixingMat.axis(2);i++)
    for (k=0; k < MixingMat.axis(1);k++) TabMMTemp(i,k) = 0;
	
    InvMixingMat = TabMMTemp;
    	
}

/****************************************************************************/

void GMCA::Update_MixingMat(fltarray &TabData,fltarray &TabSource , float &KThrd)
{
    // Definitions
    fltarray IndexS(NbrSources);
    fltarray TabMMTemp,TabSTemp,RightS;
    int i,k,NbNzS = NbrSources,NLin_Mat = MixingMat.axis(1),NbrCoef = TabSource.nx();
    
    
    
    
     // Evalve the indices of the non-zero columns of the mixing matrix :
    Where_Sources(TabSource, IndexS,NbNzS);

    // cout << "Number of active sources : " << NbNzS << endl;
    
    // Extract the non-zero columns :
     
     TabSTemp.alloc(NbrCoef,NbNzS);  /// TAILLE OK 
     int Indi = 0;
     
     for (i=0; i < NbrSources;i++)
     {
        if (IndexS(i) == 1)
        { 
                for (k=0; k < NbrCoef;k++) TabSTemp(k,Indi) = TabSource(k,i);
                Indi++;
        }
     }
     
     // Compute the pseudo-inverse:
     
     leftpi_mat(TabSTemp, RightS); 
     
     // Estimate the corresponding sources :
     
     MAT.mat_mult(TabData,RightS,TabMMTemp);  // We get the transpose of TabMMTemp
     
     // Thresholding :
     if (L1_Matrix_Constraint == True) 
     {
     	fltarray TempM;
        mat_wt1d_trans(TabMMTemp);     
     	HT_MixMat(TabMMTemp,KThrd);
     	MAT.transpose(TabMMTemp,TempM);  // Confusion transpose ou pas --- a controler
     	mat_wt1d_recons(TempM);
     	MAT.transpose(TempM,TabMMTemp);
     }
     
     // Put the columns
     Indi = 0;
               
     for (i=0; i < NbrSources;i++)
     {
        if (IndexS(i) == 1)
        { 
                for (k=0; k < NLin_Mat;k++) MixingMat(k,i) = TabMMTemp(Indi,k);// As we got the transpose of TabMMTemp
                Indi++;
        }
     }     
          
     
    
 }
 
/****************************************************************************/
void GMCA::init_TabSource(fltarray &TabSource,int &NbrSamples)
{
    int i,k;
    fltarray Temp(NbrSamples,NbrSources);
    
    for (i=0;i<NbrSamples;i++) for (k=0;k<NbrSources;k++) Temp(i,k)=0;
    
    TabSource = Temp;
}

/****************************************************************************/
void GMCA::print_mat(fltarray &Mat)
{
 int Nx = Mat.nx();
 int Ny = Mat.ny();
 int i,k;
  
 for (i=0;i<Nx;i++)
 {
  for (k=0;k<Ny;k++) cout << Mat(i,k) << endl;
  cout << "\n" << endl;  
 }   
    
}

/****************************************************************************/

void GMCA::Sort_Sources(fltarray &TabSource)
{
    void indexx(int n, double arr[], int indx[]);
    int i,k,l;
    dblarray V(NbrSources);
    dblarray VarA(NbrCannels);
    dblarray VarS(TabSource.nx());
    intarray VInd(NbrSources);
    
    
    for (i=0;i<NbrSources;i++)
    {
        for (k=0;k<NbrCannels;k++)
        {
        VarA(i) = MixingMat(k,i);
        }
       for (l=0;l<TabSource.nx();l++)
       {
        VarS(l) = TabSource(l,i);
       }
       
     V(i) = VarA.energy()*VarS.energy();   
    }
    
    // Then compute the sorted indices :

   double *P_val_l_tri = V.buffer()-1;
   int  *P_ind_tri = VInd.buffer()-1;
   indexx(NbrSources, P_val_l_tri, P_ind_tri);
   
   // Sorting the mixing matrix :
   
   fltarray MMTemp = MixingMat;
   
    for (i=0;i<NbrSources;i++)
    {
     int Indo = VInd(i)-1;
     for (k=0;k<NbrCannels;k++)
     {
      MixingMat(k,i) = MMTemp(k,Indo);
     }
    }
  }
  
/****************************************************************************/

void GMCA::RetrieveSourcesFromMixmat(fltarray &TabData,fltarray & TabSource)
{
    // Definitions
    fltarray TabMMTemp,TabSTemp,LeftA;

     // Compute the pseudo-inverse:
     
     leftpi_mat(MixingMat, LeftA);
     InvMixingMat = LeftA;

     // Estimate the corresponding sources :
     fltarray NLeftA=LeftA;
     
     MAT.transpose(NLeftA,LeftA);
     MAT.mat_mult(LeftA,TabData,TabSTemp);
     
     TabSource = TabSTemp;
}
  
  
/****************************************************************************/

void GMCA::Update_Sources(fltarray &TabData,fltarray &TabSource ,  float &KThrd)
{
    // Definitions
    fltarray IndexA(NbrSources);
    fltarray TabMMTemp,TabSTemp,LeftA;
    int i,k,NbNzS = NbrSources,NLin_Mat = MixingMat.axis(1),NCol_Mat = MixingMat.axis(2),NbrCoef = TabSource.nx();
    
     // Evalve the indices of the non-zero columns of the mixing matrix :
    Where_MixingMat(MixingMat,IndexA,NbNzS);

     // Extract the non-zero columns :
     TabMMTemp.alloc(NLin_Mat,NbNzS);
     
     int Indi = 0;
     
     for (i=0; i < NCol_Mat;i++)
     {
        if (IndexA(i) == 1)
        { 
                for (k=0; k < NLin_Mat;k++) TabMMTemp(k,Indi) = MixingMat(k,i);
                Indi++;
        }
     }
     
     // Compute the pseudo-inverse:
     
     leftpi_mat(TabMMTemp, LeftA);

     InvMixingMat = LeftA;

     // Estimate the corresponding sources :
     fltarray NLeftA=LeftA;
     
     MAT.transpose(NLeftA,LeftA);
     MAT.mat_mult(LeftA,TabData,TabSTemp);
     
     // if (GetRelax == True)
     //{
     // fltarray TabTemp = TabData - ...	
     //	MAT.mat_mult(LeftA,TabData,TabSTemp);
     //}
          
     // Thresholding :
     HT_Sources(TabSTemp, KThrd);

     // Put the sources
     Indi = 0;
          
     for (i=0; i < NCol_Mat;i++)
     {
        if (IndexA(i) == 1)
        { 
                //fltarray V(NbrCoef);
                for (k=0; k < NbrCoef;k++) TabSource(k,i) = TabSTemp(k,Indi);
                //for (k=0; k < NbrCoef;k++) V(k) = TabSTemp(k,Indi);
                Indi++;
        }
     }     
}

/****************************************************************************/

void GMCA::Where_MixingMat(fltarray &Mat, fltarray &IndexA, int &NbNzS)
{ 
    int i,k,NLin_Mat = Mat.axis(1),NCol_Mat = Mat.axis(2); 
    fltarray V(NLin_Mat);
    
    NbNzS = 0;
    
    for (i=0; i < NCol_Mat; i++)
    {    	         
    	 for (k=0;k<NLin_Mat; k++) V(k) = Mat(k,i);  
    	 
    	 float sigma  = V.sigma();
    	     	 
    	 IndexA(i) = 0;
    	 
    	 if (sigma > ErrRate) 
    	       {
    	               IndexA(i) = 1;
                       NbNzS++;
    	       }
    }
        
}

/****************************************************************************/

void GMCA::Where_Sources(fltarray &TabSource, fltarray &IndexS,int &NbNzS)
{ 
    int i,k,NbrCoef = TabSource.nx(),Nyd = TabSource.ny();
    fltarray V(NbrCoef);
        
    NbNzS = 0;
        
    for (i=0; i < Nyd; i++)
    {
    	 for (k=0;k<NbrCoef; k++) V(k) = TabSource(k,i); 
    	 
        float sigma  = V.sigma();
        
    	     	     	     	 
    	 IndexS(i) = 0;
    	 
    	 if ((sigma > ErrRate)&&(i > NbrKnownColumn-1))
    	       {
    	               IndexS(i) = 1;
                       NbNzS++;
    	       }
    	    	
    }
        
}

/****************************************************************************/

void GMCA::Normal_Mixmat()
{ 
    int i,k;
    float Norm_V;
    fltarray V(NbrCannels);
        
    for (i=0; i < NbrSources; i++)
    {
    	 Norm_V = 0;
    	 
    	 for (k=0;k<NbrCannels; k++) V(k) = MixingMat(k,i);  // Compute the norm of each column
    	 
    	 Norm_V = V.sigma();
    	     	     	 
    	 if (Norm_V > ErrRate) for (k=0;k<NbrCannels; k++) MixingMat(k,i) = MixingMat(k,i)/Norm_V; // Normalize each column
    	 if (Norm_V < ErrRate) for (k=0;k<NbrCannels; k++) MixingMat(k,i) = 0;
    }
    
}


/****************************************************************************/

float GMCA::calc_L1norm(fltarray &TabSource) 
{
    int i,k;
    float L1n = 0;
    
    for (i=0;i<TabSource.nx();i++) for (k=0;k<TabSource.ny();k++) L1n += ABS(TabSource(i,k));
    return L1n;
}

/****************************************************************************/

float GMCA::calc_L0norm(fltarray &TabSource) 
{
    int i,k;
    float L0n = 0;
    
    for (i=0;i<TabSource.nx();i++) for (k=0;k<TabSource.ny();k++) if(ABS(TabSource(i,k)) > ErrRate) L0n++;
    return L0n;
}

/****************************************************************************/

float GMCA::l2_error(fltarray & TabData, fltarray & TabSource)
{ 
    fltarray TabChannel;
    fltarray TranspMixingMat;
    fltarray Mat= MixingMat;
    MAT.transpose(Mat,TranspMixingMat);
    MAT.mat_mult(TranspMixingMat, TabSource, TabChannel);
    TabChannel -= TabData;
    return TabChannel.sigma();
}

/****************************************************************************/

float GMCA::CountResiThrd(fltarray & TabData, fltarray & TabSource, float & Thrd)
{ 
    fltarray TabChannel;
    fltarray TranspMixingMat;
    fltarray Mat= MixingMat;
    MAT.transpose(Mat,TranspMixingMat);
    MAT.mat_mult(TranspMixingMat, TabSource, TabChannel);
    TabChannel -= TabData;
    float count=0;
    int i,j;
    for (i=0; i<TabChannel.nx(); i++) for (j=0; j<TabChannel.ny(); j++) if (abs(TabChannel(i,j)) > Thrd) count++;
    return count;
}

/****************************************************************************/

float GMCA::MyKurt(fltarray & TabData)
{ 
    int Nx = TabData.nx();
    int Ny = TabData.ny(); 
    float KurtV,MeanVal,SigVal;
    int i,j,p=0;
    fltarray V(Ny);
    int NCol = Nx* Ny;
    fltarray W(NCol);
    
    for (i=0; i<Nx; i++) 
    {
    	for (j=0; j<Ny; j++) V(j) = TabData(i,j);
    	MeanVal = V.mean();
    	SigVal = V.sigma();
    	for (j=0; j<Ny; j++) 
    	{
    		W(p) = (TabData(i,j) - MeanVal)/SigVal;
    		p++;
    	}
    }

	KurtV = curtosis(W.buffer(),NCol);    
    return KurtV;
}

/****************************************************************************/

float GMCA::CalcSigmaMadResidual(fltarray & TabData, fltarray & TabSource)
{ 
    fltarray TabChannel;
    fltarray TranspMixingMat;
    fltarray Mat= MixingMat;
    float MadValue;
    fltarray V;
    int NCol_Mat;
    int p=0;
    int i,k;
    
    NCol_Mat = TabData.nx()*TabData.ny(); 
    
    V.alloc(1,NCol_Mat);
    
    MAT.transpose(Mat,TranspMixingMat);
    MAT.mat_mult(TranspMixingMat, TabSource, TabChannel);
    TabChannel -= TabData;
    
    for (i=0;i<TabData.nx();i++) 
    {
    	for (k=0;k<TabData.ny();k++)
    	{
    		V(p) = TabChannel(i,k);
    		p++;
    	}
    }
    
    MadValue = get_sigma_mad(V.buffer(), NCol_Mat);
    
    return MadValue;
}

/****************************************************************************/

float GMCA::CalcKurtosisResidual(fltarray & TabData)
{ 
    fltarray TabChannel;
    float KurtValue;

    
    // Definitions
    fltarray TabMMTemp,TabSTemp,LeftA;

    // Compute the pseudo-inverse:
     
    leftpi_mat(MixingMat, LeftA);
    InvMixingMat = LeftA;

    // Estimate the corresponding sources :
    fltarray NLeftA=LeftA;
     
    MAT.transpose(NLeftA,LeftA);
    MAT.mat_mult(LeftA,TabData,TabSTemp);
     
    fltarray TranspMixingMat;
    fltarray Mat= MixingMat;
            
    MAT.transpose(Mat,TranspMixingMat);
    MAT.mat_mult(TranspMixingMat, TabSTemp, TabChannel);
    TabChannel -= TabData;
    
    KurtValue = MyKurt(TabChannel);
    cout << KurtValue << endl;
    return KurtValue;
}

/****************************************************************************/

void GMCA::init_mixmat(fltarray &TabCannels) // Devrait être fait dans le constructeur
{ 
   int i,k;
   fltarray Temp(NbrCannels,NbrSources);
   fltarray Vec(NbrCannels);
   
   // Initialize with a PCA
   
   //pca(TabCannels,Temp);
   
   for (i=0;i < NbrCannels;i++)
   {
      for (k=0;k < NbrSources;k++)
      {
         Temp(i,k) = get_random(-1,1);
      }
   }
      
   if (NbrColumnInit > 0)   // Initialize with a known sub-matrix
   {
   	if (L1_Matrix_Constraint == True)
        {
             if (NbrColumnInit > 1)
             {
                 // Transform the known spectra
                    fltarray TempM(NbrColumnInit,NbrCannels);
                    MAT.transpose(MatColumnInit,TempM);
                    mat_wt1d_trans(TempM);
                    MAT.transpose(TempM,MatColumnInit);
              }
    
            if (NbrColumnInit == 1)
            {
            for (k=0;k<NbrCannels; k++) Vec(k) = MatColumnInit(k);
            MR_Mat.transform(Vec);
            int ind=0;
	        for (int s=0; s < MR_Mat.nbr_band (); s++)
	        {
	           for (k=0; k < MR_Mat.size_scale_np (s); k++) MatColumnInit(ind+k) = MR_Mat(s,k);
	           ind += MR_Mat.size_scale_np (s);
            }
        }
    }
	
   	
   	for (i=0;i < NbrColumnInit ;i++)
   	{
   	if (NbrColumnInit  > 1) for (k=0;k < NbrCannels;k++) Temp(k,i)=MatColumnInit(k,i);
   	if (NbrColumnInit  == 1) for (k=0;k < NbrCannels;k++) Temp(k,i)=MatColumnInit(k);
   	}
   
   	
   }
      
     
   if (NbrKnownColumn > 0)
   {
   	   //cout << "Know Columns : " << NbrKnownColumn << endl;
        if (L1_Matrix_Constraint == True)
        {
             if (NbrKnownColumn > 1)
             {
                 // Transform the known spectra
                    fltarray TempM(NbrKnownColumn,NbrCannels);
                    MAT.transpose(MatKnownColumn,TempM);
                    mat_wt1d_trans(TempM);
                    MAT.transpose(TempM,MatKnownColumn);
              }
    
            if (NbrKnownColumn == 1)
            {
            for (k=0;k<NbrCannels; k++) Vec(k) = MatKnownColumn(k);
            MR_Mat.transform(Vec);
            int ind=0;
	        for (int s=0; s < MR_Mat.nbr_band (); s++)
	        {
	           for (k=0; k < MR_Mat.size_scale_np (s); k++) MatKnownColumn(ind+k) = MR_Mat(s,k);
	           ind += MR_Mat.size_scale_np (s);
            }
        }
    }
   } 
   
   
   
   {
   	for (i=0;i < NbrKnownColumn;i++)
   	{
   	if (NbrKnownColumn > 1) for (k=0;k < NbrCannels;k++) Temp(k,i)=MatKnownColumn(k,i);
   	if (NbrKnownColumn == 1) for (k=0;k < NbrCannels;k++) Temp(k,i)=MatKnownColumn(k);
   	}
   }
   
   MixingMat = Temp;
   RecMixingMat = Temp;
}

/****************************************************************************/

void GMCA::mat_wt1d_trans(fltarray &Mat)  // Il vaudrait mieux mettre un fltarray en entree car ce sont les colonnes de TabData que nous allons transformer
{
    int i,k,s;
    fltarray Vect(NbrCannels);
    
   // cout << "1D - Spectral Wavelet Transform" << endl;
    
    for (i=0; i < Mat.nx(); i++)
    {
     for (k=0;k<NbrCannels; k++) Vect(k) = Mat(i,k);  
     
     MR_Mat.transform(Vect);
	 	 
	 int ind=0;
	 
	 for (s=0; s < MR_Mat.nbr_band (); s++)
	 {
	 for (k=0; k < MR_Mat.size_scale_np (s); k++) 
	       {
	           Mat(i,ind+k) = MR_Mat(s,k);
	           
	       }
	       
	       ind += MR_Mat.size_scale_np (s);
    }
    }
}

/****************************************************************************/

void GMCA::mat_wt1d_recons(fltarray &Mat)
{
    int i,k;
    fltarray Vect(NbrCannels);
        
    for (i=0; i < Mat.ny(); i++)
    {
	 int ind=0;
      for (int s=0; s < MR_Mat.nbr_band (); s++)
	 for (k=0; k < MR_Mat.size_scale_np (s); k++) MR_Mat(s,k) = Mat(ind++,i); 
      MR_Mat.recons(Vect);
      for (k=0;k<NbrCannels; k++) Mat(k,i) = Vect(k);  
    }
}

/****************************************************************************/

void GMCA::leftpi_mat(fltarray &Mat,fltarray &LPIMat)  // compute the inverse matrix (should be great to define a similar function to compute pseudo-inverses (left and right))
{
    
    fltarray transpM;
    fltarray cM;
    fltarray IcM;
	fltarray tempA;
		
    MAT.transpose(Mat,transpM); // X -> X^T
    MAT.mat_mult(Mat,transpM,cM); // X -> X^T*X
    //MAT.inv_mat_svd(cM, IcM); // -> (X^T*X)^(-1)
    MAT.inv_mat_iter(cM, IcM);
	MAT.transpose(Mat,transpM);		 
    MAT.mat_mult(transpM,IcM,tempA);
	LPIMat = tempA;

}

/****************************************************************************/

void GMCA::rightpi_mat(fltarray &Mat, fltarray &RPIMat)  // compute the inverse matrix (should be great to define a similar function to compute pseudo-inverses (left and right))
{
    fltarray transpM;
    fltarray cM;
    fltarray IcM;
	fltarray tempA;
	    
    MAT.transpose(Mat,transpM);
    MAT.mat_mult(Mat,transpM,cM);
    //MAT.inv_mat_svd(cM, IcM);
    MAT.inv_mat_iter(cM, IcM);
    MAT.mat_mult(transpM,IcM,tempA);
	
	RPIMat = tempA;
}

/******************** CalcMadMat ******************************/

float GMCA::CalcMadMat(fltarray &TabSource)
{
    int i,k;
    int Nxd = TabSource.nx(),Nyd = TabSource.ny();
    fltarray V(Nxd*Nyd);
    int Val = 0;
    
    for (i=0;i < Nxd;i++)
    for (k=0;k < Nyd;k++)
    {
     V(Val) = TabSource(i,k);
     Val++;    
    }
    
    float Mad = get_sigma_mad(V.buffer(), Nxd*Nyd);
    return Mad;
}

/******************** Hard Threshold the sources******************************/

void GMCA::HT_Sources(fltarray &TabSource, float &KThrd)  // Applying the matrix on the data
{
   int NbrCoef = TabSource.nx();
   int Nyd = TabSource.ny();
   int i,k;
   float Mad;
   fltarray V(NbrCoef);
   float Thrd;
   
   if (GlobThrd == True)
   {
   	Mad = CalcMadMat(TabSource); 
						      
     Thrd = KThrd*Mad;
   	
   }
                
   for (i=0; i < Nyd; i++)
   {
   	if (GlobThrd == False)
   	{
   		for (k=0; k < NbrCoef ; k++) V(k) = TabSource(k,i);	
   		
		Mad = get_sigma_mad(V.buffer(), NbrCoef); 
						      
     	Thrd = KThrd*Mad;
   	}
        								      
     for (k=0; k < NbrCoef ; k++)
     {
     	if (ABS(TabSource(k,i)) < Thrd) TabSource(k,i)=(float) 0;
     }
   }
}


/******************** Hard Threshold the spectra******************************/

void GMCA::HT_MixMat(fltarray &Mat, float &KThrd)  // Applying the matrix on the data
{
   int NLin_Mat = Mat.axis(1);
   int NCol_Mat = Mat.axis(2);
   int i,k;
   float Mad;
   fltarray V;
   
   V.alloc(1,NCol_Mat-KStart);
                 
   for (i=0; i < NLin_Mat; i++)
   {
   	for (k=0; k < NCol_Mat-KStart ; k++) V(0,k) = Mat(i,NCol_Mat-k-1); // Ok
    //if (i == 0) for (k=0; k < NCol_Mat ; k++) cout << Mat(i,k) << endl; // PAS SPARSE!!!!!!!!!
     Mad = get_sigma_mad(V.buffer(), V.n_elem());
      
     for (k=0; k < NCol_Mat-KStart ; k++)  // Should not threshold the coqrse scale
     {
		float Thrd = KThrd*Mad;
     	if (ABS(Mat(i,NCol_Mat-k-1)) < Thrd) Mat(i,NCol_Mat-k-1)=0;  

     }
   }
}

/****************************************************************************/

void GMCA::apply_mat(fltarray & TabData, fltarray &Mat, fltarray & Result)  // Applying the matrix on the data
{
   int NbrCoef = TabData.nx();
   int Nyd = TabData.ny();
   int i,k,NCol_Mat = Mat.axis(2);
   fltarray V,R;
   
   V.alloc(1,Nyd);
   R.alloc(1,NCol_Mat);
   
   if ((Result.n_elem() == 0) || (Result.naxis() != 2)) Result.alloc(NbrCoef,NCol_Mat);
   if (Result.axis(1) != NbrCoef) Result.alloc(NbrCoef,NCol_Mat);
   else if (Result.axis(2) != NCol_Mat) Result.alloc(NbrCoef,NCol_Mat);
    
   for (i=0; i < NbrCoef; i++)
   {
      for (k=0; k < Nyd ; k++) V(0,k) = TabData(i,k);
      MAT.mat_mult(Mat,V,R);
      for (k=0; k < NCol_Mat; k++) Result(i,k) = R(0,k);
   }
}

/******************** Hard Threshold the sources******************************/

void GMCA::recons_sources(fltarray &DataIn, fltarray &EstSources)  // Applying the matrix on the data
{
 int Nx = DataIn.nx();
 int Ny = DataIn.ny();
 int Nz = DataIn.nz();
 int i,k,l;
  fltarray RefData,RefSources;
  int Deb = 0;
   
 // Reform the data

 RefData.alloc(Nz,Nx*Ny);
 
 for (l = 0;l < Nz;l++)
 {
   for (i=0; i < Nx ; i++) 
   {    
        Deb = i*Ny;
        for (k=0; k < Ny ; k++) RefData(l,Deb+k) = DataIn(i,k,l);
   }
 }
 
 // Apply the mixing matrix
 MAT.mat_mult(RefData,InvMixingMat,RefSources);
 
 // Reform the sources     
 EstSources.alloc(Nx,Ny,NbrSources);
 
 for (l = 0;l < NbrSources;l++)
 {
   for (i=0; i < Nx ; i++) 
   {    
        Deb = i*Ny;
        for (k=0; k < Ny ; k++) EstSources(i,k,l) = RefSources(l,Deb+k);
   }
 }
   
}
 
  /**********************************************************************************/
 
 void GMCA::apply_mask(fltarray &Data,fltarray &Mask,fltarray &DataRec)
 {
    int i,k,s;
    for (i=0;i<Data.nx();i++)
    for (k=0;k<Data.ny();k++)
    for (s=0;s<Data.nz();s++) DataRec(i,k,s) = (1 - Mask(i,k,s))*DataRec(i,k,s) + Data(i,k,s);
 }
 
 /**********************************************************************************/
 
 void GMCA::reform_data(fltarray &TabSources,fltarray &Data)
 {
  fltarray TabCh(TabSources.nx(),NbrCannels);
  int i,k,l;
  
  for (i=0;i<TabSources.nx();i++)
  for (k=0;k<NbrCannels;k++)
  {
   float Vv = 0;
   for (l=0;l<NbrSources;l++) Vv = Vv + MixingMat(k,l)*TabSources(i,l);
   TabCh(i,k) = Vv;
  }
  
  // MAT.transpose(MixingMat,TranspMixingMat);
  // MAT.mat_mult(TranspMixingMat, TabSources, TabCh);

  recons_data(TabCh,Data);
 
 }
 
  /**********************************************************************************/
 
 void GMCA::transform_data(fltarray &DataRec,fltarray &TabCannels)
 {
    //cout << DataRec.nx() << DataRec.ny() << endl;
    transform_sources(DataRec,TabCannels);
    if (L1_Matrix_Constraint == True) mat_wt1d_trans(TabCannels);
 }
 
 /**********************************************************************************/
 
 void GMCA::recons_data(fltarray &TabCannels,fltarray &Data)
 {
  // First if 1D 
  if (L1_Matrix_Constraint == True)
  {
    
    int i,k;
    fltarray Vect(NbrCannels);
        
    for (i=0; i < TabCannels.nx(); i++)
    {
	 int ind=0;
      for (int s=0; s < MR_Mat.nbr_band (); s++)
	 for (k=0; k < MR_Mat.size_scale_np (s); k++) MR_Mat(s,k) = TabCannels(i,ind++); 
      MR_Mat.recons(Vect);
      for (k=0;k<NbrCannels; k++) TabCannels(i,k) = Vect(k);  
    }
  
  }
  // 2D reconstruction  
  
  transrecons_sources(TabCannels,Data);
    
 }
 
/****************************************************************************/