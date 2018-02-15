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


#ifndef _GMCA_H_
#define _GMCA_H_

// #include "MR2D1D.cc"
#include "DefMath.h"
#include "MatrixOper.h"  // defined in $TOOLS
#include <cmath>
#include "Array.h"
#include "NR.h"
#include "IM_IO.h"
 #include "IM_IO.h"
 #include "MR1D_Obj.h"

// #include "MR2D1D.h"

/****************************************************************************/

#define DEF_GMCA_MAX_ITER 100

enum type_gmca_threshold_decrease {GMCA_THRES_LINEAR, GMCA_THRES_EXP, GMCA_THRES_MAD, GMCA_THRES_MOM};

#define DEF_GMCA_DECREASE GMCA_THRES_LINEAR

class GMCA {
   protected:
    MatOper MAT;  // See file $Tools/MatrixOper.cc and .h
   MR_1D MR_Mat; // 1D Wavelet class for matrix transformation
    FilterAnaSynt FAS;
    
public:     
   Bool Verbose; 
    int NbrCannels;  // Number of channels 
    int NbrSources;  // Number of Sources
    int Max_GMCA_Iter; // Maximum number of iterations in the GMCA algorithm
    fltarray MixingMat; // Mixing matrice
    fltarray RecMixingMat; // When sparsity is applied on the column of the mixing matrix, RecMixingMat is the reconctructed mixing matrix
    int MatNbrScale1D; // Number of Wavelet Scale used when applying the l_1 norm constraint on the matrix
    fltarray InvMixingMat; // Inverse Matrice (i.e. the matrice to be applied to the data to recover the sources)
    int NbrKnownColumn;   // Number of column which are known in the mixing matrix (a priori)
    int NbrColumnInit;   // Number of column which are known in the mixing matrix (a priori)
    fltarray MatKnownColumn; // Matrix containing the known column
    fltarray MatColumnInit; // Initialize the mixing matrix with fixed columns
    type_gmca_threshold_decrease TypeThresholdDecrease; // Select the decreasing threshold strategy 
    Bool Inpainting;            // If true, GMCA consider missing data, and apply an inpainting technique
    fltarray Mask;              // if Inpainting==True, Mask contains the available data Mask(i) = 1 or 0 (0 for no data).
    Bool PositiveMatrix;        // if True, the mixing matrice is assume to be positive
    Bool PositiveSource;        // if True, the sources are assumed to be positive
    Bool OrthoSpectra;        // if True, the spectra are orthogonalized
    Bool L1_Matrix_Constraint;  // if True, a sparsity constraint is applied on the matrix : this is the hyperspectral case
    Bool EstimNbSources;  // if True, estimate the number of sources
    Bool DisjSpec;  // if True, the mixing matrix have columns with disjoint support for high K-Mad
    Bool GlobThrd; // if True Global Thresholding on S (should also be made for A). 
    float ErrRate; // 
    float KMin; // Last K-Mad
    fltarray SingVec; // Singular Vectors
    Bool SVConst; // Constrain the mixing matrix columns to span the same subspace as the SingVec
    int NbrCoef;
    float L2err;
    int KStart;
    bool GetRelax;
    float RelaxVal;
    
    GMCA () {NbrCannels=NbrSources=0;Max_GMCA_Iter=DEF_GMCA_MAX_ITER; NbrColumnInit=0;
             Inpainting=False;PositiveMatrix=False;PositiveSource=False;ErrRate=1e-7;KMin=3;MatNbrScale1D=0;
	     Verbose=False;DisjSpec=False;OrthoSpectra=False;}
	  
    float l2_error(fltarray & TabData,fltarray & TabSource);
    // Return the l2 norm of an estimation, i.e. || TabData - MixingMat # TabSource ||^2
    
    float CalcSigmaMadResidual(fltarray & TabData,fltarray & TabSource);
	// Return the Sigma Mad of TabData - MixingMat # TabSource in the coefficient space
	
	float CalcKurtosisResidual(fltarray & TabData);
	// Return the kurtosis of TabData - MixingMat # TabSource in the coefficient space
	
	void RetrieveSourcesFromMixmat(fltarray & TabData,fltarray & TabSource);
	// Retrieve the sources from the data and the mixing matrix
	
	void Calc_InvMixingMat(fltarray &TabSource);
	// Update the mixing matrix
	
	void pca(fltarray &Mat,fltarray &VMat);
	// Compute a PCA to initialize VMat
	
    void Normal_Mixmat();
    // Normalize the mixing matrix

    void leftpi_mat(fltarray &Mat, fltarray & LPIMat); 
    // Compute the left pseudo-inverse of a matrix
    
    void rightpi_mat(fltarray &Mat, fltarray & RPIMat); 
    // Compute the left pseudo-inverse of a matrix
	
	void Sort_Sources(fltarray &TabSource);
	// Sorting the sources with decreasing energy
	
	void Where_MixingMat(fltarray &Mat,fltarray &IndexA, int &NbNzS);
	// Look for non-zero spectra
	
	float CountResiThrd(fltarray & TabData, fltarray & TabSource, float & Thrd);
	
	float MyKurt(fltarray & TabData);
	
	void print_mixmatnorm();
	// print the std of the columns of mixmat
	
	float calc_L0norm(fltarray &TabSource); 
    // Compute the left pseudo-inverse of a matrix
    
    float calc_L1norm(fltarray &TabSource); 
    // Compute the left pseudo-inverse of a matrix
    
    float CalcMadMat(fltarray &TabSource);
	// Compute the mad of a matrix
	
	void ortho_sources_svd(fltarray &Mat,fltarray &V);
	//Orthogonalization of the sources
	
	void ortho_spectra_svd();
	//Orthogonalization of the spectra
	
	void Where_Sources(fltarray &TabSource, fltarray &IndexS, int &NbNzS);
	// Look for non-zero sources
	
	void Update_Sources(fltarray &TabData,fltarray &TabSource , float &KThrd);
	// Update the sources
	
	void Update_MixingMat(fltarray &TabData,fltarray &TabSource, float &KThrd);
	// Update the mixing matrix
		
	void HT_Sources(fltarray & TabSource, float & KThrd);  // Y faire egalement entrer le type de seuillage
	// Hard-threshold the sources
	
	void HT_MixMat(fltarray &Mat, float &KThrd);
	// Hard-threshold the mixmat
	
	void init_mixmat(fltarray &TabCannels);
	// Initialize the mixing matrix
	
	void print_mat(fltarray &Mat);
	// Print a matrix
	
	void subs_mean(fltarray &Data);
	// Substracting the mean value
	
	void force_disjspec();
	// Force the mixing matrix to have disjoint supports for KThrd > 7
	
	void force_disjsources(fltarray &TabSource);
	// Force the sources to have disjoint supports for KThrd > 7
    
    void apply_mat(fltarray & TabData,fltarray & Mat, fltarray & TabRes);
    // Apply the matrix to the data
    
    void mat_wt1d_trans(fltarray &Mat); 
    // 1D wavelet transform along the column of the matrix   MixingMat
    
    void mat_wt1d_recons(fltarray &Mat);   
    // 1D wavelet reconstruction  along the column of the matrix  MixingMat
    
    void recons_sources(fltarray &DataIn, fltarray &EstSources);
    // reconstruct the estimated sources from the mixing matrix estimate and the input data
    
    void init_TabSource(fltarray &TabSource,int & NbrSamples);
    // Initialize TabSource
    
    void ortho_hypcube(fltarray &TabSource);
    // Orthogonalize the hyperspectral source cubes 
    
    void positive_cube_constraint(fltarray &Data);
    // Constrain a cube to be positive
    
    void positive_mixmat_constraint();
    // Constrain a cube to be positive
    
    void PCAMixMat();
    // Constrain the subspace spaned by Mixing Mat
    
    void recons_data(fltarray &TabCannels,fltarray &Data);
    // To reconstruct the data
    
    void reform_data(fltarray &TabSources,fltarray &Data);
    // To reconstruct the data
        
    void transform_data(fltarray &DataRec,fltarray &TabCannels);
        
    void apply_mask(fltarray &Data,fltarray &Mask,fltarray &DataRec);
        
    void run_gmca(fltarray &TabCannels, fltarray & TabSource);
    // run the gmca method 
    
    virtual void run(fltarray &TabCannels, fltarray & TabSource, fltarray & MixingMat, fltarray & InvMixingMat) {};
    // run the gmca method calling run_gmca
    
    virtual void transrecons_sources(fltarray &TabVect,fltarray &Recdata) {};
    
    virtual void transform_sources(fltarray &Data,fltarray &TabVect) {};
    
    virtual void inpainting_run(fltarray &TabCannels, fltarray & InpData){};
    //  run the inpainting method calling run_gmca
    
    virtual ~GMCA(){};
}; 

/****************************************************************************/
 
#endif


