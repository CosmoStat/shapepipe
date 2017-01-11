
/******************************************************************************
**                   Copyright (C) 2003 by IRSN
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Philippe Querre
**
**    Date:  07/04/03
**    
**    File: IM1D_regul.h
**
*******************************************************************************
**
**    DESCRIPTION  Class for regularization
**    -----------  
**                 
******************************************************************************/

#ifndef _IM1D_REGUL_H_
#define _IM1D_REGUL_H_

/****************************************************************************/

enum type_oper_grad { OPER_LAPLACIAN,
                      OPER_MARKOV_2};

class RegulSig {
  
   inline int sgn(float Val) {return ( (Val >= 0) ? 1: -1);}
  
   /**************************************************************************/
   /**********                      ******************************************/
   /**********      markov_val2     ******************************************/
   /**********                      ******************************************/
   /**************************************************************************/   
   inline float markov_val2 (fltarray& Sig, int i) {
   
      double p = MarkovPowerParam;
      double dp=p-1;
      int N = 2;
      dblarray Tab(N);
      Tab(0) = Sig(i) - Sig(i-1,BorderType);
      Tab(1) = Sig(i) - Sig(i+1,BorderType);
      double  Val = 0;
      for (int k=0; k<N; k++) Val += p*sgn(Tab(k))*pow(ABS(Tab(k)),dp);
      return (float) Val / (float) N;
   }
  
public:  
   /**************************************************************************/
   /**********                        ****************************************/
   /**********      laplacian_val     ****************************************/
   /**********                        ****************************************/
   /**************************************************************************/   
   inline float laplacian_val (fltarray &Obj, int i) {
   
       float Val = Obj(i) - 0.5*( Obj(i+1,BorderType)
		               +  Obj(i-1,BorderType));
       return Val;
   } 

   int NbrScale; // Number of scales used in the soft thresholding
   int NbrUndecimatedScale;
                 // Number of undecimated scale. By default all scales are
		 // undecimated
   type_sb_filter TypeFilter;
                 // Type of filter. By default Haar filter is chosen, which
		 // means that SoftThresholding is equivalent to a TV constaint.
   Bool ExpDecreasingLambda;
                 // If true, the threshold descrease with the scale.
		 
 
   
   void mr1d_soft_threshold(fltarray &Ima, fltarray  &Rec, float Lambda, float NoiseLevel=0.);
   // Apply the undecimated bi-orthogonal WT to the image,
   // soft threshold the coefficients with the threshold level Lambda,
   // Reconstruct the filtered image in Rec
   void mr1d_soft_threshold(fltarray &Ima,  float Lambda, float NoiseLevel=0.)
                           {mr1d_soft_threshold(Ima,Ima, Lambda, NoiseLevel);}
        
private:       
   /**************************************************************************/
   /**********                    ********************************************/
   /**********      sig_regul     ********************************************/
   /**********                    ********************************************/
   /**************************************************************************/     
   inline void sig_regul (fltarray& Obj,  fltarray& Grad, float Regul) {
       int Nx = Obj.nx();
       int i;
       if (Grad.nx() != Nx) Grad.resize(Nx);

       if (Regul > 0) {
         switch(GradOperType) {
          
          case OPER_LAPLACIAN:
	     for (i=0; i<Nx; i++) Grad(i) = Regul*laplacian_val(Obj,i);
             break;
	  case OPER_MARKOV_2:
             for (i=0; i<Nx; i++) Grad(i) =  Regul*markov_val2(Obj,i); 
             break;
         }
      } else Grad.init();
   }
   
   
   
    
  public: 
  RegulSig() {NbrScale=4; NbrUndecimatedScale=-1; TypeFilter= F_MALLAT_7_9; // F_HAAR, F_MALLAT_7_9; 
               ExpDecreasingLambda=False;GradOperType = OPER_MARKOV_2;BorderType=I_CONT;
	       MarkovPowerParam=1.1;}
	         
    ~RegulSig () {}	
    type_oper_grad GradOperType;
    type_border BorderType;
    float MarkovPowerParam;  
    

   /**************************************************************************/
   /**********                    ********************************************/
   /**********      obj_regul     ********************************************/
   /**********                    ********************************************/
   /**************************************************************************/   
   inline void obj_regul (fltarray& Obj, fltarray& Grad, float Regul) {
   
       int Nx = Obj.nx();
       fltarray ImaAux (Nx,"aux");
       sig_regul (Obj,ImaAux,Regul);
       Grad -= ImaAux;
   }
};



#endif
