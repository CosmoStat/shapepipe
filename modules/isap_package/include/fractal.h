#ifndef _FRACTAL_H
#define _FRACTAL_H

#include "GlobalInc.h"
#include "IM_Obj.h"
#include "IM_IO.h"
#include "MR1D_Obj.h"


#define LENGTH_SEARCH_CHILD 80
#define LENGTH_MAX_SEARCH_TWO_CHILD 7

// Bollean Tab of position of Max (1 if Max is detected, 0 else)
//intarray* spo_DetectMax= new intarray (96,729+2*729/3);
//intarray* spo_AllMax= new intarray (96,729+2*729/3);


extern intarray* spo_DetectMax;
extern intarray* spo_AllMax;
extern int si_compt;

#define MAX_CHAIN_NUMBER 100



class to_GenDevilStairCase {

  fltarray* po_Level;   // level of DevileStairCase
  int i_NbPoint;        // nb points of DevilStairCase only (without Border)
  fltarray* po_Prob;    // coeff of all intervals of DSC
  int i_NbPBorder;      // if != 0, we add i_NbPBorder points 0 before
                        // DSC and i_NbPBorder points 1 after DSC.
public:
  to_GenDevilStairCase (int pi_NbPoint, int pi_NbPBorder=0,
                        fltarray* ppo_Prob=NULL);
  fltarray& get_level ();
private:
  void compute ();
  void random_compute ();
};



class to_LineElt {
public:
  int i_scale;
  int i_ind;
  to_LineElt* po_suivant1;
  to_LineElt* po_suivant2;
  to_LineElt (int pi_scale, int pi_ind) {
	i_scale = pi_scale;
    i_ind = pi_ind;
	po_suivant1 = po_suivant2 = (to_LineElt*)NULL;
  }
  //Bool search_max_line (to_LineElt* ppo_father);
};

int draw_line (to_LineElt* ppo_Elt, Ifloat& pro_MatRec, int pi_NumLine);
Bool search_max_line (to_LineElt* ppo_father);
Bool is_best_father(to_LineElt* ppo_Childreen, to_LineElt* ppo_father);

enum te_type_chain {
   E_CHAIN_OK = 0,
   E_LENGTH_CHAIN_TOO_SHORT = 1
};

char* str_type_chaine (te_type_chain pe_TypeChain);

class to_Skeleton {
private:
  MR_1D ro_Mr1dMaxData;                 // Max Wavelet coef
  fltarray* po_IntSkel;                 // skel of max
  int i_NbScale;                        // number of scale of ro_Mr1dMaxData
  int i_NbPoint;                        // number of point of ro_Mr1dMaxData
  int i_NbMaxChain;                     // number max of chain
  int i_NbChain;                        // number of chain
  to_LineElt** pto_Elt;                 // tab of begin of chain
  int *pi_NbElt;                        // number elt of chain
  te_type_chain *pi_TypeChain;          // type chain : ok , remove ..
  intarray* po_IndTabMax;               // Ind (scale, pos) of all max
  
public:
  Bool Verbose;
  to_Skeleton (MR_1D& pro_MaxData, int pi_NbMaxChain=100);
  void init (int pi_NbMaxChain);
  fltarray& compute ();
  fltarray& get_skel ();
  int number_of_chain () {return i_NbChain;};
  void draw_line (to_LineElt* ppo_Elt, int pi_NumLine, Bool pe_info,fltarray& po_Skel);
  fltarray& draw_lines (Bool pe_info=False);
  void take_max_along_inf_scale ();
  void remove_thin_skel (int pi_Length=0);
  void remove_low_level (float pf_DynRange);
  void remove_all_max_at_min_dist (float pf_EpsDistMinBetweenMax);
  
private:
  void search_max_line (to_LineElt* ppo_father);
  int left_search_max_line (to_LineElt* ppo_father, int pi_LengthSearch=3);
  int right_search_max_line (to_LineElt* ppo_father, int pi_LengthSearch=3);
  void create_childreen (to_LineElt* ppo_father, int pi_IndMax);
  int number_of_scale () {return ro_Mr1dMaxData.nbr_scale();};
  int number_of_point () {return ro_Mr1dMaxData.size_ima_np();};
  float & operator() (int s, int i) {return ro_Mr1dMaxData(s,i);};
};



class to_ThermoDynRepartFunction {
private:
  fltarray o_Max;                // WTMM of input signal
  fltarray o_SupportMax;         // skel of o_Max
  int i_NbPoint;
  int i_NbScale;
  int i_NbExp;
  int i_NbAlpha;                 
  fltarray* po_q;                // list of exponents 
  fltarray* po_Z;                // thermo partition

public:
  to_ThermoDynRepartFunction (fltarray& pro_MaxWavTransf, 
                              fltarray& pro_SupportMaxWavTransf);
  void init_default_q();
  void compute_thermo_partition (fltarray* ppo_q=NULL, char* ppc_NameOut=NULL);
  void write_thermo_partition (int pi_IndMinscale, int pi_IndMaxScale);
};



class to_ThermoDynAnalysis {
private:
             
  fltarray* po_s;                // list of scale 
  fltarray* po_q;                // list of exponents 
  fltarray* po_Z;                // thermo partition
  fltarray* po_Tau;              // 
  fltarray* po_SigTau;              // 
  int i_IndMinScale;             // ind min scale for tau calcul
  int i_IndMaxScale;             // ind mac scale for tau calcul
  fltarray* po_Alpha;            // abs used for frac spectrum po_FracSpectrum
  fltarray* po_FracSpectrum;     // frac spectrum
  fltarray* po_ValQMin;          // q min used for frac spectrum

public:
  int i_NbScale;
  int i_NbExp;
  int i_NbAlpha;      

  to_ThermoDynAnalysis (char* ppc_FileName);
  void init_default_alpha();
  void compute_tau (int pi_IndMinscale, int pi_IndMaxScale);
  void display_tau ();
  void compute_frac_spectrum (fltarray* ppo_Alpha=NULL);
  void display_frac_spectrum ();
  void alpha_for_q_equals_zero ();

};



#endif
