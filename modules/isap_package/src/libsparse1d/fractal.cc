
#include "fractal.h"
#include "IM_IO.h"
#include "IM1D_IO.h"

#define EPS 1.e-08
#define FRACFLTMAX MAXFLOAT
#define LENGTH_SEARCH 3
#define NOT_FIND -1


#define EPS_FRACTAL 1e-07


char* str_type_chaine (te_type_chain pe_TypeChain) {
   switch (pe_TypeChain) {
   case E_CHAIN_OK : return (""); break;
   case E_LENGTH_CHAIN_TOO_SHORT : return ("Chain remove : too short"); break;
   }
   return ("Pb in function str_typr_chaine");
};

int xi_compt;

to_GenDevilStairCase::to_GenDevilStairCase (int pi_NbPoint, 
                                            int pi_NbPBorder,
                                            fltarray* ppo_Prob) {
  
  
  // coef not defined => standard DevilStaircase with p1=.5, p2=0 and p3=.5
  if (ppo_Prob == (fltarray*)NULL) {
    int ai_NbInterval = 3;
    po_Prob = new fltarray (ai_NbInterval);
    (*po_Prob)(0)=.5; (*po_Prob)(1)=0; (*po_Prob)(2)=.5;
  } else {
    // init coef 
    po_Prob = new fltarray (ppo_Prob->n_elem());
    (*po_Prob) = (*ppo_Prob);
  }

  i_NbPoint = pi_NbPoint;
  i_NbPBorder = pi_NbPBorder;
  po_Level = new fltarray (i_NbPoint);
  for (int i=0;i<i_NbPoint;i++) (*po_Level)(i) = 1./2.;
  random_compute ();
}



fltarray& to_GenDevilStairCase::get_level () {

  if (i_NbPBorder == 0) return (*po_Level);
  else {

    fltarray* apo_RetTab = new fltarray (i_NbPoint+2*i_NbPBorder);
    int i,j;
    for (i=0,j=0; i<i_NbPBorder; i++,j++) 
      (*apo_RetTab)(i) = 0.;
      //(*apo_RetTab)(i) = 0. + 0.01*sin(8*3.1415*j/i_NbPoint);
    for (i=i_NbPBorder; i<i_NbPBorder+i_NbPoint; i++,j++)
      //(*apo_RetTab)(i) = (*po_Level)(i-i_NbPBorder)+
      //                   0.01*sin(8*3.1415*j/i_NbPoint); 
      (*apo_RetTab)(i) = (*po_Level)(i-i_NbPBorder);
    for (i=i_NbPBorder+i_NbPoint; i<i_NbPBorder+i_NbPoint+i_NbPBorder; i++,j++)
      (*apo_RetTab)(i) = 1.;
      //(*apo_RetTab)(i) = 1. + 0.01*sin(8*3.1415*j/i_NbPoint);
    return (*apo_RetTab);
  }
}


void to_GenDevilStairCase::compute () {

  int ai_Cmpt=0;
  fltarray ao_Temp (i_NbPoint);
  //f_p1=.69; f_p2=.46; f_p3=-.46; f_p4=.31;
  int ai_NbInterval = po_Prob->n_elem();

  while (True) {
    ai_Cmpt++;

    if (i_NbPoint/pow((float)ai_NbInterval,(float)ai_Cmpt) < ai_NbInterval) break;
    int ai_NbPointInInterval = i_NbPoint/ai_NbInterval;

    // loop on all intervals
    for (int i=0;i<ai_NbInterval;i++) {
 
      float af_Sigma=0; int j;
      for (j=0;j<i;j++) af_Sigma += (*po_Prob)(j);
      
      // loop on all points of current interval
      for (j=i*ai_NbPointInInterval;j<(i+1)*ai_NbPointInInterval;j++) {
     
	ao_Temp(j) = af_Sigma + (*po_Prob)(i) *
                                (*po_Level)(ai_NbInterval*j-(i*i_NbPoint));
      }
    }
    *po_Level = ao_Temp;
  }
}


void to_GenDevilStairCase::random_compute () {

  int ai_Cmpt=0;
  fltarray ao_Temp (i_NbPoint);
  //f_p1=.69; f_p2=.46; f_p3=-.46; f_p4=.31;
  int ai_NbInterval = po_Prob->n_elem();
  
  
  while (i_NbPoint/pow((float)ai_NbInterval,(float)ai_Cmpt) >= ai_NbInterval) ai_Cmpt++;
  fltarray* po_RandProb = new fltarray (ai_NbInterval, ai_Cmpt);
  //intarray* po_Position = new intarray (ai_NbInterval);

  /*for (int i=0; i<ai_Cmpt; i++) {
    for (int k=0; k<ai_NbInterval; k++) (*po_Position)(k)=k;
    for (int j=0; j<ai_NbInterval-1; j++) {
      // get a pos in [0,ai_NbInterval-j-1] ([0,3], [0,2], [0,1])
      int ai_pos = (ai_NbInterval-j-1)*rand()/(pow(2,15)-1);
      // set jth prob with prob at position pos(rand)
      (*po_RandProb)(j, i) = (*po_Prob)( (*po_Position)(ai_pos) );
      // reinit po_Position
      for (int k=ai_pos; k<ai_NbInterval-1; k++) 
      	(*po_Position)(k) = (*po_Position)(k+1);
    }
    // get the last pos in po_Position
    (*po_RandProb)(ai_NbInterval-1, i)=(*po_Prob)( (*po_Position)(0));
  }*/
  
  for (int i=0; i<ai_Cmpt; i++) {
    for (int j=0; j<ai_NbInterval; j++) {
      (*po_RandProb)(j, i) = (*po_Prob)(j);
    }
  }
  

  // po_RandProb->display(ai_NbInterval*ai_Cmpt);
  
  ai_Cmpt = 0;
  while (True) {
    ai_Cmpt++;

    int ai_NbPointInInterval = i_NbPoint/ai_NbInterval;

    // loop on all intervals
    for (int i=0;i<ai_NbInterval;i++) {
 
      float af_Sigma=0; int j;
      for (j=0;j<i;j++) af_Sigma += (*po_RandProb)(j,ai_Cmpt-1);

      
      // loop on all points of current interval
      for (j=i*ai_NbPointInInterval;j<(i+1)*ai_NbPointInInterval;j++) {
     
	ao_Temp(j) = af_Sigma + (*po_RandProb)(i,ai_Cmpt-1) *
                                (*po_Level)(ai_NbInterval*j-(i*i_NbPoint));
      }
    }
    *po_Level = ao_Temp;

    if (i_NbPoint/pow((float)ai_NbInterval,(float)ai_Cmpt) < ai_NbInterval) break;

  }
}






to_Skeleton::to_Skeleton (MR_1D& pro_MaxData, int pi_NbMaxChain) {
  Verbose = False;
  ro_Mr1dMaxData = pro_MaxData;
  i_NbScale = ro_Mr1dMaxData.nbr_scale();
  i_NbPoint = ro_Mr1dMaxData.size_ima_np();
  i_NbChain = 0;
  po_IntSkel = new fltarray (i_NbPoint, i_NbScale);
  po_IntSkel->init(0);
  init (pi_NbMaxChain);
}

void to_Skeleton::init (int pi_NbMaxChain) {
  i_NbMaxChain = pi_NbMaxChain;
  pto_Elt = (to_LineElt**) new to_LineElt* [i_NbMaxChain];
  pi_NbElt = new int [i_NbMaxChain];
  pi_TypeChain = new te_type_chain [i_NbMaxChain];
  for (int i=0;i<i_NbMaxChain;i++) pi_TypeChain[i]=E_CHAIN_OK;
  po_IndTabMax = new intarray (i_NbScale, i_NbPoint);
  // init max if (s,i) max then po_IndTabMax=1 else po_IndTabMax=0
  for (int s=0;s<i_NbScale-1;s++) {
    for (int i=0;i<i_NbPoint;i++) {
      (*po_IndTabMax)(s,i) = 0;
      if (ro_Mr1dMaxData(s,i) > EPS || ro_Mr1dMaxData(s,i) < -EPS)
	    (*po_IndTabMax)(s,i) = 1;
    }
  }
} 
	
fltarray& to_Skeleton::get_skel () {return *po_IntSkel;}

fltarray& to_Skeleton::compute (){

  // number of current max chain
  int ai_IndCurrentLine = 0;
  
  // for all scales and all points
  for (int s=i_NbScale-1; s>=0; s--) {
	
    // search for points all max at scale s
    for (int i=0; i<i_NbPoint; i++) {

      // is there a max at point (s,i)
      if ((*po_IndTabMax)(s,i) > EPS || (*po_IndTabMax)(s,i) < -EPS) { 

        // first element of new maxima line
        to_LineElt* apo_FirstEltMaximaLine = new to_LineElt (s,i);

        // remove Max(s,i) of Boolean Tab, init po_IntSkel
        (*po_IndTabMax)(s,i) = 0;
        (*po_IntSkel)(i,s) = 1;

        // one elt in current chain
        xi_compt = 1;

        // mem of begin line
        pto_Elt[ai_IndCurrentLine] = apo_FirstEltMaximaLine;
        if (ai_IndCurrentLine > i_NbMaxChain) {
          cout << "Too much chain" << endl;
          exit (-1);
        }

        // search for chidreen
        search_max_line (apo_FirstEltMaximaLine);

        // mem number of elt of current chain
        pi_NbElt[ai_IndCurrentLine] = xi_compt;

        ai_IndCurrentLine++; // for next maxima line
      }
    }
  }

  // end number of chain
  i_NbChain = ai_IndCurrentLine;

  return *po_IntSkel;
}


void to_Skeleton::search_max_line (to_LineElt* ppo_father) {

  // search firstmax on rigth and left
  int ai_RightInd = right_search_max_line (ppo_father, LENGTH_SEARCH);
  int ai_LeftInd = left_search_max_line (ppo_father, LENGTH_SEARCH);

  // test if not end of Max line
  if (ai_RightInd != NOT_FIND || ai_LeftInd != NOT_FIND) {

    // take best max
    int ai_IndMax = (ai_RightInd > ai_LeftInd ? ai_LeftInd : ai_RightInd);
    if (ai_IndMax == NOT_FIND) 
      ai_IndMax = (ai_RightInd > ai_LeftInd ? ai_RightInd : ai_LeftInd);
    if (ai_IndMax == ai_LeftInd) ai_IndMax = -ai_LeftInd;

	xi_compt++; // add one elt to current chain
    create_childreen (ppo_father, ai_IndMax);
  }
}


void to_Skeleton::create_childreen (to_LineElt* ppo_father,
                                    int pi_IndMax) {

  // create elt 
  to_LineElt* apo_Childreen = new to_LineElt (ppo_father->i_scale-1,
                                              ppo_father->i_ind+pi_IndMax);

  // chain father <-> child
  ppo_father->po_suivant1 = apo_Childreen;

  // remove Max(s,i) of Boolean Tab, init po_IntSkel
  (*po_IndTabMax)(ppo_father->i_scale-1, ppo_father->i_ind+pi_IndMax) = 0;
  (*po_IntSkel)(ppo_father->i_ind+pi_IndMax, ppo_father->i_scale-1) = 1;
 
  // search for rest of maxima line
  search_max_line (apo_Childreen);
}


int to_Skeleton::left_search_max_line (to_LineElt* ppo_father, 
                                       int pi_LengthSearch) {

  int ai_FatherScale = ppo_father->i_scale;
  int ai_FatherInd = ppo_father->i_ind;

  for (int i=0; i<pi_LengthSearch; i++) {

    if (    (ai_FatherInd - i > 0)
	     && (  (*po_IndTabMax)(ai_FatherScale-1,ai_FatherInd-i) >  EPS
             ||(*po_IndTabMax)(ai_FatherScale-1,ai_FatherInd-i) < -EPS)) { 

      // Max detected at sacle-1 and ind-i
      return (i);
    }
  }
  return (NOT_FIND);
}
 
   
int to_Skeleton::right_search_max_line (to_LineElt* ppo_father, 
                                        int pi_LengthSearch) {

  int ai_FatherScale = ppo_father->i_scale;
  int ai_FatherInd = ppo_father->i_ind;

  for (int i=0; i<pi_LengthSearch; i++) {

    if (    (ai_FatherInd + i < i_NbPoint)
         && (  (*po_IndTabMax)(ai_FatherScale-1,ai_FatherInd+i) >  EPS 
             ||(*po_IndTabMax)(ai_FatherScale-1,ai_FatherInd+i) < -EPS)) {

      // Max detected at sacle-1 and ind-i
      return (i);
    }
  }
  return (NOT_FIND);
}


void to_Skeleton::draw_line (to_LineElt* ppo_Elt, int pi_NumLine, Bool pe_info,
                             fltarray& po_Skel){

  // print info
  //if (    (xi_compt==0) 
  //    &&  ((*po_IntSkel)( ppo_Elt->i_ind, ppo_Elt->i_scale) != 0)) {
  if (xi_compt==0) {
    if (pe_info) 
	cout << "pere: (" << ppo_Elt->i_scale << "," << ppo_Elt->i_ind
	     << ")->" << ro_Mr1dMaxData (ppo_Elt->i_scale-1, ppo_Elt->i_ind)
             <<", ";
    xi_compt++;
  }

  // Level of line pi_NumLine is set to 10*(pi_NumLine+1)
  // if type chain is OK
  // if curent point is OK
  if (pi_TypeChain[pi_NumLine] == E_CHAIN_OK)
     po_Skel(ppo_Elt->i_ind, ppo_Elt->i_scale) = 10*(pi_NumLine+1);

  if (ppo_Elt->po_suivant1 != NULL) {
	
	// if elt not remove by -s dyn
	if ((*po_IntSkel)(ppo_Elt->po_suivant1->i_ind, 
                      ppo_Elt->po_suivant1->i_scale) != 0) {

      // number of elt of chain pi_NumChai
      xi_compt++;

      if (pe_info) cout << "fils"<<xi_compt<<": (" 
                        << ppo_Elt->po_suivant1->i_scale << "," 
                        << ppo_Elt->po_suivant1->i_ind << ") ->"
                        << ro_Mr1dMaxData (ppo_Elt->po_suivant1->i_scale,
			                   ppo_Elt->po_suivant1->i_ind) 
                        << ", ";
	}
  } 

  // if next Elt exist
  if (ppo_Elt->po_suivant1 != NULL) {
    to_LineElt* apo_Elt = ppo_Elt->po_suivant1;
    draw_line (apo_Elt, pi_NumLine, pe_info, po_Skel);
  }
}


fltarray& to_Skeleton::draw_lines (Bool pe_info) {

  fltarray* apo_skel = new fltarray (i_NbPoint, i_NbScale);
  for (int i=0;i<i_NbPoint;i++)
     for (int j=0;j<i_NbScale;j++)
        (*apo_skel)(i,j) = (*po_IntSkel)(i,j);
	
  for (int i=0; i<number_of_chain(); i++) {
        xi_compt = 0;
        //if (pe_info) 
	   cout << "Line " << i+1 << 
	           " : (nb elt:" << pi_NbElt[i] << 
	           "), " << str_type_chaine(pi_TypeChain[i]) << endl;
        draw_line (pto_Elt[i], i, pe_info, *apo_skel);   
	if (pe_info) cout << endl;			   
  }
  return *apo_skel;
}


void to_Skeleton::take_max_along_inf_scale () {

  for (int i=0; i<number_of_chain(); i++) {

    fltarray ao_Level (pi_NbElt[i]);
    int ai_Ind=0; 
    to_LineElt* ppo_Elt = pto_Elt[i];
    ao_Level(ai_Ind) =   ro_Mr1dMaxData (ppo_Elt->i_scale, ppo_Elt->i_ind);

    // get all level of maxima line in ao_Level
    while (ppo_Elt->po_suivant1 != NULL) {

    ai_Ind++;
    ppo_Elt = ppo_Elt->po_suivant1;
    ao_Level(ai_Ind) =  ro_Mr1dMaxData (ppo_Elt->i_scale, ppo_Elt->i_ind);
    }

    // sort 
    for (int j=pi_NbElt[i]-1; j>0; j--) 
      if (ao_Level(j) > ao_Level(j-1)) ao_Level(j-1)=ao_Level(j);
		
    // set all level of maxima line in ao_Level
    ai_Ind=0; ppo_Elt = pto_Elt[i];
    ro_Mr1dMaxData (ppo_Elt->i_scale,ppo_Elt->i_ind) = ao_Level(0);
    while (ppo_Elt->po_suivant1 != NULL) {

      ai_Ind++;
      ppo_Elt = ppo_Elt->po_suivant1;
      ro_Mr1dMaxData (ppo_Elt->i_scale,ppo_Elt->i_ind) = ao_Level(ai_Ind);
    }

  }
}     


void to_Skeleton::remove_thin_skel (int pi_Length) {

  if (pi_Length > 0 ) {

    to_LineElt* apo_Elt;
    for (int i=0; i<number_of_chain(); i++) {
      if (pi_NbElt[i] < pi_Length) {

	pi_TypeChain[i] = E_LENGTH_CHAIN_TOO_SHORT;

        // remove this chain => all level have zero now
        apo_Elt = pto_Elt[i];
        (*po_IntSkel)(apo_Elt->i_ind, apo_Elt->i_scale) = 0;
        while (apo_Elt->po_suivant1 != NULL) {
          apo_Elt = apo_Elt->po_suivant1;
          (*po_IntSkel)(apo_Elt->i_ind, apo_Elt->i_scale) = 0;
        }
      }
    }
  }
}      


void to_Skeleton::remove_low_level (float pf_DynRange) {

  if (pf_DynRange > 0) {

    int ai_NbRemove = 0;
    fltarray ao_LevelMaxAtScale (i_NbScale);
    for (int s=0; s<i_NbScale; s++) {

      // search max at all scale
      int i;
      ao_LevelMaxAtScale(s) = 0;
      for (i=0; i<i_NbPoint;i++) {
        if (ABS(ro_Mr1dMaxData(s,i)) > ao_LevelMaxAtScale(s))
        ao_LevelMaxAtScale(s) = ABS(ro_Mr1dMaxData(s,i));
      }

      // remove level < dynrange * max at scale s
      for (i=0; i<i_NbPoint;i++) {
        if (    (*po_IntSkel)(i,s) != 0 
             && ABS(ro_Mr1dMaxData(s,i)) < pf_DynRange/100. * ao_LevelMaxAtScale(s)) {

          cout << "remove max (Dyn range) at (scale:" << s << ", pos:" << i << ")" << endl;
          (*po_IntSkel)(i,s) = 0;
          ai_NbRemove++;
        }
      }
    }
    if (Verbose == True) cout << "number of points remove : " << ai_NbRemove << endl;
  }
}



void to_Skeleton::remove_all_max_at_min_dist (float pf_EpsDistMinBetweenMax) {

  // for all scales and all points
  for (int s=i_NbScale-1; s>=0; s--) {
	
    // search for points all max at scale s
    for (int i=0; i<i_NbPoint; i++) {

      // is there a max at point (s,i)
      if ((*po_IndTabMax)(s,i) > EPS || (*po_IndTabMax)(s,i) < -EPS) { 
      
        //is there another max at pf_EpsDistMinBetweenMax*scale
        int ai_NbPoint = (int) (pf_EpsDistMinBetweenMax * s);
	int bmin = (i-ai_NbPoint < 0 ? 0 : i-ai_NbPoint);
	int bmax = (i+ai_NbPoint > i_NbPoint-1 ? i_NbPoint-1 : i-ai_NbPoint);
	for (int k=bmin; k<=bmax; k++) {
          if (   ABS((*po_IndTabMax)(s,k)) > EPS
	      && ABS((*po_IndTabMax)(s,k)) > ABS((*po_IndTabMax)(s,k))) {
            (*po_IndTabMax)(s,k) = 0;
	    cout << "remove max (min dist between max) at (scale:" << s << ", pos:" << i << ")" << endl;
          } else if ( ABS((*po_IndTabMax)(s,k)) > EPS
	           && ABS((*po_IndTabMax)(s,k)) < ABS((*po_IndTabMax)(s,k))) {
	    (*po_IndTabMax)(s,k) = 0;
	    cout << "remove max (min dist between max) at (scale:" << s << ", pos:" << k << ")" << endl;
	  }
        }
      }
    }
  }
}









to_ThermoDynRepartFunction::to_ThermoDynRepartFunction (
                               fltarray& pro_MaxWavTransf, 
			                   fltarray& pro_SupportMaxWavTransf) {

  i_NbPoint = pro_MaxWavTransf.nx();
  i_NbScale = pro_MaxWavTransf.ny();
  o_Max = pro_MaxWavTransf;
  o_SupportMax = pro_SupportMaxWavTransf;
}

void to_ThermoDynRepartFunction::init_default_q () {
  po_q = new fltarray (5);
  for (int i=0;i<5;i++) (*po_q)(i) = -10.0+5*i;
}


void to_ThermoDynRepartFunction::compute_thermo_partition (fltarray* ppo_q,
                                                           char* ppc_NameOut) {
  int s;
  extern type_1d_format IO_1D_Format;
  if (ppo_q == (fltarray*)NULL) init_default_q();
  else po_q = ppo_q;
  i_NbExp = po_q->n_elem();
  po_Z = new fltarray (i_NbExp, i_NbScale);
  //po_Tau = new fltarray (i_NbExp);

  // for all scales
  for (s=0; s<i_NbScale-1; s++) {

    // for all exponents
    for (int e=0; e<i_NbExp; e++) {

      (*po_Z)(e,s) = 0;
      // calc thermodyn partition at scale s and exponent q(e)
      for (int i=0; i<i_NbPoint; i++) {

	    //if ( o_Max(i,s) != 0  ) {
        if (o_SupportMax(i,s) != 0) {

	      (*po_Z)(e,s) += pow (ABS(o_Max(i,s)), (*po_q)(e));
	    }
      }
    }  
  }

  // write Z(q,s), po_q, and sacle on disk for mf_analyse
  char atc_FileName[80];
  sprintf (atc_FileName, "Z_%s", ppc_NameOut);
  if (IO_1D_Format == F1D_FITS)
      fits_write_fltarr(atc_FileName, *po_Z);
  else io_write2d_ascii(atc_FileName, *po_Z);
  sprintf (atc_FileName, "q_%s", ppc_NameOut);
  io_1d_write_data(atc_FileName, *po_q);
  fltarray ao_Scale (i_NbScale);
  for (s=0; s<i_NbScale; s++) 
	ao_Scale(s) = s/12. + log(1./sqrt(3.))/log(2.);
  sprintf (atc_FileName, "s_%s", ppc_NameOut);
  io_1d_write_data(atc_FileName, ao_Scale);
}

void to_ThermoDynRepartFunction::write_thermo_partition (int pi_IndMinscale, 
                                                         int pi_IndMaxScale) {

  char atc_FileName[80];

  // X axe
  fltarray ao_AxeXLogScale (pi_IndMaxScale-pi_IndMinscale);
  for (int s=pi_IndMinscale; s<pi_IndMaxScale; s++) 
    ao_AxeXLogScale(s-pi_IndMinscale) = s/12. + log(1./sqrt(3.))/log(2.);
  io_1d_write_data("X_LogZfunctLoga", ao_AxeXLogScale);

  // for all exponents Y axe
  for (int e=0; e<i_NbExp; e++) {

    fltarray ao_AxeYLogZ (pi_IndMaxScale-pi_IndMinscale);
    for (int s=pi_IndMinscale; s<pi_IndMaxScale; s++) 
      ao_AxeYLogZ(s-pi_IndMinscale) = 
//	                          log ((double)(*po_Z)(e,s))/log(2)/(*po_q)(e);
  	                          log ((double)(*po_Z)(e,s))/log(2.);
    sprintf (atc_FileName, "Y_LogZfunctLoga_%d", e);
    io_1d_write_data(atc_FileName, ao_AxeYLogZ);
  }
}

 






to_ThermoDynAnalysis::to_ThermoDynAnalysis (char* ppc_FileName) {

  // read Z(q,s), po_q, and sacle on disk for mf_analyse
  char* apc_FileName;
  apc_FileName = new char[strlen(ppc_FileName)+8];
  sprintf (apc_FileName, "Z_%s", ppc_FileName);
  strcpy (apc_FileName+strlen(ppc_FileName)+7,"\0");
  po_Z = new fltarray();
  type_1d_format  DatFormat = io_detect_1dformat(apc_FileName);
  if (DatFormat == F1D_FITS) fits_read_fltarr(apc_FileName, *po_Z);
  else io_read2d_ascii(apc_FileName, *po_Z);
  
  char atc_FileName[256];
  sprintf (atc_FileName, "q_%s", ppc_FileName);
  strcpy (atc_FileName+strlen(ppc_FileName)+8,"\0");
  po_q = new fltarray();
  // fits_read_fltarr(atc_FileName, *po_q);
  io_1d_read_data(atc_FileName, *po_q);

  i_NbScale = po_Z->ny();
  i_NbExp = po_q->n_elem();

  fltarray ao_Scale (i_NbScale);
  sprintf (atc_FileName, "s_%s", ppc_FileName);
  strcpy (atc_FileName+strlen(ppc_FileName)+8,"\0");
  po_s = new fltarray();
  // fits_read_fltarr(atc_FileName, *po_s);
  io_1d_read_data(atc_FileName, *po_s);

  po_Tau = new fltarray (i_NbExp);
  po_SigTau = new fltarray (i_NbExp);
}



void to_ThermoDynAnalysis::init_default_alpha () {
  po_q = new fltarray (21);
  for (int i=0;i<21;i++) (*po_q)(i) = 1./20.*i;
}

void to_ThermoDynAnalysis::compute_tau (int pi_IndMinscale, 
                                        int pi_IndMaxScale) {

  i_IndMinScale = pi_IndMinscale;
  i_IndMaxScale = pi_IndMaxScale;
  int ai_NbScale = pi_IndMaxScale-pi_IndMinscale;
  // for all exponents
  for (int e=0; e<i_NbExp; e++) {

    // Calcul de Tau = f(q) (pour chaque courbe log(Z(q,a)) = f (log(a))
    // on alcule la pente Tau(q) en fonction de q

    float sx=0, sxx=0, sy=0, sxy=0, sigma_in=0;
    float delta;
    for (int s=i_IndMinScale; s<i_IndMaxScale; s++) {

      if ((*po_Z)(e,s) > 1e-99) { // pb if not max at scale s

        sx += (s/12. + log(1./sqrt(3.))/log(2.));
        sxx += (s/12. + log(1./sqrt(3.))/log(2.)) *
               (s/12. + log(1./sqrt(3.))/log(2.));
        sy += (log((double)(*po_Z)(e,s))/log(2.));
        sxy += (s/12. + log(1./sqrt(3.))/log(2.)) *
               (log((double)(*po_Z)(e,s))/log(2.));
      }
    }
    delta = ai_NbScale*sxx-sx*sx;
    float a = (ai_NbScale*sxy - sx*sy) / delta;
    float b = (sxx*sy - sx*sxy) / delta;
    for (int s=i_IndMinScale; s<i_IndMaxScale; s++) {
       float y     = log((double)(*po_Z)(e,s))/log(2.);
       float y_est = a * (s/12. + log(1./sqrt(3.))/log(2.)) + b;
       float diff  = (y-y_est)*(y-y_est);
       sigma_in += diff;
       //sigma_in += sigma_in + pow ( (log((double)(*po_Z)(e,s))/log(2.) 
       //            - ( a * (s/12. + log(1./sqrt(3.))/log(2.)) + b)), 2.);
    }
    sigma_in = sigma_in /(i_IndMaxScale-i_IndMinScale);
                   
    float sigma = sqrt(sigma_in * ai_NbScale / delta); //(eq 15.2.9 Numerical receipes)

    //cout << "q=" << (*po_q)(e) <<", a="<<a<< endl;
 
    (*po_Tau)(e) = a;(*po_SigTau)(e) = 3*sigma;
  }
}
 
void to_ThermoDynAnalysis::display_tau () {
  
  char atc_FileName[80];
  fltarray ao_XTau (i_NbExp);
  fltarray ao_YTauTheo (i_NbExp);

  for (int e=0; e<i_NbExp; e++) {
    ao_XTau(e) = (*po_q)(e);
    ao_YTauTheo(e) = ((*po_q)(e)-1)*log(2.)/log(3.);
  }    
  sprintf (atc_FileName, "X_Tau");
  io_1d_write_data(atc_FileName, ao_XTau);

  sprintf (atc_FileName, "Y_Tau");
  io_1d_write_data(atc_FileName, *po_Tau);
  
  sprintf (atc_FileName, "Sig_Tau");
  io_1d_write_data(atc_FileName, *po_SigTau);
  
  //sprintf (atc_FileName, "Y_TauTheo");
  //fits_write_fltarr(atc_FileName, ao_YTauTheo);
}




void to_ThermoDynAnalysis::compute_frac_spectrum (fltarray* ppo_Alpha) {

  float af_Temp;

  if (ppo_Alpha==NULL) init_default_alpha();
  else po_Alpha = ppo_Alpha;
  i_NbAlpha = po_Alpha->n_elem();
  
  // is Tau(q)=f(q) a line ?
  float sx=0, sxx=0, sy=0, sxy=0, sigma_in=0, delta;  
  for (int e=0; e<i_NbExp; e++) {
     sx += (*po_q)(e);
     sxx += (*po_q)(e)*(*po_q)(e);
     sy += (*po_Tau)(e);
     sxy += (*po_q)(e) * (*po_Tau)(e);
  }  
  delta = i_NbExp*sxx-sx*sx;
  float a = (i_NbExp*sxy - sx*sy) / delta;
  float b = (sxx*sy - sx*sxy) / delta; 
  for (int e=0; e<i_NbExp; e++) 
    sigma_in += ((*po_Tau)(e) - ( a * (*po_q)(e) + b )) * ((*po_Tau)(e) - ( a * (*po_q)(e) + b )) ;
  sigma_in /= i_NbExp;
  float sigma = sqrt(sigma_in * i_NbExp / delta); //(eq 15.2.9 Numerical receipes)    
  float sigma_rel =   sigma/a;
                             
  if (ABS(sigma_rel) < 5e-3) {        // !! on peut tester aussi sur b..., 
                                        // changer le seuil...
    po_FracSpectrum = new fltarray (1);
    delete(po_Alpha);
    po_Alpha = new fltarray (1);
    (*po_FracSpectrum)(0) = -b;
    (*po_Alpha)(0) = a;
  
  } else {
    po_FracSpectrum = new fltarray (i_NbAlpha);
    po_ValQMin = new fltarray (i_NbAlpha);

    for (int alpha=0; alpha<i_NbAlpha; alpha++) {

      (*po_FracSpectrum)(alpha) = FRACFLTMAX;
      for (int e=0; e<i_NbExp; e++) {

        if ( (af_Temp = (*po_q)(e) * ((*po_Alpha)(alpha)) - (*po_Tau)(e)) 
	     < (*po_FracSpectrum)(alpha) ) {
	  (*po_FracSpectrum)(alpha) = af_Temp;
	  (*po_ValQMin)(alpha) = (*po_q)(e);
        }
      }
      //cout << "alpha="<<(*po_Alpha)(alpha)<<", q="<<(*po_ValQMin)(alpha)
      //     <<", D(alpha)="<<(*po_FracSpectrum)(alpha)<<endl;
    }
  }
}


void to_ThermoDynAnalysis::display_frac_spectrum () {

  char atc_FileName[80];

  sprintf (atc_FileName, "X_FracSpectrum");
  io_1d_write_data(atc_FileName, *po_Alpha);

  sprintf (atc_FileName, "Y_FracSpectrum");
  io_1d_write_data(atc_FileName, *po_FracSpectrum);
  //sprintf (atc_FileName, "ValQMin");
  //fits_write_fltarr(atc_FileName, *po_ValQMin);
}


void to_ThermoDynAnalysis::alpha_for_q_equals_zero () {

  
  float xf_MinAlpha = 0;
  float xf_MaxAlpha = 1;
 
  while (xf_MaxAlpha-xf_MinAlpha>1e-04) {
  
    int ai_IndMin=0;
    int ai_IndMax=0;
    for (int alpha=0; alpha<i_NbAlpha-1; alpha++) {
    
      if (ABS((*po_ValQMin)(alpha)) <= EPS_FRACTAL) {
        ai_IndMin = ai_IndMax = alpha;
        break;
      }
      
      if (   (*po_ValQMin)(alpha) > 0 && (*po_ValQMin)(alpha+1) < 0
          || (*po_ValQMin)(alpha) < 0 && (*po_ValQMin)(alpha+1) > 0) {
        ai_IndMin = alpha; ai_IndMax = alpha+1;
	break;
      }
	
    }
      
    cout << "Min    : " << (*po_Alpha)(ai_IndMin) << ", q=" <<
                           (*po_ValQMin)(ai_IndMin) << endl;
    cout << "Max    : " << (*po_Alpha)(ai_IndMax) << ", q=" <<
                           (*po_ValQMin)(ai_IndMax) << endl;
			 
    if (ai_IndMin == ai_IndMax) break;
			 
    xf_MinAlpha = ((*po_Alpha)(ai_IndMin) > (*po_Alpha)(ai_IndMax) ?
                   (*po_Alpha)(ai_IndMax) : (*po_Alpha)(ai_IndMin));
    xf_MaxAlpha = ((*po_Alpha)(ai_IndMin) > (*po_Alpha)(ai_IndMax) ?
                   (*po_Alpha)(ai_IndMin) : (*po_Alpha)(ai_IndMax));		       			 
    
    int ai_NbAlpha = 100;
    fltarray ao_Alpha (ai_NbAlpha+1);
    for (int i=0; i<ai_NbAlpha+1; i++) 
      ao_Alpha(i) = xf_MinAlpha + (xf_MaxAlpha-xf_MinAlpha)/ai_NbAlpha*i;
    
    compute_frac_spectrum (&ao_Alpha);
  }
     
}











