/******************************************************************************
**                   Copyright (C) 1994 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 3.1
**
**    Author: Jean-Luc Starck
**
**    Date:  96/05/02
**    
**    File:  IM_Maxima.cc
**
*******************************************************************************
**
**    DESCRIPTION  maxima detection  
**    ----------- 
**                 
**  void im_detect(Ifloat &Imag_in, Ifloat &Imag_Detect, 
**                 float **x_ef, float **y_ef, float **f_ef, int *nmax)
**
**  Imag_in: input image
**  Imag_detect: image which contains only the maxima
**  x_ef: x position array
**  y_ef: y position array
**  nmax: number of maxima
**
*********************************************************************/


#include "IM_Obj.h"

void im_detect(Ifloat &Imag_in, Ifloat &Imag_Detect, 
               float **x_ef, float **y_ef,  float **f_ef, int *nmax)
{

  float    *x_e, *y_e;
  float	   *f_e;
  int   *suivi0_x, *suivi0_y, *suivi_x, *suivi_y;
  int   lig, col, lig1, lig2, lig3;
  int   i, j, taille_allocation ;
  int Nl = Imag_in.nl();
  int Nc = Imag_in.nc();

 taille_allocation = Nl*Nc;

 suivi0_x = new int [Nc];
 suivi0_y = new int [Nc];
 suivi_x = new int [Nc];
 suivi_y = new int [Nc];

 x_e = new float [taille_allocation];
 y_e = new float [taille_allocation];
 f_e = new float [taille_allocation];

  *nmax = 0;

  for (i = 0; i < Nl; i++) 
  for (j = 0; j < Nc; j++) Imag_Detect(i,j) = 0.;

  for (i = 0; i < Nc; i++) 
  {
    suivi0_x[i] = 1;
    suivi0_y[i] = 1;
  }

 
  for (lig = 0; lig < Nl; lig++)
  {
    lig1 = lig - 1;
    lig2 = lig;
    lig3 = lig + 1;

    if (lig1 < 0)        lig1 = lig;
    if (lig3 >=  Nl) lig3 = lig;

    for (col = 0; col < Nc; col++)
    {
      if (col == 0) 
      {
	suivi_x[col] = 1;
	if (lig == 1) 
        {
	  suivi_y[col] = 1;
	  goto test;
	}
 /*a*/  if (Imag_in(lig2,col) < Imag_in(lig1,col))
        {
	  suivi_y[col] = 0;
	  goto cont_col;
	}
 /*b*/	if (Imag_in(lig2,col) > Imag_in(lig1,col)) 
        {
	  suivi_y[col] = 1;
	  goto test;
	}
 /*c*/	if (Imag_in(lig2,col) == Imag_in(lig1,col)) 
        {
	  if (suivi0_y[col] == 0) 
          {
	    suivi_y[col] = 0;
	    goto cont_col;
	  } 
          else 
          {
	    suivi_y[col] = 1;
	    goto test;
	  }
	}
      }

      if (col != 0) 
      {
	
  /*1*/	 if (Imag_in(lig2,col) < Imag_in(lig2,col-1)) 
         {
	  suivi_x[col] = 0;
	  if (lig == 0) 
          {
	    suivi_y[col] = 1;
	    goto cont_col;
	  }

   /*a*/  if (Imag_in(lig2,col) < Imag_in(lig1,col)) 
          {
	    suivi_y[col] = 0;
	    goto cont_col;
	  }

   /*b*/  if (Imag_in(lig2,col) > Imag_in(lig1,col)) 
          {
	    suivi_y[col] = 1;
	    goto cont_col;
	  }

   /*c*/  if (Imag_in(lig2,col) == Imag_in(lig1,col)) 
          {
	    if (suivi0_y[col] == 1) suivi_y[col] = 1;
	    else suivi_y[col] = 0;
	    goto cont_col;
	  }
	}

/*2*/	if (Imag_in(lig2,col) == Imag_in(lig2,col-1)) 
        {
	  if (suivi_x[col - 1] == 1)  suivi_x[col] = 1;
	  else  suivi_x[col] = 0;

	  if (lig == 0) 
          {
	    suivi_y[col] = 1;
	    if (suivi_x[col] == 1)  goto test;
	    else  goto cont_col;
	  }
   /*a*/  if (Imag_in(lig2,col) < Imag_in(lig1,col)) 
          {
	    suivi_y[col] = 0;
	    goto cont_col;
	  }
   /*b*/  if (Imag_in(lig2,col) > Imag_in(lig1,col)) 
          {
	    suivi_y[col] = 1;
	    if (suivi_x[col] == 0)
	      goto cont_col;
	    else
	      goto test;
	  }
  /*c*/	  if (Imag_in(lig2,col) == Imag_in(lig1,col)) 
          {
	    if (suivi0_y[col] == 0) 
            {
	      suivi_y[col] = 0;
	      goto cont_col;
	    } 
            else 
            {
	      suivi_y[col] = 1;
	      if (suivi_x[col] == 1)
		goto test;
	      else
		goto cont_col;
	    }
	  }

	}

/*3*/	if (Imag_in(lig2,col) > Imag_in(lig2,col-1))
        {
	  suivi_x[col] = 1;
	  if (lig == 0) 
          {
	    suivi_y[col] = 1;
	    goto test;
	  }

   /*a*/ if (Imag_in(lig2,col) < Imag_in(lig1,col)) 
         {
	    suivi_y[col] = 0;
	    goto cont_col;
	  }

  /*b*/	  if (Imag_in(lig2,col) > Imag_in(lig1,col))
          {
	    suivi_y[col] = 1;
	    goto test;
	  }

  /*c*/	  if (Imag_in(lig2,col) == Imag_in(lig1,col)) 
          {
	    if (suivi0_y[col] == 0) 
            {
	      suivi_y[col] = 0;
	      goto cont_col;
	    } 
            else 
            {
	      suivi_y[col] = 1;
	      goto test;
	    }
	  }

	}

      }

test:

  /* ------------------------------------------------------------------ */
  /*----         debut de la procedure de test                      --- */
  /* ------------------------------------------------------------------ */

      if (col != 0) 
      {
	if (col != (Nc - 1)) 
        {
         /* ----------------------------------------------------------- */
	 /* ---    test sur les elements horizontaux et verticaux   --- */
         /* ----------------------------------------------------------- */

	  if (Imag_in(lig2,col) < Imag_in(lig2,col-1)) 
	    goto cont_col;
	  if (Imag_in(lig2,col) < Imag_in(lig1,col)) 
            goto cont_col;
	  if (Imag_in(lig2,col) <= Imag_in(lig2,col+1)) 
	    goto cont_col;
	  if (Imag_in(lig2,col) <= Imag_in(lig3,col)) 
	    goto cont_col;

         /* ----------------------------------------------------------- */
	 /* ---            test sur les elements diagonaux          --- */
         /* ----------------------------------------------------------- */

	  if ( (Imag_in(lig2,col) < Imag_in(lig1,col-1)) && 
               (Imag_in(lig2,col) < Imag_in(lig1,col+1)) ) 
	    goto cont_col;
	  if (Imag_in(lig2,col) < Imag_in(lig1,col-1))
	    goto cont_col;
	  if (Imag_in(lig2,col) < Imag_in(lig1,col+1))
	    goto cont_col;
	  if ( (Imag_in(lig2,col) <= Imag_in(lig3,col-1)) && 
               (Imag_in(lig2,col) <= Imag_in(lig3,col+1)) )
	    goto cont_col;
	  if (Imag_in(lig2,col) <= Imag_in(lig3,col-1)) 
	    goto cont_col;
	  if (Imag_in(lig2,col) <= Imag_in(lig3,col+1)) 
	    goto cont_col;

	} 
        else 
       {

         /* ----------------------------------------------------------- */
	 /* ---    test sur les elements horizontaux et verticaux   --- */
         /* ----------------------------------------------------------- */

	  if (Imag_in(lig2,col) < Imag_in(lig2,col-1)) 
	    goto cont_col;
	  if (Imag_in(lig2,col) < Imag_in(lig1,col)) 
	    goto cont_col;
	  if (Imag_in(lig2,col) <= Imag_in(lig3,col)) 
	    goto cont_col;

         /* ----------------------------------------------------------- */
	 /* ---           test sur les elements diagonaux           --- */
         /* ----------------------------------------------------------- */
	  if (Imag_in(lig2,col) < Imag_in(lig1,col-1)) 
	    goto cont_col;
	  if (Imag_in(lig2,col) <= Imag_in(lig3,col-1)) 
	    goto cont_col;
          }
        } 
        else 
        {

        /* ----------------------------------------------------------- */
	/* ---    test sur les elements horizontaux et verticaux   --- */
        /* ----------------------------------------------------------- */
	if (Imag_in(lig2,col) < Imag_in(lig1,col)) 
	  goto cont_col;
	if (Imag_in(lig2,col) <= Imag_in(lig2,col+1))
	  goto cont_col;
	if (Imag_in(lig2,col) <= Imag_in(lig3,col)) 
	  goto cont_col;

        /* ----------------------------------------------------------- */
	/* ---           test sur les elements diagonaux           --- */
        /* ----------------------------------------------------------- */
	if (Imag_in(lig2,col) < Imag_in(lig1,col+1)) 
	  goto cont_col;
	if (Imag_in(lig2,col) <= Imag_in(lig3,col+1)) 
	  goto cont_col;

      }

      Imag_Detect(lig2,col) = Imag_in(lig2,col);
      x_e[*nmax] = col;
      y_e[*nmax] = lig;
      f_e[*nmax] = Imag_in(lig2,col);

/*
      printf(" y_e[%3d] : %4.0f \t x_e[%3d] : %4.0f \t f_e[%3d] : %4e",
              *nmax, y_e[*nmax], *nmax, x_e[*nmax], *nmax, f_e[*nmax]);
      printf("  \n");*/

      (*nmax)++;


cont_col: ;

    }  /* end col */

    for (i = 0; i < Nc; i++) suivi0_y[i] = suivi_y[i];

  }  /* end ligne */

 delete [] suivi0_x;
 delete [] suivi0_y;
 delete [] suivi_x;
 delete [] suivi_y;

 *x_ef = new float [*nmax];
 *y_ef = new float [*nmax];
 *f_ef = new float [*nmax];
 
 for (i=0; i< (*nmax); i++)
   {
     (*x_ef)[i] = x_e[i] ;
     (*y_ef)[i] = y_e[i] ;
     (*f_ef)[i] = f_e[i] ;
   }

delete [] x_e;
delete [] y_e;
delete [] f_e;
}

