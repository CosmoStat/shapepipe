
// pyramidal median transform with integer
#include "GlobalInc.h"
// #include "NR.h"
#include "MR_Obj.h"
#include "OptMedian.h"

#define MAX_NBR_SCALE 10

/* coefficients array size */
#define SIZE_T_COEF_INTER 10 

#define C1 .00245
#define C2 -.00915
#define C3 .03414
#define C4 -0.12720
#define C5 0.60048

#define CS ((C1+C2+C3+C4+C5)*2.)
#define C1n (C1/CS)
#define C2n (C2/CS)
#define C3n (C3/CS)
#define C4n (C4/CS)
#define C5n (C5/CS)
#define ISWAP(a,b) temp=(a);(a)=(b);(b)=temp;

static float T_Coeff_Interp[SIZE_T_COEF_INTER] =
                              {C1n,C2n,C3n,C4n,C5n,C5n,C4n,C3n,C2n,C1n};

#define BSPLINE 1

void mri_pos_ind_medpyr (int *Tab_Nl, int *Tab_Col, int *Tab_Pos, 
                         int Nl, int Nc, int Nbr_Plan)
{
    int i;
    Tab_Nl [0] = Nl;
    Tab_Col [0] = Nc;
    Tab_Pos [0] = 0;

    for (i = 1; i < Nbr_Plan; i++)
    {
        Tab_Nl [i] = (Tab_Nl [i-1] - 1) / 2 + 1;
        Tab_Col [i] = (Tab_Col [i-1] - 1) / 2 + 1;
        Tab_Pos [i] = Tab_Pos [i-1] + Tab_Nl [i-1] * Tab_Col [i-1];
    }
}

/****************************************************************************/

int mri_size_medpyr (int Nl, int Nc, int Nbr_Plan)
{
    int Size,i;
    int Nl1,Nc1,Nl2=Nl,Nc2=Nc,Pos1,Pos2=0;

    Nl1 = Nl;
    Nc1 = Nc;
    Pos1 = 0;
    for (i = 1; i < Nbr_Plan; i++)
    {
        Nl2 = (Nl1 - 1) / 2 + 1;
        Nc2 = (Nc1 - 1) / 2 + 1;
        Pos2  = Pos1 + Nl1  * Nc1;
        Nl1 = Nl2;
        Nc1 = Nc2;
        Pos1 = Pos2;
    }
    
    Size = Pos2 + Nl2 * Nc2;
    return (Size);
}


/*****************************************************************************/

/* interpolates an image by a factor 2 */
static void im_increase_size_2 (float *Pict_in, int Nl2, int Nc2, float *Pict_out, int Nl, int Nc)
{
    int i,j,k,Pos_I;
    int ind_dep_pyr,Ind_Pict;
    int pi, pj, indk, lastnc;

    /* Pixels at node interpolation */
    for (i = 0; i < Nl2; i ++)
    {
        pi = 2*i*Nc;
        if (2*i < Nl)
        {
          Pos_I = i * Nc2;
          for (j = 0; j < Nc2; j ++)
          {
            pj = 2*j;
            if (pj < Nc) Pict_out [pi + pj] = Pict_in [Pos_I + j];
          }
        }
    }

    /* pixel node in column */
    for (i = 0; i < Nl2; i ++)
    {
        pi = 2*i*Nc;
        if (2*i < Nl)
        {    
          Pos_I = i * Nc2;
          for (j = 0; j < Nc2; j ++)
          {
            pj = 2*j+1;
            if (pj < Nc)
            {
               Ind_Pict = pi + pj;
               Pict_out [Ind_Pict] = 0.;
               for (k = 0; k < SIZE_T_COEF_INTER; k ++)
               {
                indk = j - 4 + k;
                if (indk < 0) indk = 0;
                else if (indk >= Nc2) indk = Nc2-1;
                Pict_out [Ind_Pict] += T_Coeff_Interp[k]*Pict_in [Pos_I+indk];
               }
            }
          }
        }
    }
    
    /* pixel node in line */
    for (i = 0; i < Nl2; i ++)
    {
        pi = (2*i+1)*Nc;
        if (2*i+1 < Nl)
        {
            for (j = 0; j < Nc2; j ++)
            {
               pj = 2*j;
               if (pj < Nc)
               {
                  Ind_Pict = pi + pj;
                  Pict_out [Ind_Pict] = 0.;
                  for (k = 0; k < SIZE_T_COEF_INTER; k ++)
                  {
                     indk = i - 4 + k;
                     if (indk < 0) indk = 0;
                     else if (indk >= Nl2) indk = Nl2-1;
                     Pict_out [Ind_Pict] += 
                              T_Coeff_Interp [k]*Pict_in [indk * Nc2 + j];
                  }
               }
            }
        }
    }

    for (i = 0; i < Nl2; i ++)
    {
        pi = (2*i+1)*Nc;
        if (2*i+1 < Nl)
        {
            /* pixel node in line and column */
           if (Nc % 2 == 0) lastnc = Nc - 2;
           else lastnc = Nc - 1;
           for (j = 0; j < Nc2; j ++)
           {
               pj = 2*j+1;
               if (pj < Nc)
               {
                  Ind_Pict = pi  + pj;
                  Pict_out [Ind_Pict] = 0.;
                  for (k = 0; k < SIZE_T_COEF_INTER; k ++)
                  {
                     indk = 2*(j - 4 + k);
                     ind_dep_pyr = pi  + indk;
                     if (indk < 0) ind_dep_pyr = pi;
                     else if (indk >= Nc) ind_dep_pyr = pi+lastnc;
                     Pict_out [Ind_Pict] +=
                                 T_Coeff_Interp [k]*Pict_out [ind_dep_pyr];
                  }
               }
           }
       }
    }
}


/***************************************************************************/

/* interpolates an image by a factor 2 */
static void imI_increase_size_2 (int *Pict_in, int Nl2, int Nc2, 
                                 float *Pict_out, int Nl, int Nc)
{
    int i,j,k,Pos_I;
    int ind_dep_pyr,Ind_Pict;
    int pi, pj, indk, lastnc;

    /* Pixels at node interpolation */
    for (i = 0; i < Nl2; i ++)
    {
        pi = 2*i*Nc;
        if (2*i < Nl)
        {
          Pos_I = i * Nc2;
          for (j = 0; j < Nc2; j ++)
          {
            pj = 2*j;
            if (pj < Nc) Pict_out [pi + pj] = Pict_in [Pos_I + j];
          }
        }
    }

    /* pixel node in column */
    for (i = 0; i < Nl2; i ++)
    {
        pi = 2*i*Nc;
        if (2*i < Nl)
        {    
          Pos_I = i * Nc2;
          for (j = 0; j < Nc2; j ++)
          {
            pj = 2*j+1;
            if (pj < Nc)
            {
               Ind_Pict = pi + pj;
               Pict_out [Ind_Pict] = 0.;
               for (k = 0; k < SIZE_T_COEF_INTER; k ++)
               {
                indk = j - 4 + k;
                if (indk < 0) indk = 0;
                else if (indk >= Nc2) indk = Nc2-1;
                Pict_out [Ind_Pict] += T_Coeff_Interp[k]*Pict_in [Pos_I+indk];
               }
            }
          }
        }
    }
    
    /* pixel node in line */
    for (i = 0; i < Nl2; i ++)
    {
        pi = (2*i+1)*Nc;
        if (2*i+1 < Nl)
        {
            for (j = 0; j < Nc2; j ++)
            {
               pj = 2*j;
               if (pj < Nc)
               {
                  Ind_Pict = pi + pj;
                  Pict_out [Ind_Pict] = 0.;
                  for (k = 0; k < SIZE_T_COEF_INTER; k ++)
                  {
                     indk = i - 4 + k;
                     if (indk < 0) indk = 0;
                     else if (indk >= Nl2) indk = Nl2-1;
                     Pict_out [Ind_Pict] += 
                              T_Coeff_Interp [k]*Pict_in [indk * Nc2 + j];
                  }
               }
            }
        }
    }

    for (i = 0; i < Nl2; i ++)
    {
        pi = (2*i+1)*Nc;
        if (2*i+1 < Nl)
        {
            /* pixel node in line and column */
           if (Nc % 2 == 0) lastnc = Nc - 2;
           else lastnc = Nc - 1;
           for (j = 0; j < Nc2; j ++)
           {
               pj = 2*j+1;
               if (pj < Nc)
               {
                  Ind_Pict = pi  + pj;
                  Pict_out [Ind_Pict] = 0.;
                  for (k = 0; k < SIZE_T_COEF_INTER; k ++)
                  {
                     indk = 2*(j - 4 + k);
                     ind_dep_pyr = pi  + indk;
                     if (indk < 0) ind_dep_pyr = pi;
                     else if (indk >= Nc) ind_dep_pyr = pi+lastnc;
                     Pict_out [Ind_Pict] +=
                                 T_Coeff_Interp [k]*Pict_out [ind_dep_pyr];
                  }
               }
           }
       }
    }
}

/***************************************************************************/
/*
static void im_bilinear_interp (float *INimage, int INlin, int INcol, float *OUTimage, int OUTlin, int OUTcol)
{		
    register int i, j, ii, jj, i1, j1, i2, j2;
    int *indice_lin, *indice_col;
    float delta, delta_x, delta_y, delta_yx;

    if ((INlin > OUTlin) || (INcol > OUTcol))
    {
        fprintf (stderr, "Error in im_bilinear_interp: bad image size\n");
        exit (-1);
    }

    indice_lin = (int *) malloc ((unsigned) INlin * sizeof (int));
    indice_col = (int *) malloc ((unsigned) INcol * sizeof (int));

    for (i = 0; i < INlin; i++)  indice_lin [i] = (i * (OUTlin-1)) / (INlin-1);
    for (i = 0; i < INcol; i++)  indice_col [i] = (i  *(OUTcol-1)) / (INcol-1);

    for (i=0; i<INlin; i++)
    {
        ii = indice_lin[i];
 	for (j=0; j<INcol; j++)
        {
            jj = indice_col [j];
            OUTimage [ii * OUTcol + jj] = INimage [i * INcol + j];
        }
    }


    for (i=0; i<INlin; i++)
    {
        ii = indice_lin[i];
        for (j=0; j<(INcol-1); j++)
        {
            j1 = indice_col [j];
            j2 = indice_col [j+1];
            delta_x = (float) j2-j1;
            delta_y = OUTimage [ii* OUTcol +j2]-OUTimage [ii* OUTcol +j1];
            delta_yx = delta_y/delta_x;
            for (jj=j1; jj<j2; jj++)
            {
                delta = (float) jj-j1;
                OUTimage[ii* OUTcol +jj] = OUTimage [ii* OUTcol +j1]
                                                  + (delta * delta_yx);
            }
        }
     }

     for (jj = 0; jj < OUTcol; jj++)
     {
         for (i = 0; i < (INlin-1); i++)
         {
             i1 = indice_lin[i];
             i2 = indice_lin[i+1];
             delta_x = (float) i2-i1;
             delta_y = OUTimage [i2 * OUTcol + jj] - OUTimage [i1 * OUTcol +jj];
             delta_yx = delta_y / delta_x;
             for (ii = i1; ii < i2; ii++)
             {
                 delta = (float) ii-i1;
                 OUTimage [ii* OUTcol +jj] = OUTimage[i1 * OUTcol +jj]
                                                  + (delta * delta_yx);
             }
         }
    }
    free ((char *) indice_lin);
    free ((char *) indice_col);
}
 */
/*****************************************************************************/

void imI_bilinear_interp (int *INimage, int INlin, int INcol, 
                          int *OUTimage, int OUTlin, int OUTcol)
{		
    register int i, j, ii, jj, i1, j1, i2, j2;
    int *indice_lin, *indice_col;
    float delta, delta_x, delta_y, delta_yx;
    int Val;

    if ((INlin > OUTlin) || (INcol > OUTcol))
    {
        fprintf (stderr, "Error in im_bilinear_interp: bad image size\n");
        exit (-1);
    }

    indice_lin = (int *) malloc ((unsigned) INlin * sizeof (int));
    indice_col = (int *) malloc ((unsigned) INcol * sizeof (int));

    for (i = 0; i < INlin; i++)  indice_lin [i] = (i * (OUTlin-1)) / (INlin-1);
    for (i = 0; i < INcol; i++)  indice_col [i] = (i  *(OUTcol-1)) / (INcol-1);

    for (i=0; i<INlin; i++)
    {
        ii = indice_lin[i];
 	for (j=0; j<INcol; j++)
        {
            jj = indice_col [j];
            OUTimage [ii * OUTcol + jj] = INimage [i * INcol + j];
        }
    }


    for (i=0; i<INlin; i++)
    {
        ii = indice_lin[i];
        for (j=0; j<(INcol-1); j++)
        {
            j1 = indice_col [j];
            j2 = indice_col [j+1];
            delta_x = (float) j2-j1;
            delta_y = OUTimage [ii* OUTcol +j2]-OUTimage [ii* OUTcol +j1];
            delta_yx = delta_y/delta_x;
            for (jj=j1; jj<j2; jj++)
            {
                delta = (float) (jj-j1)* delta_yx;
                Val = (delta >= 0) ? (int)(delta+0.5) : (int)(delta-0.5);
                OUTimage[ii* OUTcol +jj] = OUTimage [ii* OUTcol +j1] + Val;
            }
        }
     }

     for (jj = 0; jj < OUTcol; jj++)
     {
         for (i = 0; i < (INlin-1); i++)
         {
             i1 = indice_lin[i];
             i2 = indice_lin[i+1];
             delta_x = (float) i2-i1;
             delta_y = OUTimage [i2 * OUTcol + jj] - OUTimage [i1 * OUTcol +jj];
             delta_yx = delta_y / delta_x;
             for (ii = i1; ii < i2; ii++)
             {
                 delta = (float) (ii-i1)*delta_yx;
                 Val = (delta >= 0) ? (int)(delta+0.5) : (int)(delta-0.5);
                 OUTimage [ii* OUTcol +jj] = OUTimage[i1 * OUTcol +jj] + Val;
             }
         }
    }
    free ((char *) indice_lin);
    free ((char *) indice_col);
}
 
/****************************************************************************/
/*
static int iselect(int k, int n, int *arr)
{
  unsigned long i,ir,j,l,mid;
  int a, temp;
  
  l=0;
  ir = n-1;
  for(;;){
	 if(ir<=l+1){
		if (ir == l+1 && arr[ir] < arr[l]){
		  ISWAP(arr[l],arr[ir])
		}
		return arr[k];
	 }else {
		mid =(l+ir) >> 1;
		ISWAP(arr[mid],arr[l+1])
		if(arr[l+1] > arr[ir]){
		  ISWAP(arr[l+1],arr[ir])
		}
		if(arr[l] > arr[ir]){
		  ISWAP(arr[l],arr[ir])
		}
		if(arr[l+1] > arr[l]){
		  ISWAP(arr[l+1],arr[l])
		}
		i=l+1;
		j=ir;
		a=arr[l];
		for(;;){
		  do i++; while (arr[i] < a);
		  do j--; while (arr[j] > a);
		  if (j<i) break;
		  ISWAP(arr[i],arr[j])
		}
		arr[l]=arr[j];
		arr[j]=a;
		if (j>=k) ir=j-1;
		if (j<=k) l=i;
	 }
  }
}
*/
/****************************************************************************/

void mri_pyrmed (int *Ima, int *Pyramid, int *Tab_Nl, 
                         int *Tab_Col, int *Tab_Pos, int Nbr_Etap,
                         int MedianWindowSize)
{
    int Nl,Nc,Nls,Ncs,s,i,j,k,l,Size,Window_Size=MedianWindowSize;
    int *Pict, *Plan, *Plan_Next;
    unsigned long ind_fen;
    int *fenetre;
    int Window2 = (Window_Size - 1) / 2;
    int si, sj,ei,ej;
#if BSPLINE
    float *Pictf;
#else
    int *Pictf;
#endif

    fenetre = new int [Window_Size*Window_Size+1];
 
    Size = Tab_Nl [0]*Tab_Col [0];   
    Pict = i_alloc (Size);
#if BSPLINE
    Pictf = f_alloc(Size);
#else
    Pictf = i_alloc(Size);
#endif

    for (i=0;i<Size;i++) Pyramid[i] = Ima[i];


    for (s = 0; s < Nbr_Etap; s++)
    {
         Plan = Pyramid + Tab_Pos[s];
         Plan_Next = Pyramid + Tab_Pos[s+1];

         Nl = Tab_Nl [s];
         Nc = Tab_Col [s];
         Nls = Tab_Nl [s+1];
         Ncs = Tab_Col [s+1];

        for (i = 0; i < Nls; i++) 
        {
           si = 2*i - Window2;
           ei = 2*i + Window2;
           if (si < 0) si = 0;
           if (ei >= Nl) ei = Nl-1;
           for (j = 0; j < Ncs; j++) 
           {
               ind_fen = 0;
               sj = 2*j - Window2;
               ej = 2*j + Window2;
               if (sj < 0) sj = 0;
               if (ej >= Nc) ej = Nc-1;

               for (k = si; k <= ei; k ++)
               for (l = sj; l <= ej; l ++)
                          fenetre[ind_fen++] = Plan[k*Nc + l];
               // sort(ind_fen-1, fenetre);
               // Plan_Next[i*Ncs + j] = fenetre[ind_fen/2];
	       if (ind_fen == 9) Plan_Next[i*Ncs + j] = opt_med9(fenetre);  
              else  Plan_Next[i*Ncs + j] = get_median(fenetre, ind_fen);
            }
        }

        /* computes  the multiresolution coefficients */
#if BSPLINE
        imI_increase_size_2 (Plan_Next, Nls, Ncs,  Pictf, Nl, Nc);
        for (i = 0; i < Nl*Nc; i++) Plan [i] -= (int) (Pictf [i]+0.5);
#else
        imI_bilinear_interp (Plan_Next, Nls, Ncs, Pictf, Nl, Nc);
        for (i = 0; i < Nl*Nc; i++) Plan [i] -= Pictf [i];
#endif
    }
    i_free (Pict);
    free_buffer ((char *) Pictf);
    delete [] fenetre;
}

/****************************************************************************/

void mri_pyrmedian_transform (int *Data, int Nl, int Nc, int **Median, int Nbr_Plan, int MedianWindowSize)
{
    int Size;
    int Tab_Nl[MAX_NBR_SCALE], Tab_Col[MAX_NBR_SCALE], Tab_Pos[MAX_NBR_SCALE];

    Size = mri_size_medpyr (Nl, Nc, Nbr_Plan);
    *Median = i_alloc (Size);
    mri_pos_ind_medpyr (Tab_Nl,Tab_Col,Tab_Pos, Nl, Nc, Nbr_Plan);
    mri_pyrmed (Data, *Median, Tab_Nl, Tab_Col, Tab_Pos, 
                Nbr_Plan-1, MedianWindowSize);
}

/****************************************************************************/

void mri_pyrmedian_rec (int *Pyramid, int *PictInt, 
                              int Nl, int Nc, int Nbr_Plan, float *Tab_Level)
{
    int s,i,Nls, Ncs, *Plan;
    int Tab_Nl[MAX_NBR_SCALE], Tab_Col[MAX_NBR_SCALE], Tab_Pos[MAX_NBR_SCALE];
    int Nbr_Etap = Nbr_Plan - 1;
    float *PictIn;
    float *PictOut;

    PictIn = f_alloc (Nl*Nc);
    PictOut = f_alloc (Nl*Nc);
    mri_pos_ind_medpyr (Tab_Nl,Tab_Col,Tab_Pos, Nl, Nc, Nbr_Plan);

    /* copy the last scale to PictOut */
    s = Nbr_Etap;
    Nls = Tab_Nl [s];
    Ncs = Tab_Col[s];
    Plan = Pyramid + Tab_Pos[s];
    for (i = 0; i < Nls*Ncs; i++) PictIn[i] = (float) (Plan[i] * Tab_Level[s]);

    for (s = Nbr_Etap - 1; s >= 0; s--)
    {
        Nls = Tab_Nl [s];
        Ncs = Tab_Col[s];
        Plan = Pyramid + Tab_Pos[s];

        /* pyramid interpolation plane */
#if BSPLINE
       im_increase_size_2 (PictIn, Tab_Nl [s+1], Tab_Col[s+1], 
                            PictOut, Nls, Ncs); 
#else
        im_bilinear_interp (PictIn, Tab_Nl [s+1], Tab_Col[s+1], 
                            PictOut, Nls, Ncs);
#endif

        /* addition of the pyramid plane to the image */
        if (Tab_Level[s] > 1)
           for (i = 0; i < Nls*Ncs; i++) 
                       PictIn[i] = PictOut[i] + Plan[i]*Tab_Level[s];
        else
           for (i = 0; i < Nls*Ncs; i++) 
                       PictIn[i] = PictOut[i] + Plan[i];
    }

    for (i = 0; i < Nl*Nc; i++) PictInt[i] = iround(PictIn[i]);

    f_free (PictIn);
    f_free (PictOut);
}

/***************************************************************************/
  
void morphoi_noise_rec (int *Pyramid, int *PictInt, 
                              int Nl, int Nc, int Nbr_Plan)
{
    int s,i,Nls, Ncs, *Plan;
    int Tab_Nl[MAX_NBR_SCALE], Tab_Col[MAX_NBR_SCALE], Tab_Pos[MAX_NBR_SCALE];
    int Nbr_Etap = Nbr_Plan - 1;
    float *PictIn;
    float *PictOut;
    
    PictIn = f_alloc (Nl*Nc);
    PictOut = f_alloc (Nl*Nc);
    mri_pos_ind_medpyr (Tab_Nl,Tab_Col,Tab_Pos, Nl, Nc, Nbr_Plan);
 
    /* copy the last scale to PictOut */
    s = Nbr_Etap;
    Nls = Tab_Nl [s];
    Ncs = Tab_Col[s];
    Plan = Pyramid + Tab_Pos[s];
   
    for (i = 0; i < Nls*Ncs; i++) PictIn[i] = (float) Plan[i];
 
    for (s = Nbr_Etap - 1; s >= 0; s--)
    {
      Nls=Tab_Nl[s];
      Ncs=Tab_Col[s];
      
      /* pyramid interpolation plane */
#if BSPLINE
       im_increase_size_2 (PictIn, Tab_Nl [s+1], Tab_Col[s+1], 
                            PictOut, Nls, Ncs); 
#else
        im_bilinear_interp (PictIn, Tab_Nl [s+1], Tab_Col[s+1], 
                            PictOut, Nls, Ncs);
#endif
      if (s > 0) for (i = 0; i < Nls*Ncs; i++) PictIn[i] = PictOut[i];
   }
 
   
   for (i = 0; i < Nl*Nc; i++) PictInt[i] = iround(PictOut[i]);
 
    f_free (PictIn);
    f_free (PictOut);
}
/************************************************************************/
