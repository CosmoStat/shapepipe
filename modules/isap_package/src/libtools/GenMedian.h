

#include "DefMath.h"

#define PIX_SORT(a,b) { if ((a)>(b)) PIX_SWAP((a),(b)); }
#define PIX_SWAP(a,b) { elem_type temp=(a);(a)=(b);(b)=temp; }

/*----------------------------------------------------------------------------
 * Function :   opt_med3()
 * In       :   pointer to array of 3 pixel values
 * Out      :   a elem_type
 * Job      :   optimized search of the median of 3 pixel values
 * Notice   :   found on sci.image.processing
 *              cannot go faster unless assumptions are made
 *              on the nature of the input signal.
 *--------------------------------------------------------------------------*/

elem_type opt_med3(elem_type  *p)
{
    PIX_SORT(p[0],p[1]) ; PIX_SORT(p[1],p[2]) ; PIX_SORT(p[0],p[1]) ;
    return(p[1]) ;
}

/*----------------------------------------------------------------------------
 * Function :   opt_med5()
 * In       :   pointer to array of 5 pixel values
 * Out      :   a elem_type
 * Job      :   optimized search of the median of 5 pixel values
 * Notice   :   found on sci.image.processing
 *              cannot go faster unless assumptions are made
 *              on the nature of the input signal.
 *--------------------------------------------------------------------------*/

elem_type opt_med5(elem_type  *p)
{
    PIX_SORT(p[0],p[1]) ; PIX_SORT(p[3],p[4]) ; PIX_SORT(p[0],p[3]) ;
    PIX_SORT(p[1],p[4]) ; PIX_SORT(p[1],p[2]) ; PIX_SORT(p[2],p[3]) ;
    PIX_SORT(p[1],p[2]) ; return(p[2]) ;
}

/*----------------------------------------------------------------------------
 * Function :   opt_med7()
 * In       :   pointer to array of 7 pixel values
 * Out      :   a elem_type
 * Job      :   optimized search of the median of 7 pixel values
 * Notice   :   found on sci.image.processing
 *              cannot go faster unless assumptions are made
 *              on the nature of the input signal.
 *--------------------------------------------------------------------------*/

elem_type opt_med7(elem_type  *p)
{
    PIX_SORT(p[0], p[5]) ; PIX_SORT(p[0], p[3]) ; PIX_SORT(p[1], p[6]) ;
    PIX_SORT(p[2], p[4]) ; PIX_SORT(p[0], p[1]) ; PIX_SORT(p[3], p[5]) ;
    PIX_SORT(p[2], p[6]) ; PIX_SORT(p[2], p[3]) ; PIX_SORT(p[3], p[6]) ;
    PIX_SORT(p[4], p[5]) ; PIX_SORT(p[1], p[4]) ; PIX_SORT(p[1], p[3]) ;
    PIX_SORT(p[3], p[4]) ; 
    return (p[3]) ;
}

/*----------------------------------------------------------------------------
 * Function :   opt_med9()
 * In       :   pointer to an array of 9 elem_types
 * Out      :   a elem_type
 * Job      :   optimized search of the median of 9 elem_types
 * Notice   :   in theory, cannot go faster without assumptions on the
 *              signal.
 *              Formula from:
 *              XILINX XCELL magazine, vol. 23 by John L. Smith
 *
 *              The input array is modified in the process
 *              The result array is guaranteed to contain the median
 *              value
 *              in middle position, but other elements are NOT sorted.
 *--------------------------------------------------------------------------*/

elem_type opt_med9(elem_type  *p)
{
    PIX_SORT(p[1], p[2]) ; PIX_SORT(p[4], p[5]) ; PIX_SORT(p[7], p[8]) ;
    PIX_SORT(p[0], p[1]) ; PIX_SORT(p[3], p[4]) ; PIX_SORT(p[6], p[7]) ;
    PIX_SORT(p[1], p[2]) ; PIX_SORT(p[4], p[5]) ; PIX_SORT(p[7], p[8]) ;
    PIX_SORT(p[0], p[3]) ; PIX_SORT(p[5], p[8]) ; PIX_SORT(p[4], p[7]) ;
    PIX_SORT(p[3], p[6]) ; PIX_SORT(p[1], p[4]) ; PIX_SORT(p[2], p[5]) ;
    PIX_SORT(p[4], p[7]) ; PIX_SORT(p[4], p[2]) ; PIX_SORT(p[6], p[4]) ;
    PIX_SORT(p[4], p[2]) ; 
    return(p[4]) ;
}

#undef PIX_SORT
#undef PIX_SWAP



#define ELEM_SWAP(a,b) { register elem_type t=(a);(a)=(b);(b)=t; }

/*---------------------------------------------------------------------------
   Function :   kth_smallest()
   In       :   array of elements, # of elements in the array, rank k
   Out      :   one element
   Job      :   find the kth smallest element in the array
   Notice   :   use the median() macro defined below to get the median. 

                Reference:

                  Author: Wirth, Niklaus 
                   Title: Algorithms + data structures = programs 
               Publisher: Englewood Cliffs: Prentice-Hall, 1976 
    Physical description: 366 p. 
                  Series: Prentice-Hall Series in Automatic Computation 

 ---------------------------------------------------------------------------*/


elem_type kth_smallest(elem_type a[], int n, int k)
{
	register int i,j,l,m ;
	register elem_type x ;

    l=0 ; m=n-1 ;
    while (l<m) {
        x=a[k] ;
        i=l ;
        j=m ;
        do {
            while (a[i]<x) i++ ;
            while (x<a[j]) j-- ;
            if (i<=j) {
                ELEM_SWAP(a[i],a[j]) ;
                i++ ; j-- ;
            }
        } while (i<=j) ;
        if (j<k) l=i ;
        if (k<i) m=j ;
    }
    return a[k] ;
}

elem_type get_median(elem_type a[], int n)
{
   return kth_smallest(a,n,(((n)&1)?((n)/2):(((n)/2)-1)));
}

elem_type abs_kth_smallest(elem_type a[], int n, int k)
{
	register int i,j,l,m ;
	register elem_type x ;

    l=0 ; m=n-1 ;
    while (l<m) {
        x=ABS(a[k]);
        i=l ;
        j=m ;
        do {
            while (ABS(a[i])<x) i++ ;
            while (x<ABS(a[j])) j-- ;
            if (i<=j) {
                ELEM_SWAP(a[i],a[j]) ;
                i++ ; j-- ;
            }
        } while (i<=j) ;
        if (j<k) l=i ;
        if (k<i) m=j ;
    }
    return ABS(a[k]) ;
}

elem_type get_abs_median(elem_type a[], int n)
{
   return abs_kth_smallest(a,n,(((n)&1)?((n)/2):(((n)/2)-1)));
}

#undef ELEM_SWAP


/******************************** hmedian ***********************************/

elem_type hmedian(elem_type  *ra, int n)
/* Median using Heapsort algorithm  (based on Num.Rec algo.).
   Warning: changes the order of data! 
*/
{
   int		l, j, ir, i;
   elem_type	rra;

  if (n<2)
    return *ra;
  ra--;
  for (l = ((ir=n)>>1)+1;;)
    {
    if (l>1)
      rra = ra[--l];
    else
      {
      rra = ra[ir];
      ra[ir] = ra[1];
      if (--ir == 1)
        {
        ra[1] = rra;
        return n&1? ra[n/2+1] : (elem_type)((ra[n/2]+ra[n/2+1])/2.0);
        }
      }
    for (j = (i=l)<<1; j <= ir;)
      {
      if (j < ir && ra[j] < ra[j+1])
        ++j;
      if (rra < ra[j])
        {
        ra[i] = ra[j];
        j += (i=j);
        }
      else
        j = ir + 1;
      }
    ra[i] = rra;
    }

/* (the 'return' is inside the loop!!) */
  }

