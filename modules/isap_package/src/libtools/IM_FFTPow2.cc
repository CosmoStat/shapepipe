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
**    File:  IM_FftFloat.cc
**
*******************************************************************************
**
**    DESCRIPTION 
**    -----------    
**
**    this module contains the routines converning the Fourier transform
**		One should note that the inverse transforms usually set the input
**		fourier coefficients to the ifft result even in fft2d(image, imageFFT, -1)
**
**************************************************************************
**  
** ft_cf_any_power_of_2(dat,direction,length)
** complex_float *dat;
** int direction;
** int length;
**
** Takes the array pointed to by dat which is assumed to be of side-length
** length and performs a fourier transform of its contents
**
** External function calls
** None
** 
** INPUT dat = complex_floats images
**       direction = 1 ==>   directe Fourier transform
**       direction = -1 ==>  inverse Fourier transform
**       length = number of lines X number of columns
**
**************************************************************************
**
** void fft2d (const Icomplex &Im_in, const Icomplex &Im_out, int Dir)
** performs a fourier transform of Im_in
**   Im_out = FFT(Im_in) if Dir != -1
**   Im_out = INV_FFT(Im_in) if Dir = -1
**
** External function calls
** None
** 
** INPUT Im_in = input images
**       direction = 1 ==>   directe Fourier transform
**       direction = -1 ==>  inverse Fourier transform
** OUTPUT Im_out  = output images
**
**************************************************************************
**
** void fft1d (const Ifloat &Signal, const Icomplex_f & Signal_cf, int Dir)
** performs a fourier transform of Signal
**   Signal_cf = FFT(Signal) if Dir != -1
**       direction = 1 ==>   directe Fourier transform
**       direction = -1 ==>  inverse Fourier transform
**       length = number of elements
**
**************************************************************************
** 
** void fft1d (const Icomplex_f &Signal, const Icomplex_f & Signal_cf, int Dir)
** 
** performs a fourier transform of complexe Signal
**   Signal_cf = FFT(Signal) if Dir != -1
**       direction = 1 ==>   directe Fourier transform
**       direction = -1 ==>  inverse Fourier transform
**       length = number of elements
**
**************************************************************************/ 

// static char sccsid[] = "@(#)IM_FftFloat.cc 3.1 96/05/02 CEA 1994 @(#)";

#include "IM_Obj.h"
#include <math.h>


#define INT_POW(x,y,z) { int l,xx,yy; xx = (x) ; yy = (y);  for (l=0,(z)=1;l<yy;++ l,z *= xx); } 
#ifndef PI
#define PI 3.1415926536 
#endif

typedef struct  {
                        float   cos;
                        float   sin;
                } trig_table_struct;
typedef struct  {
                        int                     table_len;
                        int                     length;
                        int                     len_exp;
                        trig_table_struct       *table;
                } trig_table_type;


typedef struct  {
                        int     source;
                        int     destin;
                 } rev_map_struct;
typedef struct  {
                        int             map_len;
                        int             length;
                        rev_map_struct  *map;
                } rev_map_type;

typedef	 struct {
      		  float 	re;
       		  float 	im;
      		 }       complex_float;


/******************************22222222222*********************************/

static complex_float *b_s_row_to_vector(complex_float *dat, int index, int length)
{
	complex_float	transfer;
	register	int		loop;
	complex_float	*vector,*vrem;
	
	vrem = vector = (dat + index*length);
	length /= 2;
 	for ( loop = 0 ; loop < length ; ++loop, ++vector){
		transfer = *(vector);
		*(vector) = *(vector+ length);
		*(vector+length) = transfer;  
	}
	return (vrem);
}

/*************************2222222222222222222****************************/

static void b_s_row_from_vector(complex_float *vector, int length)
{
 	/* this routine assumes that the correct address of thr row vector */
	/* within the data array is correctly specified by the pointer vector*/
	/* ie vector has not been changed since it was set by b_s_row_to_vector*/

	complex_float	transfer;
	register	int		loop;

	length /= 2;
 	for ( loop = 0 ; loop < length ; ++loop, ++vector){
		transfer = *(vector);
		*(vector) = *(vector+ length);
		*(vector+length) = transfer;
	}
}

/******************************22222222222*********************************/

static	void b_s_col_to_vector(complex_float *dat, complex_float *vector, int index, int length)
{
	register	int	ofst;
	register	int	loop;
	register	int	len_b_2;
	
	len_b_2 = length/2;
	ofst = len_b_2*length;
	dat += index;
	for ( loop = 0 ; loop < len_b_2 ; ++loop, ++vector, dat += length){
		*(vector) = *(dat+ ofst);
		*(vector+len_b_2) = *dat;
	}
}

/******************************22222222222*********************************/

static	void b_s_col_from_vector(complex_float *dat, complex_float *vector, int index, int length)
{
	register	int	ofst;
	register	int	loop;
	register	int	len_b_2;

	len_b_2 = length/2;
	ofst = len_b_2*length;
	dat += index;
	for ( loop = 0 ; loop < len_b_2 ; ++loop , ++vector , dat += length){
		*(dat + ofst) = *vector;
		*dat = *(vector+len_b_2);
	}
}

/******************************22222222222*********************************/

static void bitreverse(complex_float *vector, rev_map_type rev_map)
{
	register	rev_map_struct	*map;
	complex_float	transfer;
	register	int		loop;

	map = rev_map.map;	
	for ( loop = 0 ; loop < rev_map.map_len ; ++loop, ++map){
		 transfer = *(vector+map->destin);
		*(vector+map->destin) = *(vector+map->source);
		*(vector+map->source) = transfer;
 	}
}

/******************************22222222222*********************************/

static void ft_1d(complex_float *vector, trig_table_type trig_table)
{
	register	float	c,s;
	register int	butterfly = 1 ;	
	int 	rank;
	int	l2,end_of_loop2;
	register	complex_float 	*top_btrfly_ptr,*mid_btrfly_ptr;	
	register 	int 	trig_table_pos;
	register	float	re_transfer_buffer,im_transfer_buffer;
	trig_table_struct	*table;

	table = trig_table.table;
	for (rank = 1 ;  rank  <= trig_table.len_exp ; ++rank){
		top_btrfly_ptr = vector +0;
		mid_btrfly_ptr =vector+ butterfly;
		end_of_loop2 = trig_table.length >> rank;/* 2^(len_exp - rank) */
		for ( l2 = 1 ; l2 <= end_of_loop2 ; ++l2){
			for ( trig_table_pos = 0 ;
			      trig_table_pos < trig_table.table_len ;
			      trig_table_pos += end_of_loop2 ){
				c = (table+trig_table_pos)->cos;
				s = (table+trig_table_pos)->sin;
				re_transfer_buffer = (*mid_btrfly_ptr).re*c - 
						     (*mid_btrfly_ptr).im*s;
				im_transfer_buffer = (*mid_btrfly_ptr).im*c +
						     (*mid_btrfly_ptr).re*s;
				(*mid_btrfly_ptr).re = (*top_btrfly_ptr).re -
						    re_transfer_buffer;
				(*mid_btrfly_ptr).im = (*top_btrfly_ptr).im -
						    im_transfer_buffer;
				(*top_btrfly_ptr).re += re_transfer_buffer; 
				(*top_btrfly_ptr).im += im_transfer_buffer;
				top_btrfly_ptr += 1;
				mid_btrfly_ptr += 1;

			}
			top_btrfly_ptr += butterfly;
			mid_btrfly_ptr += butterfly;
		}
		butterfly *= 2;
	}
}
/******************************22222222222*********************************/

static void set_trig_table(trig_table_type trig_table, int direction)
{
        register	int	l;
	double 	arg;
	trig_table_struct	*table;

	table = trig_table.table;
	for ( l = 0 ; (l < trig_table.table_len) ; ++l, ++table){
		arg =direction*PI*l/trig_table.table_len;
		table->cos = (float) cos(arg);
		table->sin = (float) sin(arg);
	}
}

/******************************22222222222*********************************/

static void set_rev_map(rev_map_type rev_map)
{
    register int source,destin,s_mask,d_mask ;
    rev_map_struct *map;
    int adrs = 0;
	
    map = rev_map.map;
    for ( source = 0 ; source < rev_map.length -1 ; ++source )
    {
       destin = 0;
       d_mask = (int) rev_map.length/2;
       for ( s_mask = 1 ; (s_mask < rev_map.length) ; )
       {
           if ((s_mask & source) > 0)
           {
	       destin = destin | d_mask;
           }
           d_mask = d_mask >> 1;
           s_mask = s_mask << 1;
       }
       if ( destin > source)
       {
           (map+adrs)->source = source;
           (map+adrs)->destin = destin;
           ++adrs;
       }
   } 	  	
}

/******************************22222222222*********************************/


static void norm2d (complex_float *ary, int length)
{
    int l1,l2;

    for( l1 = 0 ; l1 < length ; ++l1)
    {
       for( l2 = 0 ; l2 < length ; ++l2, ++ary)
       {
          ary->re /= (float)(length * length);
          ary->im /= (float)(length * length);
       }
    }
}


/******************************22222222222*********************************/

void ft_cf_any_power_of_2(complex_float *dat, int direction, int length, bool normalize)

/* Takes the array pointed to by dat which is assumed to be of side-length
length and performs a fourier transform of its contents

External function calls
	None

*/
{
 	complex_float	*vector;
        complex_float	*vect_space_holder;
	int		index;
	rev_map_type	rev_map;
	trig_table_type	trig_table;
	int		len_exp;
	int		temp;

/* First we must check length to make sure it is a power of 2 */
	len_exp = (int)(0.3+log((double)(length))/(log(2.0)));
	INT_POW(2,len_exp,temp);
	if ( length != temp ){
		printf ("error: is not a power of two\n");
		exit (-1);
	}
	else{
/* Now we must reserve space for the various working arrays */
		INT_POW(2,len_exp-1,rev_map.map_len);
		INT_POW(2,(len_exp-1)/2,temp);
		rev_map.map_len -= temp;
		rev_map.length = length;
		rev_map.map = new rev_map_struct [rev_map.map_len];
		INT_POW(2,len_exp-1,trig_table.table_len);
		trig_table.length = length;
		trig_table.len_exp = len_exp;
		trig_table.table = new trig_table_struct [trig_table.table_len];
		vect_space_holder = new complex_float [length];
		set_trig_table(trig_table,direction);
		set_rev_map(rev_map);

/* row vectors actually occur in the data array so all we need to do*/
/* is referance their start position		*/
		for (index = 0 ; index < length ; ++index){
			vector = b_s_row_to_vector(dat,index,length); 
			bitreverse(vector,rev_map);
			ft_1d(vector,trig_table);   
			b_s_row_from_vector(vector,length); 
		}

/* column vectors are dispersed in the data array and thus have to */
/* collected together the variable vect_space_holder reserves enough*/
/* space in thi routine to hold the colledted vector		*/

		vector = &(vect_space_holder[0]);
		 for (index = 0 ; index < length ; ++index){
			b_s_col_to_vector(dat,vector,index,length);
			bitreverse(vector,rev_map);
			ft_1d(vector,trig_table);  
			b_s_col_from_vector(dat,vector,index,length);
		}

/*Now free-up space used by the working arrays */
		delete [] vect_space_holder;
		delete [] trig_table.table;
		delete [] rev_map.map;
	}
	if (direction == -1) norm2d (dat,length);
	if(normalize)
	{
		float val;
		if(direction==1) val = 1./float(length);
		else val = float(length);
		for(int i=0;i<length*length;i++)
		{
			dat[i].re*=val;
			dat[i].im*=val;
		}
	}
}

/*************************11111111111111118*******************************/

void  fft2d (Icomplex_f & Im_InOut, int Dir, bool normalize)
{
    complex_float *FFT_Dat;
    int Nl = Im_InOut.nl();

    FFT_Dat = (complex_float *)  &(Im_InOut(0,0));
    if (Dir != -1) ft_cf_any_power_of_2(FFT_Dat, 1, Nl, normalize);
    else ft_cf_any_power_of_2(FFT_Dat, -1, Nl, normalize);
}

/*************************11111111111111118*******************************/

void  fft2d (const Icomplex_f &Im_in, const Icomplex_f &Im_out, int Dir, bool normalize)
{
    complex_float *FFT_Dat;
    int Nl = Im_in.nl();
    int Nc = Im_in.nc();

    FFT_Dat = (complex_float *)  &(Im_out(0,0));

    if (Dir != -1)
	{
	    for (int i = 0; i < Nl; i++)
	    for (int j = 0; j < Nc; j++)
	    {
	        FFT_Dat[i*Nc+j].re = Im_in(i,j).real();
	        FFT_Dat[i*Nc+j].im = Im_in(i,j).imag();
	    }
		ft_cf_any_power_of_2(FFT_Dat, 1, Nl, normalize);
	}
    else
	{
		ft_cf_any_power_of_2(FFT_Dat, -1, Nl, normalize);
	    for (int i = 0; i < Nl; i++)
	    for (int j = 0; j < Nc; j++)
	        Im_in(i,j) = complex_f(FFT_Dat[i*Nc+j].re , FFT_Dat[i*Nc+j].im ) ;
	}
}

/******************************22222222222*********************************/

void fft2d (const Iint &Im_in, const Icomplex_f &Im_out, int Dir, bool normalize)
{
    complex_float *FFT_Dat;
    int Nl = Im_in.nl();
    int Nc = Im_in.nc();
    int i,j;

    FFT_Dat = (complex_float *)  &(Im_out(0,0));

    if (Dir != -1)
	{
    	for ( i = 0; i < Nl; i++)
    	for ( j = 0; j < Nc; j++)
    	{
        	FFT_Dat[i*Nc+j].re = Im_in(i,j);
        	FFT_Dat[i*Nc+j].im = 0.;
    	}
		ft_cf_any_power_of_2(FFT_Dat, 1, Nl, normalize);
	}
    else
	{
		ft_cf_any_power_of_2(FFT_Dat, -1, Nl, normalize);
    	for ( i = 0; i < Nl; i++)
    	for ( j = 0; j < Nc; j++)
        	Im_in(i,j) = FFT_Dat[i*Nc+j].re ;
	}
}

/******************************22222222222*********************************/

void fft2d (const Ifloat &Im_in, const Icomplex_f &Im_out, int Dir, bool normalize)
{
    complex_float *FFT_Dat;
    int Nl = Im_in.nl();
    int Nc = Im_in.nc();
    int i,j;

    FFT_Dat = (complex_float *)  &(Im_out(0,0));

    if (Dir != -1)
	{
    	for ( i = 0; i < Nl; i++)
    	for ( j = 0; j < Nc; j++)
    	{
        	FFT_Dat[i*Nc+j].re = Im_in(i,j);
        	FFT_Dat[i*Nc+j].im = 0.;
    	}
		ft_cf_any_power_of_2(FFT_Dat, 1, Nl, normalize);
	}
    else
	{
		ft_cf_any_power_of_2(FFT_Dat, -1, Nl, normalize);
    	for ( i = 0; i < Nl; i++)
    	for ( j = 0; j < Nc; j++)
        	Im_in(i,j) = FFT_Dat[i*Nc+j].re ;
	}
}

/******************************22222222222*********************************/

void fft2d_conv(const Ifloat &Im_in1, const Ifloat &Im_in2, Ifloat &Im_out)
{
    int Nl = Im_in1.nl();
    int Nc = Im_in1.nc();
    Icomplex_f  Im_in1_cf(Nl, Nc, "Im_in1_cf");
    Icomplex_f  Im_in2_cf(Nl, Nc, "Im_in2_cf");

    fft2d (Im_in1, Im_in1_cf, 1);
    fft2d (Im_in2, Im_in2_cf, 1);
    Im_in2_cf *= Im_in1_cf;
    fft2d (Im_in2_cf, Im_in1_cf, -1);

    for (int i = 0; i < Nl; i++)
    for (int j = 0; j < Nc; j++)
        Im_out(i,j) = Im_in1_cf(i,j).real();
}

/******************************22222222222*********************************/

static void norm1d (complex_float *ary,int length)
{
        int i;

        for (i = 0 ;i < length ; ++i, ++ary)
        {
            ary->re /= (float)(length);
            ary->im /= (float)(length);
        }
}

/*************************11111111111111118*******************************/

static void ft_cf_1d (complex_float *dat, int direction, int length)
/* Takes the array pointed to by dat which is assumed to be of side-length
length and performs a fourier transform of its contents

External function calls
	None

Return codes
0 :	No problems encountered.
1 :	Length was not a power of 2 FATAL_ERROR

*/

{
	int		error = 0;
 	complex_float	*vector;
        complex_float	*vect_space_holder;
	int		index;
	rev_map_type	rev_map;
	trig_table_type	trig_table;
	int		len_exp;
	int		temp;

/* First we must check length to make sure it is a power of 2 */
	len_exp = (int)(0.3+log((double)(length))/(log(2.0)));
	INT_POW(2,len_exp,temp);
	if ( length != temp ){
		error = 1;
		fprintf(stderr,"ft_cf_any_power_of_2... array size= %d,not a power of 2\n",length);
	}
	else{
/* Now we must reserve space for the various working arrays */
		INT_POW(2,len_exp-1,rev_map.map_len);
		INT_POW(2,(len_exp-1)/2,temp);
		rev_map.map_len -= temp;
		rev_map.length = length;
		rev_map.map = new rev_map_struct [rev_map.map_len]; 
		INT_POW(2,len_exp-1,trig_table.table_len);
		trig_table.length = length;
		trig_table.len_exp = len_exp;
		trig_table.table = new trig_table_struct[trig_table.table_len];
		vect_space_holder = new complex_float[length];
		set_trig_table(trig_table,direction);
		set_rev_map(rev_map);

/* row vectors actually occur in the data array so all we need to do*/
/* is referance their start position		*/
		index = 0;
		vector = b_s_row_to_vector(dat,index,length); 
		bitreverse(vector,rev_map);
		ft_1d(vector,trig_table);   
		b_s_row_from_vector(vector,length);  

/*Now free-up space used by the working arrays */
		delete [] vect_space_holder;
		delete [] trig_table.table;
		delete [] rev_map.map;
	}
	if (direction == -1) norm1d (dat,length);
}


/******************************22222222222*********************************/

void fft1d (const Ifloat &Signal, const Icomplex_f & Signal_cf, int Dir)
{
    int i;
    int N = Signal.nc();
    complex_float *FFT_Dat;

 /*   FFT_Dat = (complex_float *)  &(Signal_cf(0,0));*/

    FFT_Dat = new complex_float [N];
    
	if(Dir==1)
	{
		for (i=0; i<N; i++)
    	{
    	    FFT_Dat[i].re = Signal(i);
    	    FFT_Dat[i].im = 0.;
    	}
	}
	else
	{
		for (i=0; i<N; i++)
    	{
    	    FFT_Dat[i].re = Signal_cf(i).real();
    	    FFT_Dat[i].im = Signal_cf(i).imag();
    	}
	}

   	ft_cf_1d (FFT_Dat, Dir, N);
	
	if(Dir==1)
	    for (i=0; i<N; i++)
            Signal_cf(i) = complex_f (FFT_Dat[i].re, FFT_Dat[i].im);
	else
		for (i=0; i<N; i++)
            Signal(i) = FFT_Dat[i].re;
		
    delete [] FFT_Dat;
}

/***************************************************************************/

void fft1d (const Icomplex_f &Signal, const Icomplex_f & Signal_cf, int Dir)
{
    int i;
    int N = Signal.nc();
    complex_float *FFT_Dat;

    FFT_Dat = new complex_float [N];
    if(Dir==1)
	{
		for (i=0; i<N; i++)
    	{
        	FFT_Dat[i].re = Signal(i).real();
        	FFT_Dat[i].im = Signal(i).imag();
    	}
	}
	else
	{
		for (i=0; i<N; i++)
    	{
        	FFT_Dat[i].re = Signal_cf(i).real();
        	FFT_Dat[i].im = Signal_cf(i).imag();
    	}
	}
    ft_cf_1d (FFT_Dat, Dir, N);
    if(Dir==1)
	    for (i=0; i<N; i++)
            Signal_cf(i) = complex_f (FFT_Dat[i].re, FFT_Dat[i].im);
	else
	    for (i=0; i<N; i++)
            Signal(i) = complex_f (FFT_Dat[i].re, FFT_Dat[i].im);
    delete [] FFT_Dat;
}

/***************************************************************************/
