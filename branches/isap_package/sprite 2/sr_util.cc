#include "sr_util.h"

int compute_centroid(Ifloat *img,double sig,double *cent,int niter)
{
  int Nx = (*img).nx();
  int Ny = (*img).ny();
  cent[0]=0;
  cent[1]=0;

  // Loop to compute the centroid 
  int i,j,k;
  Ifloat w(Nx,Ny);
  Ifloat wimg(Nx,Ny);
  for (i=0;i<Nx;i++) 
    for (j=0;j<Ny;j++) w(i,j)=1;

  for (k=0; k < niter; k++)
    {
      double tot=0; 
      for (i=0;i<Nx;i++) 
	for (j=0;j<Ny;j++) 
	  {
	    if (k>0)
	      {
		w(i,j)=exp(-(pow(i+1-cent[0],2)+pow(j+1-cent[1],2))/(2*pow(sig,2)));
	      }
	    wimg(i,j)=(*img)(i,j)*w(i,j);
	   
	    tot+=wimg(i,j);
	  }
      cent[0]=0;
      cent[1]=0;
      for (i=0;i<Nx;i++)
	{
	  for (j=0;j<Ny;j++)
	    {
	      cent[0]+=wimg(i,j)*(i+1);
	      cent[1]+=wimg(i,j)*(j+1);
	    }
	} 
      cent[0]/=tot;
      cent[1]/=tot;
      
    }
  
  // Memory free
  w.free();
  wimg.free();
  if (isnan(cent[0]) || isinf(cent[0]) || isnan(cent[1]) || isinf(cent[1]))
    {
      return 0;
    }
  else return 1; 
}

int compute_centroid_arr(fltarray *data,double sig,double **cent,int niter)
{
  int Nz = (*data).nz(), Nx = (*data).nx(), Ny = (*data).ny();
  int k,i,j;
  Ifloat temp_img(Nx,Ny);
  double *cent_temp = (double*)malloc(2*sizeof(double));
  int flag = 1;
  for (k=0;k<Nz;k++)
    {
      for (i=0;i<Nx;i++)
	for(j=0;j<Ny;j++) temp_img(i,j)=(*data)(i,j,k);
      
      int flag_temp = compute_centroid(&temp_img,sig,cent_temp,niter);
      if (flag_temp==0) flag=0;
      cent[0][k]=cent_temp[0];
      cent[1][k]=cent_temp[1];
    }
  // Memory free 
  temp_img.free();
  free(cent_temp);
  return flag;
}

void thresholding(Ifloat *data_in,Ifloat *data_out,Ifloat *thresh,bool thresh_type)
{
  // If thresh_type = 0 => hard thresholding ; otherwise => soft thresholding ; default is hard
  int Nx = (*data_in).nx(),Ny = (*data_in).ny();
  int i,j;
  for (i=0;i<Nx;i++)
    for(j=0;j<Ny;j++)
      {
	if (abs((*data_in)(i,j)) < (*thresh)(i,j)) (*data_out)(i,j) = 0;
	else (*data_out)(i,j) = (*data_in)(i,j) - (*thresh)(i,j)*thresh_type*sign_num((*data_in)(i,j));
      }
}

void thresholding(fltarray *data_in,fltarray *data_out,fltarray* thresh,bool thresh_type)
{
  // If thresh_type = 0 => hard thresholding ; otherwise => soft thresholding ; default is hard
  int Nx = (*data_in).nx(),Ny = (*data_in).ny(),Nz=(*data_in).nz();
  int i,j,k;
  for (i=0;i<Nx;i++)
    for(j=0;j<Ny;j++)
      for(k=0;k<Nz;k++)
	{
	 
	  if (abs((*data_in)(i,j,k)) < (*thresh)(i,j,k)) (*data_out)(i,j,k) = 0;
	  else (*data_out)(i,j,k) = (*data_in)(i,j,k) - (*thresh)(i,j,k)*thresh_type*sign_num((*data_in)(i,j,k));
	}
}

void circ_thresh(Ifloat *data_in,Ifloat *data_out,double r,double *cent)
{
  int Nx=data_in->nx(),Ny=data_in->ny(),i,j;
  for (i=0;i<Nx;i++)
    for(j=0;j<Ny;j++)
      {
	(*data_out)(i,j)=(*data_in)(i,j);
	if ((pow(i-cent[0],2)+pow(j-cent[1],2)) > pow(r,2)) (*data_out)(i,j)=0;
      }
}

void circ_thresh(fltarray *data_in,fltarray *data_out,double r,double **cent)
{
  int Nx=(*data_in).nx(),Ny=(*data_in).ny(),Nz=(*data_in).nz(),i,j,k;
  for (k=0;k<Nz;k++)
    for (i=0;i<Nx;i++)
      for(j=0;j<Ny;j++)
	{
	  (*data_out)(i,j,k)=(*data_in)(i,j,k);
	  if ((pow(i-cent[0][k],2)+pow(j-cent[1][k],2)) > pow(r,2)) (*data_out)(i,j,k)=0;
	}
}

void flux_est(Ifloat *data, double r,double *cent,double*flux)
{
  int Nx=data->nx(),Ny=data->ny();
  Ifloat data_temp(Nx,Ny);
  circ_thresh(data,&data_temp,r,cent);
  *flux = data_temp.total();
  data_temp.free();
}

void flux_est(fltarray *data, double r,double **cent,double*flux,int ref_im_ind)
{
  int Nx=data->nx(),Ny=data->ny(),Nz=data->nz(),i,j,k;
  fltarray data_temp(Nx,Ny,Nz);
  Ifloat im_temp(Nx,Ny);
  circ_thresh(data,&data_temp,r,cent);
  for (i=0;i<Nx;i++)
    for(j=0;j<Ny;j++) im_temp(i,j)=data_temp(i,j,ref_im_ind);
  flux[ref_im_ind] = im_temp.total();
  for (k=0;k<Nz;k++) 
    {
      for (i=0;i<Nx;i++)
	for(j=0;j<Ny;j++) im_temp(i,j)=data_temp(i,j,k);
      if (k!=ref_im_ind)
	{
	  flux[k] = im_temp.total();
	  flux[k] = flux[k]/flux[ref_im_ind];
	}  
    }
  flux[ref_im_ind] = 1;
  im_temp.free();
  data_temp.free();
}

int wl_trans(Ifloat *Dat,mr_opt opt,MultiResol* MR_Data)
{
  type_transform transf = opt.transf;
  int nb_sc = opt.nb_sc;
  int nb_usc = opt.nb_usc;
  sb_type_norm Norm = opt.Norm;
  type_undec_filter U_Filter = opt.U_Filter;
  type_sb_filter SB_Filter = opt.SB_Filter;
  type_border Bord = opt.Bord;
  Bool Verbose = False;
  FilterAnaSynt *FAS = opt.FAS;
  MR_Data->alloc ((*Dat).nl(), (*Dat).nc(), nb_sc, transf, FAS, Norm, nb_usc,U_Filter);
  MR_Data->Border = Bord;
  MR_Data->Verbose = Verbose;
  MR_Data->transform ((*Dat));
  
  if (MR_Data != NULL) return 1;
  else return 0;
}
int sign_num (double a)
{
  if (a==0) return 0;
  else return a/abs(a);
}

mr_opt mr_opt_init(int Nx,type_transform transf,int nb_sc,int nb_usc,sb_type_norm Norm,type_undec_filter U_Filter,type_sb_filter SB_Filter,type_border Bord,Bool Verbose )
{
  FilterAnaSynt *FAS ;
  mr_opt opt ;
  opt.transf=transf ;
  if (nb_sc<0)  
    opt.nb_sc = floor(log(Nx)/log(2))-1 ;
  else 
    opt.nb_sc = nb_sc;
  if (nb_usc<0)  
    opt.nb_usc = opt.nb_sc ;
  else 
    opt.nb_usc = nb_usc;
  opt.Norm = Norm ;
  opt.U_Filter = U_Filter ;
  opt.SB_Filter = SB_Filter ; 
  opt.Bord = Bord ;
  opt.FAS = new FilterAnaSynt();
  if ((transf == TO_MALLAT) || (transf == TO_UNDECIMATED_MALLAT))
    {
      (opt.FAS)->Verbose = Verbose;
      (opt.FAS)->alloc(SB_Filter);
    }
  return opt;
}
void wl_thresholding(MultiResol*wav_coeff,fltarray *thresh,bool thresh_type)
{
  int Nx = (*wav_coeff).size_band_nc(0),Ny=(*wav_coeff).size_band_nl(0),Nz=(*wav_coeff).nbr_band()-1,i,j,k; // The coarse scale is not thresholded
  fltarray temp_coeff(Nx,Ny,Nz), coeff_thresh(Nx,Ny,Nz);
  for (i=0;i<Nx;i++)
    for (j=0;j<Ny;j++)
      for(k=0;k<Nz;k++) temp_coeff(i,j,k) = (*wav_coeff)(k,j,i);
  thresholding(&temp_coeff,&coeff_thresh,thresh,thresh_type);
  for (i=0;i<Nx;i++)
    for (j=0;j<Ny;j++)
      for(k=0;k<Nz;k++)  (*wav_coeff)(k,j,i)=coeff_thresh(i,j,k);
  temp_coeff.free();
  coeff_thresh.free();
}

int gnoise_est(Ifloat* data,Ifloat* noise_arr) 
{
  int Nx = data->nx(),Ny=data->ny(),i,j;
 
  double sig = get_noise (*data,METHOD_MAD);
  for (i=0;i<Nx;i++)
    for (j=0;j<Ny;j++) (*noise_arr)(i,j)=sig;
  if (noise_arr != NULL) return 1;
  else return 0;
}

int noise_est(fltarray* data,fltarray* noise_arr)
{
  int Nx = data->nx(),Ny=data->ny(),Nz=data->nz(),i,j,k; 
  Ifloat temp(Nx,Ny);
  Ifloat wl_temp(Nx,Ny);
  mr_opt den_opt = mr_opt_init(Nx,TO_UNDECIMATED_MALLAT,2,2,NORM_L2);
  for (k=0;k<Nz;k++)
    {
      for (i=0;i<Nx;i++)
	for (j=0;j<Ny;j++) temp(i,j) = (*data)(i,j,k);
      MultiResol MR_Data;
      int flag = wl_trans(&temp,den_opt,&MR_Data);
      for (i=0;i<Nx;i++)
	for (j=0;j<Ny;j++) wl_temp(i,j) = MR_Data(0,j,i);
      double sig = get_noise(wl_temp,METHOD_MAD);
      for (i=0;i<Nx;i++)
	for (j=0;j<Ny;j++) (*noise_arr)(i,j,k)=sig;
      
      MR_Data.free();
    }
  
  temp.free();
  wl_temp.free();
  if (noise_arr != NULL) return 1;
  else return 0;
}

int gnoise_est(fltarray* data,fltarray* noise_arr) 
{
  int Nx = data->nx(),Ny=data->ny(),Nz=data->nz(),i,j,k; 
  Ifloat temp(Nx,Ny);
  for (k=0;k<Nz;k++)
    {
      for (i=0;i<Nx;i++)
	for (j=0;j<Ny;j++) temp(i,j) = (*data)(i,j,k);
      double sig = get_noise(temp,METHOD_MAD);
      for (i=0;i<Nx;i++)
	for (j=0;j<Ny;j++) (*noise_arr)(i,j,k)=sig;
    }
  temp.free();
  if (noise_arr != NULL) return 1;
  else return 0;
}

int wl_gnoise_est(MultiResol* wav_coeff,fltarray* noise_arr) 
{
  int Nx = (*wav_coeff).size_band_nc(0),Ny=(*wav_coeff).size_band_nl(0),Nz=(*wav_coeff).nbr_band()-1,i,j,k; 
  Ifloat temp(Nx,Ny);
  for (k=0;k<Nz;k++)
    {
      for (i=0;i<Nx;i++)
	for (j=0;j<Ny;j++) temp(i,j) = (*wav_coeff)(k,j,i);
      double sig = get_noise (temp,METHOD_MAD);
      
      if(sig < 0)
	cout << "Noise estimation failed"<< endl;
      for (i=0;i<Nx;i++)
	for (j=0;j<Ny;j++) (*noise_arr)(i,j,k)=sig;
    }
  temp.free();
  return 1;
}

void wl_filter(Ifloat *img,Ifloat*img_filt,mr_opt opt,double nsig,bool thresh_type,fltarray* noise_map,fltarray* coeff_thresh,Ifloat *coarse_scale)
{
  MultiResol wav_coeff;

  int flag_trans = wl_trans(img,opt,&wav_coeff);
  int Nx = wav_coeff.size_band_nc(0),Ny=wav_coeff.size_band_nl(0),Nz=wav_coeff.nbr_band()-1,i,j,k; // The coarse scale is not thresholded
  fltarray noise_arr(Nx,Ny,Nz);
  if (noise_map!=NULL)
    {
      for (i=0;i<Nx;i++)
	for (j=0;j<Ny;j++)
	  for(k=0;k<Nz;k++)  noise_arr(i,j,k)=(*noise_map)(i,j,k)*nsig;
      }
  else
    {
      int flag_noise = wl_gnoise_est(&wav_coeff,&noise_arr);
      for (i=0;i<Nx;i++)
	for (j=0;j<Ny;j++)
	  for(k=0;k<Nz;k++)  noise_arr(i,j,k)=noise_arr(i,j,k)*nsig;
    }
  

  wl_thresholding(&wav_coeff,&noise_arr,thresh_type);
   
  if (coarse_scale!=NULL)
    
    for (i=0;i<Nx;i++)
      for (j=0;j<Ny;j++)
	{
	  
	  (*coarse_scale)(i,j) = wav_coeff(Nz,j,i);
	}
  wav_coeff.recons(*img_filt,opt.Bord);
  
  if (coeff_thresh!=NULL)
     for (i=0;i<Nx;i++)
       for (j=0;j<Ny;j++) 
	 for (k=0;k<Nz;k++) (*coeff_thresh)(i,j,k) = wav_coeff(k,j,i);
  
  wav_coeff.free();
  noise_arr.free();
  
}


void wl_arr_filter(fltarray *data,fltarray*data_filt,mr_opt opt,double nsig,bool thresh_type)
{
  int Nz = data->nz(), Nx = data->nx(), Ny = data->ny(),i,j,k;
  for (k=0;k<Nz;k++)
    {
      Ifloat temp_im(Nx,Ny);
      Ifloat temp_im_filt(Nx,Ny);
      for (i=0;i<Nx;i++)
	for(j=0;j<Ny;j++) temp_im(i,j)=(*data)(i,j,k);
     
      wl_filter(&temp_im,&temp_im_filt,opt,nsig,thresh_type);
       
      for (i=0;i<Nx;i++)
	for(j=0;j<Ny;j++) (*data_filt)(i,j,k)=temp_im_filt(i,j);
      temp_im.free();
      temp_im_filt.free();	  
    }  	 
}


void coadd_all_frames(fltarray Data, Ifloat &Ima,Bool GetMedian)
{
  int Nx, Ny, Nz;
  Nx = Data.nx();
  Ny = Data.ny();
  Nz = Data.nz();
  int x,y,z;
  double Val;
  if (GetMedian == False)
    {
    
      for (x=0; x < Nx; x++)
	for (y=0; y < Ny; y++)
	  {	
	    Val = 0.;
	    for (z=0; z < Nz; z++) Val += Data(x,y,z);
	    if (Nz > 0) Ima(x,y) = (float) (Val / (double) Nz);
	    else Ima(x,y) = 0;
	  }
    }
  else
    {
      
      float *V = new float [Nz];
      for (x=0; x < Nx; x++)
	for (y=0; y < Ny; y++)
	  {
	    if (Nz > 0) 
	      {
		for (z=0; z < Nz; z++)  V[z] = Data(x,y,z);
	   
		Ima(x,y) = get_median(V, Nz);
	      }
	    else Ima(x,y) = 0.;
	 
	  }

      delete [] V;

    }
   
   
}

/*********************************************************************/

void coadd_all_frames(fltarray Dat, Ifloat &Ima, float Zoom,
		      Bool GetMedian,int MaxDist, 
		      int Surface, type_interp_corr TypeInterp,
		      Bool OptNoOffset,Bool MeanSub,Bool Verbose)
{

  int PosX = -1;
  int PosY = -1;
  int ValRet, i, nz = Dat.nz();
  int Nx1 = (int)(Zoom * Dat.nx());
  int Ny1 = (int)(Zoom * Dat.ny());
  int k,l;
  fltarray DatOut(Nx1, Ny1, Dat.nz());   
  xcorr_def Xcor;
   Xcor.xc_np = nz;
   Xcor.xc_x = new double [nz];
   Xcor.xc_y = new double [nz];
   Xcor.xc_level = new double [nz];
   for (i=0; i < nz; i++) Xcor.xc_x[i] = Xcor.xc_y[i] = Xcor.xc_level[i] = 0.;
   Xcor.xc_dx = MaxDist;
   Xcor.xc_dy = MaxDist;
   if ((Surface < 1) || (Surface > Dat.nx()/2))
   {
      Xcor.xc_hx = MAX(0, Dat.nx()/2 -1 - MaxDist);
      Xcor.xc_hy = MAX(0, Dat.ny()/2 -1 - MaxDist);
   }
   else
   {
      Xcor.xc_hx = Surface;
      Xcor.xc_hy = Surface;
   }
   if (PosX < 0) PosX =  Dat.nx() / 2;
   if (PosY < 0) PosY =  Dat.ny() / 2;
   int D1 = MIN(PosX-MaxDist,PosY-MaxDist);
   int D2 = MIN(Dat.nx()-PosX-1-MaxDist,Dat.ny()-PosY-1-MaxDist);
   int DM = MIN(D1,D2);
   if (DM < 2)
   {
      cout << "Error: search position cannot be at the border of the image ... " << endl;
      exit(-1);
   } 
   if (Surface >= DM) 
   {
      Surface = DM-1;
      Xcor.xc_hx = Surface;
      Xcor.xc_hy = Surface;
      cout << "Warning: surface must be decreased: new value = " << Surface << endl;
   } 
       
   //Xcor.xc_hx = 10;
   //Xcor.xc_hy = 10;
   Xcor.xc_method = 0;
   char *interp_method= (char*)StringCorrInterp(TypeInterp);
   
   float **TabFrame = new float * [nz];
   float * TabMean = new float [nz];
   for (i = 0; i < nz; i++) TabFrame[i] = Dat.buffer() + i*Dat.nx()*Dat.ny();
   
   
   if (OptNoOffset == False) // otherwise we calculte the offset frame by correlation
   {  
      Ifloat ImaFrame;
      Xcor.xc_init = 0;
      
      // Test position and surface
      
      for (i = 0; i < nz; i++) 
      {
          Xcor.xc_x[i] = Xcor.xc_y[i] = 0.;
	  if (MeanSub == True)
	  {
	     ImaFrame.alloc(TabFrame[i], Dat.ny(), Dat.nx());
	     TabMean[i] = (float) average(ImaFrame);
 	     for (k=0; k < ImaFrame.nl(); k++)
	     for (l=0; l < ImaFrame.nc(); l++)  ImaFrame(k,l) -= TabMean[i];
	  }
      }
      float *Pattern = TabFrame[0];
      
      ValRet = cube_get_offset(TabFrame, Pattern, 
                               Dat.nx(), Dat.ny(), Dat.nz(), 
			       PosX, PosY, &Xcor);
      for (i = 0; i < nz; i++) 
      {
 	  ImaFrame.alloc(TabFrame[i], Dat.ny(), Dat.nx());
	  if (MeanSub == True)
	  {
 	     for (k=0; k < ImaFrame.nl(); k++)
	     for (l=0; l < ImaFrame.nc(); l++)  ImaFrame(k,l) += TabMean[i];
	  }
      }
      // cout <<"ValRet = " << ValRet << endl; 
      if (Verbose == True)
        for (i = 0; i < nz; i++)
           cout << "Offset Frame: " << i+1 << ", dx = " << Xcor.xc_x[i] << " dy = " << Xcor.xc_y[i] << endl;     

      
   }
   
  
   
   if ((OptNoOffset == False) || (Zoom != 1))
   {
       if (Verbose == True) cout << "Resample cube ... " << endl;
       
      
       int ValRet, i;
    
    
    float **TabFrameOut;
    
    TabFrameOut = new float * [nz];   
    for (i = 0; i < nz; i++) 
        TabFrameOut[i] = DatOut.buffer() + i*DatOut.nx()*DatOut.ny();
    
    ValRet = cube_resample(TabFrame, TabFrameOut, &Xcor, Dat.nx(), Dat.ny(), Dat.nz(),
                           DatOut.nx(), DatOut.ny(), interp_method); 
    
    coadd_all_frames(DatOut, Ima, GetMedian);
    free(TabFrameOut);
       
   }
else 
    {
      coadd_all_frames(Dat, Ima);
    }
  DatOut.free();
  free(Xcor.xc_x);
  free(Xcor.xc_y);
  free(Xcor.xc_level);
  free(TabFrame);
  free(TabMean);
   
}

void shift_est(fltarray *data, double sig_gfit, double *sig_vec,double **cent,double **shift,int ref_im_ind,double nsig)
{
  int Nx  = data->nx();
  int Ny  = data->ny();
  int Nz  = data->nz();
  fltarray thresh(Nx,Ny,Nz);
  fltarray data_filt(Nx,Ny,Nz);
  int i,j,k;
  for(k=0;k<Nz;k++)
    {
      Ifloat temp(Nx,Ny);
      for (i=0;i<Nx;i++)
	for (j=0;j<Ny;j++) temp(i,j)=(*data)(i,j,k);
      int psnr=floor(temp.max()/sig_vec[k]);
      if (psnr<nsig)
	{
	  for (i=0;i<Nx;i++)
	    for (j=0;j<Ny;j++) thresh(i,j,k)=(psnr-1)*sig_vec[k]; 
	}
      else 
	{
	  for (i=0;i<Nx;i++)
	    for (j=0;j<Ny;j++) thresh(i,j,k)=nsig*sig_vec[k]; 
	}
      temp.free();
    }
  thresholding(data,&data_filt,&thresh);
  int flag_cen = compute_centroid_arr(&data_filt,sig_gfit,cent);
  thresh.free();    
  for(k=0;k<Nz;k++)
    {
      shift[0][k] = cent[0][k]-cent[0][ref_im_ind];
      shift[1][k] = cent[1][k]-cent[1][ref_im_ind];
    }
  data_filt.free();
}

double sinc(double x)
{
  double out;
  if (x==0) out = 1;
  else out = sin(x)/x;
  return out;
}

void lanczos(double u,double*mask,int mask_rad)
{
  int i;
  for (i=0;i<2*mask_rad+1;i++)  mask[i] = sinc(Pi*(u-(i-mask_rad)))*sinc(Pi*(u-(i-mask_rad))/mask_rad);
	      
}

void lanczos(double u[],double**mask,int mask_rad)
{
  int i,j;
  for (i=0;i<2*mask_rad+1;i++) 
    for (j=0;j<2*mask_rad+1;j++) mask[i][j] = sinc(Pi*(u[0]-(i-mask_rad)))*sinc(Pi*(u[0]-(i-mask_rad))/mask_rad)* sinc(Pi*(u[1]-(j-mask_rad)))*sinc(Pi*(u[1]-(j-mask_rad))/mask_rad);
}

void lanczos(double u[],Ifloat &mask,int mask_rad)
{
  int i,j;
  for (i=0;i<2*mask_rad+1;i++) 
    for (j=0;j<2*mask_rad+1;j++) mask(i,j) = sinc(Pi*(u[0]-(i-mask_rad)))*sinc(Pi*(u[0]-(i-mask_rad))/mask_rad)* sinc(Pi*(u[1]-(j-mask_rad)))*sinc(Pi*(u[1]-(j-mask_rad))/mask_rad);
}

void decim(Ifloat *img,Ifloat *img_filt,Ifloat *img_dec,int D,Ifloat* mask,Bool filt_en,Bool fft_en) 
{
  int Nx = img->nx(),Ny = img->ny(),i,j;
  
  if ((filt_en==True) & (D>1))
    if (fft_en==True) psf_convol((*img), (*mask), (*img_filt),True,False);
  if (D==1) (*img_dec)=(*img);
  else
    {
      if (filt_en==True)
	{
	  for(i=0;i<floor(Nx/D);i++)
	    for(j=0;j<floor(Ny/D);j++) (*img_dec)(i,j)=(*img_filt)(D*i,D*j);
	}
      else
	{
	  for(i=0;i<floor(Nx/D);i++)
	    for(j=0;j<floor(Ny/D);j++) (*img_dec)(i,j)=(*img)(D*i,D*j);	  
	}
    }	
}

void rotate(Ifloat *input,Ifloat *output,int dir) // Same as in IDL
{
  dir = dir%8;
  int Nx = input->nx(), Ny = input->ny(),x,y,x1,y1;
  switch(dir){
  case 0 : 
    for (x=0;x<Nx;x++)
      for (y=0;y<Ny;y++)
	{
	  x1=x;
	  y1=y;
	  (*output)(x1,y1)=(*input)(x,y);
	}
    break;
  case 1 : 
    for (x=0;x<Nx;x++)
      for (y=0;y<Ny;y++)
	{
	  x1=Ny-y-1;
	  y1=x;
	  (*output)(x1,y1)=(*input)(x,y);
	}
    break;
  case 2 : 
    for (x=0;x<Nx;x++)
      for (y=0;y<Ny;y++)
	{
	  x1=Nx-x-1;
	  y1=Ny-y-1;
	  (*output)(x1,y1)=(*input)(x,y);
	}
    break;
  case 3 : 
    for (x=0;x<Nx;x++)
      for (y=0;y<Ny;y++)
	{
	  x1=y;
	  y1=Nx-x-1;
	  (*output)(x1,y1)=(*input)(x,y);
	}
    break;
  case 4 : 
    for (x=0;x<Nx;x++)
      for (y=0;y<Ny;y++)
	{
	  x1=y;
	  y1=x;
	  (*output)(x1,y1)=(*input)(x,y);
	}
    break;
  case 5 : 
    for (x=0;x<Nx;x++)
      for (y=0;y<Ny;y++)
	{
	  x1=Nx-x-1;
	  y1=y;
	  (*output)(x1,y1)=(*input)(x,y);
	}
    break;
  case 6 : 
    for (x=0;x<Nx;x++)
      for (y=0;y<Ny;y++)
	{
	  x1=Ny-y-1;
	  y1=Nx-x-1;
	  (*output)(x1,y1)=(*input)(x,y);
	}
    break;
  case 7 : 
    for (x=0;x<Nx;x++)
      for (y=0;y<Ny;y++)
	{
	  x1=x;
	  y1=Ny-y-1;
	  (*output)(x1,y1)=(*input)(x,y);
	}
    break;
  }
    
}

void transpose_decim(Ifloat*im_in,Ifloat*im_out,int D)
{
  int Nxlr = im_in->nx(),Nylr = im_in->ny(),Nxhr = im_out->nx(),Nyhr = im_out->ny(),i,j;
  im_in->reform (Nxlr*Nylr);
  im_out->reform (Nxhr*Nyhr);
  for (i=0;i<Nxlr*Nylr;i++)
    {
      int Li = i/Nxlr;
      int Ci = i - Li*Nxlr;
      (*im_out)(Li*Nxlr*pow((double) D,(double)2)+Ci*D)=(*im_in)(i);
    }
  im_in->reform (Nxlr,Nylr);
  im_out->reform(Nxhr,Nyhr);
}



void power_meth(void(*opname)(double*,double*),int size_input,int *nb_iter,double *spec_rad,int nb_iter_max,double tol)
{
  double L_old=0;
  double L=1;
  *nb_iter=0;
  double * input = (double*)malloc(size_input*sizeof(double));
  //double * output = (double*)malloc(size_output*sizeof(double));
  int i;
  
  for (i=0;i<size_input;i++)
    {
      input[i]=rand();
    }
  double l0 = norm2(input,size_input);
  scale (input,size_input,1.0/l0);
    while ((abs(L_old-L)/L <tol) & (*nb_iter <nb_iter_max))
      {
	scale(input,size_input,1.0/L);
	(*opname)(input,input);
	L_old = L;
	L = norm2(input,size_input);
	(*nb_iter)++;
      }
  *spec_rad = L*(1+tol); // Better have a larger value
  if (*nb_iter==nb_iter_max) cout << "Warning: max number of iteration reached in spectral radius estimation"<< endl;
  free(input);
}

double norm2 (double * x,int siz)
{
  int i;
  double out=0;
  for (i=0;i<siz;i++) out+=pow(x[i],2);
  out = sqrt(out);
  return out;
			    
}
		        
void scale (double * x,double siz,double a)
{
  int i;
  for (i=0;i<siz;i++) x[i]*=a;
}

void ineq_cons(Ifloat &A, Ifloat B,double b) // projects A onto  E = {X/X(i,j)>=b*B(i,j)}
{
  int i,j,Nx=A.nx(),Ny=B.ny();
  for (i=0;i<Nx;i++)
    for (j=0;j<Ny;j++) 
      if (A(i,j)<b*B(i,j)) A(i,j) = b*B(i,j);
}


void reverse(double*u,int length,double*ur)
{
  int i;
  for(i=0;i<length;i++)
    ur[i]=u[length-1-i];
}

void holes(double*h,int sc_nb,int length,double*hout)
{
  int i;
  for (i=0;i<length;i++)
    hout[i*(int)pow((double) 2,(double) sc_nb)]=h[i];
}

void convol1D(double*in,double*out,double*h,int length,int filt_length)
// Edges zero convention
{
  int i,k;
  for (i=0;i<length;i++)
    out[i]=0;
  for (k=0;k<min(i,filt_length);k++)
       out[i] = out[i]+in[k]*h[filt_length-1-k];
}


void sep_filt2d(Ifloat* im_in,Ifloat*im_out,double *h,double *g,int lh,int lg)
{
  int i,j,k,Nx,Ny;
  Nx = im_in->nx();
  Ny = im_in->ny();
  Ifloat im_temp(Nx,Ny);
  // Columns filtering
  for(i=0;i<Nx;i++)      
      for(j=0;j<Ny;j++)
	for (k=0;k<min(j,lg);k++)
	  im_temp(i,j) = im_temp(i,j) + (*im_in)(i,k)*g[lg-k-1];
  // Lines filtering
  for(j=0;j<Ny;i++)      
      for(i=0;i<Nx;i++)
	for (k=0;k<min(i,lh);k++)
	  (*im_out)(i,j) = (*im_out)(i,j) + im_temp(k,j)*h[lh-k-1];
  im_temp.free();
}

void transp_sep_filt2d(Ifloat* im_in,Ifloat*im_out,double *h,double *g,int lh,int lg)
{
  double *hr = new double[lh];
  double *gr = new double[lg];
  reverse(h,lh,hr);
  reverse(g,lg,gr);
  sep_filt2d(im_in,im_out,hr,gr,lh,lg);
}



/*void randomn(double *x,double sig,int length,double mean)

{
  int i;  
  random_device rd;
  default_random_engine generator(rd());
  normal_distribution<double> distribution(mean,sig);
  for(i=0;i<length;i++)
    {  
      x[i] = distribution(generator);     
    }
    }*/

void randomngsl(double *x,double sig,int length,double mean)
{
  int i;  
  gsl_rng *rng;
  const gsl_rng_type * T;
  T = gsl_rng_default;
  rng = gsl_rng_alloc (T);
  int s = renewSeed();
  gsl_rng_set(rng, s);
  
  for(i=0;i<length;i++)
    {  
      x[i] =gsl_ran_gaussian(rng,sig)+mean;     
    }
}

/*void randomn(Ifloat *mat, double sig,double mean)
{
  random_device rd;
  default_random_engine generator(rd());
  normal_distribution<double> distribution(mean,sig);
  int i,j,Nx,Ny;
  Nx = mat->nx();
  Ny = mat->ny();
  for (i=0;i<Nx;i++)
    for (j=0;j<Ny;j++)
      {
      	double x;
	(*mat)(i,j) = distribution(generator);
      }
      }*/

void randomngsl(Ifloat *mat, double sig,double mean)
{
   gsl_rng *rng;
  const gsl_rng_type * T;
  T = gsl_rng_default;
  rng = gsl_rng_alloc (T);
  int s = renewSeed();
  gsl_rng_set(rng, s);
  int i,j,Nx,Ny;
  Nx = mat->nx();
  Ny = mat->ny();
  for (i=0;i<Nx;i++)
    for (j=0;j<Ny;j++)
      {
      	double x;
	(*mat)(i,j) = gsl_ran_gaussian(rng,sig)+mean;
      }
}

/*void randomn(fltarray *mat,double *sig,double mean)
{
  random_device rd;
  default_random_engine generator(rd());
  int i,j,k,Nx,Ny,Nz;
  Nx = mat->nx();
  Ny = mat->ny();
  Nz = mat->nz();
  for (k=0;k<Nz;k++)
    {
    normal_distribution<double> distribution(mean,sig[k]);
    for (i=0;i<Nx;i++)
      for (j=0;j<Ny;j++)
	{
	  
	  (*mat)(i,j,k) = distribution(generator);
	  
	}
    }
    }*/

void randomngsl(fltarray *mat,double *sig,double mean)
{
   gsl_rng *rng;
  const gsl_rng_type * T;
  T = gsl_rng_default;
  rng = gsl_rng_alloc (T);
  int s = renewSeed();
  gsl_rng_set(rng, s);
  int i,j,k,Nx,Ny,Nz;
  Nx = mat->nx();
  Ny = mat->ny();
  Nz = mat->nz();
  for (k=0;k<Nz;k++)
    {
      for (i=0;i<Nx;i++)
	for (j=0;j<Ny;j++)
	  {
	    
	    (*mat)(i,j,k) = gsl_ran_gaussian(rng,sig[k])+mean;
	    
	    
	  }
    }
}

/*void randomn(fltarray *mat,double *sig,double*mean)
{
  random_device rd;
  default_random_engine generator(rd());
  int i,j,k,Nx,Ny,Nz;
  Nx = mat->nx();
  Ny = mat->ny();
  Nz = mat->nz();
  for (k=0;k<Nz;k++)
    {
    normal_distribution<double> distribution(mean[k],sig[k]);
    for (i=0;i<Nx;i++)
      for (j=0;j<Ny;j++)
	{
	  
	  (*mat)(i,j,k) = distribution(generator);
	}
    }
    }*/

void randomngsl(fltarray *mat,double *sig,double*mean)
{
   gsl_rng *rng;
  const gsl_rng_type * T;
  T = gsl_rng_default;
  rng = gsl_rng_alloc (T);
  int s = renewSeed();
  gsl_rng_set(rng, s);
  int i,j,k,Nx,Ny,Nz;
  Nx = mat->nx();
  Ny = mat->ny();
  Nz = mat->nz();
  for (k=0;k<Nz;k++)
    {
    for (i=0;i<Nx;i++)
      for (j=0;j<Ny;j++)
	{
	  
	  (*mat)(i,j,k) = gsl_ran_gaussian(rng,sig[k])+mean[k];
	}
    }
}

double var_iter(double x,int n,double *mean,double*M)
{
  double delta = x - *mean;
  *mean = *mean+delta/(n+1);
  *M = *M+delta*(x-*mean);
  double var = 0;
  if (n>0)
    var = *M/n;
  return var;
}

void mr_support_filt(Ifloat img,Ifloat &img_filt, mr_opt opt,double nsig,int nb_iter,double*mse,Bool Pos_coeff,double lambda,Bool coarse_cons,Bool pos_cons, Bool drop_coarse,Bool iso_cons, double sig)
{
  int Nx = img.nx(), Ny = img.ny(),i,j,k,n;
  if (!coarse_cons) 
    lambda =1;
  if (drop_coarse) 
    lambda = 0;
  Ifloat resi(Nx,Ny);
  Ifloat zeros_mat(Nx,Ny);
  Ifloat mask(Nx,Ny);
  double mu = 1.0;
  if(iso_cons== True)
    {
      double *cent = new double[2];
      int flag = compute_centroid(&img,sig,cent,10);
      for (i=0;i<Nx;i++) 
	for (j=0;j<Ny;j++) 
	  {
	    mask(i,j)=exp((pow(i+1-cent[0],2)+pow(j+1-cent[1],2))/(2*pow(sig,2)));
	  }
      free(cent);
      mu = mu/pow(mask.max(),2);
      cout << "mu = "<<mu<<endl;
      
    }
  for (i=0;i<Nx;i++)
      for (j=0;j<Ny;j++)
	img_filt(i,j) = 0; // The algorithm is initialized at 0
  // Support setting
  
  MultiResol mr_data;
  MultiResol mr_data_filt;
  MultiResol mr_res;
  int flag = wl_trans(&img,opt,&mr_data);
  flag = wl_trans(&img,opt,&mr_res);
  int nb_band = mr_data.nbr_band();
  fltarray supp(Nx,Ny,nb_band);
  for (i=0;i<Nx;i++)
    for (j=0;j<Ny;j++)
      for (k=0;k<nb_band;k++)
	supp(i,j,k) = 1;
  Ifloat suppk(Nx,Ny);
  //Ifloat band_smooth(Nx,Ny);
  //Ifloat band_temp(Nx,Ny);
  
  for (k=0;k<nb_band-1;k++)
    {
      Ifloat bandk = mr_data.extract_band (k);
      //smooth_bspline (bandk, band_smooth);
      //for (i=0;i<Nx;i++)
      //for (j=0;j<Ny;j++)
      //band_temp(i,j)=bandk(i,j)-band_smooth(i,j);
      double sig = get_noise(bandk,METHOD_MAD);
      //cout << "sig = "<<sig<<endl;
      Bool abs_en = (Bool)(!Pos_coeff);
      
      check_ineq(bandk,nsig*sig,abs_en,suppk);
      for (i=0;i<Nx;i++)
	for (j=0;j<Ny;j++)
	  supp(i,j,k) = suppk(j,i);
      bandk.free();
    }
  
  if (drop_coarse==True or coarse_cons==True)
    {
      for (i=0;i<Nx;i++)
	for (j=0;j<Ny;j++)
	  supp(i,j,nb_band-1) =0;
    }
  
  fltarray supp_cons(Nx,Ny,nb_band);
  for (k=0;k<nb_band-1;k++)
    for (i=0;i<Nx;i++)
      for (j=0;j<Ny;j++)
	supp_cons(i,j,k) = supp(i,j,k);   
  for (i=0;i<Nx;i++)
    for (j=0;j<Ny;j++)
      supp_cons(i,j,nb_band-1) =lambda;

  for (k=0;k<nb_band;k++)
    for (i=0;i<Nx;i++)
      for (j=0;j<Ny;j++)
	mr_data(k,j,i) = mr_data(k,j,i)*supp(i,j,k);  
  
  for (n=0;n<nb_iter;n++)
    {
      if (mse!=NULL)
	 mse[n] = 0;
      wl_trans(&img_filt,opt,&mr_data_filt);
      for (k=0;k<nb_band;k++)
	for (i=0;i<Nx;i++)
	  for (j=0;j<Ny;j++)
	    {
	      mr_res(k,j,i) = mr_data(k,j,i) - mr_data_filt(k,j,i)*supp_cons(i,j,k);
	      if (mse!=NULL)
		mse[n] = mse[n]+ pow(mr_res(k,j,i),2);
	    }
      mr_res.recons(resi,opt.Bord);
      for (i=0;i<Nx;i++)
	for (j=0;j<Ny;j++)
	  img_filt(i,j) = img_filt(i,j)+mu*(resi(i,j)-iso_cons*img_filt(i,j)*pow(mask(i,j),2));
      if (pos_cons==True)
	ineq_cons(img_filt, zeros_mat); 
    }
  resi.free();
  zeros_mat.free();
  supp.free();
  suppk.free();  
  supp_cons.free();
  mr_data.free();
  mr_data_filt.free();
  mr_res.free();
  mask.free();
  //band_smooth.free();
  //band_temp.free();
}

void check_ineq(Ifloat img,double thresh,Bool abs_en,Ifloat &flag)
{
   int Nx = img.nx(), Ny = img.ny(),i,j;
   for (i=0;i<Nx;i++)
    for (j=0;j<Ny;j++)
      {
	flag(i,j)=0;
	if(abs_en==True)
	  {
	    if (abs(img(i,j))>thresh)
	      flag(i,j)=1;
	  }
	else
	  if (img(i,j)>thresh)
	    flag(i,j)=1;
      }
}

void decim_conv_mat(Ifloat mask,int im_siz[],int decim_fact,Ifloat &output,double flux,double sig) 
{
  int i,j,l,m,Nx=mask.nx(),Ny = mask.ny();
  int mask_rad = (Nx-1)/2; // The kernel mask is supposed to have odd sides lengths
  Ifloat mask_extend(im_siz[0],im_siz[1]);
  //Ifloat mask_decim(floor(im_siz[0]/decim_fact),floor(im_siz[1]/decim_fact)); 
  Ifloat mask_rot(Nx,Ny);
  rotate(&mask,&mask_rot,2);
  int shift[2];
  int D = decim_fact;
  for (i=0;i<im_siz[0];i++)
    for (j=0;j<im_siz[1];j++)
    {
      shift[0] = i;
      shift[1] = j;
      //cout << "Debugg1" << endl;
      convmask_shift_extend(mask_rot,im_siz,shift,mask_extend);
      //cout << "Debugg2" << endl;
      for(l=0;l<floor(im_siz[0]/decim_fact);l++)
	    for(m=0;m<floor(im_siz[1]/decim_fact);m++) 
	      if ((D*l<shift[0]+mask_rad+1) or (D*l>shift[0]-mask_rad-1) or (D*m<shift[1]+mask_rad+1) or (D*m>shift[1]-mask_rad-1))
		{
		  //mask_decim(l,m)=mask_extend(D*l,D*m);	 
		  output(j*im_siz[0]+i,m*floor(im_siz[0]/decim_fact)+l) = output(j*im_siz[0]+i,m*floor(im_siz[0]/decim_fact)+l)+mask_extend(D*l,D*m)*flux/sig; // ! output has to be correctly initiliazed 
		}
      
      }
  mask_extend.free();
  mask_rot.free();
  // mask_decim.free();
}

void convmask_shift_extend(Ifloat mask,int im_siz[],int shift[],Ifloat &output)
{
  
  int Nx = mask.nx(),Ny = mask.ny(),i,j;
  int mask_rad = (Nx-1)/2; // The kernel mask is supposed to have odd sides lengths
  for (i = 0;i<Nx;i++)
    for (j = 0;j<Ny;j++)
      {
	if((i+shift[0]-mask_rad>-1) or (i+shift[0]-mask_rad<im_siz[0]) or (j+shift[1]-mask_rad>-1) or (i+shift[1]-mask_rad<im_siz[1]))
	  output(i+shift[0]-mask_rad,j+shift[1]-mask_rad) = mask(Nx-i,Ny-j);	
      }
    
}

void noise_map_comp(double **shifts, double *sig, double *flux,int nb_im,int lancrad, int im_siz[],mr_opt opt, int decim_fact,fltarray &output)
// Computes a noise standard deviation map
{
  int k,i,j,l,m;
  MultiResol mr_data; 
  int nb_band = mr_data.computeNbBand(im_siz[0],im_siz[1],opt.nb_sc,opt.transf,opt.nb_usc);
  Ifloat decim_conv(im_siz[0]*im_siz[1],im_siz[0]*im_siz[1]/pow(decim_fact,2));
  Ifloat im_temp(im_siz[0],im_siz[1]);
  Ifloat mask(2*lancrad+1,2*lancrad+1);
  double utemp[2];
 
  for(k=0;k<nb_im;k++)
    {
      utemp[0] = shifts[0][k];
      utemp[1] = shifts[1][k];
      lanczos(utemp,mask,lancrad);
      decim_conv_mat(mask,im_siz,decim_fact,decim_conv,flux[k],sig[k]);  
    }
   
  for(l=0;l<im_siz[0]*im_siz[1]/pow(decim_fact,2);l++)
    {
      for(i=0;i<im_siz[0];i++)
	for(j=0;j<im_siz[1];j++)
	  {
	    im_temp(i,j) = decim_conv(j*im_siz[0]+i,l);
	  }
      int flag = wl_trans(&im_temp,opt,&mr_data);
      for (m=0;m<nb_band;m++)
	for(i=0;i<im_siz[0];i++)
	  for(j=0;j<im_siz[1];j++)
	    output(i,j,m) = sqrt(pow(output(i,j,m),2)+pow(mr_data(m,j,i),2));
      }
  decim_conv.free();
  mask.free();
  im_temp.free();
  mr_data.free();
}

int renewSeed()
{
  FILE *file = fopen("/dev/urandom","r");   
  int seed;
  fread(&seed, sizeof(int), 1, file);
  fclose(file);
  //printf("Random generator seed = %d\n", seed);
  return seed;
}

/****************************************************************************/
