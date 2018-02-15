#include "sprite.h"

void Sprite::mod_gen(int ref_im_ind,int denoising_type,float zoom,double nsig,bool thresh_type)
{
  int Nz = this->data->nz(), Nx = this->data->nx(), Ny = this->data->ny(),i,j,k;
  fltarray data_cp(Nx,Ny,Nz);
  for (k=0;k<Nz;k++)
    for (i=0;i<Nx;i++)
      for(j=0;j<Ny;j++) data_cp(i,j,k) = (*data)(i,j,(k+ref_im_ind)%Nz);
  fltarray data_filt(Nx,Ny,Nz);
  if (denoising_type==LR_DENOISING)
    wl_arr_filter(&data_cp,&data_filt,this->mod_opt,nsig,thresh_type);
  
  if (zoom ==-1) zoom = this->upfact;
  if (zoom != this->upfact) cout << "Warning: Model upsampling factor different from object intrinsec upsampling factor!" << endl;
  Ifloat model_filt(Nx*zoom,Ny*zoom);
  
  if (denoising_type==LR_DENOISING)
    coadd_all_frames(data_filt, model, zoom);
  else
    coadd_all_frames(data_cp, model, zoom);
 
  if (denoising_type==MR_SUPP_DENOISING)
    {
      
      mr_support_filt(model,model_filt, this->mod_opt);
      for (i=0;i<Nx*zoom;i++)
	for(j=0;j<Ny*zoom;j++) model(i,j)=model_filt(i,j);
      
    }
   
  
  if (pos_cons==True)
    {
      Ifloat zeros_mat(model.nx(),model.ny());
      ineq_cons(model, zeros_mat);
      zeros_mat.free();
    }
  data_filt.free();
  data_cp.free();
  model_filt.free();
  
}      

Sprite::Sprite(int upfact,fltarray *data,int upfact_0,mr_opt mod_opt,mr_opt sr_opt,int nb_iter,double r,double sig_gfit,int mod_denoising_type,double nsig,string output_dir,string output_file_name,string res_file_name,int noise_est_meth,int min_id,Bool eq_flux,Bool eq_noise,int win,double tol,bool Verbose,char * shift_filename,Bool Sparse_en,char * dist_weights_filename,double spectral_rad,Bool fft_en,int ref_im_ind,Bool pos_cons,Bool deconv_en,Ifloat*conv_mask,int lanc_rad)
{
  int i,j,k;
  this->upfact = upfact;
  this->data = data;
  Nx  = data->nx();
  Ny  = data->ny();
  Nz  = data->nz();
  this->mod_opt = mod_opt;
  this->ref_im_ind = ref_im_ind;
  this->sr_opt = sr_opt;
  this->conv_mask = conv_mask;
  this->deconv_en = deconv_en;
  this->lanc_rad = lanc_rad;
  this->sig_gfit = sig_gfit;
  this->spectral_rad = spectral_rad;
  this->nb_iter = nb_iter;
  this->output_dir = output_dir;
  this->output_file_name = output_file_name; 
  this->resi_file_name = res_file_name; 
  this->min_id = min_id;
  this->noise_est_meth = noise_est_meth;
  this->mod_denoising_type = mod_denoising_type;
  this->eq_flux = eq_flux;
  this->eq_noise = eq_noise;
  this->Verbose = Verbose;
  this->Sparse_en = Sparse_en;
  cv_win = win;
  cv_tol = tol;
  
  if (deconv_en==True & conv_mask==NULL)
    {
      cerr << "Deconvolution enabled and no deconv mask provided, quitting..."<< endl;
      exit(0);
    }
  this->fft_en = fft_en; 
  this->pos_cons = pos_cons;
  // -------------------------- Various allocations ------------------------------ //
  iter_diff_buff = new double [win];
  nb_iter_cv = 0;
  sr_rec.alloc(Nx*upfact,Ny*upfact);
  sr_rec_lr.alloc(Nx,Ny,Nz);
  resi.alloc(Nx,Ny,Nz);
  resi_out.alloc(Nx,Ny,Nz);
  gradient.alloc(Nx*upfact,Ny*upfact);
  coarse_scale.alloc(Nx*upfact,Ny*upfact);
  op_mat_output.alloc(Nx,Ny,Nz);
  model.alloc(upfact*Nx,upfact*Ny);
  cout << "Computing a first guess... " << endl;
  mod_gen(ref_im_ind,mod_denoising_type); 
  
    mse_steep = (double*)malloc(nb_iter*sizeof(double));
  this->nsig = (double*)malloc(nb_iter*sizeof(double));
  for (i=0;i<nb_iter;i++) (this->nsig)[i] = nsig;
  if(min_id==0)
    {
      MultiResol mr_test; 
      int nb_band = mr_test.computeNbBand(Nx*upfact, Ny*upfact,sr_opt.nb_sc,sr_opt.transf,sr_opt.nb_usc);
      mr_test.free();
      noise_est_trace.alloc(nb_iter,nb_band-1);
      wvl_noise_map.alloc(Nx*upfact,Ny*upfact,nb_band-1);
      
    }
   // ----------------------- Weights related to PSF spatial variations -------------- //
  dist_weights = (double*)malloc(Nz*sizeof(double));
  for(k=0;k<Nz;k++)
      dist_weights[k]=1;
  
    
  
  if(dist_weights_filename!=NULL)
    {
      
      FILE *pFile;
      float weight;
      pFile = fopen (dist_weights_filename,"r");
      
      for (k=0;k<Nz;k++)
	 {
	   if(pFile==NULL)
	     cout << "Warning: Distance weights file not found, weights will be equal to 1. " << endl;
	   else
	     {
	       fscanf(pFile,"%f\n",&weight);
	       dist_weights[k]=weight;
	       
	   }
	}
      fclose(pFile);
    }
  // -------------------------------------------------------------------------------- //
  // ------------------------ Noise map and thresh setting ----------------------- //
  data_nrm.alloc(Nx,Ny,Nz);
  noise_sig_map.alloc(Nx,Ny,Nz);
  fltarray thresh(Nx,Ny,Nz);
  double *sig_vect = (double*)malloc(Nz*sizeof(double));
  
  int flag_sig = noise_est(data,&noise_sig_map);
  
  sig_max = noise_sig_map.max();
  
  for(k=0;k<Nz;k++)
    { 
      
      sig_vect[k]=noise_sig_map(0,0,k);
      
      for (i=0;i<Nx;i++)
	for (j=0;j<Ny;j++) 
	  {
	    double w = sig_max;
	    if (eq_noise == False)
	      w = noise_sig_map(i,j,k);
	    data_nrm(i,j,k) = dist_weights[k]*(*data)(i,j,k)/w; 
	  }
    }
  // --------------------------------------------------------------------------------- //

  // ------------------------------- Shift omputing ---------------------------------- //
  shift = (double**)malloc(2*sizeof(double*));
  shift[0] = (double*)malloc(Nz*sizeof(double));
  shift[1] = (double*)malloc(Nz*sizeof(double));
  shift_cent = (double**)malloc(2*sizeof(double*));
  shift_cent[0] = (double*)malloc(Nz*sizeof(double));
  shift_cent[1] = (double*)malloc(Nz*sizeof(double));
  
  
  cout << "Subpixel shifts estimation... "<< endl;
  cent = (double**)malloc(2*sizeof(double*));
  cent[0] = (double*)malloc(Nz*sizeof(double));
  cent[1] = (double*)malloc(Nz*sizeof(double));
  
  shift_est(data, sig_gfit/(upfact*upfact_0),sig_vect,cent,shift_cent,ref_im_ind);
  for(k=0;k<Nz;k++)
    {
      shift[0][k] = upfact*shift_cent[0][k];
      shift[1][k] = upfact*shift_cent[1][k];
      if (this->Verbose) printf("Shift Est[%d] = %f, %f ; %f, %f\n",k,shift[0][k],shift[1][k],shift_cent[0][k],shift_cent[1][k]);
    }
  
  if(shift_filename!=NULL)
    {
      
      FILE *pFile;
      float shift_x,shift_y;
      pFile = fopen (shift_filename,"r");
      if(this->Verbose) printf("shift = %s\n",shift_filename);
      for (k=0;k<Nz;k++)
	 {
	   fscanf(pFile,"%f %f\n",&shift_x,&shift_y);
	   shift[0][k]=shift_x*upfact;
	   shift[1][k]=shift_y*upfact;
	   if(this->Verbose) printf("Read shift [%d]=%f %f\n",k,shift_x,shift_y);
	}
      fclose(pFile);
    }
 // --------------------------------------------------------------------------------- //
  
  // ----------------------- Registration kernels setting --------------------------- //
  lanczos_ker.alloc(2*lanc_rad+1,2*lanc_rad+1,Nz);
  double **temp_ker = (double **)malloc((2*lanc_rad+1)*sizeof(double*));
  for (k=0;k<2*lanc_rad+1;k++) temp_ker[k] = (double*)malloc((2*lanc_rad+1)*sizeof(double));
  double u[2];  
  for(k=0;k<Nz;k++)
    {
      u[0] = shift[0][k];
      u[1] = shift[1][k];
      lanczos(u,temp_ker,lanc_rad);
      for (i=0;i<2*lanc_rad+1;i++) 
	for (j=0;j<2*lanc_rad+1;j++) 
	  {
	  lanczos_ker(i,j,k) = temp_ker[i][j];
      
       }
    }
  // ----------------------------------------------------------------------------------//
  

  // ------------------------------ Flux computing ----------------------------------- //
  cout << "Flux estimation... "<< endl;
  flux = (double*)malloc(Nz*sizeof(double));
  if (eq_flux == False)
    flux_est(data, r,cent,flux,ref_im_ind);
  else
    for(k=0;k<Nz;k++)
      flux[k]=1;

  // --------------------------------------------------------------------------------- //

  free(sig_vect);
  for (k=0;k<2*lanc_rad+1;k++) free(temp_ker[k]);
  free(temp_ker);
 
}

void Sprite::sprite_op_mat(double *input,double *output)
{
  int i,j,k;
  if (input!=NULL)
    {
      
      for (i=0;i<Nx*upfact;i++)
	for (j=0;j<Ny*upfact;j++) sr_rec(i,j)=input[i+j*Nx*upfact];
       }
  
  Ifloat sr_rec_temp(Nx*upfact,Ny*upfact);
  Ifloat sr_rec_filt(Nx*upfact,Ny*upfact);
  Ifloat sr_dec(Nx,Ny);
  Ifloat lanc_ker_temp(2*lanc_rad+1,2*lanc_rad+1);
  for (k=0;k<Nz;k++)
    {
      if (deconv_en==True) 
	if(fft_en==True) psf_convol(sr_rec, *conv_mask, sr_rec_temp,True,False);
     
      for (i=0;i<2*lanc_rad+1;i++) 
	for (j=0;j<2*lanc_rad+1;j++) 
	  {
	    lanc_ker_temp(i,j) = lanczos_ker(i,j,k);	    
	  }
      if (deconv_en==True) decim(&sr_rec_temp,&sr_rec_filt,&sr_dec,upfact,&lanc_ker_temp,True,fft_en); 
      else decim(&sr_rec,&sr_rec_filt,&sr_dec,upfact,&lanc_ker_temp,True,fft_en); 
      
      for (i=0;i<Nx;i++)
	for (j=0;j<Ny;j++) 
	  {
	    double w = sig_max;
	    if (eq_noise == False)
	      w = noise_sig_map(i,j,k);
	    op_mat_output(i,j,k)=sr_dec(i,j)*flux[k]*dist_weights[k]/w;
	    if (output!=NULL) output[i+j*Nx+k*Nx*Ny]=op_mat_output(i,j,k);
	  }
      
    }
  sr_rec_temp.free();
  lanc_ker_temp.free();
  sr_rec_filt.free();
  sr_dec.free();
}

void Sprite::sprite_op_mat(Ifloat input,fltarray &output)
{
  int i,j,k;
  Ifloat sr_rec_temp(Nx*upfact,Ny*upfact);
  Ifloat sr_rec_filt(Nx*upfact,Ny*upfact);
  Ifloat sr_dec(Nx,Ny);
  Ifloat lanc_ker_temp(2*lanc_rad+1,2*lanc_rad+1);
  for (k=0;k<Nz;k++)
    {
      
      if (deconv_en==True) 
	if(fft_en==True) psf_convol(input, *conv_mask, sr_rec_temp,True,False);
     
      for (i=0;i<2*lanc_rad+1;i++) 
	for (j=0;j<2*lanc_rad+1;j++) 
	  {
	    lanc_ker_temp(i,j) = lanczos_ker(i,j,k);	    
	  }
      if (deconv_en==True) decim(&sr_rec_temp,&sr_rec_filt,&sr_dec,upfact,&lanc_ker_temp,True,fft_en); 
      else decim(&input,&sr_rec_filt,&sr_dec,upfact,&lanc_ker_temp,True,fft_en); 
      
      for (i=0;i<Nx;i++)
	for (j=0;j<Ny;j++) 
	  {
	    double w = sig_max;
	    if (eq_noise == False)
	      w = noise_sig_map(i,j,k);
	    output(i,j,k)=sr_dec(i,j)*flux[k]*dist_weights[k]/w;
	  }
      
    }
  
  sr_rec_temp.free();
  lanc_ker_temp.free();
  sr_rec_filt.free();
  sr_dec.free();
}


void Sprite::sprite_transp_op_mat(double *input,double *output)
{
  
  Ifloat lr_temp(Nx,Ny);
  Ifloat hr_temp(Nx*upfact,Ny*upfact);
  Ifloat hr_temp_filt(Nx*upfact,Ny*upfact);
  Ifloat hr_temp_filt_2(Nx*upfact,Ny*upfact);
  Ifloat lanc_ker_temp(2*lanc_rad+1,2*lanc_rad+1);
  Ifloat lanc_ker_temp_rot(2*lanc_rad+1,2*lanc_rad+1);
  int i,j,k;
  for (i=0;i<Nx*upfact;i++)
    for (j=0;j<Ny*upfact;j++) gradient(i,j)=0;
  for (k=0;k<Nz;k++)
    {
      for (i=0;i<Nx;i++)
	for (j=0;j<Ny;j++) 
	  {
	    double w = sig_max;
	    if (eq_noise == False)
	      w = noise_sig_map(i,j,k);
	    if (input != NULL) lr_temp(i,j)=input[i+j*Nx+k*Nx*Ny]/w;
	  else lr_temp(i,j) = resi(i,j,k)/w;
	  }
      transpose_decim(&lr_temp,&hr_temp,upfact);  
      for (i=0;i<2*lanc_rad+1;i++) 
	for (j=0;j<2*lanc_rad+1;j++) lanc_ker_temp(i,j) = lanczos_ker(i,j,k);
      rotate(&lanc_ker_temp,&lanc_ker_temp_rot,2);
      if(fft_en==True) psf_convol(hr_temp, lanc_ker_temp_rot, hr_temp_filt,True,False);
      for (i=0;i<Nx*upfact;i++)
	for (j=0;j<Ny*upfact;j++) 
	  gradient(i,j)=gradient(i,j)+flux[k]*dist_weights[k]*hr_temp_filt(i,j);
    }      
  
  if (deconv_en==True) 
    {
      int siz1=conv_mask->nx(),siz2=conv_mask->ny();
      Ifloat conv_mask_rot(siz1,siz2);
      rotate(conv_mask,&conv_mask_rot,2);
      if(fft_en==True) psf_convol(gradient, conv_mask_rot, hr_temp_filt_2,True,False);
      for (i=0;i<Nx*upfact;i++)
	for (j=0;j<Ny*upfact;j++) 
	  gradient(i,j)= hr_temp_filt_2(i,j);
      conv_mask_rot.free();
    }
  if (output != NULL)
    for (i=0;i<Nx*upfact;i++)
      for (j=0;j<Ny*upfact;j++) output[i+j*Nx*upfact] = gradient(i,j);
  
  lr_temp.free();
  hr_temp.free();
  hr_temp_filt.free();
  hr_temp_filt_2.free();
  lanc_ker_temp.free();
  lanc_ker_temp_rot.free();
  
}

void Sprite::sprite_transp_op_mat(fltarray input,Ifloat &output)
{
  Ifloat lr_temp(Nx,Ny);
  Ifloat hr_temp(Nx*upfact,Ny*upfact);
  Ifloat hr_temp_filt(Nx*upfact,Ny*upfact);
  Ifloat hr_temp_filt_2(Nx*upfact,Ny*upfact);
  Ifloat lanc_ker_temp(2*lanc_rad+1,2*lanc_rad+1);
  Ifloat lanc_ker_temp_rot(2*lanc_rad+1,2*lanc_rad+1);
  int i,j,k;
  for (i=0;i<Nx*upfact;i++)
    for (j=0;j<Ny*upfact;j++) output(i,j)=0;
  for (k=0;k<Nz;k++)
    {
      for (i=0;i<Nx;i++)
	for (j=0;j<Ny;j++)
	  {
	    double w = sig_max;
	    if (eq_noise == False)
	      w = noise_sig_map(i,j,k);
	    lr_temp(i,j) = input(i,j,k)/w;
	  }
      transpose_decim(&lr_temp,&hr_temp,upfact);  
      for (i=0;i<2*lanc_rad+1;i++) 
	for (j=0;j<2*lanc_rad+1;j++) lanc_ker_temp(i,j) = lanczos_ker(i,j,k);
      rotate(&lanc_ker_temp,&lanc_ker_temp_rot,2);
      if(fft_en==True) psf_convol(hr_temp, lanc_ker_temp_rot, hr_temp_filt,True,False);
      for (i=0;i<Nx*upfact;i++)
	for (j=0;j<Ny*upfact;j++) 
	  output(i,j)=output(i,j)+flux[k]*dist_weights[k]*hr_temp_filt(i,j);
    }      
  
  if (deconv_en==True) 
    {
      int siz1=conv_mask->nx(),siz2=conv_mask->ny();
      Ifloat conv_mask_rot(siz1,siz2);
      rotate(conv_mask,&conv_mask_rot,2);
      if(fft_en==True) psf_convol(output, conv_mask_rot, hr_temp_filt_2,True,False);
      for (i=0;i<Nx*upfact;i++)
	for (j=0;j<Ny*upfact;j++) 
	  output(i,j)= hr_temp_filt_2(i,j);
      conv_mask_rot.free();
    }
  
  lr_temp.free();
  hr_temp.free();
  hr_temp_filt.free();
  hr_temp_filt_2.free();
  lanc_ker_temp.free();
  lanc_ker_temp_rot.free();
}
void Sprite::sprite_global_op_mat(double *input,double *output)
{
  double *op_mat_output = (double*)malloc(Nx*Ny*Nz*sizeof(double*));
  
  sprite_op_mat(input,op_mat_output);
  int i;
  
  sprite_transp_op_mat(op_mat_output,output);
  free(op_mat_output);
}

void Sprite::spectral_rad_est(int *nb_iter,int nb_iter_max,double tol)
{
  double L_old=0;
  double L=1;
  *nb_iter=0;
  int size_input = Nx*Ny*pow((double) upfact,(double)2);
  double * input = (double*)malloc(size_input*sizeof(double));
  int i;
  
  for (i=0;i<size_input;i++)
    {
      input[i]=rand();
    }
  double l0 = norm2(input,size_input);
  scale (input,size_input,1.0/l0);
 
  while ((fabs(L_old-L)/L >tol) & (*nb_iter <nb_iter_max))
    {
      scale(input,size_input,1.0/L);
      
      sprite_global_op_mat(input,input);
      
      L_old = L;
      L = norm2(input,size_input);
      
      (*nb_iter)++;
    }
   
  spectral_rad = L*(1+tol); // Better have a larger value
  
  if (*nb_iter==nb_iter_max) cout << "Warning: max number of iteration reached in spectral radius estimation"<< endl;
  free(input);
}
void Sprite::min1(Ifloat &delta,Bool thresh_type,int nb_rw) 

{
  if (spectral_rad==0)
    {
      cerr << "No spectral radius provided, quitting!" << endl;
      exit(0);
    }
  // Memory allocation for weights
  MultiResol mr_test;
  int i,j,k,l,m;
  
  int nb_band = mr_test.computeNbBand(Nx*upfact, Ny*upfact,sr_opt.nb_sc,sr_opt.transf,sr_opt.nb_usc);
  mr_test.free();
  fltarray weights(Nx*upfact,Ny*upfact,nb_band);
  
  for (i=0;i<Nx*upfact;i++)
    for (j=0;j<Ny*upfact;j++) 
      for (k=0;k<nb_band;k++) weights(i,j,k) = 1;
  fltarray coeff_thresh(Nx*upfact,Ny*upfact,nb_band-1);
  // Variables init 
  double t=1;
  Ifloat delta_old(Nx*upfact,Ny*upfact);
  Ifloat delta_old_2(Nx*upfact,Ny*upfact);
  Ifloat grad_weighted(Nx*upfact,Ny*upfact);
  fltarray weighted_wvl_noise_map(Nx*upfact,Ny*upfact,nb_band-1);
  MultiResol MR_grad;
  int ptr_b,ptr_e;
  
  for (l=0;l<nb_rw;l++)
    {
      for (k=0;k<nb_band-1;k++)
	for (i=0;i<Nx*upfact;i++)
	  for (j=0;j<Ny*upfact;j++) 
	    weighted_wvl_noise_map(i,j,k)=weights(i,j,k)*wvl_noise_map(i,j,k);
      double mean_iter_diff = 0;
      nb_iter_cv = 0;
      for (m=0;m<nb_iter;m++)
	{
	  nb_iter_cv+=1;
	  double iter_diff=0;
	  // --------------- Gradient step ------------------ //
  	  for (i=0;i<Nx*upfact;i++)
	    for (j=0;j<Ny*upfact;j++) sr_rec(i,j) = delta(i,j)+model(i,j);
	 
	  sprite_op_mat();
	  
	  for (i=0;i<Nx;i++)
	    for (j=0;j<Ny;j++) 
	      for (k=0;k<Nz;k++) 
		{
		  resi(i,j,k) = data_nrm(i,j,k)-op_mat_output(i,j,k);
		  double w = sig_max;
		  if (eq_noise == False)
		    {
		      w = noise_sig_map(i,j,k);
		      resi_out(i,j,k) = resi(i,j,k)*w/dist_weights[k];
		      
		      sr_rec_lr(i,j,k) = (*data)(i,j,k)-resi_out(i,j,k);
		    }
		}
	  mse_steep[m]=resi.energy();
	  //cout << mse_steep[m] << endl;
	  sprite_transp_op_mat();
	  for (i=0;i<Nx*upfact;i++)
	    for (j=0;j<Ny*upfact;j++) 
	      {
		delta_old_2(i,j) = delta(i,j);
		grad_weighted(i,j) = gradient(i,j)/spectral_rad;
		delta(i,j) = delta(i,j)+grad_weighted(i,j);
	      }
	  // ------------------------------------------------ //
	  
	  // ------ Noise estimation using the gradient ----- //
	  if(noise_est_meth==0)
	    {
	      int flag_wl = wl_trans(&grad_weighted,sr_opt,&MR_grad);
	      int flag_sig = wl_gnoise_est(&MR_grad,&wvl_noise_map);
	      MR_grad.free();
	      for (k=0;k<nb_band-1;k++)
		for (i=0;i<Nx*upfact;i++)
		  for (j=0;j<Ny*upfact;j++) 
		    weighted_wvl_noise_map(i,j,k)=weights(i,j,k)*wvl_noise_map(i,j,k);
	    }
	  if (l==0)
	    for (k=0;k<nb_band;k++)
	      {
		/*if (k==0)
		  cout << wvl_noise_map(0,0,k) << endl;*/
		noise_est_trace(m,k)=wvl_noise_map(floor(Nx*upfact/2),floor(Ny*upfact/2),k);
	      }
	  
	    
	  // ------------------------------------------------ //

	  // ------------------ Sparse constraint ----------- //
	  if(Sparse_en==True)
	    wl_filter(&delta,&delta,sr_opt,nsig[m],thresh_type,&weighted_wvl_noise_map,&coeff_thresh,&coarse_scale);
	  
	  // ------------------------------------------------ //
	  
	  // ---------- Gradient step acceleration ---------- //
	  double t_temp = (1+sqrt(4*pow(t,2)+1))/2;
	  double lambda = 1 + (t-1)/t_temp;
	  t = t_temp;
	  for (i=0;i<Nx*upfact;i++)
	    for (j=0;j<Ny*upfact;j++) 
	      {
		double temp = delta(i,j);
		delta(i,j) = delta_old(i,j)+lambda*(delta(i,j)-delta_old(i,j));
		delta_old(i,j) = temp;
	      }
	  if (pos_cons==True)
	    {
	      ineq_cons(delta,model,-1);
	    }
	  // ---------------- Stopping criteria ----------------- //
      double delta_old_l1=0;
      for (i=0;i<Nx*upfact;i++)
	for (j=0;j<Ny*upfact;j++)
	  { 
	    iter_diff += fabs(delta(i,j)-delta_old_2(i,j)); 
	    delta_old_l1+=fabs(delta_old_2(i,j));
	  }
      ptr_e = m%cv_win;
      iter_diff_buff[ptr_e]=iter_diff;
      mean_iter_diff+=iter_diff/cv_win;
      if (m>=cv_win)
	{
	  ptr_b = (ptr_e+1)%cv_win;
	  mean_iter_diff-=iter_diff_buff[ptr_b]/cv_win;
	  if(mean_iter_diff<cv_tol*delta_old_l1)
	    break;
	}
	}
      
      // --------------------- Reweighting ------------------ //
      for (i=0;i<Nx*upfact;i++)
	for (j=0;j<Ny*upfact;j++) 
	  for (k=0;k<nb_band-1;k++) weights(i,j,k) = 1/(1+fabs(coeff_thresh(i,j,k))/(3*wvl_noise_map(i,j,k)));
      if(Verbose==True)
	if(nb_iter_cv==nb_iter)
	  cout << "Warning: max number of iterations reached in the minimizer"<< endl;
    }
  for (i=0;i<Nx*upfact;i++)
    for (j=0;j<Ny*upfact;j++) 
      {
	sr_rec(i,j) = delta(i,j)+model(i,j);
	
      }
  
  delta_old.free();
  grad_weighted.free();
  coeff_thresh.free();
  weights.free();
  weighted_wvl_noise_map.free();
}
void Sprite::noise_simulator(int nb_monte)
{
  MultiResol mr_test;
  int nb_band = mr_test.computeNbBand(Nx*upfact, Ny*upfact,sr_opt.nb_sc,sr_opt.transf,sr_opt.nb_usc);
  fltarray mean_est(Nx*upfact, Ny*upfact,nb_band-1);
  fltarray M(Nx*upfact, Ny*upfact,nb_band-1);
  double *sig = new double[Nz];
  int k,i,l,m,n;
  for (k=0;k<Nz;k++)
    {
      sig[k] = 1/spectral_rad;
    }
  fltarray noise_sim(Nx,Ny,Nz);
  Ifloat op_out(Nx*upfact,Ny*upfact);
  MultiResol MR_noise;
  for (i=0;i<nb_monte;i++)
    {
      randomngsl(&noise_sim,sig);
      
      
      sprite_transp_op_mat(noise_sim,op_out);
      int flag_wl = wl_trans(&op_out,sr_opt,&MR_noise);
      for (l=0;l<Nx*upfact;l++)
	for (m=0;m<Ny*upfact;m++) 
	  for (n=0;n<nb_band-1;n++) 
	    {
	      double mean_temp = mean_est(l,m,n);
	      double M_temp = M(l,m,n);
	      double var = var_iter(MR_noise(n,m,l),i,&mean_temp,&M_temp);
	      mean_est(l,m,n) = mean_temp;
	      M(l,m,n) = M_temp;
	      wvl_noise_map(l,m,n) = sqrt(var);
	     } 
      
   }

  noise_sim.free();
  mean_est.free();
  M.free();
  op_out.free();
  MR_noise.free();
  mr_test.free();
}

void Sprite::noise_computer(void)
{
  double *sig = new double[Nz];
  int i;
  for (i=0;i<Nz;i++)
    {
      double w = sig_max;
      if (eq_noise == False)
	w = noise_sig_map(0,0,i);
      sig[i] = w/spectral_rad;
    }
  int im_siz[2];
  im_siz[0] = Nx*upfact;
  im_siz[1] = Ny*upfact;
  noise_map_comp(shift,sig,flux,Nz,lanc_rad,im_siz,sr_opt,upfact,wvl_noise_map);
}


void Sprite::write_outputs(Bool writeRes,Bool writeParam)
{
  int i ;
  size_t size = output_dir.size() + 1024;
  char* sol_file = new char[size];
  char* sol_file_lr = new char[size];
  struct stat s;
  int err = stat(output_dir.c_str(), &s);
  MatOper MAT;
  fltarray TempArray,tTempArray;
  if(-1 == err) {
    if(ENOENT == errno) {
      int e;
      e = mkdir(output_dir.c_str(), S_IRWXU);
      if(e!=0) {
	printf("Cannot create directory\n");
	exit(1);
      }
    } else {
      printf("Cannot create directory: file existing with same name\n");
      exit(1);
    }
  }   
  sprintf(sol_file,"%s/%s",output_dir.c_str(),output_file_name.c_str());
  sprintf(sol_file_lr,"%slr_%s",output_dir.c_str(),output_file_name.c_str());
  TempArray.alloc(sr_rec.buffer(),sr_rec.ny(),sr_rec.nx());
//  io_write_ima_float(sol_file,sr_rec);
  MAT.transpose(TempArray,tTempArray);
  fits_write_fltarr(sol_file,tTempArray);
  TempArray.alloc(sr_rec_lr.buffer(),sr_rec_lr.ny(),sr_rec_lr.nx());
  MAT.transpose(TempArray,tTempArray);
  fits_write_fltarr(sol_file_lr,tTempArray);
  delete[] sol_file;
  delete[] sol_file_lr;
  if (writeRes==True)
    {
      char* res_file = new char[size];
      sprintf(res_file,"%s/%s",output_dir.c_str(),resi_file_name.c_str());
      
      io_3d_write_data(res_file,resi);
      delete[] res_file;
    }
  if(writeParam==True)
    {
      struct stat s;
      int err = stat(output_dir.c_str(), &s);
      if(-1 == err) {
	if(ENOENT == errno) {
	  int e;
	  e = mkdir(output_dir.c_str(), S_IRWXU);
	  if(e!=0) {
	    printf("Cannot create directory\n");
	    exit(1);
	  }
	} else {
	  printf("Cannot create directory: file existing with same name\n");
	  exit(1);
	}
      } 
      char* mod_file = new char[size];
      sprintf(mod_file,"%s/%s",output_dir.c_str(),"model.fits");
      char* coarse_file = new char[size];
      sprintf(coarse_file,"%s/%s",output_dir.c_str(),"coarse_scale.fits");
      char* wvlt_noise_file = new char[size];
      sprintf(wvlt_noise_file,"%s/%s",output_dir.c_str(),"wvl_sig_noise.fits");
      char* noise_trace_file = new char[size];
      sprintf(noise_trace_file,"%s/%s",output_dir.c_str(),"noise_est_trace.fits");
      char* cent_file = new char[size];
      sprintf(cent_file,"%s/%s",output_dir.c_str(),"centroids.txt");
      char* shift_file = new char[size];
      sprintf(shift_file,"%s/%s",output_dir.c_str(),"shift.txt");
      char* flux_file = new char[size];
      sprintf(flux_file,"%s/%s",output_dir.c_str(),"flux.txt");
      char* sig_est_file = new char[size];
      sprintf(sig_est_file,"%s/%s",output_dir.c_str(),"sig_est.txt");
      char* mse_file = new char[size];
      sprintf(mse_file,"%s/%s",output_dir.c_str(),"mse.txt");
      char* user_parameters_file = new char[size];
      sprintf(user_parameters_file,"%s/%s",output_dir.c_str(),"param.txt");
  
      FILE *pFile1,*pFile2,*pFile3,*pFile4,*pFile5,*pFile6;
      pFile1 = fopen (cent_file,"w");
      pFile2 = fopen (shift_file,"w");
      pFile3 = fopen (flux_file,"w");
      pFile4 = fopen (mse_file,"w");
      pFile5 = fopen (user_parameters_file,"w");
      pFile6 = fopen (sig_est_file,"w");
      for (i=0;i<Nz;i++)
	{
      
	  fprintf(pFile1,"%f %f\n",cent[0][i],cent[1][i]);
	  fprintf(pFile2,"%f %f\n",shift[0][i],shift[1][i]);
	  fprintf(pFile3,"%f\n",flux[i]);
	  fprintf(pFile6,"%lf\n",noise_sig_map(0,0,i));
     
	}
      for (i=0;i<nb_iter;i++)
	{
	  fprintf(pFile4,"%lf ",mse_steep[i]);
	}
      fprintf(pFile5,"Estimated spectral radius: %f\n",spectral_rad);
      fprintf(pFile5,"nsigma: %f\n",nsig[0]);
      fprintf(pFile5,"Type of transform: %d\n",sr_opt.transf);
      fprintf(pFile5,"Type of filter: %d\n",sr_opt.SB_Filter);
      fprintf(pFile5,"Noise estimation method: %d\n",noise_est_meth);
      fprintf(pFile5,"Model denoising method: %d\n",mod_denoising_type);
  
      TempArray.alloc(model.buffer(),model.ny(),model.nx());
      MAT.transpose(TempArray,tTempArray);
      fits_write_fltarr(mod_file,tTempArray);
//      io_write_ima_float(mod_file,model);
      TempArray.alloc(coarse_scale.buffer(),coarse_scale.ny(),coarse_scale.nx());
      MAT.transpose(TempArray,tTempArray);
      fits_write_fltarr(coarse_file,tTempArray);
//      io_write_ima_float(coarse_file,coarse_scale);
      if(min_id==0)
	{
	  io_write_ima_float(noise_trace_file,noise_est_trace);
	  io_3d_write_data(wvlt_noise_file,wvl_noise_map);
	}
      fclose (pFile1);
      fclose (pFile2);
      fclose (pFile3);
      fclose (pFile4);
      fclose (pFile5);
      fclose (pFile6);
      delete[] mod_file;
      delete[] coarse_file;
      delete[] wvlt_noise_file;
      delete[] noise_trace_file;
      delete[] cent_file;
      delete[] shift_file;
      delete[] flux_file;
      delete[] sig_est_file;
      delete[] mse_file;
      delete[] user_parameters_file;
    }
}

// ------------- Main procedure --------------- //
void Sprite::do_it(Bool wr_en,Bool writeRes,Bool writeParam,int noise_sim_iter,Bool thresh_type,int nb_rw)
{
   int nb_iter_spec;
   cout << "Spectral radius estimation..."<< endl;
   spectral_rad_est(&nb_iter_spec);
   Ifloat delta(Nx*upfact,Ny*upfact);
   if(noise_est_meth==1)
     {
       cout << "Noise simulation..."<< endl;
       if(noise_sim_iter<1)
	 noise_sim_iter = nb_rw*nb_iter;
       noise_simulator(noise_sim_iter);
     }
   /*if(noise_est_meth==2)
     {
       cout << "Calcultating noise levels..."<< endl;
       noise_computer();
       }*/
   cout << "Minimizer..."<< endl;
   if(min_id==0)
      min1(delta,thresh_type,nb_rw);
   if (wr_en==True) 
     {
       cout << "Outputs saving..."<< endl;
       write_outputs(writeRes,writeParam);
     }
   delta.free();
}

Ifloat Sprite::get_sol()
{
  return sr_rec;
}
Ifloat Sprite::get_mod()
{
  return model;
}
double** Sprite::get_cent()
{
  return cent;
}
double** Sprite::get_shift()
{
  return shift;
}
double* Sprite::get_flux()
{
  return flux;
}
void Sprite::set_spectral_rad(double spectral_rad)
{
  this->spectral_rad = spectral_rad; 
}


Sprite::~Sprite()
{
  model.free();
  sr_rec.free();
  gradient.free();
  lanczos_ker.free();
  data_nrm.free();
  resi.free();
  resi_out.free();
  sr_rec_lr.free();
  op_mat_output.free();
  noise_sig_map.free();
  noise_est_trace.free();
  wvl_noise_map.free();
  free(cent[0]);
  free(cent[1]);
  free(cent);
  free(shift[0]);
  free(shift[1]);
  free(shift);
  free(shift_cent[0]);
  free(shift_cent[1]);
  free(shift_cent);
  free(flux);
  free(dist_weights);
  free(mse_steep);
  free(nsig);
}
