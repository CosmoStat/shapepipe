
#include "IM_Obj.h"
#include "IM_IO.h"
#include "MR1D_Obj.h"
#include "IM_Noise.h"
#include "MR1D_NoiseModel.h"







/*void building_imag_ascii(char *Name_Imag_In,Iint &Event_Image,Ifloat &Image)
{
 int c;
 FILE *input;
 int Nl,Nc,pixel_x,pixel_y,x_cur,y_cur,x_pos,y_pos,p,q;
 float x_lu, y_lu,distance_x,distance_y,energy_x,energy_y;
	

 // Read the size of the input image 
 input = fopen(Name_Imag_In,"r");
 if (input == NULL) {
   cout << "Error in allocation of file " <<  Name_Imag_In << " ... or file doesn't exist" << endl;
   exit(-1);
 }
 fscanf(input,"%d\t%d\n",&Nl,&Nc);
 // create the resulting image 
 Event_Image.resize(Nl,Nc);
 Event_Image.init(0);
 Image.resize(Nl,Nc);
 Image.init(0.0);

 // Read the input file 
 while ( (c=getc(input)) != EOF )	
 {
 	ungetc(c,input);
	if (fscanf(input,"%f\t%f\n",&x_lu,&y_lu) == 2)
	{
	if ((x_lu < 0) || (x_lu >= Nl-0.5) || (y_lu < 0) || (y_lu >= Nc-0.5))
	{
	   cerr << "Error: incorrect event coordinates: " << x_lu << " " << y_lu << endl;
	   cerr << "       x coordinates must be between 0 and " << Nl-0.5 << endl;
	   cerr << "       y coordinates must be between 0 and " << Nc-0.5 << endl;
	   exit (-1);
	}

	// determine the coordinates of the left low corner pixel 
	// containing (x_lu,y_lu).                                
	pixel_x = (int) floor(x_lu);
	pixel_y = (int) floor(y_lu);
	
	// counting the number of events 
	(ABS(x_lu - pixel_x)>0.5) ? (p=pixel_x + 1) : (p=pixel_x);
	(ABS(y_lu - pixel_y)>0.5) ? (q=pixel_y + 1) : (q=pixel_y);
	Event_Image(p,q) ++;
	
	// examine the neighbourhood of real pixel for computing 
	// Bspline filtering                                     
	for (x_cur=-2;x_cur<=2;x_cur++) 
	{
	   x_pos = x_cur + pixel_x;
	   distance_x = x_pos-x_lu;
	   if (distance_x <= 2.0)
 	     energy_x = (CUBE_FABS(distance_x - 2) + CUBE_FABS(distance_x + 2)
	          - 4 * (CUBE_FABS(distance_x - 1) + CUBE_FABS(distance_x + 1))
		  + 6 * CUBE_FABS(distance_x)) / 12.;
	   else  energy_x = 0.;
	   for (y_cur=-2;y_cur<=2;y_cur++) 
	   {
	        y_pos = y_cur + pixel_y;
		distance_y = y_pos-y_lu;
 		if (distance_y <= 2.0)
		{
 	          energy_y = (CUBE_FABS(distance_y - 2) + CUBE_FABS(distance_y + 2)
	          - 4 * (CUBE_FABS(distance_y - 1) + CUBE_FABS(distance_y + 1))
		  + 6 * CUBE_FABS(distance_y)) / 12.;  	           
 		}   
		else energy_y = 0.;
		if ((x_pos >= 0) && (y_pos >= 0) && (x_pos < Nl) && (y_pos < Nc)) 
		                       Image(x_pos,y_pos) += energy_x * energy_y ;
	   }   
	}   	
      }
    }

}
   

*/

#define     CUBE_FABS(x)   (fabs(x) * fabs(x) * fabs(x))
/************************************************************************/

void building_signal_signal(fltarray &Signal, intarray &Event_Image) {	

   int i;
   
   // init_random();
   int Nx = Signal.nx();
   Event_Image.alloc(Nx);
 
   //fits_write_fltarr ("Bef.fits", Signal);

   /* the initial signal already containts the number of events */
   for (i=0;i<Nx;i++)
      Event_Image(i) = (int) (Signal(i) + 0.5);
      
   Signal.init();
 
   //io_write_ima_float ("Event.fits", Event_Image);
   /* Filter with Bspline */
   for (i=0;i<Nx;i++) {
    
      for (int e=0; e <  Event_Image(i); e++) {
       
         float Dxr=0.;
         // Dxr = get_random(-0.5,0.5);
         for (int x_cur=-2;x_cur<=2;x_cur++) {
	   
	    float distance_x = x_cur + Dxr;
 	    float energy_x = ( CUBE_FABS(distance_x - 2) + CUBE_FABS(distance_x + 2)
	                - 4 * (CUBE_FABS(distance_x - 1) + CUBE_FABS(distance_x + 1))
                        + 6 * CUBE_FABS(distance_x)) / 12.;


 	    if ((i+x_cur >= 0) && (i+x_cur < Nx)) {
	       Signal(i+x_cur)  += energy_x;
	    }  
         }
      }
   }
 
  //fits_write_fltarr ("After.fits", Signal);

}











void event1d_one_scale (intarray&   EventSignal, 
                        int         s, 
		        intarray&   EventCount, 
		        type_border Border,
		        Bool        WriteAllInfo) {

   // support at scale s 
   int Size = (int) (pow((double)2., (double)(s+2)) + 0.5);
   
   // test border
   if (Border != I_MIRROR  && Border != I_CONT)
      cout << "Only Border I_MIRROR and I_CONT are implemented " 
           << "in event count (func event_one_scale)!" << endl;
	   
   int NbPoint = EventSignal.n_elem();
   for (int i=0; i<NbPoint; i++) {
       
      int Total = 0;
      if (Border == I_MIRROR) {
         for (int k =-Size; k <= Size; k++) {
      
            int ind = (i+k < 0 ? abs(i+k) : i+k);
	    ind = (ind > (NbPoint-1) ? (NbPoint-1)-(ind - (NbPoint-1)) : ind);
            Total += EventSignal(ind);
         }
      }
      if (Border == I_CONT) {
         for (int k =-Size; k <= Size; k++) {
      
            int ind = (i+k < 0 ? 0 : i+k);
	    ind = (ind > (NbPoint-1) ? (NbPoint-1) : ind);
            Total += EventSignal(ind);
         }
      }     
      

      EventCount(i) = Total;
   }
    
   if (WriteAllInfo) {     
      cout << "Scale:" << s << ", Size:" << 2*Size+1 << endl;
      char FileName[256];
      sprintf (FileName, "EventCount_%d", s);
      fltarray prov(EventSignal.n_elem());
      for (int i=0;i<EventSignal.n_elem();i++) prov(i)=(float)EventCount(i);
      fits_write_fltarr(FileName,prov); 
   }  
}







void event1d_set_support(MR_1D&            Mr1d_Data, 
                         int               CurrentScale, 
                         intarray&         EventSignal, 
		         type_border       Border,
                         const Ifloat&     Abaque, 
		         MR1DNoiseModel&   Mr1d_NoiseModel,
		         Bool              WriteAllInfo) {

   // compute 4^(CurrentScale-1)
   //for (int s=0,float alpha=1.0;s<CurrentScale;s++) alpha *= 1;
   float Alpha=1.;
   if (CurrentScale != 0) Alpha = pow((double)2., (double)(CurrentScale-1));
   int NbPoint = Mr1d_Data.size_ima_np();
   
   // number of event
   intarray EventCount (NbPoint);
   event1d_one_scale (EventSignal, CurrentScale, EventCount, Border,WriteAllInfo);
  
   if (WriteAllInfo) { 
      char FileName[256];
      sprintf (FileName, "Mr1d_%d", CurrentScale);
      fltarray prov(EventSignal.n_elem());
      Mr1d_Data.scale (prov,CurrentScale);
      fits_write_fltarr(FileName,prov);
   }      
   
   // for all points !!???!!! same nb points at all scale !!???!!!
   //
   for (int i=0;i<NbPoint;i++) {

      // reading the number of events
      //
      int NEventReal = EventCount(i);

      // compute "power" : 2^(power-1) < n_event_real <= 2^(power) 
      //
      int Power = 0;
      int NEvent = 1;
      while (NEvent<NEventReal) {Power++; NEvent *= 2;}
      
      
      // interpolation of thresholds using abaque
      //
      float SeuilMin, SeuilMax;
      if (Power > Abaque.n_elem()) {
        
	 // nb scale in abaque to small, take the max value in abaque     
         Power = Abaque.n_elem();
	 SeuilMin = Abaque(Power,0);
	 SeuilMax = Abaque(Power,1);     
     
      } else if (NEventReal == NEvent) {
     
         // NEventReal is a power of two = NEvent, read value in Abaque tab
         SeuilMin = Abaque(Power,0);
         SeuilMax = Abaque(Power,1);
	 
      } else {
      
         // interpolate in Abaque tab between level power-1 and power
	 if (Power==0) {
	    SeuilMin = Abaque(0,0);
	    SeuilMax = Abaque(0,1);
	 } else {
	    // linear interp 
	    SeuilMin = -2*(Abaque(Power,0) - Abaque(Power-1,0))
                          *(NEvent - NEventReal) / NEvent 
			  + Abaque(Power,0);
            SeuilMax = -2*(Abaque(Power,1) - Abaque(Power-1,1))
                          *(NEvent - NEventReal) / NEvent 
			  + Abaque(Power,1);
         }
      }	
      
      // sigma of wavelet IN (function psy of bspline_histo_1D)
      //
      float Sigma=0.121070;
      //Sigma=0.0405078;    
      if (Power > 0) {
         SeuilMin *= sqrt((float) NEventReal) * Sigma / Alpha;
	 SeuilMax *= sqrt((float) NEventReal) * Sigma / Alpha;
      } else {
	 SeuilMin *= Sigma  / Alpha;
	 SeuilMax *= Sigma  / Alpha;
      }
           

      // detection of signal : comparison with thresholds
      //
      int s = CurrentScale;
      Mr1d_NoiseModel.support(s,i) = VAL_SupNull;
      if ((Mr1d_Data(s,i)<=SeuilMin) || (Mr1d_Data(s,i)>=SeuilMax)) {
	     
         Mr1d_NoiseModel.support(s,i) = VAL_SupOK;
	 	 
	 // 
	 if (NEventReal<Mr1d_NoiseModel.MinEventNumber) 
            Mr1d_NoiseModel.support(s,i) = VAL_SupMinEv;
	 
	 //
	 if ((Mr1d_NoiseModel.OnlyPositivDetect) && (Mr1d_Data(s,i)<0))
	    Mr1d_NoiseModel.support(s,i) = VAL_SupNull;
	 
	 //    
	 if (Mr1d_NoiseModel.FirstDectectScale > s) 
	    Mr1d_NoiseModel.support(s,i) = VAL_SupFirstScale;
      }
      
      // Threshold level
      //
      float Sm = (ABS(SeuilMin) > ABS(SeuilMax)) ? ABS(SeuilMin) : ABS(SeuilMax);
      Mr1d_NoiseModel.sigma(s,i) = Sm / Mr1d_NoiseModel.NSigma[s];
   }
}

void mr1d_psupport(MR_1D&          Mr1d_Data, 
                   MR1DNoiseModel& Mr1d_NoiseModel, 
                   type_border     Border,
		   Bool            WriteAllInfo) {
   
   if  (WriteAllInfo) cout << "Compute threshold..." << endl;			
   int NScale = Mr1d_Data.nbr_scale();
   for (int s=0; s< NScale-1; s++) {
      Mr1d_NoiseModel.CFewEventPoisson1d->find_threshold (
                        (Mr1d_NoiseModel.TabEps)[s],
			WriteAllInfo);
      event1d_set_support(Mr1d_Data, s, Mr1d_NoiseModel.EventSignal,  Border, 
                          Mr1d_NoiseModel.CFewEventPoisson1d->_Threshold, 
			  Mr1d_NoiseModel, WriteAllInfo);
   }
}


void mr1d_psupport(intarray&       EventSignal, 
                   MR_1D&          Mr1d_Data, 
                   Ifloat&         Abaque, 
		   MR1DNoiseModel& Mr1d_NoiseModel, 
                   type_border     Border, 
		   Bool            WriteAllInfo) {
   
   if  (WriteAllInfo) cout << "Compute threshold..." << endl; 
   int NScale = Mr1d_Data.nbr_scale();

   // writing on disk the coefficients of wavelet 
   // !!!!!!! if (WriteAll == True)  Mr1d_NoiseModel.write("Wavelet_Support");

   // detecting signal : comparison with the ABAQUE of thresholds */
   if  (WriteAllInfo) cout << "Detect signal ... " << endl;
   for (int s=0;s< NScale-1; s++) 
      event1d_set_support(Mr1d_Data, s, EventSignal, Border, Abaque, 
                          Mr1d_NoiseModel, WriteAllInfo);
}

