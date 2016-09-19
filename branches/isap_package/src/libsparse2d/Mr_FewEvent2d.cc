
#include "IM_Obj.h"
#include "IM_IO.h"
#include "MR_Obj.h"
#include "IM_Noise.h"
#include "MR_NoiseModel.h"

#include "MR_Psupport.h" ;to remove

void event2d_one_scale (Iint&   EventImage, 
                        int         s, 
		        Iint&       EventCount, 
		        type_border Border,
		        Bool        WriteAllInfo) {

   // support at scale s 
   int Size = (int) (pow((double)2., (double)(s+2)) + 0.5);
   	   
   int Nl = EventImage.nl();
   int Nc = EventImage.nc();
   for (int i=0; i<Nl; i++) {
       
      int Total = 0;
      int j,k;
      for (k =-Size; k <= Size; k++) 
         for (j =-Size; j <= Size; j++) 
	    Total += EventImage(i+k, j, Border);
       
      EventCount(i,0) = Total;
      for (j =1; j < Nc; j++) {
         for (k = -Size; k <= Size; k++) {
            Total -= (int) EventImage(i+k,j-Size-1,Border);
            Total += (int) EventImage(i+k,j+Size, Border);
         }
         EventCount(i,j) = Total;
      }	    

   }
}



void event2d_set_support(MultiResol&       Mr2d_Data, 
                         int               CurrentScale, 
                         Iint&             EventSignal, 
		         type_border       Border,
                         const Ifloat&     Abaque, 
		         MRNoiseModel&     Mr2d_NoiseModel,
		         Bool              WriteAllInfo) {
			 
   std::cout << "event2d_set_support" << std::endl;

   // compute 4^(CurrentScale-1)
   //for (int s=0,float alpha=1.0;s<CurrentScale;s++) alpha *= 1;
   float Alpha=1.;
   if (CurrentScale != 0) Alpha = pow((double)4., (double)(CurrentScale));
   int Nl = Mr2d_Data.size_ima_nl();
   int Nc = Mr2d_Data.size_ima_nc();
   
   // number of event
   Iint EventCount (Nl,Nc,"");
   event2d_one_scale( EventSignal, CurrentScale, EventCount, Border,WriteAllInfo);
   
   // for all points !!???!!! same nb points at all scale !!???!!!
   //
   for (int i=0;i<Nl;i++)
      for (int j=0;j<Nc;j++) {
      
      // reading the number of events
      //
      int NEventReal = EventCount(i,j);

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
	    float OldSeuilMin = -2*(Abaque(Power,0) - Abaque(Power-1,0))
                          *(NEvent - NEventReal) / NEvent 
			  + Abaque(Power,0);
            float OldSeuilMax = -2*(Abaque(Power,1) - Abaque(Power-1,1))
                          *(NEvent - NEventReal) / NEvent 
			  + Abaque(Power,1);
			  
            double weight = Power - log((double)NEventReal)/log(2.);
            float seuilMaxAfter = Abaque(Power,1);
	    float seuilMaxBefore = Abaque(Power-1,1);
            SeuilMax =   weight * ( seuilMaxBefore ) + 
	                   + ( 1.-weight ) * ( seuilMaxAfter );
            float seuilMinAfter = Abaque(Power,0);
	    float seuilMinBefore = Abaque(Power-1,0);
            SeuilMin =    weight * ( seuilMinBefore )  + 
	                   + ( 1.-weight ) * ( seuilMinAfter );
			   
            if( WriteAllInfo ) {
	       std::cout << "      seuil min in [" << seuilMinAfter
	                 << "," << seuilMinBefore << "]" 
			 << ", seuilmin = " << SeuilMin  
			 << ", old seuilmin = " << OldSeuilMin << std::endl;
	       std::cout << "      seuil max in [" << seuilMaxAfter
	                 << "," << seuilMaxBefore << "]" 
			 << ", seuilmax = " << SeuilMax 
			 << ", old seuilmax = " << OldSeuilMax << std::endl;
	    }

         }
      }	
      
      // reduction of the wavelet coef
      //float Sigma=0.040717; //
      float Sigma=0.0405078;   
      int s = CurrentScale;
      float redWavCoef;    
      if( Power>0 )
         redWavCoef = Mr2d_Data(s,i,j) 
                       * Alpha / sqrt((float) NEventReal) / Sigma;
      else 
         redWavCoef = Mr2d_Data(s,i,j) * Alpha / Sigma;
      
      
      // detection of signal : comparison w red with thresholds
      //
      Mr2d_NoiseModel.support(s,i,j) = VAL_SupNull;
      if( (redWavCoef<=SeuilMin) || (redWavCoef>=SeuilMax) ) {    
      
         Mr2d_NoiseModel.support(s,i,j) = VAL_SupOK;
	 	 
	 // 
	 if (NEventReal<Mr2d_NoiseModel.MinEventNumber) 
            Mr2d_NoiseModel.support(s,i,j) = VAL_SupMinEv;
	 
	 //
	 if ((Mr2d_NoiseModel.OnlyPositivDetect) && (Mr2d_Data(s,i,j)<0))
	    Mr2d_NoiseModel.support(s,i,j) = VAL_SupNull;
	 
	 //    
	 if (Mr2d_NoiseModel.FirstDectectScale > s) 
	    Mr2d_NoiseModel.support(s,i,j) = VAL_SupFirstScale;
      }
    

      // Threshold level
      //      
      float Sm = (ABS(SeuilMin) > ABS(SeuilMax)) ? ABS(SeuilMin) : ABS(SeuilMax);
      if( Power > 0 ) 
         Mr2d_NoiseModel.sigma(s,i,j) = Sm 
	                              * sqrt((float) NEventReal) * Sigma / Alpha
				      / Mr2d_NoiseModel.NSigma[s];
      else 
         Mr2d_NoiseModel.sigma(s,i,j) = Sm * Sigma / Alpha 
	                              / Mr2d_NoiseModel.NSigma[s];
	 
	 


   }
}


void mr2d_psupport(MultiResol&     Mr2d_Data, 
                   MRNoiseModel&   Mr2d_NoiseModel, 
                   type_border     Border,
		   Bool            WriteAllInfo) {
   
   if  (WriteAllInfo) cout << "Compute threshold..." << endl;			
   int NScale = Mr2d_Data.nbr_scale();
   for (int s=0; s< NScale-1; s++) {
      Mr2d_NoiseModel.CFewEvent2d->find_threshold (
                        (Mr2d_NoiseModel.TabEps)[s],
			WriteAllInfo);
      event2d_set_support(Mr2d_Data, s, Mr2d_NoiseModel.Event_Image,  Border, 
                          Mr2d_NoiseModel.CFewEvent2d->mThreshold, 
			  Mr2d_NoiseModel, WriteAllInfo);
   }
}


void mr2d_psupport(Iint&           EventSignal, 
                   MultiResol&     Mr2d_Data, 
                   Ifloat&         Abaque, 
		   MRNoiseModel&   Mr2d_NoiseModel, 
                   type_border     Border, 
		   Bool            WriteAllInfo) {
   
   if  (WriteAllInfo) cout << "Compute threshold..." << endl; 
   int NScale = Mr2d_Data.nbr_scale();

   // detecting signal : comparison with the ABAQUE of thresholds */
   if  (WriteAllInfo) cout << "Detect signal ... " << endl;
   for (int s=0;s< NScale-1; s++) 
      event2d_set_support(Mr2d_Data, s, EventSignal, Border, Abaque, 
                          Mr2d_NoiseModel, WriteAllInfo);
}


