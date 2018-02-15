
//void building_signal_ascii (char *Name_Imag_In, Iint &Event_Image, 
//                                 Ifloat &Image);


void building_signal_signal (fltarray &Signal, intarray &EventSignal);

void event1d_one_scale (intarray&   EventSignal, 
                        int         s, 
		        intarray&   EventCount, 
		        type_border Border,
		        Bool        WriteAllInfo=False);

void event1d_set_support(MR_1D&            Mr1d_Data, 
                         int               CurrentScale, 
                         intarray&         EventSignal, 
		         type_border       Border,
                         const Ifloat&     Abaque, 
		         MR1DNoiseModel&   Mr1d_NoiseModel,
		         Bool              WriteAllInfo=False);


void mr1d_psupport(MR_1D&          Mr1d_Data, 
                   MR1DNoiseModel& Mr1d_NoiseModel, 
                   type_border     Border,
		   Bool            WriteAllInfo);

void mr1d_psupport(intarray&       EventSignal, 
                   MR_1D&          Mr1d_Data, 
                   Ifloat&         Abaque, 
		   MR1DNoiseModel& Mr1d_NoiseModel, 
                   type_border     Border, 
		   Bool            WriteAllInfo);
