
void event2d_one_scale (Iint&   EventImage, 
                        int         s, 
		        Iint&       EventCount, 
		        type_border Border,
		        Bool        WriteAllInfo);

void event2d_set_support(MultiResol&       Mr2d_Data, 
                         int               CurrentScale, 
                         Iint&             EventSignal, 
		         type_border       Border,
                         const Ifloat&     Abaque, 
		         MRNoiseModel&     Mr2d_NoiseModel,
		         Bool              WriteAllInfo);

void mr2d_psupport(MultiResol&     Mr2d_Data, 
                   MRNoiseModel&   Mr2d_NoiseModel, 
                   type_border     Border,
		   Bool            WriteAllInfo);

void mr2d_psupport(Iint&           EventSignal, 
                   MultiResol&     Mr2d_Data, 
                   Ifloat&         Abaque, 
		   MRNoiseModel&   Mr2d_NoiseModel, 
                   type_border     Border, 
		   Bool            WriteAllInfo);
