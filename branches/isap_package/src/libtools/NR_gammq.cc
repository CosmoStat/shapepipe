#include "NR.h"

float gammq(float a,float x)
{
	float gamser,gammcf,gln;

	if (x < 0.0 || a <= 0.0) nrerror("Invalid arguments in routine GAMMQ");
	if (x < (a+1.0)) {
		gser(&gamser,a,x,&gln);
		return 1.0-gamser;
	} else {
		gcf(&gammcf,a,x,&gln);
		return gammcf;
	}
}
