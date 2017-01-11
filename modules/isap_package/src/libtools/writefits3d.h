#ifndef _IOFITS_H_
#define _IOFITS_H_


#include "IM_IO.h"
#include "IM_Obj.h"
#include "IM_IOTools.h"

void writefltarr(char * filename, fltarray &Mat);
void writefltarr(char * filename, dblarray &Mat);
void writecfarr(char * filename, cfarray &Mat);



#endif
