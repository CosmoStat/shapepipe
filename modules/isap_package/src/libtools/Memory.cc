/*******************************************************************************
**
**    UNIT
**
**    Version: 3.2
**
**    Author: Jean-Luc Starck
**
**    Date:  96/06/13 
**    
**    File:  Memory.cc
**
*******************************************************************************
**
**    DECRIPTION  This module contains allocation and free procedures
**    ---------- 
**
******************************************************************************/

#include "GlobalInc.h"

/****************************************************************************/

TempCMem<int> MemInt;
TempCMem<float> MemFloat;
TempCMem<double> MemDouble;
TempCMem<unsigned long> MemULong;//FCS ADDED
TempCMem<complex_f> MemCF;
TempCMem<complex_d> MemCD;
Bool UseVMS = False;

void memory_abort () {
    cerr << "Error: cannot allocate memory ... " << endl;
    exit (0);
}

char *alloc_buffer(size_t  Nelem) {
   char *Ptr;
   extern char *mem_alloc_buffer(size_t  Nelem);
   if (UseVMS == True) Ptr =  mem_alloc_buffer(Nelem);
   else Ptr = (char *) malloc (Nelem);
   if (Ptr == NULL) memory_abort ();       
   return Ptr;
}

void free_buffer(char *Ptr) {
   extern void mem_free_buffer(char *Ptr);
   if (Ptr != NULL) {
      if (UseVMS == True) mem_free_buffer(Ptr);
      else free (Ptr);
      Ptr = NULL;
   }
}
