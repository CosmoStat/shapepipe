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
**    File:  Memory.h
**
*******************************************************************************
**
**    DESCRIPTION  Memory definitions
**    ----------- 
**                 
**    PARAMETRES    
**    ----------    
** 
**    RESULTS      
**    -------  
**
**
******************************************************************************/

#ifndef _MEMORY_H_
#define _MEMORY_H_

void memory_abort ();
char *alloc_buffer(size_t  Nelem) ;
void free_buffer(char *Ptr);

#include "GlobalInc.h"
#include "TempMemory.h"

#endif
