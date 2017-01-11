/******************************************************************************
**                   Copyright (C) 1998 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 2.1
**
**    Author: Jean-Luc Starck
**
**    Date:  98/05/12 
**    
**    File:  SoftInfo.h
**
*******************************************************************************
**
**    DESCRIPTION  
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

#ifndef _SOFT_H_
#define _SOFT_H_

#include <string.h>

#define ADRESS "(DAPNIA CEA-Saclay France)"
#define MR1_RELEASE 4.0
#define MR1_NAME "MR/1"
 
#define MR2_RELEASE 1.2
#define MR2_NAME "MR/2"

#define MR3_RELEASE 2.0
#define MR3_NAME "MR/3"

#define MR4_RELEASE 1.0
#define MR4_NAME "MR/4"



class softinfo {
     float Release;
     char Name[256];
     char Banner[256];
     void soft_init()
     {
#ifdef KBUFF    
        setbuf(stdout, NULL);
        setbuf(stdin, NULL);
        setbuf(stderr, NULL);
#endif
        mr1();
    }
    public:
     void iso()
     {
        Release = 1.0;
	strcpy(Name, "ISO");
	strcpy(Banner,  "ISO (DAPNIA CEA-Saclay France)");
     }
     void mr1()
     {
        Release = MR1_RELEASE;
	strcpy(Name, MR1_NAME);
	sprintf(Banner, "%s V%2.1f %s", Name, MR1_RELEASE, ADRESS);
     }
     void mr2()
     {
        Release = MR2_RELEASE;
	strcpy(Name, MR2_NAME);
 	sprintf(Banner, "%s V%2.1f %s", Name, MR2_RELEASE, ADRESS);
     }
     void mr3()
     {
        Release = MR3_RELEASE;
	strcpy(Name, MR3_NAME);
 	sprintf(Banner, "%s V%2.1f %s", Name, MR3_RELEASE, ADRESS);
     }
     void mr4()
     {
        Release = MR4_RELEASE;
	strcpy(Name, MR4_NAME);
 	sprintf(Banner, "%s V%2.1f %s", Name, MR4_RELEASE, ADRESS);
     }
     softinfo()  { soft_init();}
     float release() { return Release;}
     char *name() { return Name;} 
     char *banner() { return Banner;}
};
// extern softinfo Soft;

#endif



