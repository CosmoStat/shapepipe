/******************************************************************************
**                   Copyright (C) 1994 CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 3.1
**
**    Author: J.L. Starck
**
**    Date:  96/05/02 
**    
**    File:  GetOpt.cc
**
*******************************************************************************
**
**    DECRIPTION  this is a public domain version of getopt  
**    ---------- 
**
*****************************************************************************/

// static char sccsid[] = "@(#)getopt.cc 3.1 96/05/02 CEA 1994 @(#)";

#include <stdio.h>
#include <string.h>

#define ERR(ms,cc) if(OptErr) {(void) fprintf(stderr,"%s%s%c\n",argv[0],ms,cc);}

extern int  strcmp();
extern char *strchr();

int  OptErr = 1;
int  OptInd = 1;
int  OptOpt;
char *OptArg;

int GetOpt(int argc, char **argv, char *opts)
{
static int sp = 1;
int c;
char *cp;

   if(sp == 1)
   {
      if(OptInd >= argc || argv[OptInd][0] != '-' || argv[OptInd][1] == '\0')
			return(EOF);
      else if(strcmp(argv[OptInd], "--") == 0) 
           {
	      OptInd++;
	      return(EOF);
           }
    }

    OptOpt = c = argv[OptInd][sp];
    if(c == ':' || (cp=strchr(opts, c)) == NULL) 
    {
	ERR(": unknown option, -", c);
	if(argv[OptInd][++sp] == '\0') 
        {
	   OptInd++;
	   sp = 1;
	}
	return('?');
    }
	
    if(*++cp == ':') 
    {
	if(argv[OptInd][sp+1] != '\0')
			OptArg = &argv[OptInd++][sp+1];
	else if(++OptInd >= argc) 
             {
	        ERR(": argument missing for -", c);
	        sp = 1;
	        return('?');
	     }
             else OptArg = argv[OptInd++];
	sp = 1;
     } 
     else
     {
	if(argv[OptInd][++sp] == '\0') 
        {
	   sp = 1;
	   OptInd++;
	}
	OptArg = NULL;
     }
    return(c);
}
