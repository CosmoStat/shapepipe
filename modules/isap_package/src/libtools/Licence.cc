/******************************************************************************
**                   Copyright (C) 1998 CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.01
**
**    Author: J.L. Starck
**
**    Date:  7/12/98 
**    
**    File:  Licence.cc
**
*******************************************************************************
**
**    DECRIPTION    License file
**    ---------- 
**
*****************************************************************************/

#include "GlobalInc.h"
#include "Licence.h"

#define DemoDat    "30-dec-05"
#define DemoUser   "jstarck"
#define DemoHost   "rgmiller"

DemoLic MRDEMO;

#ifdef USELM

extern "C"
{
int CnxSrv(char *host, int *sock, int log);
int TestCEACode(int product);
}
#endif

void lic_test_date()
{
#ifdef DEMODAT
   char *ms;
   int  i,j, m, a;
   char timeLimit[80];	
   time_t  tloc;
   struct tm *tm;
   static char *Months[12] = { "jan", "feb", "mar", "apr", "may", "jun",
		"jul", "aug", "sep", "oct", "nov", "dec" };
   strcpy(timeLimit,  DemoDat);
   // cout << "limit date is : " << timeLimit << endl;
   j = atoi(timeLimit);
   ms = strchr(timeLimit, '-');
   if (ms) 
   {
	ms = strdup(strchr(timeLimit, '-') + 1);
	ms[3] = 0;
   } 
   else  ms = "aaa";
     
   if (strrchr(timeLimit, '-'))  a = atoi(strrchr(timeLimit, '-') + 1) + 100;
   else  a = 0;
   m = 0;
   for (i = 0; i < 12; i++) if (!strcmp(ms, Months[i]))  m = i;
   time(&tloc);
   tm = localtime(&tloc);
   if (tm->tm_year > a ||  (tm->tm_year == a &&
          (tm->tm_mon > m || (tm->tm_mon == m && tm->tm_mday > j)))) 
   {
	cerr << "Error: your license has expired ... " << endl;
	cout << "Expiration date is : " << timeLimit << endl;
	exit(-1);
   }
#endif 
}

void lic_test_user()
{
#ifdef DEMOUSER
   char *logname;
   if (getenv("LOGNAME"))  logname = getenv("LOGNAME");
   else if (getenv("USER")) logname = getenv("USER");
   else logname = "unknown";
   if (strcmp(logname, DemoUser))
   {
      cerr << "Error: license not available for user " << logname << endl;
      exit(-1);
   }
#endif 
}

void lic_test_host()
{
#ifdef DEMOHOST
    char    Host[255];
    Host[0] = 0;
    gethostname(Host, 255);
	
   if (strcmp(Host, DemoHost))
   {
      cerr << "Error: license not available for host " << Host << endl;
      exit(-1);
   }
#endif 
}

#ifdef USELM
static void err_connect(int Ret)
{
   if (Ret > 1)
     cerr << "Error: no license available. " << endl;
   else 
     cerr << "Error: cannot connect to the license manager ... " << endl;
   exit(-1);
}
#endif  

#ifdef USELM
static int lm_lic_check(int TypeLic)
{
   int sock;
#ifdef NOSOC
    int V = TestCEACode(TypeLic);
    // cout << "lm_lic_check : NOSOC " << V << endl;
    return (V);
#else
    int V = CnxSrv(NULL, &sock, TypeLic);
    // cout << "lm_lic_check : LM " << V << endl;
    return (V);
#endif
}
#endif
  
 
void lm_check(int TypeLic)
{
   Bool DebugLM = False;
#ifdef USELM 
     int ret=1;
     if ((TypeLic > 0) && (TypeLic < NBR_LIC))
     {
     switch (TypeLic)
     {
        case LIC_ALL:
	  ret = lm_lic_check(LIC_ALL);
	  // if (ret != 1) err_connect(ret);
	  break;
	case LIC_MRA:
	  ret = lm_lic_check(LIC_MRA);
	  // if (ret != 1) err_connect(ret);
	  break;
	case LIC_MR1:
	  ret = lm_lic_check(LIC_MR1);
	  if (ret != 1) ret = lm_lic_check(LIC_MRA);
	  // if (ret != 1) err_connect(ret);
	  break;
	case LIC_MR2:
	  ret = lm_lic_check(LIC_MRA);
	  if (ret != 1) ret = lm_lic_check(LIC_MR2);
	  // if (ret != 1) err_connect(ret);
	  break;
	case LIC_MR3:
	  ret = lm_lic_check(LIC_MRA);
	  if (ret != 1) ret = lm_lic_check(LIC_MR3);
	  // if (ret != 1) err_connect(ret);
	  break;
	case LIC_MR4:
	  ret = lm_lic_check(LIC_MRA);
	  if (ret != 1) ret = lm_lic_check(LIC_MR4);
	  // if (ret != 1) err_connect(ret);
	  break;
	case LIC_POL:
	  ret = lm_lic_check(LIC_POL);
	  if (ret != 1) ret = lm_lic_check(LIC_POL);
	  if (ret != 1) ret = lm_lic_check(LIC_MRA);
 	  break;
// 	case LIC_CMB:
// 	   ret = lm_lic_check(LIC_SAP);
// 	  if (ret != 1) ret = lm_lic_check(LIC_CMB);
// 	  // if (ret != 1) err_connect(ret);
// 	  break;
	case LIC_M1D:
	  ret = lm_lic_check(LIC_M1D);
	  if (ret != 1) ret = lm_lic_check(LIC_MR1);
	  if (ret != 1) ret = lm_lic_check(LIC_MRA);
 	  break;   
        default:
          cerr << "Error: unknown licence type ... " << endl;
	  exit(-1);
	  break;
   }
   if (ret != 1) 
   {
       if (DebugLM == True) cout << "  SET ACTIVE " << endl;
       MRDEMO.Active = True;
   }
   else if (DebugLM == True) cout << "  SET NO ACTIVE " << endl;
   }
#else
 if (DebugLM == True) cout << "  NO LM  " << endl;
 int a = TypeLic; // just to suppress the warning
 if (a) a = 0;
#endif

#ifdef DEMODAT
lic_test_date();
#endif
#ifdef DEMOUSER
lic_test_user();
#endif
#ifdef DEMOHOST
lic_test_host();
#endif
}

/****************************************************************************/

void DemoLic::test(int Nx)
{
    if ((Active == True) && (Limit1D != Nx))
    {
        printf("Error: need a signal of size = %d\n", Limit1D);
	exit(-1); 
    }
    if ((Active == True) && (Verbose == True)) cout << " DEMO MODE " << Nx << " OK " << endl;
}
// **************************

void DemoLic::test(int Nx, int Ny)
{
    if ((Active == True) && ((Limit2D != Nx) || (Limit2D != Ny)))
    {
	printf("Error: need an image of size = %d X %d\n", Limit2D, Limit2D);
	exit(-1); 
    }
    if ((Active == True) && (Verbose == True)) cout << " DEMO MODE " << endl;
}

// **************************

void DemoLic::test(int Nx, int Ny, int Nz)
{ 
    if (Nz != 3)
    {
        if ((Active == True) && ((Limit3D != Nx) || (Limit3D != Ny) || (Limit3D != Nz)))
	{
	    printf("Error: need a cube of size = %d X %d X %d\n", Limit3D, Limit3D, Limit3D);
	    exit(-1); 
	 }
    }
    else if ((Active == True) && ((Limit2D != Nx) || (Limit2D != Ny)))
    {
	printf("Error: need an image of size = %d X %d\n", Limit2D, Limit2D);
	exit(-1); 
    }
    if ((Active == True) && (Verbose == True)) cout << " DEMO MODE " << endl;
} 
		    
/****************************************************************************/
