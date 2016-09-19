/******************************************************************************
**                   Copyright (C) 1995 CEA
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
**    File:  Licence.h
**
*******************************************************************************
**
**    DECRIPTION    License file
**    ---------- 
**
*****************************************************************************/
 
#ifndef __LIC__
#define __LIC__

#include <ctype.h>
#include <time.h>
#ifdef SOL3
#include <sys/systeminfo.h>
#endif

#ifdef SYSINFO
#include <sys/systeminfo.h>
#endif

#ifdef RTU
#include <ctype.h>
#endif

#define NBR_LIC        10
#define LIC_NO         -1   /* No Licence   */
#define LIC_ALL         0   /* All products */
#define LIC_MRA         1   /* all multiresolution products */
#define LIC_MR1         2   /*  MR1 product */
#define LIC_MR2         3   /*  MR2 product */
#define LIC_MR3         4   /*  MR3 product */
#define LIC_MR4         5   /*  MR4 product */
#define LIC_POL         6   /*  ISO product */
#define LIC_XMM         7   /*  XMM  product */
#define LIC_CMB         8   /*  Planck  product */
#define LIC_M1D         9   /*  1D software */

void soft_init();
void lic_test_date();
void lic_test_user();
void lic_test_host();
void lm_check(int TypeLic);

#define DEMO1D_LIMIT_SIZE 512
#define DEMO2D_LIMIT_SIZE 256
#define DEMO3D_LIMIT_SIZE 30
#define DEMOCOL_LIMIT_SIZE DEMO2D_LIMIT_SIZE

class DemoLic
{
  public:
    Bool Verbose;
    int Limit1D;
    int Limit2D;
    int Limit3D;
    int LimitCol;
    Bool Active;
    DemoLic() {Limit1D=DEMO1D_LIMIT_SIZE;
	            Limit2D=DEMO2D_LIMIT_SIZE;
	            Limit3D=DEMO3D_LIMIT_SIZE;
	            LimitCol=DEMO2D_LIMIT_SIZE;Active=False;
		    Verbose=False;
	      }
    void test(int Nx);
    void test(int Nx, int Ny);
    void test(int Nx, int Ny, int Nz);
    ~DemoLic() {}
};


#endif

