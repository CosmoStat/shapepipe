/******************************************************************************
**                   Copyright (C) 2000 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  24/02/00
**    
**    File:  IM_Lut.h
**
**
*******************************************************************************
**
**    DESCRIPTION  LUT information
**    ----------- 
**                 
******************************************************************************/

#ifndef _IM_LUT_H_
#define _IM_LUT_H_

#include"IM_Color.h"

#define IO_STD_SIZE_TAB_LUT 256



class C_IO_RGB {
   int SizeLut;
   byte  *r, *g, *b;
  public:
    C_IO_RGB() {r = NULL; g = NULL; b = NULL; SizeLut = 0;}
    int lut_size() { return SizeLut; }
    inline byte red(int i) {return r[i];}
    inline byte green(int i) {return g[i];}
    inline byte blue(int i) {return b[i];}
    byte *red () {return r;}
    byte *green() {return g;}
    byte *blue() {return b;}
    void alloc(int InSizeLut) 
    {
       int i;
       if (SizeLut != 0)
       {
          if (r != NULL) delete [] r;
	  if (g != NULL) delete [] g;
	  if (b != NULL) delete [] b;
       }
       SizeLut = InSizeLut;
       r = new byte[MAX(SizeLut,256)];
       g = new byte[MAX(SizeLut,256)];
       b = new byte[MAX(SizeLut,256)];       
       for (i=0;i < SizeLut; i++)
       {
          r[i] = g [i] = b[i] = (byte) i;
       }
       for (i= SizeLut;i < 256; i++)
       {
          r[i] = g [i] = b[i] = (byte) 0;
       }
    }
    void alloc_color() 
    {
       alloc(256);
       for (int i=0;i < SizeLut; i++)
       {
          r[i] = 16+(i & 224);
          g[i] = 16+((i & 28)<<3);
          b[i] = 32+((i & 3)<<6);
       }
    }
    void alloc(int InSizeLut, byte *ri, byte *gi, byte *bi)
    {
       alloc(InSizeLut);
       for (int i=0;i < SizeLut; i++)
       {
          r[i] = ri[i];
	  g[i] = gi[i];
	  b[i] = bi[i];
       }
    }
    ~C_IO_RGB() 
    { 
       if (SizeLut != 0)
       {
          if (r != NULL) delete []r;
	  if (g != NULL) delete []g;
	  if (b != NULL) delete []b;
       }
       SizeLut = 0;
    }
    inline void pseudo_to_color(byte Val, byte &Red, byte &Green, byte &Blue)
    {
       Red = red(Val);
       Green = green(Val);
       Blue = blue(Val);
    }
    inline byte color_to_pseudo(byte Red, byte Green, byte Blue)
    {
        byte ValRet = ((Red>>5)<<5 | (Green>>5)<<2 | (Blue>>6));
        return ValRet;
    }
};
extern C_IO_RGB IO_RGB;

#endif







