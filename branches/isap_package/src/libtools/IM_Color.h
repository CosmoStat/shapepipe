/******************************************************************************
**                   Copyright (C) 2001 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  24/04/01
**    
**    File:  IM_Color.h
**
**
*******************************************************************************
**
**    DESCRIPTION  Color conversion
**    ----------- 
**                 
******************************************************************************/

#ifndef _IM_COL_H_
#define _IM_COL_H_

void XYZ_to_LUV (float  X, float  Y, float  Z,
                 float  & L, float  & U, float  & V);
void LUV_to_XYZ (float  L, float  U, float  V, 
                 float  & X, float  & Y, float &  Z);		 
		 
void RGB_to_XYZ (float  R, float G, float  B,
                 float & X, float & Y, float  &Z);
void XYZ_to_RGB(float  X, float Y, float  Z,
                float & R, float & G, float  &B);		    		 

void RGB_to_LUV (float  R, float  G, float  B,
                 float  & L, float  & U, float  & V);
void LUV_to_RGB (float L, float  U, float  V,
                 float  & R, float  & G, float  & B);
		    
void RGB_to_HSV(float r, float g, float b, float & h, float & s, float & v );
void HSV_to_RGB( float h, float s, float v, float & r, float & g, float & b);

void XYZ_to_LAB (float X, float Y, float  Z, float & L, float &A, float & B);

void rgb_to_luv(fltarray &Data);
void luv_to_rgb(fltarray &Data);
void hsv_to_rgb(fltarray &Data);
void rgb_to_hsv(fltarray &Data);

inline void  RGB_to_YUV (float R, float G, float B, float & Y, float &  Cr, float & Cb)
{
   Y  =  (0.257 * R) + (0.504 * G) + (0.098 * B) + 16;
   Cr =  (0.439 * R) - (0.368 * G) - (0.071 * B) + 128;
   Cb = -(0.148 * R) - (0.291 * G) + (0.439 * B) + 128;
}

inline void YUV_to_RGB (float Y, float Cr, float Cb, float & R, float & G, float & B)
{
   B = 1.164*(Y - 16)                   + 2.018*(Cb - 128);
   G = 1.164*(Y - 16) - 0.813*(Cr - 128) - 0.391*(Cb - 128);
   R = 1.164*(Y - 16) + 1.596*(Cr - 128);
}

void col_rescale(fltarray &Data, float MinVal=0., float MaxVal=255.);

/***************************************************/

inline float saturation(float Val, float VatSatur=255.)
{
    if (Val < 0) return 0;
    else if (Val > VatSatur) return VatSatur;
    else return Val;
}
 
/***************************************************/
inline void rgb_to_yuv(fltarray &Data)
{
   int i,j;
   for (i = 0; i < Data.nx(); i++)
   for (j = 0; j < Data.ny(); j++)
   {
       float Y,  Cr, Cb;
       RGB_to_YUV(Data(i,j,0), Data(i,j,1),Data(i,j,2), Y,  Cr, Cb);
       Data(i,j,0) = Y; // saturation(Y, 255);
       Data(i,j,1) = Cr; // saturation(Cr, 255);
       Data(i,j,2) = Cb; // saturation(Cb, 255);
   }
}
 
inline void yuv_to_rgb(fltarray &Data)
{
   int i,j;
   for (i = 0; i < Data.nx(); i++)
   for (j = 0; j < Data.ny(); j++)
   {
       float R,G,B;
       YUV_to_RGB(Data(i,j,0), Data(i,j,1),Data(i,j,2), R,G,B);
       Data(i,j,0) = saturation(R, 255);
       Data(i,j,1) = saturation(G, 255);
       Data(i,j,2) = saturation(B, 255);
   }
}
 
/***************************************************/

#endif







