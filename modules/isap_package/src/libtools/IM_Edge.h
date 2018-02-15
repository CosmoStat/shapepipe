/******************************************************************************
**                   Copyright (C) 1998 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  96/06/13 
**    
**    File:  IM_Edge.h
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

#ifndef __IM_EDGE__
#define __IM_EDGE__

#include "IM_Obj.h"

#define NBR_EDGE 15

enum type_edge {ED_CANNY, ED_PIXDIFF, ED_SEPPIXDIFF, ED_SOBEL, ED_PREWITT,
                ED_ROBERTS, ED_FREI_CHEN, ED_LAP1, ED_LAP2, ED_LAP3, ED_LOG,
		ED_COMPASS_PREWITT, ED_COMPASS_KIRSCH, ED_COMPASS_ROBIN3L,
		ED_COMPASS_ROBIN5L};
		
#define DEF_EDGE_METHOD ED_SOBEL


inline const char * StringEdge (type_edge  type)
{
    switch (type)
    {
    case ED_CANNY:
      return ("Canny (derivative of a Gaussian)");break;
    case ED_PIXDIFF: 
      return ("Pixel difference");break;
    case ED_SEPPIXDIFF: 
      return ("Separated pixel difference");break;
    case ED_SOBEL: 
      return ("Sobel");break;
    case ED_PREWITT: 
      return ("Prewitt");break;
    case ED_ROBERTS:
      return ("Roberts");break;
    case ED_FREI_CHEN: 
      return ("Frei Chen");break;
   case  ED_LAP1: 
      return ("Laplacian: filter 1");break;
   case  ED_LAP2: 
      return ("Laplacian: filter 2");break;
   case  ED_LAP3: 
      return ("Laplacian: filter 3");break;
   case ED_LOG:
      return ("Marr and Hildrith: Laplacian of Gaussian (LoG)");break;   
   case ED_COMPASS_PREWITT : 
      return ("Prewitt compass gradient");break;
   case ED_COMPASS_KIRSCH: 
      return ("Kirsch");break;
   case ED_COMPASS_ROBIN3L: 
      return ("Robinson 3-level");break;
   case ED_COMPASS_ROBIN5L: 
      return ("Robinson 5-level");break;
    default:
      return ("Error: bad edge detection method");break;
    }
}
inline void edge_usage()
{
    fprintf(OUTMAN, "         [-M edge_detection_method]\n");
    for (int i = 0; i < NBR_EDGE; i++)
    fprintf(OUTMAN, "              %d: %s \n",i+1,
                                            StringEdge ((type_edge)i));
    fprintf(OUTMAN, "             default is Sobel method\n");
}



enum type_gradient {FIRST_DERIV, SECOND_DERIV, TEMPLATE_GRAD};
inline int angle_dir (float Phase)
/*
   return the angle direction as a number between 1 and 8
   Phase must be expressed in radian.

   return angle direction     
    Angle = 1 ==> 0 deg = 0 rd
          = 2 ==> 45     = PI/4
          = 3 ==> 90     = PI /2
          = 4 ==> 135    = 3 PI /4
          = 5 ==> 180    = PI
          = 6 ==> 225    = 5 PI / 4
          = 7 ==> 270    = 3 PI /2
          = 8 ==> 315    = 7 PI / 4
*/
{
   int Angle;
   float Pha;
   
   // Pha = floatting angle between 1 and 8
   Pha =  ( Phase >= 0) ? Phase : 2.*PI + Phase;
   Pha = Pha / PI * 4. + 1;
      
   // Angle takes the closer integer value
   Angle = (int)(Pha + 0.5);
   if (Angle >= 9) Angle = 1;
   return Angle;
} 


class EDGE {
    type_edge Method; // Method for the edge detection
    int NbrFilter;    // number of filters used my this method
    int SizeFilter;   // size of the filters
    fltarray Filter;  // 3D array which contains all used filters
    float NormFilter; // Normalisation value
    type_gradient TypeGrad; // Type of gradient
    void init_param(type_edge Edge);
    

 
    void template_gradient(const Ifloat &Im1, Ifloat &Mod, Ifloat &Phase);
    // calculate the modulus and phase in case of template matching detection
    // methods. The phase is return as integer between 0 and 8
    
    void LoG(const Ifloat &Data, Ifloat &Lap, float Scale);
    // calculate the laplacian of a Gaussian
    
    void canny(const Ifloat &Data, Ifloat &Gradx, Ifloat &Grady, float Scale);
    // convolution with the derivative of a Gaussian
    
    public:
      float Sigma; // gaussian standard deviation: used only with 
                   // LoG and Canny methods
      type_border Bord; // border managment
      Bool SelectMaxima; // if SelectMaxima is set to True then,
                         // select zero-crossing  or 
                         // maxima in the gradient direction
      EDGE(type_edge Edge) 
      {
         init_param(Edge);
      }
      
      void detection(Ifloat &Data, Ifloat &Gradx);
      // edge detection
      // if SelectMaxima == True, select zero-crossing  or 
      // maxima in the gradient direction

      void edge_from_1deriv(Ifloat & Mod, Ifloat & Phase, Bool PhaseInt=False);
      // Mod and Phase are the modulus and the phase map of the gradient
      // for all (i,j)  threshold Mod(i,j) 
      //        if Mod(i,j) is not a maximum in the gradient direction
      // if PhaseInt == False then phase must be expressed in radian
      // else it must be expressed by an integer value in [1,8]
    
    void edge_from_2deriv(Ifloat & Mod);
    // threshold all values where no zero crossing are detected
          
      void kill_isol(Ifloat &Edge);
      // suppress isolated pixel edges
      
      void pix_dir_angle(int Angle, int i, int j, 
                   int & i1, int & j1, int & i2, int & j2, int DistPix=1);
     /*
       return the coordinated of the two points related to a pixel i,j
       in a direction given in Angle
       Angle  = in: Contour type (value between 1 and 8)
       i,j = in: pixel position
       i1,j1 = out: position of the first pixel in the contour
       i2,j2 = out: position of the second pixel in the contour
       DistPix = in: distance between two consecutive pixel
     */    
     
      void  pix_dir_angle(float Phase, int i, int j, 
                   int & i1, int & j1, int & i2, int &j2, int DistPix=1)
      {
         pix_dir_angle(angle_dir(Phase), i, j,i1, j1, i2, j2, DistPix);
      }
     ~EDGE(){};
};




#endif
