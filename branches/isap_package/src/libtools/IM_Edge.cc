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
**    File:  IM_Maxima.cc
**
*******************************************************************************
**
**    DESCRIPTION  edge detection  
**    ----------- 
**                 
**   
**
*********************************************************************/

// static char sccsid[] = "@(#)IM_Maxima.cc 3.1 96/05/02 CEA 1994 @(#)";


#include "IM_Obj.h"
#include "IM_Edge.h"
// #include "IM_IO.h"

/****************************************************************************/

static float PixDiffr[3][3] = { {0,0,0}, {0,1,-1}, {0,0,0}};
static float PixDiffc[3][3] = { {0,-1,0}, {0,1,0}, {0,0,0}};
#define NORM_PIXDIFF 1.

static float PixSepDiffr[3][3] = { {0,0,0}, {1,0,-1}, {0,0,0}};
static float PixSepDiffc[3][3] = { {0,-1,0}, {0,0,0}, {0,1,0}};
#define NORM_PIXSEPDIFF 1.

static float PixRobertsr[3][3] = { {0,0,-1}, {0,1,0}, {0,0,0}};
static float PixRobertsc[3][3] = { {-1,0,0}, {0,1,0}, {0,0,0}};
#define NORM_ROBERTS 1.

static float Sobelr[3][3] = { {1,0,-1}, {2,0,-2}, {1,0,-1}};
static float Sobelc[3][3] = { {-1,-2,-1}, {0,0,0}, {1,2,1}};
#define NORM_SOBEL 4.

static float Prewittr[3][3] = { {1,0,-1}, {1,0,-1}, {1,0,-1}};
static float Prewittc[3][3] = { {-1,-1,-1}, {0,0,0}, {1,1,1}};
#define NORM_PREWITT 1./3.

#define S2 sqrt(2.)
static float FreiChenr[3][3] = { {1,0,-1}, {S2,0,-S2}, {1,0,-1}};
static float FreiChenc[3][3] = { {-1,-S2,-1}, {0,0,0}, {1,S2,1}};
#define NORM_FREICHEN 1./(2+S2)

static float Laplacian1[3][3] = { {0,-1,0}, {-1,4,-1}, {0,-1,0}};
#define NORM_LAP1 1./4.

static float Laplacian2[3][3] = { {-1,-1,-1}, {-1,8,-1}, {-1,-1,-1}};
#define NORM_LAP2 1./8.

static float Laplacian3[3][3] = { {1,-2,1}, {-2,4,-2}, {1,-2,1}};
#define NORM_LAP3 1./8.


static float CompassPrewitt1[3][3] = { {1,1,-1}, {1,-2,-1}, {1,1,-1}};
static float CompassPrewitt2[3][3] = { {1,-1,-1}, {1,-2,-1}, {1,1,1}};
static float CompassPrewitt3[3][3] = { {-1,-1,-1}, {1,-2,1}, {1,1,1}};
static float CompassPrewitt4[3][3] = { {-1,-1,1}, {-1,-2,1}, {1,1,1}};
static float CompassPrewitt5[3][3] = { {-1,1,1}, {-1,-2,1}, {-1,1,1}};
static float CompassPrewitt6[3][3] = { {1,1,1}, {-1,-2,1}, {-1,-1,1}};
static float CompassPrewitt7[3][3] = { {1,1,1}, {1,-2,1}, {-1,-1,-1}};
static float CompassPrewitt8[3][3] = { {1,1,1}, {1,-2,-1}, {1,-1,-1}};
#define NORM_COMPASSPREWITT 1./5.

static float CompassKirsch1[3][3] = { {5,-3,-3}, {5,0,-3}, {5,-3,-3}};
static float CompassKirsch2[3][3] = { {-3,-3,-3}, {5,0,-3}, {5,5,-3}};
static float CompassKirsch3[3][3] = { {-3,-3,-3}, {-3,0,-3}, {5,5,5}};
static float CompassKirsch4[3][3] = { {-3,-3,-3}, {-3,0,5}, {-3,5,5}};
static float CompassKirsch5[3][3] = { {-3,-3,5}, {-3,0,5}, {-3,-3,5}};
static float CompassKirsch6[3][3] = { {-3,5,5}, {-3,0,5}, {-3,-3,-3}};
static float CompassKirsch7[3][3] = { {5,5,5}, {-3,0,-3}, {-3,-3,-3}};
static float CompassKirsch8[3][3] = { {5,5,-3}, {5,0,-3}, {-3,-3,-3}};
#define NORM_COMPASSKIRSCH 1./15.

static float CompassRobin3L1[3][3] = { {1,0,-1}, {1,0,-1}, {1,0,-1}};
static float CompassRobin3L2[3][3] = { {0,-1,-1}, {1,0,-1}, {1,1,0}};
static float CompassRobin3L3[3][3] = { {-1,-1,-1}, {0,0,0}, {1,1,1}};
static float CompassRobin3L4[3][3] = { {-1,-1,0}, {-1,0,1}, {0,1,1}};
static float CompassRobin3L5[3][3] = { {-1,0,1}, {-1,0,1}, {-1,0,1}};
static float CompassRobin3L6[3][3] = { {0,1,1}, {-1,0,1}, {-1,-1,0}};
static float CompassRobin3L7[3][3] = { {1,1,1}, {0,0,0}, {-1,-1,-1}};
static float CompassRobin3L8[3][3] = { {1,1,0}, {1,0,-1}, {0,-1,-1}};
#define NORM_COMPASSROBIN3L 1./3.

static float CompassRobin5L1[3][3] = { {1,0,-1}, {2,0,-2}, {1,0,-1}};
static float CompassRobin5L2[3][3] = { {0,-1,-2}, {1,0,-1}, {2,1,0}};
static float CompassRobin5L3[3][3] = { {-1,-2,-1}, {0,0,0}, {1,2,1}};
static float CompassRobin5L4[3][3] = { {-2,-1,0}, {-1,0,1}, {0,1,2}};
static float CompassRobin5L5[3][3] = { {-1,0,1}, {-2,0,2}, {-1,0,1}};
static float CompassRobin5L6[3][3] = { {0,1,2}, {-1,0,1}, {-2,-1,0}};
static float CompassRobin5L7[3][3] = { {1,2,1}, {0,0,0}, {-1,-2,-1}};
static float CompassRobin5L8[3][3] = { {2,1,0}, {1,0,-1}, {0,-1,-2}};
#define NORM_COMPASSROBIN5L 1./4.


/****************************************************************************/
 
void EDGE::template_gradient(const Ifloat &Im1, Ifloat &Mod, Ifloat &Phase)
{
    int f,i,j,k,l;
    int ii,jj;
    float Modulus;
    float Max=0.;
    int Angle;
    
    for (i=0; i < Im1.nl(); i++)
    for (j=0; j < Im1.nc(); j++)
    {
      ii = i +  SizeFilter  / 2;
      jj = j +  SizeFilter  / 2;
      Max = 0.;
      Angle = 0;
      for (f=0; f < NbrFilter; f++)
      {
         Modulus  = 0;
         for (k = 0; k < SizeFilter; k++)
         for (l = 0; l < SizeFilter; l++)
            Modulus  += Im1(ii - k, jj - l, Bord) * Filter(k, l, f);
         if (ABS(Modulus) > Max) 
         {
            Max =  ABS(Modulus);
	    Angle = 1 + f;
         }
       }
       Mod(i,j) = Max / NormFilter;
       Phase(i,j) = Angle;
    }
}

/****************************************************************************/

void EDGE::LoG(const Ifloat &Data, Ifloat &Lap, float Scale)
{
   int i,j,k,l;
   int Nl = Data.nl();
   int Nc = Data.nc();
   float x,s1 = Scale;
   float s2 = s1 * 1.6;
   float Val;
   Ifloat G1x(Nl,Nc, "buff");
   Ifloat G2x(Nl,Nc, "buff");
   
   cout << "s1 = " << s1 << " s2 = " << s2 << endl;
   
   // convolution following x-direction
    for (i = 0; i < Nl; i++)
    for (j = 0; j < Nc; j++)
    {
        k = (int) (5. * s1);
        G1x (i,j) = 0.;
        for (l = j-k; l < j+k; l++)
        {
            x = (float)(j - l)/ s1;
            G1x (i,j) +=  exp(-x*x * .5) * Data(i,l,Bord);
        }
 	G1x (i,j) /=  s1;
	
	k = (int) (5. * s2);
        G2x (i,j) = 0.;
        for (l = j-k; l < j+k; l++)
        {
            x = (float)(j - l)/ s2;
            G2x (i,j) +=  exp(-x*x * .5) * Data(i,l,Bord);
        }
 	G2x (i,j) /= s2;
    }
    //INFO(G1x, "G1x");    
    //io_write_ima_float("lap1x.fits", G1x);

    //INFO(G2x, "G2x");
    //io_write_ima_float("lap1y.fits", G2x);
    
    // convolution following y-direction
    for (i = 0; i < Nl; i++)
    for (j = 0; j < Nc; j++)
    {
        k = (int) (5. * s1);
        Lap(i,j) = 0.;
        for (l = i-k; l < i+k; l++)
        {
            x = (float)(i - l)/ s1;
            Lap(i,j) +=  exp(-x*x * .5) * G1x(l,j,Bord);
        }
 	Lap(i,j) /= s1;
	
	k = (int) (5. * s2);
        Val = 0.;
        for (l = i-k; l < i+k; l++)
        {
           x = (float)(i - l)/ s2;
           Val += exp(-x*x*.5) *  G2x(l,j,Bord);
        }
 	Lap(i,j) -= Val / s2;
    }
    //INFO(Lap, "Lap");
    //io_write_ima_float("lap.fits",Lap);
}

/****************************************************************************/


void EDGE::canny(const Ifloat &Data, Ifloat &Gradx, Ifloat &Grady, float Scale)
{
   int i,j,k,l;
   int Nl = Data.nl();
   int Nc = Data.nc();
   float x,Scale2;
   Ifloat Gy(Nl,Nc, "buff");
   Ifloat Gx(Nl,Nc, "buff");
   
    k = (int) (5. * Scale);
    for (i = 0; i < Nl; i++)
    for (j = 0; j < Nc; j++)
    {
       Gy(i,j) = 0.;
       for (l = i-k; l < i+k; l++)
       {
           x = (float)(i - l)/ Scale;
           Gy(i,j) +=  exp(-x*x * .5) *  Data(l,j,Bord);
        }
        Gx (i,j) = 0.;
        for (l = j-k; l < j+k; l++)
        {
            x = (float)(j - l)/ Scale;
            Gx (i,j) +=  exp(-x*x * .5) * Data(i,l,Bord);
        }
      }
    
    Scale2 = Scale*Scale;
    for (i = 0; i < Nl; i++)
    for (j = 0; j < Nc; j++)
    {
        Gradx (i,j) = 0.;
        for (l = j-k; l < j+k; l++)
        {
            x = (float)(j - l)/ Scale;
            Gradx (i,j) -=  x * exp(-x*x * .5) * Gy(i,l,Bord);
        }
 	Gradx (i,j) /= Scale;
	
        Grady(i,j) = 0.;
        for (l = i-k; l < i+k; l++)
        {
           x = (float)(i - l)/ Scale;
           Grady(i,j) -=  x * exp(-x*x * .5) *  Data(l,j,Bord);
        }
  	Grady(i,j) /= Scale;
    }     
}

/****************************************************************************/

void EDGE::pix_dir_angle(int Angle, int i, int j, 
                   int & i1, int & j1, int & i2, int & j2, int DistPix)
/*
   return the coordinated of the two points related to a pixel i,j
   in a direction given in Angle
   Tcont = in: Contour type (value between 1 and 8)
   i,j = in: pixel position
   i1,j1 = out: position of the first pixel in the contour
   i2,j2 = out: position of the second pixel in the contour
   DistPix = in: distance between two consecutive pixel
*/
{
   int Di1, Di2, Dj1, Dj2;
   Di1 =  Dj1 =  Di2 = Dj2 = 0;
   switch(Angle)
   {
        case 0: break;
        case 3: 
	case 7: Di1-=1; Di2+=1; break;
	case 4:
	case 8: Di1+=1; Dj1-=1; Di2-=1; Dj2+=1; break;
	case 1:
	case 5: Dj1-=1; Dj2+=1; break;
	case 2:
	case 6: Di1-=1; Dj1-=1; Di2+=1; Dj2+=1; break;
        default:
	   cerr << "Error: bad orientation ... " << endl;
	   exit(-1);
	   break;  
   }  
   i1 = i + Di1 * DistPix;
   i2 = i + Di2 * DistPix;
   j1 = j + Dj1 * DistPix;
   j2 = j + Dj2 * DistPix; 
}   

/*********************************************************************/  
  
void EDGE::edge_from_1deriv(Ifloat & Mod, Ifloat & Phase, Bool PhaseInt)
// threshold Mod(i,j) if Mod(i,j) is not a maximum in the gradient direction
{
   int Nl = Mod.nl();
   int Nc = Mod.nc();
   int i,j,i1,j1,i2,j2;
   int K=1, Angle;
   Ifloat MM(Nl,Nc,"mod");
   
   MM = Mod;
   for (i=0;i< Nl;i++)
   for (j=0;j< Nc;j++)
   {
      // the modulus value must be a maximum in the 
      // gradient direction  
      if (PhaseInt == False) Angle = angle_dir(Phase(i,j));
      else Angle = (int) Phase(i,j);
      pix_dir_angle(Angle, i, j,i1, j1, i2, j2, K);
       
      if ((ABS(Mod(i1,j1,Bord)) >=  ABS(Mod(i,j))) 
           || (ABS(Mod(i2,j2,Bord))  >=  ABS(Mod(i,j)))) 
      {
           MM(i,j) = 0;
	   Phase(i,j) = 0.;
      }
      
//       else 
//       {
//           Phase(i,j) = Phase(i,j) / PI * 180.;
// 	  if (Phase(i,j) < 0) Phase(i,j) += 360.;
//       }
   }
   Mod = MM;
   
}
/*********************************************************************/

void EDGE::kill_isol(Ifloat &Edge)
{
   int i,j;
   int Nl = Edge.nl();
   int Nc = Edge.nc();  
   
   for (i=0;i < Nl; i++)
   for (j=0;j < Nc; j++) 
   {
        // supress isolated pixel
        if (isolated_pixel(Edge,i,j,Bord) == True)  Edge(i,j) = 0.;       
   }
}
   
/*********************************************************************/

void EDGE::edge_from_2deriv(Ifloat & Mod)
// threshold all values where no zero crossing are detected
{
   int Nl = Mod.nl();
   int Nc = Mod.nc();  
   int i,j;
   
    for (i = 0; i < Nl; i++)
    for (j = 0; j < Nc; j++)
    {
        // test zero crossing
        if ((  Mod(i,j) * Mod(i+1, j, Bord) > 0) &&
            (  Mod(i,j) * Mod(i+1, j+1, Bord) > 0) &&
            (  Mod(i,j) * Mod(i, j+1, Bord) > 0))  Mod(i,j) = 0.;
    }
}

/*********************************************************************/

                         
void EDGE::detection(Ifloat &Data, Ifloat & Gradx)
{
   int Nl = Data.nl();
   int Nc = Data.nc();
   Ifloat Grady(Nl,Nc,"Gx");
   int i,j;
   Ifloat Gx(SizeFilter, SizeFilter, "kernel");
   Ifloat Gy(SizeFilter, SizeFilter, "kernel");
   float Mod,Phase;
   extern void convol(const Ifloat &Im1, const Ifloat &Im2, Ifloat &Result,
            type_border Type_Border);

   switch (TypeGrad)
   {
      case FIRST_DERIV:
         if (Method == ED_CANNY) canny(Data,  Gradx,  Grady,  Sigma);
	 else
	 {
            for (i =0; i < SizeFilter; i++)
            for (j =0; j < SizeFilter; j++)
            {
               Gx(i,j) = Filter(i,j,0) / NormFilter;
               Gy(i,j) = Filter(i,j,1) / NormFilter;
            }
            convol(Data, Gx, Gradx, Bord);
            convol(Data, Gy, Grady, Bord);
	 }
 	 for (i =0; i < Nl*Nc; i++) 
         {
           Mod = sqrt(Gradx(i)* Gradx(i) + Grady(i)*Grady(i));
           // compute the angle corresponding to edge
           ARG (Gradx(i), Grady(i), Phase);  
           if (Method == ED_ROBERTS) Phase += PI / 4.;
           Gradx(i) = Mod;
           Grady(i) = (int) angle_dir(Phase);
	   if (Mod < FLOAT_EPSILON) Grady(i) = 0.;
         }
	 if (SelectMaxima == True) edge_from_1deriv(Gradx, Grady, True);
	 break;
      case SECOND_DERIV:
         if (Method == ED_LOG) LoG(Data,  Gradx , Sigma);
         else
	 {
	    for (i =0; i < SizeFilter; i++)
            for (j =0; j < SizeFilter; j++)
              Gx(i,j) = Filter(i,j,0) / NormFilter;
            convol(Data, Gx, Gradx, Bord);
	 }
         if (SelectMaxima == True) edge_from_2deriv(Gradx);
	 break;
      case TEMPLATE_GRAD:
         template_gradient(Data, Gradx, Grady);
         if (SelectMaxima == True) edge_from_1deriv(Gradx, Grady,True);
         break;
    }
}
     
/*********************************************************************/

void EDGE::init_param(type_edge Edge)
{
   int i,j;
   Method = Edge;
   SizeFilter= 3;
   Sigma = 0.;
   Bord = I_MIRROR;
   SelectMaxima = True;
   
   switch(Method)
   {
          case ED_PIXDIFF:  
	      NbrFilter = 2;
	      Filter.alloc(SizeFilter, SizeFilter, NbrFilter);
	      for (i=0; i < SizeFilter; i++)
	      for (j=0; j < SizeFilter; j++)
	      {
	        Filter(i,j,0) = PixDiffr[i][j];
		Filter(i,j,1) = PixDiffc[i][j];
	      }
	      NormFilter =  NORM_PIXDIFF;
	      TypeGrad = FIRST_DERIV;
	      break;
           case ED_SEPPIXDIFF:
	      NbrFilter = 2;
	      Filter.alloc(SizeFilter, SizeFilter, NbrFilter);
	      for (i=0; i < SizeFilter; i++)
	      for (j=0; j < SizeFilter; j++)
	      {
	        Filter(i,j,0) =  PixSepDiffr[i][j];
		Filter(i,j,1) =  PixSepDiffc[i][j];
	      }
	      NormFilter = NORM_PIXSEPDIFF;
	      TypeGrad = FIRST_DERIV;
 	      break;
          case ED_SOBEL: 
 	      NbrFilter = 2;
	      Filter.alloc(SizeFilter, SizeFilter, NbrFilter);
	      for (i=0; i < SizeFilter; i++)
	      for (j=0; j < SizeFilter; j++)
	      {
	        Filter(i,j,0) =  Sobelr[i][j];
		Filter(i,j,1) =  Sobelc[i][j];
	      }
	      NormFilter =  NORM_SOBEL;
	      TypeGrad = FIRST_DERIV;
 	      break;
          case ED_PREWITT: 
  	      NbrFilter = 2;
	      Filter.alloc(SizeFilter, SizeFilter, NbrFilter);
	      for (i=0; i < SizeFilter; i++)
	      for (j=0; j < SizeFilter; j++)
	      {
	        Filter(i,j,0) = Prewittr[i][j];
		Filter(i,j,1) = Prewittc[i][j];
	      }
	      NormFilter = NORM_PREWITT;
	      TypeGrad = FIRST_DERIV;
	      break;        
          case ED_ROBERTS:
	      NbrFilter = 2;
	      Filter.alloc(SizeFilter, SizeFilter, NbrFilter);
	      for (i=0; i < SizeFilter; i++)
	      for (j=0; j < SizeFilter; j++)
	      {
	        Filter(i,j,0) =  PixRobertsr[i][j];
		Filter(i,j,1) =  PixRobertsc[i][j];
	      }
	      NormFilter =  NORM_ROBERTS;
	      TypeGrad = FIRST_DERIV;
	      break;
           case ED_FREI_CHEN: 
 	      NbrFilter = 2;
	      Filter.alloc(SizeFilter, SizeFilter, NbrFilter);
	      for (i=0; i < SizeFilter; i++)
	      for (j=0; j < SizeFilter; j++)
	      {
	        Filter(i,j,0) =   FreiChenr[i][j];
		Filter(i,j,1) =   FreiChenc[i][j];
	      }
	      NormFilter =  NORM_FREICHEN;
	      TypeGrad = FIRST_DERIV;
	      break;     	     
          case  ED_LAP1: 
	     NbrFilter = 1;
	      Filter.alloc(SizeFilter, SizeFilter, NbrFilter);
	      for (i=0; i < SizeFilter; i++)
	      for (j=0; j < SizeFilter; j++)
	      {
	        Filter(i,j,0) =   Laplacian1 [i][j];
 	      }
	      NormFilter =  NORM_LAP1 ;
 	      TypeGrad = SECOND_DERIV;
	      break;
          case  ED_LAP2: 
	     NbrFilter = 1;
	      Filter.alloc(SizeFilter, SizeFilter, NbrFilter);
	      for (i=0; i < SizeFilter; i++)
	      for (j=0; j < SizeFilter; j++)
	      {
	        Filter(i,j,0) =   Laplacian2 [i][j];
 	      }
	      NormFilter =  NORM_LAP2 ;
	      TypeGrad  =  SECOND_DERIV;
	      break;
          case  ED_LAP3: 
	      NbrFilter = 1;
	      Filter.alloc(SizeFilter, SizeFilter, NbrFilter);
	      for (i=0; i < SizeFilter; i++)
	      for (j=0; j < SizeFilter; j++)
	      {
	        Filter(i,j,0) =   Laplacian3 [i][j];
 	      }
	      NormFilter =  NORM_LAP3 ;
 	      TypeGrad  =  SECOND_DERIV;
	      break;
          case ED_COMPASS_PREWITT:
	      NbrFilter = 8;
	      Filter.alloc(SizeFilter, SizeFilter, NbrFilter);
	      for (i=0; i < SizeFilter; i++)
	      for (j=0; j < SizeFilter; j++)
	      {
	        Filter(i,j,0) = CompassPrewitt1[i][j];
		Filter(i,j,1) = CompassPrewitt2[i][j];
		Filter(i,j,2) = CompassPrewitt3[i][j];
		Filter(i,j,3) = CompassPrewitt4[i][j];
		Filter(i,j,4) = CompassPrewitt5[i][j];
		Filter(i,j,5) = CompassPrewitt6[i][j];
		Filter(i,j,6) = CompassPrewitt7[i][j];
		Filter(i,j,7) = CompassPrewitt8[i][j];
 	      }
	      NormFilter = NORM_COMPASSPREWITT;
	      TypeGrad  =  TEMPLATE_GRAD;
	      break; 
          case ED_COMPASS_KIRSCH: 
	      NbrFilter = 8;
	      Filter.alloc(SizeFilter, SizeFilter, NbrFilter);
	      for (i=0; i < SizeFilter; i++)
	      for (j=0; j < SizeFilter; j++)
	      {
	        Filter(i,j,0) = CompassKirsch1[i][j];
		Filter(i,j,1) = CompassKirsch2[i][j];
		Filter(i,j,2) = CompassKirsch3[i][j];
		Filter(i,j,3) = CompassKirsch4[i][j];
		Filter(i,j,4) = CompassKirsch5[i][j];
		Filter(i,j,5) = CompassKirsch6[i][j];
		Filter(i,j,6) = CompassKirsch7[i][j];
		Filter(i,j,7) = CompassKirsch8[i][j];
 	      }
	      NormFilter = NORM_COMPASSKIRSCH;
	      TypeGrad  =  TEMPLATE_GRAD;
	      break; 
          case ED_COMPASS_ROBIN3L: 
	      NbrFilter = 8;
	      Filter.alloc(SizeFilter, SizeFilter, NbrFilter);
	      for (i=0; i < SizeFilter; i++)
	      for (j=0; j < SizeFilter; j++)
	      {
	        Filter(i,j,0) = CompassRobin3L1[i][j];
		Filter(i,j,1) = CompassRobin3L2[i][j];
		Filter(i,j,2) = CompassRobin3L3[i][j];
		Filter(i,j,3) = CompassRobin3L4[i][j];
		Filter(i,j,4) = CompassRobin3L5[i][j];
		Filter(i,j,5) = CompassRobin3L6[i][j];
		Filter(i,j,6) = CompassRobin3L7[i][j];
		Filter(i,j,7) = CompassRobin3L8[i][j];
 	      }
	      NormFilter = NORM_COMPASSROBIN3L;
	      TypeGrad  =  TEMPLATE_GRAD;
	      break; 
          case ED_COMPASS_ROBIN5L:        
	      NbrFilter = 8;
	      Filter.alloc(SizeFilter, SizeFilter, NbrFilter);
	      for (i=0; i < SizeFilter; i++)
	      for (j=0; j < SizeFilter; j++)
	      {
	        Filter(i,j,0) = CompassRobin5L1[i][j];
		Filter(i,j,1) = CompassRobin5L2[i][j];
		Filter(i,j,2) = CompassRobin5L3[i][j];
		Filter(i,j,3) = CompassRobin5L4[i][j];
		Filter(i,j,4) = CompassRobin5L5[i][j];
		Filter(i,j,5) = CompassRobin5L6[i][j];
		Filter(i,j,6) = CompassRobin5L7[i][j];
		Filter(i,j,7) = CompassRobin5L8[i][j];
 	      }
	      NormFilter = NORM_COMPASSROBIN5L;
	      TypeGrad  =  TEMPLATE_GRAD;
	      break; 
	  case ED_LOG:
	      Sigma = 1. /sqrt(3.);
	      NbrFilter = 0;
	      TypeGrad  = SECOND_DERIV;
	      break;
	  case ED_CANNY:
	      Sigma = 1. /sqrt(3.);
	      NbrFilter = 0;
	      TypeGrad  = FIRST_DERIV;
	      break;
    }
}
/**********************************************************/
/**********************************************************/

void im_threshold_kirsh_edge(Ifloat &Edge, Iint &Angle, float Noise, float Nsigma, int Step)
{      
   int Nl = Edge.nl();
   int Nc = Edge.nc();
   int i,j,i1,j1,i2,j2;
   Ifloat ImaEdge(Nl,Nc,"ImaEdge");
   int K = POW2(Step);
   float Level = Noise*Nsigma/sqrt(3.);
   type_border Border=I_MIRROR;

   ImaEdge = Edge;
   for (i=0;i<Nl;i++)
   for (j=0;j<Nc;j++)
   {
      Edge(i,j) = ImaEdge(i,j);
      if (ABS(ImaEdge(i,j)) > Level)
      {
         i1 = i2 = i;
         j1 = j2 = j;
         switch(Angle(i,j))
         {
            case 1: j1-=K; j2+=K; break;
            case 2: i1-=K; i2+=K; break;
            case 3: i1-=K; j1-=K; i2+=K; j2+=K; break;
            case 4: i1-=K; j1+=K; i2+=K; j2-=K;break;
            case 5: i1-=K; j2-=K;break;
            case 6: i1-=K; j2+=K;break;
            case 7: j1-=K; i2+=K; break;
            case 8: j1+=K; i2+=K; break;
            default: 
              cerr << "Error: bad orientation ... " << endl;
              exit(-10);
         }
         if (     (ABS(ImaEdge(i1,j1,Border)) < Level) 
               || (ABS(ImaEdge(i2,j2,Border)) < Level))
            Edge(i,j) = 0.;
      }
      else Edge(i,j) = 0.;
   } 
}

/**********************************************************/
/*
   Angle = 1 ==> 0 deg
         = 2 ==> 45 
         = 3 ==> 90
         = 4 ==> 135
         = 5 ==> 180
         = 6 ==> 225
         = 7 ==> 270
         = 8 ==> 315
*/
void im_edge_kirsh(Ifloat &Ima, Iint &Angle, Ifloat &ImaEdge, int Step)
{
   int Nl = Ima.nl();
   int Nc = Ima.nc();
   int i,j;
   fltarray Tab(8);
   int Ind_Angle, Ind_Angle_Min;
   type_border Border=I_MIRROR;
   int K = POW2(Step);
   float Min;
   cout << "K = " << K << endl;

   for (i=0;i<Nl;i++)
   for (j=0;j<Nc;j++)
   {
      Tab(0) = (Ima(i,j-K,Border) + Ima(i,j) + Ima(i,j+K,Border)) / 3.;
      Tab(1) = (Ima(i-K,j,Border) + Ima(i,j) + Ima(i+K,j,Border)) / 3.;
      Tab(2) = (Ima(i-K,j-K,Border) + Ima(i,j) + Ima(i+K,j+K,Border)) / 3.;
      Tab(3) = (Ima(i-K,j+K,Border) + Ima(i,j) + Ima(i+K,j-K,Border)) / 3.;
      Tab(4) = (Ima(i-K,j,Border) + Ima(i,j) + Ima(i,j-K,Border)) / 3.;
      Tab(5) = (Ima(i-K,j,Border) + Ima(i,j) + Ima(i,j+K,Border)) / 3.;
      Tab(6) = (Ima(i,j-K,Border) + Ima(i,j) + Ima(i+K,j,Border)) / 3.;
      Tab(7) = (Ima(i,j+K,Border) + Ima(i,j) + Ima(i+K,j,Border)) / 3.;
      //ImaEdge(i,j) = max(Tab, Ind_Angle);
      ImaEdge(i,j) = Tab.max(Ind_Angle);
      //Min = min (Tab, Ind_Angle_Min);
      Min = Tab.min (Ind_Angle_Min);
      if (ABS(Min) > ABS(ImaEdge(i,j)))
      {
         ImaEdge(i,j) = Min;
         Angle(i,j) = Ind_Angle_Min;
      }
      Angle(i,j) = Ind_Angle + 1;
   }
}

/*********************************************************************/
 

