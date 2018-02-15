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
**    File:  IM_Segment.cc
**
*******************************************************************************
**
**    DESCRIPTION  image segmentation  
**    ----------- 
**                 
**  void im_segment (Ifloat &Imag, Ifloat &Segment, int &NLabel, float Level)
**
**
**  Imag : input image we want to segment
**  Segment: output segmented image
**  NLabel: Number of labels or regions
**  Level: segmentation level
**
*******************************************************************************
** 
** void im_histo (Ifloat &Pict, int **Tab_Histo, int &Nbr_Val)
** 
** computes the histogram of an image
** 
** Pict = image
** Tab_Histo = OUT: histogram
** Nbr_Val = OUT: number of values in the histogram = histogram size
**
*********************************************************************/

// static char sccsid[] = "@(#)IM_Segment.cc 3.1 96/05/02 CEA 1994 @(#)";


#include "IM_Obj.h"


#define		PILE_SIZE_BLOC		100000

class Pile {
	protected:
		int		bloc;
	public:
		int		num;
		int		*data;

		Pile (void);
		Pile (int);
		~Pile (void);
		void		add_data (int);
		void		alloue (void);
		void		trie (int);
		void		remplace (int, int);
		int		max (void);
	};


/*********************************************************************/
/*

void im_segment (Ifloat &Data, Ifloat &Segment, int &NLabel, float Level)
{
    int Nl = Data.nl();
    int Nc = Data.nc();
    Ifloat Imag(Nl, Nc, "Segment");
    int NbRealLabel;
    float Error;

    Imag = Data;
    NLabel = 0;
    NbRealLabel = 0;
    Segment.init();
    for (int i=0; i < Nl; i++)
    for (int j=0; j < Nc; j++)
    {
       if (Imag(i,j) > Level)
       {
          if ((i-1 >= 0) && (j-1>=0) && (Segment(i-1,j-1) != 0))
              Segment(i,j) = Segment(i-1,j-1);
          else if ((i-1 >= 0) &&  (Segment(i-1,j) != 0))
              Segment(i,j) = Segment(i-1,j);
          else if ((i-1 >= 0) && (j+1 < Nc) && (Segment(i-1,j+1) != 0))
              Segment(i,j) = Segment(i-1,j+1);
          else if ((j-1 >= 0) && (Segment(i,j-1) != 0))
              Segment(i,j) = Segment(i,j-1);
          Error = 0;
          if (Segment(i,j) != 0)
          {
             if ((i-1 >= 0) && (j-1>=0) && (Segment(i-1,j-1) != 0))
                 if (Segment(i,j) != Segment(i-1,j-1)) Error = Segment(i-1,j-1);
             else if ((i-1 >= 0) &&  (Segment(i-1,j) != 0))
                 if (Segment(i,j) != Segment(i-1,j)) Error = Segment(i-1,j);
             else if ((i-1 >= 0) && (j+1 < Nc) && (Segment(i-1,j+1) != 0))
                 if (Segment(i,j) != Segment(i-1,j+1)) Error = Segment(i-1,j+1);
             else if ((j-1 >= 0) && (Segment(i,j-1) != 0))
                 if (Segment(i,j) != Segment(i,j-1)) Error = Segment(i,j-1);

             if (Error > 0)
             {
                NbRealLabel --;
             }
             else
             {
                NLabel++;
                NbRealLabel++;
                Segment(i,j) = NLabel;
             }
          }
       }
    }
}

*/

/*********************************************************************/

void im_histo (Ifloat &Pict, int **Tab_Histo, int &Nbr_Val)
{
    int i,j,ind;
    float Min,Max;
    int Nl = Pict.nl();
    int Nc = Pict.nc();
 
    Min = min(Pict);
    Max = max (Pict);

    /* Calcul de l'entropie */
    Nbr_Val = (int) (Max - Min) + 1;

    *Tab_Histo = new int [Nbr_Val];    
    for (i = 0; i < Nbr_Val; i++) (*Tab_Histo) [i] = 0;
    
    for (i = 0; i < Nl; i++)
    for (j = 0; j < Nc; j++)
    {
        ind = iround(Pict(i,j) - Min);
        (*Tab_Histo) [ind] ++;
    }
 }

/*********************************************************************/

void	im_segment (Ifloat & Im, Ifloat &Segment, int &nb_seg, float S, 
        Bool CleanBord, int FirstLabel)
{
    int Nl = Im.nl();
    int Nc = Im.nc();
    int l,j,i,n;
    int Nvoisin = 4;
    int TabNeighbour[Nvoisin];
    Pile P;
    n = 1;
    P.add_data (0);

    if ((Segment.nl() != Nl) || (Segment.nc() != Nc)) Segment.resize(Nl,Nc);
    Segment.init();

    // image segmentation
    for (i=0; i<Nl; i++)
    for (j=0; j<Nc; j++)
    {
       if (ABS(Im(i,j)) > S)  
       {
          Bool NewLabel= True;
          int UseLabel = 0;

          TabNeighbour[0] = (int) Segment(i-1,j,I_ZERO);
          TabNeighbour[1] = (int) Segment(i,j-1,I_ZERO);
          TabNeighbour[2] = (int) Segment(i-1,j-1,I_ZERO);
          TabNeighbour[3] = (int) Segment(i-1,j+1,I_ZERO);

          // Find the label to apply to the pixel
          for (l=0; l < Nvoisin; l++) 
          { 
              if (TabNeighbour[l] != 0) 
              {
                 int L = TabNeighbour[l];
 		 int ColNeighbour = P.data[L];
                 NewLabel = False;
                 if (L > n)
                 {
                      cout << "Error: Segmentation problem ... " << endl;
                      cout << " n = " << n << " TabNeighbour[l] = " << L << endl;
                      exit(-1);
                 }
                 if (UseLabel == 0) UseLabel = ColNeighbour;
                 else if (UseLabel > ColNeighbour) UseLabel = ColNeighbour;
             }
          }

          // New label
          if (NewLabel == True)
	  {
             Segment(i,j) = n;
	     P.add_data (n++);
	  }
	  else 
          {
	      Segment(i,j) = UseLabel;
              // neighbourhood label must have the same label
              // and therefore may need to be changed
              for (l=0; l < Nvoisin; l++)
              {
                  if (TabNeighbour[l] != 0)
		  {
		     int L = TabNeighbour[l];
		     int ColNeighbour = P.data[L]; 
		     if (ColNeighbour != UseLabel)
		     {
                        if (L > n)
                        {
                        cout << "Error: Segmentation problem .... " << endl;
                        cout << " n = " << n << " TabNeighbour[l] = " << L << endl;
                        exit(-1);
                        }
                        P.remplace (P.data[L], UseLabel);
		     }
                  }
              }
           }
       }
   }

   if (CleanBord == True)
   {
       for (i=0 ; i< Nc; i++) 
       {
           if (Segment(0,i) != 0)    P.remplace((int) Segment(0,i),0);
	   if (Segment(Nl-1,i) != 0) P.remplace((int) Segment(Nl-1,i), 0);
       }
       for (i=0 ; i< Nl ; i++)
       {
 	  if (Segment(i,0) != 0)     P.remplace((int)Segment(i,0),0);
	  if (Segment(i,Nc-1) != 0)  P.remplace((int)Segment(i, Nc-1),0);
       }
    }  

    P.trie (n);
    int FL = FirstLabel-1;
    for (i=0 ; i < Nl; i++)
    for (j=0 ; j < Nc; j++)
    {
        int L = (int) Segment(i,j);
        if (L > n)
        {
            cout << "Error: Segmentation problem ..... " << endl;
            cout << " n = " << n << " TabNeighbour[l] = " << L << endl;
            exit(-1);
        }
	Segment(i,j) = (float)(P.data[L]);
	if (Segment(i,j) != 0) Segment(i,j) += FL;
    }
    nb_seg = P.max();
    // for (i=1; i <= n; i++)
    //  cout << "i: " << i << " Col = " << P.data[i] << endl;
}

//----------------------------------------------------------
//	Pile
//----------------------------------------------------------
Pile::Pile(void)
{
	data = NULL;
	num = 0;
	bloc = 0;
}

//----------------------------------------------------------
//	Pile
//----------------------------------------------------------
Pile::Pile (int n)
{
	data = NULL;
	bloc = n;
	if (data != NULL) data = (int *)realloc (data, bloc * sizeof (int));
	else data = (int *) malloc (bloc * sizeof (int));

	num = 0;
}

//----------------------------------------------------------
//	Pile :: destructeur
//----------------------------------------------------------
Pile::~Pile (void)
{
	if (data) free (data);
}

//----------------------------------------------------------
//	Pile :: max retourne le plus grand ellement
//----------------------------------------------------------
int	Pile::max (void)
{
	int	max = -100000;
	for (int i=0 ; i<num ; i++)
		max = (max > data[i]) ? max : data[i];
	return max;
}

//----------------------------------------------------------
//	Pile :: allocation de la memoire de data
//----------------------------------------------------------
void	Pile::alloue (void)
{
	bloc += PILE_SIZE_BLOC;
	if (data != NULL) data = (int *)realloc (data, bloc * sizeof (int));
	else data = (int *) malloc (bloc * sizeof (int));
}

//----------------------------------------------------------
//	Pile :: ajoute une donnee
//----------------------------------------------------------
void	Pile::add_data (int d)
{
	if (num >= bloc)
		alloue ();
	data[num++] = d;
}

//----------------------------------------------------------
//	Pile :: trie les donnes
//----------------------------------------------------------

void	Pile::trie (int max)
{
	int	i, c;
	Pile	histo (max);
	for (i=0 ; i<num ; i++) histo.data[i] = 0;
	for (i=0 ; i<num ; i++) histo.data[data[i]]++;
	c = 1;
	for (i=1 ; i<num ; i++)
		if ((histo.data[i] != 0) && (i != 0))
			remplace (i, c++);
}

//----------------------------------------------------------
//	Pile :: remplace d par f
//----------------------------------------------------------
void	Pile::remplace (int d, int f)
{
	for (int i=0 ; i<num ; i++)
		if (data[i] == d)
			data[i] = f;
}
