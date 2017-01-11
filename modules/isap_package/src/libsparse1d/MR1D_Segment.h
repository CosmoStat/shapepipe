/******************************************************************************
**                   Copyright (C) 1996 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 3.3
**
**    Author: Jean-Luc Starck
**
**    Date:  96/06/13 
**    
**    File:  MR1D_Segment.h
**
**    Modification history :
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

#ifndef _SEGMENT1D_H_
#define _SEGMENT1D_H_

#if _NOISE1D
#include "NoiseInterface.h"
#else
#include "MR1D_NoiseModel.h"
#endif


// structure (a structure is defined for each connected group 
// of pixels detected at a given scale)
// A object is contained several structures
typedef struct
{
   int NumObj;        // structure number
   int NumConnectObj; // object number
   int Scale;        // scale of the structure
   int NumScale;     // scale where the maximum of the object is found  
   int Depx;         // starting pixel of the structure
   int Maxx;         // position of the maximum of the structure
   float ValMax;     // maximum value of the structure
   int Sizex;        // size of the structure
   float Noise;      // noise at this position
   Bool SubObj;      // equal to True, if the structure is the 
                     // starting point of a sub-object
} obj1d_struct;

#define MAX_REC_ITER_OBJ1D 20

class MR1D_Seg {
   MR_1D Segment;          // Multiresolution object
                           // Segment(s,i) gives the object number
                           // of the wavelet coefficient (s,i)
                           // if Segment(s,i) = 0, nothing is detected
                           // at this scale
   obj1d_struct *TabStruct; // Array of structure  
                            // TabStruct[0] is not used
                            // TabStruct[1..NbrStruct] gives information
                            // about each detected structure
   int TabNbrStruct[MAX_SCALE_1D]; // number of structure per scale
   int Nbr_Plan;           // number of scale
   int Np;                 // number of pixels
   int NbrStruct;           // number of structures
   int NbrObj;              // number of objects: one object can contain
                            // several structure: NbrObj <= NbrStruct
   float *TabMax;           // TabMax[1..NbrObj] gives the maximum value
                            // of each object
   int *TabScaleMax;        // TabScaleMax[1..NbrObj] gives the scale
                            // where the maximum is found 
                            // scales are ranged between 0 and Nbr_Plan-2
   int *TabPosMax;          // TabPosMax[1..NbrObj] gives the position
                            // where the maximum is found
   int *TabNbrStructPerObj;  // TabNbrStructPerObj[1..NbrObj] gives the
                            // number of structure per object. 
   int **TabObjIndStruct;   // TabObjIndStruct[Num_Obj][1..NStruct_per_obj]
                            // gives the index in table TabNbrStruct.
   void set_obj_info();     // set the tables related to the objects
   void recons_adjoint_iter(MR_1D & MR_Data, int NumObj, 
                             fltarray & SignalRec, Bool Pos, Bool Neg);
   public:
      Bool Verbose;        // Print processing results
      int FirstScale;      // First scale to be used for the reconstruction
      int LastScale;       // Last scale to be used for the reconstruction
      Bool PosDetect;      // Reconstruct only positive structrure
      Bool NegDetect;      // Reconstruct only negative structrure
      int MaxRecIter;      // Maximum number of iterations
      MR1D_Seg();
      MR1D_Seg(MR_1D & MR_Data);
      void init(MR_1D & MR_Data);
      void print_info ();
      void print_obj_info ();
      int operator() (int s, int i) const {return (int) Segment(s,i);};
      int nbr_obj() { return NbrObj; };
      int nbr_struct() { return NbrStruct; };
      int nbr_struct_in_obj(int NumObj);
      int scale_of_obj(int NumObj);
      float max_of_obj(int NumObj);
      int posmax_of_obj(int NumObj);
      int firstpos_of_obj(int NumObj);
      int lastpos_of_obj(int NumObj);

      obj1d_struct & struct_in_obj(int NumObj, int NumStruct);

      void write_data_seg (char *FileName);

#if _NOISE1D
      void segment(MR_1D & MR_Data, INoiseModel1d& pro_NoiseModel);
#else
      void segment(MR_1D & MR_Data, MR1DNoiseModel & NoiseModel);
#endif
                  // MR_Data: input = multiresolution transform of the data
                  // MR1DNoiseModel: input = noise model
                  // set the field of the structure

      void recons(MR_1D & MR_Data, int NumObj, fltarray & SignalRec,
                  Bool Pos=False, Bool Neg=False);
                  // MR_Data: input = multiresolution transform of the data
                  // object number to be reconstruct
                  // SignalRec: input = noise model
                  // Pos: input: the object must be positive
                  // Neg: input: the object must be negative

      void recons_direct(MR_1D & MR_Data, int NumObj, fltarray & SignalRec,
                  Bool Pos=False, Bool Neg=False, Bool UseAdjointRec=False);
                  // MR_Data: input = multiresolution transform of the data
                  // object number to be reconstruct
                  // SignalRec: input = noise model
                  // Pos: input: the object must be positive
                  // Neg: input: the object must be negative

      void recons(MR_1D & MR_Data, fltarray & Tab_SignalRec);
                  // MR_Data: input = multiresolution transform of the data
                  // Tab_SignalRec: output: Tab_SignalRec(0..Np-1, obj)
                  //                is the reconstructed object of number obj-1
                  // obj = 0..NbrObj-1

      void recons(MR_1D & MR_Data, fltarray & Tab_SignalRec, fltarray & Tab_Add);
                  // MR_Data: input = multiresolution transform of the data
                  // Tab_SignalRec: output: Tab_SignalRec(0..Np-1, obj)
                  //                is the reconstructed object of number obj-1
                  // obj = 0..NbrObj-1
                  // Tab_Add : output: sum of all reconstructed objects
#if _NOISE1D
      void recons (MR_1D & MR_Data, fltarray & Tab_SignalRec, fltarray & Tab_Add, fltarray & Data, INoiseModel1d& pro_NoiseModel);
#else
      void recons (MR_1D & MR_Data, fltarray & Tab_SignalRec, fltarray & Tab_Add, fltarray & Data, MR1DNoiseModel &NoiseModel);
#endif
                  // MR_Data: input = multiresolution transform of the data
                  // Tab_SignalRec: output: Tab_SignalRec(0..Np-1, obj)
                  //                is the reconstructed object of number obj-1
                  // obj = 0..NbrObj-1
                  // Tab_Add : output: sum of all reconstructed objects
                  // Data: input data
                  // NoiseModel: input noise model
                  // This version of recons takes into account that data
                  // can have been transformed in case of Poisson noise
                  // or Gaussian + Poisson noise. It is not the case 
                  // in the previous

      void free();
      void clean(); // reset detection structures without deallocation
      ~MR1D_Seg();
};




#endif
