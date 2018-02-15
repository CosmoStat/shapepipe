
#ifndef _TEMPMEMORY_H
#define _TEMPMEMORY_H

#include "GlobalInc.h"
#include <omp.h>
 
template <class PARAM_TYPE> class TempCMem;
template <class PARAM_TYPE> class TempBuffMem;

extern char *alloc_buffer(size_t  Nelem);
extern void free_buffer(char *Ptr);
extern void memory_abort ();

#ifdef LARGE_BUFF
#define  VMS_DIR             "CEA_VM_DIR"
#define  VMS_SIZE            "CEA_VM_SIZE"

void vms_init(int UserSize, char * UserName, Bool Verbose);
#endif

#undef MEM_NOT_MANADGE
#define MEM_NOT_MANADGE 1 //FCS Modified

#undef DEBUG_MEM
#define DEBUG_MEM 0

#define MAX_IMA_IN_MEM 1500


//******************************************************************************
// external var 
//*****************************************************************************/
extern TempCMem<int> MemInt;
extern TempCMem<float> MemFloat;
extern TempCMem<double> MemDouble;
extern TempCMem<unsigned long> MemULong;//FCS ADDED
extern TempCMem<complex_f> MemCF;
extern TempCMem<complex_d> MemCD;
extern Bool UseVMS;



//*****************************************************************************/
// Template TempCMem class
//*****************************************************************************/
template <class PARAM_TYPE> 
class TempCMem {

private:
   TempBuffMem<PARAM_TYPE> TabBuffMem[MAX_IMA_IN_MEM];
   
public:
   TempCMem (){}
   ~TempCMem (){}
   PARAM_TYPE* alloc (int Nelem);
   void free (PARAM_TYPE *Ptr_Data);
};

 
//*****************************************************************************/
// Template TempBuffMem class
//*****************************************************************************/
template <class PARAM_TYPE> 
class TempBuffMem {

private:
   PARAM_TYPE* Ptr;
   int Size;
   Bool Use;
        
public:
   TempBuffMem () {Size = 0; Use = False;}
   PARAM_TYPE* alloc (int Nelem) ;
   void give_back () {Use = False;}
   void take () {Use = True;}
   int size () {return Size;}
   Bool use () {return Use;}
   PARAM_TYPE* buffer () { return Ptr;}
   ~TempBuffMem ();
};

//------------------------------------------------------------------------------
// class TempCMem<>::alloc ()
//------------------------------------------------------------------------------
template <class PARAM_TYPE> 
PARAM_TYPE* TempCMem<PARAM_TYPE>::alloc (int Nelem) {
     

   PARAM_TYPE* Pf=NULL;

	//FCS MODIFIED TO TAKE THE BUFFER IN THE WHILE LOOP AND BE CRITICAL FOR OPENMP

#if MEM_NOT_MANADGE    
	#ifdef _OPENMP
	#pragma omp critical ( mem_alloc )
	{
   	#endif        
   	
   	Pf = (PARAM_TYPE *) malloc (Nelem*sizeof(PARAM_TYPE));
  	if (Pf == NULL) memory_abort ();    	
  		
	#ifdef _OPENMP	
		#if  DEBUG_MEM
			int thr_num=omp_get_thread_num();
			cout << "Thr: " << thr_num << " Alloc Index Memory Ptr = " << Pf   <<" Nelem= "<<Nelem << endl;
		#endif
   	}
  	#endif
    return Pf;
#else

   int i=0;
   Bool Find=False;

   #ifdef _OPENMP
	#pragma omp critical ( mem_alloc )
	{
   #endif
   while (   (TabBuffMem[i].size() != 0) 
          && (!Find) && (i < MAX_IMA_IN_MEM)) {
      if ((TabBuffMem[i].use()) || (TabBuffMem[i].size() != Nelem)) i++;
      else {
      	Find = True;
   	TabBuffMem[i].take ();
	Pf = TabBuffMem[i].buffer();
	break;
      }
   } 
   if (i < MAX_IMA_IN_MEM) {
   	Pf = TabBuffMem[i].alloc (Nelem);
	#if DEBUG_MEM
		 #ifdef _OPENMP
			int thr_num=omp_get_thread_num();
			cout << "Thr: " << thr_num << " Alloc Index Memory = " << i << "  Ptr = " << Pf << "  Of size = " << TabBuffMem[i].size()    << endl;
		 #else
			cout << "Alloc Create: " << i << "  Ptr = " << Pf << endl;
	 	#endif
	#endif
   }   else {
     cerr << "Error: CMemInt cannot allocate memory ... " << endl;
      system ("pstat -s");
      i=0;
      while ((TabBuffMem[i].size() != 0) &&  (i < MAX_IMA_IN_MEM)) {
         cout << "Buffer " << i << " Size = " << TabBuffMem[i].size();
         if (TabBuffMem[i].use()) cout << " USE " << endl;
         else cout << " NOT USE " << endl;
         i++;
      }
      exit (0);
   }
   #ifdef _OPENMP	
	}
   #endif
        return Pf;
      

#endif

};


//------------------------------------------------------------------------------
//  class TempCMem<>::free (PARAM_TYPE *Ptr_Data)
//------------------------------------------------------------------------------
template <class PARAM_TYPE>
void TempCMem<PARAM_TYPE>::free (PARAM_TYPE *Ptr_Data) {

 #if MEM_NOT_MANADGE
	#ifdef _OPENMP
		#pragma omp critical ( mem_alloc )
   		{
   		#if  DEBUG_MEM
  		
			int thr_num=omp_get_thread_num();
			cout << "Thr: " << thr_num << "  Dealloc Index Memory Ptr = " << Ptr_Data  <<  endl;
		#endif
    #endif
   free_buffer((char *) Ptr_Data);
   
   #ifdef _OPENMP	
	}
   #endif
   
#else

   int i=0;
   Bool Find=False;

   #ifdef _OPENMP

   #pragma omp critical ( mem_alloc )
   {
    #endif
   while (   (TabBuffMem[i].size() != 0) 
          && (!Find) && (i < MAX_IMA_IN_MEM)) {
      if (Ptr_Data == TabBuffMem[i].buffer()) Find = True;
      else i ++;
   }
   #if DEBUG_MEM
		#ifdef _OPENMP
			int thr_num=omp_get_thread_num();
			cout << "Thr: " << thr_num << "  Dealloc Index Memory = " << i << "  Ptr = " << Ptr_Data << "  Of size = " << TabBuffMem[i].size() <<  endl;
		#else  
   			cout << "DeAlloc: " << i << "  Ptr = " << Ptr_Data << endl;
   		#endif
	#endif
	   if (Find) TabBuffMem[i].give_back();
	   else {  
    	  cerr << "Error: CMemInt cannot deallocate the memory ... "<< Ptr_Data << endl;
    	  exit (0);  
   		}
   	#ifdef _OPENMP	
	}
   #endif
   		
#endif
};


		 
//------------------------------------------------------------------------------
// class TempBuffMem<>::alloc (int Nelem)
//------------------------------------------------------------------------------
template <class PARAM_TYPE>		 
PARAM_TYPE* TempBuffMem<PARAM_TYPE>::alloc (int Nelem) {  
   
   Ptr = (PARAM_TYPE*) alloc_buffer((size_t) (Nelem * sizeof(PARAM_TYPE)));
   Size = Nelem; 
   Use = True;
   return Ptr;
};


//------------------------------------------------------------------------------
// class ~TempBuffMem ()
//------------------------------------------------------------------------------
template <class PARAM_TYPE>
TempBuffMem<PARAM_TYPE>::~TempBuffMem() {
   if (Size != 0) free_buffer ((char *) Ptr);
   Size = 0; 
   Use = False;
};	


//******************************************************************************
//  template free function 
//*****************************************************************************/

template <class PARAM_TYPE> 
inline PARAM_TYPE* temp_alloc(int Nelem, PARAM_TYPE& Dummy) {
   PARAM_TYPE* Ptr;
   Ptr = (PARAM_TYPE*) alloc_buffer((size_t) (Nelem*sizeof(PARAM_TYPE)));
   return Ptr;
};
inline float * f_alloc(int Nelem) {float Dummy;return temp_alloc(Nelem,Dummy);};
inline int * i_alloc(int Nelem) {int Dummy;return temp_alloc(Nelem,Dummy);};
inline unsigned int * ui_alloc(int Nelem) {unsigned int Dummy;return temp_alloc(Nelem,Dummy);};
inline short * s_alloc(int Nelem) {short Dummy; return temp_alloc(Nelem,Dummy);};
inline unsigned  short * us_alloc(int Nelem) {unsigned  short Dummy;return temp_alloc(Nelem,Dummy);};
inline char * c_alloc(int Nelem) {char Dummy; return temp_alloc(Nelem,Dummy);};
inline unsigned char * uc_alloc(int Nelem) {unsigned char Dummy; return temp_alloc(Nelem,Dummy);};

template <class PARAM_TYPE>
inline void temp_free(PARAM_TYPE* Ptr) {
   free_buffer((char *) Ptr);
};
inline void f_free(float *ptr) {temp_free(ptr);};
inline void i_free(int *ptr) {temp_free(ptr);};
inline void s_free(short *ptr) {temp_free(ptr);};
inline void c_free(char *ptr) {temp_free(ptr);};
inline void ui_free(unsigned int *ptr) {temp_free(ptr);};
inline void us_free(unsigned short *ptr) {temp_free(ptr);};
inline void uc_free(unsigned char *ptr) {temp_free(ptr);};


//******************************************************************************
// some alloc and free .... 
//*****************************************************************************/
template <class PARAM_TYPE>
inline PARAM_TYPE* vector_alloc(int Nelem, PARAM_TYPE& Dummy) {
   PARAM_TYPE* Vector;
   Vector = new PARAM_TYPE[Nelem];
   if (Vector == NULL) memory_abort();
   return Vector;
};
inline double *d_vector_alloc(int Nbr_Elem) {
   double Dummy;return vector_alloc(Nbr_Elem,Dummy);};
inline float *f_vector_alloc(int Nbr_Elem) {
   float Dummy;return vector_alloc(Nbr_Elem,Dummy);};
inline int *i_vector_alloc(int Nbr_Elem) {
   int Dummy;return vector_alloc(Nbr_Elem,Dummy);};
inline complex_f *cf_vector_alloc(int Nbr_Elem) {
   complex_f Dummy;return vector_alloc(Nbr_Elem,Dummy);};
  
   
template <class PARAM_TYPE>
inline void matrix_free(PARAM_TYPE **matrix, int nbr_lin) {
   for (int i=0; i<nbr_lin; i++)  delete [] matrix[i];
   delete [] matrix;
} 
inline void i_matrix_free(int **matrix, int nbr_lin) {matrix_free (matrix, nbr_lin);}
inline void f_matrix_free(float **matrix, int nbr_lin) {matrix_free (matrix, nbr_lin);}
inline void cf_matrix_free(complex_f **matrix, int nbr_lin) {matrix_free (matrix, nbr_lin);}

template <class PARAM_TYPE>
inline PARAM_TYPE** matrix_alloc(int nbr_lin, int nbr_col,PARAM_TYPE Dummy) {
   auto PARAM_TYPE** matrix;
   register int i;

   matrix = new  PARAM_TYPE* [nbr_lin];
   if (matrix == NULL) memory_abort();

   for (i=0; i<nbr_lin; i++) {
      matrix[i] = new PARAM_TYPE [nbr_col];
      if (matrix[i] == NULL) memory_abort();
   }
   return(matrix);
}
inline int** i_matrix_alloc(int nbr_lin, int nbr_col) {
   int Dummy=0; return matrix_alloc(nbr_lin,nbr_col,Dummy);
}
inline float** f_matrix_alloc(int nbr_lin, int nbr_col) {
   float Dummy=0; return matrix_alloc(nbr_lin,nbr_col,Dummy);
}
inline complex_f** cf_matrix_alloc(int nbr_lin, int nbr_col) {
   complex_f Dummy; return matrix_alloc(nbr_lin,nbr_col,Dummy);
}


/**********************************************************/
/**********************************************************/
/**********************************************************/

inline void MemMg_free (float* po_Buffer) {MemFloat.free (po_Buffer);}
inline float* MemMg_alloc (int Size, float Dummy) {return (MemFloat.alloc (Size));}

inline void MemMg_free (double* po_Buffer) {MemDouble.free (po_Buffer);}
inline double* MemMg_alloc (int Size, double Dummy) {return (MemDouble.alloc (Size));}

inline void MemMg_free (int* po_Buffer) {MemInt.free (po_Buffer);}
inline int* MemMg_alloc (int Size, int Dummy) {return (MemInt.alloc (Size));}

//FCS ADDED
inline void MemMg_free (unsigned long* po_Buffer) {MemULong.free (po_Buffer);}
inline unsigned long* MemMg_alloc (int Size, unsigned long Dummy) {return (MemULong.alloc (Size));}
//END FCS ADDED

inline void MemMg_free (complex_f* po_Buffer) {MemCF.free (po_Buffer);}
inline complex_f* MemMg_alloc (int Size, complex_f Dummy) {return (MemCF.alloc (Size));}

inline void MemMg_free (complex_d* po_Buffer) {MemCD.free (po_Buffer);}
inline complex_d* MemMg_alloc (int Size, complex_d Dummy) {return (MemCD.alloc (Size));}


	
#endif
