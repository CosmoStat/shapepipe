#ifndef _TEMPARRAY_H
#define _TEMPARRAY_H


#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>

#include "Border.h"

#define MAX_NBR_AXIS 3
#define TA_MIN_SIZE_FOR_MEM_ALLOC_CALL 50000

#undef _USEMEM
#define _USEMEM 0

#ifdef _USEMEM
#include "Memory.h"
#endif

#define CHECK_DIM 0

using std::string;
 
//******************************************************************************
// Template pratial specialization 
//*****************************************************************************/

class NewArray {
public:
   bool operator() () {return true;}
};

class Old2dArray {
public:
   bool operator() () {return false;}
};

#define fltarray to_array<float,true>
#define dblarray to_array<double,true>
#define intarray to_array<int,true>
#define bytearray to_array<byte,true>
#define cfarray to_array<complex_f,true>
#define cdarray to_array<complex_d,true>

#define Ifloat to_array<float,false>
#define Idouble to_array<double,false>
#define Iint to_array<int,false>
#define Icomplex_f to_array<complex_f,false>
#define Icomplex_d to_array<complex_d,false>
 


//******************************************************************************
// Template Array class
//*****************************************************************************/

template <class PARAM_TYPE, bool ARRAY_TYPE>
class to_array {

private:
  PARAM_TYPE*  po_Buffer;    //  vector of i_NbElem elt
  int          i_NbElem;     // number of elt
  int          i_NbAxis;      // number of axis
  int          pto_TabNaxis[MAX_NBR_AXIS];   // number of point per axis
  //char         tc_NameArray[SIZE_NAME];
  string       o_NameArray;
  bool         e_UseClassMemAlloc;
  bool         e_GetBuffer;
  bool        is_ima;
public:
  PARAM_TYPE  JunkVar;
  void set2ima() {if  (i_NbAxis == 2) is_ima=true;}
  void set2tab() {if  (i_NbAxis == 2) is_ima=false;}

  to_array();   
  to_array (int pi_Nx, const char *Name);
  to_array (int pi_Nx, int pi_Ny, const char *Name);
  to_array (int pi_Nx, int pi_Ny, int pi_Nz, const char *Name);
  to_array (int pi_Nx, int pi_Ny=0, int pi_Nz=0);
  ~to_array ();
  void free();
  PARAM_TYPE* buffer(); 
  PARAM_TYPE* const buffer() const;
  void init (const to_array<PARAM_TYPE,ARRAY_TYPE>& pro_Mat);
  void init ();
  void init (PARAM_TYPE Val);
  void alloc (int pi_Nx, const char* Name=0);  
  void alloc (int pi_Nx, int pi_Ny, const char* Name=0);  
  void alloc (int pi_Nx, int pi_Ny, int pi_Nz, const char* Name=0);
  void alloc (PARAM_TYPE *BuffData, int Nbr_Line, int Nbr_Col, const char *Name=0, bool MemManag=false);
  void alloc (PARAM_TYPE *BuffData, int Nx, int Ny, int Nz, const char *Name=0, bool MemManag=false);
  void reform (const int pi_Nx, const int pi_Ny=0, const int pi_Nz=0);
  void resize (const int pi_Nx, const int pi_Ny=0, const int pi_Nz=0);
  
  inline PARAM_TYPE& operator() (int x)  const;
  inline PARAM_TYPE  operator() (int x, type_border bord)  const;
  inline PARAM_TYPE& operator() (int x, int y) const;  
  inline PARAM_TYPE  operator() (int x, int y, type_border bord) const;    
  inline PARAM_TYPE& operator() (int x, int y, int z) const;
  inline PARAM_TYPE  operator() (int x, int y, int z, type_border bord) const;
  
  inline PARAM_TYPE&  setx (int x, type_border bord)  ;
  inline PARAM_TYPE&  setxy (int x, int y, type_border bord) ;
  inline PARAM_TYPE&  setxyz (int x, int y, int z, type_border bord) ;

  /*
  inline PARAM_TYPE& operator[] (int x)  const;
  inline PARAM_TYPE&  operator[]  (int x, type_border bord)  const;
  inline PARAM_TYPE& operator[]  (int x, int y) const;  
  inline PARAM_TYPE&  operator[]  (int x, int y, type_border bord) const;    
  inline PARAM_TYPE& operator[]  (int x, int y, int z) const;
  inline PARAM_TYPE & operator[]  (int x, int y, int z, type_border bord) const;
*/
  const to_array<PARAM_TYPE,ARRAY_TYPE>& operator = (const to_array<PARAM_TYPE,ARRAY_TYPE>& pro_Mat);
  const to_array<PARAM_TYPE,ARRAY_TYPE>& operator += (const to_array<PARAM_TYPE,ARRAY_TYPE>& pro_Mat);
  const to_array<PARAM_TYPE,ARRAY_TYPE>& operator *= (const to_array<PARAM_TYPE,ARRAY_TYPE>& pro_Mat);
  const to_array<PARAM_TYPE,ARRAY_TYPE>& operator *= (const double coef);
  const to_array<PARAM_TYPE,ARRAY_TYPE>& operator -= (const to_array<PARAM_TYPE,ARRAY_TYPE>& pro_Mat);
  const to_array<PARAM_TYPE,ARRAY_TYPE>& operator /= (const to_array<PARAM_TYPE,ARRAY_TYPE>& pro_Mat);
  const to_array<PARAM_TYPE,ARRAY_TYPE>& operator ^ (const double pf_coef);
  
  void info(string Name="");
  void display (int pi_NbElem=0);
  void rampgen ();
   
  void sup_threshold (float ThresholLevel);
  void inf_threshold (float ThresholLevel);
  
  int n_elem() const {return i_NbElem;} 
  int naxis() const { return i_NbAxis;}
  int axis(int pi_NumAxis) const { return pto_TabNaxis[pi_NumAxis-1];}
  int nx() const { return axis(1);}
  int ny() const { return axis(2);}
  int nz() const { return axis(3);} 
  //string get_name () const {return o_NameArray;}
  bool get_isbuffer() const { return e_GetBuffer;}
  void set_isbuffer(bool b) { e_GetBuffer = b;}
  bool get_buf() const { return e_GetBuffer;}
  bool get_memalloc() const { return e_UseClassMemAlloc;}
  
  int nc() const { return axis(1);}
  int nl() const { return axis(2);}
  
  PARAM_TYPE min ();
  PARAM_TYPE max (); 
  PARAM_TYPE maxfabs (); 
  PARAM_TYPE min (int& pri_ind);
  PARAM_TYPE max (int& pri_ind);
  PARAM_TYPE maxfabs (int& pri_ind);
  double total () const;
  double energy () const;
  double sigma () const;
  double mean () const;
  void sigma_clip (float& pf_Mean, float& pf_Sigma, int pi_Nit=3) const;
  float sigma_clip (int pi_Nit=3) const;
  
  to_array<PARAM_TYPE,ARRAY_TYPE> (const to_array<PARAM_TYPE,ARRAY_TYPE>& pro_Obj);

private:
  void set_attrib();
};




//------------------------------------------------------------------------------
// to_array ()
//------------------------------------------------------------------------------
template <class PARAM_TYPE, bool ARRAY_TYPE>
to_array<PARAM_TYPE,ARRAY_TYPE>::to_array () {
  set_attrib();
}

//------------------------------------------------------------------------------
// to_array (int pi_Nx, const char* Name)
//------------------------------------------------------------------------------  
template <class PARAM_TYPE, bool ARRAY_TYPE>
to_array<PARAM_TYPE,ARRAY_TYPE>::to_array (int pi_Nx, const char* Name) {
  set_attrib();
  alloc(pi_Nx, 0, 0, Name);
}

//------------------------------------------------------------------------------
// to_array (int pi_Nx, int pi_Ny, const char* Name)
//------------------------------------------------------------------------------  
template <class PARAM_TYPE, bool ARRAY_TYPE>
to_array<PARAM_TYPE,ARRAY_TYPE>::to_array (int pi_Nx, int pi_Ny, const char* Name) {
  set_attrib();  
  alloc(pi_Nx, pi_Ny, 0, Name);
}

//------------------------------------------------------------------------------
// to_array (int pi_Nx, int pi_Ny, int pi_Nz, const char* Name)
//------------------------------------------------------------------------------  
template <class PARAM_TYPE, bool ARRAY_TYPE>
to_array<PARAM_TYPE,ARRAY_TYPE>::to_array (int pi_Nx, int pi_Ny, int pi_Nz, const char* Name) {
  set_attrib();   
  alloc(pi_Nx, pi_Ny, pi_Nz, Name);
}

//------------------------------------------------------------------------------
// to_array (int pi_Nx, int pi_Ny, int pi_Nz)
//------------------------------------------------------------------------------  
template <class PARAM_TYPE, bool ARRAY_TYPE>
to_array<PARAM_TYPE,ARRAY_TYPE>::to_array (int pi_Nx, int pi_Ny, int pi_Nz) {
  set_attrib();   
  alloc(pi_Nx, pi_Ny, pi_Nz);
}

//------------------------------------------------------------------------------
// ~to_array ()
//------------------------------------------------------------------------------  
template <class PARAM_TYPE, bool ARRAY_TYPE>
to_array<PARAM_TYPE,ARRAY_TYPE>::~to_array () {free();}

//------------------------------------------------------------------------------
// free ()
//------------------------------------------------------------------------------
template <class PARAM_TYPE, bool ARRAY_TYPE>
void to_array<PARAM_TYPE,ARRAY_TYPE>::free() {
  if (e_UseClassMemAlloc == true){
#ifdef _USEMEM
    MemMg_free (po_Buffer);
#endif
  } else {
    if (i_NbElem != 0 && e_GetBuffer == false) delete[] po_Buffer;
  } 
  i_NbElem=0;i_NbAxis=0;o_NameArray="";//tc_NameArray[0]='\0';
  e_UseClassMemAlloc=false;
  for (int i=0;i<MAX_NBR_AXIS;i++) pto_TabNaxis[i]=0;
}

//------------------------------------------------------------------------------
// buffer ()
//------------------------------------------------------------------------------
template <class PARAM_TYPE, bool ARRAY_TYPE>
PARAM_TYPE* to_array<PARAM_TYPE,ARRAY_TYPE>::buffer() {
   return po_Buffer;
}

//------------------------------------------------------------------------------
// buffer ()
//------------------------------------------------------------------------------
template <class PARAM_TYPE, bool ARRAY_TYPE>
PARAM_TYPE* const to_array<PARAM_TYPE,ARRAY_TYPE>::buffer() const {
   return po_Buffer;
}

//------------------------------------------------------------------------------
// init ()
//------------------------------------------------------------------------------
template <class PARAM_TYPE, bool ARRAY_TYPE>
void to_array<PARAM_TYPE,ARRAY_TYPE>::init (const to_array<PARAM_TYPE,ARRAY_TYPE>& pro_Mat) { 
   if (n_elem() != 0) free();
   if (ARRAY_TYPE==true) {
      alloc(pro_Mat.nx(), pro_Mat.ny(), pro_Mat.nz());
   } else {
      alloc(pro_Mat.ny(), pro_Mat.nx(), pro_Mat.nz());
   }
}

//------------------------------------------------------------------------------
// init (PARAM_TYPE Val=0)
//------------------------------------------------------------------------------
template <class PARAM_TYPE, bool ARRAY_TYPE>
// !!!!! convert function from int to PARAM_TYPE must exist
void to_array<PARAM_TYPE,ARRAY_TYPE>::init (PARAM_TYPE Val) { 
   for (int i=0;i<n_elem();i++) po_Buffer[i] = Val;
}

template <class PARAM_TYPE, bool ARRAY_TYPE>
// !!!!! convert function from int to PARAM_TYPE must exist
void to_array<PARAM_TYPE,ARRAY_TYPE>::init () 
{ 
   PARAM_TYPE Val=0;
   for (int i=0;i<n_elem();i++) po_Buffer[i] = Val;
}
//------------------------------------------------------------------------------
// alloc (int pi_Nx, const char* Name)
//------------------------------------------------------------------------------
template <class PARAM_TYPE, bool ARRAY_TYPE>
void to_array<PARAM_TYPE,ARRAY_TYPE>::alloc (int pi_Nx, const char* Name) {
  alloc (pi_Nx, 0, 0, Name);
}
 
//------------------------------------------------------------------------------
// alloc (int pi_Nx, int pi_Ny, const char* Name)
//------------------------------------------------------------------------------
template <class PARAM_TYPE, bool ARRAY_TYPE>
void to_array<PARAM_TYPE,ARRAY_TYPE>::alloc (int pi_Nx, int pi_Ny, const char* Name) {
  alloc (pi_Nx, pi_Ny, 0, Name);
}
 
//------------------------------------------------------------------------------
// alloc (int pi_Nx, int pi_Ny, int pi_Nz, const char* Name)
//------------------------------------------------------------------------------
template <class PARAM_TYPE, bool ARRAY_TYPE>
void to_array<PARAM_TYPE,ARRAY_TYPE>::alloc (int pi_Nx, int pi_Ny, int pi_Nz, const char* Name) {

  if (i_NbElem != 0) free();
  
  if (pi_Nz != 0) i_NbElem = pi_Nz*pi_Ny*pi_Nx;
  else if (pi_Ny != 0) i_NbElem = pi_Ny*pi_Nx;
  else i_NbElem = pi_Nx;
  
  if (i_NbElem > TA_MIN_SIZE_FOR_MEM_ALLOC_CALL) {
#ifdef _USEMEM
    PARAM_TYPE Dummy=0;
    po_Buffer = MemMg_alloc (i_NbElem,Dummy);
    e_UseClassMemAlloc = true;
#else    
    e_UseClassMemAlloc = false; 	
    po_Buffer = new PARAM_TYPE [i_NbElem];
    if (po_Buffer == 0) cout << " Not enought memory " << endl;
#endif
  } else if (i_NbElem != 0) {
    e_UseClassMemAlloc = false;
    po_Buffer = new PARAM_TYPE [i_NbElem];
    if (po_Buffer == 0) cout << " Not enought memory " << endl;
  } else {
    e_UseClassMemAlloc = false;
    po_Buffer = (PARAM_TYPE*)NULL; 
    e_GetBuffer=false;
  }
  e_GetBuffer = false;
  pto_TabNaxis[2] = (pi_Nz != 0) ? pi_Nz : 0;
  if (ARRAY_TYPE==true) {  
    pto_TabNaxis[1] = (pi_Ny != 0) ? pi_Ny : 0;
    pto_TabNaxis[0] = (pi_Nx != 0) ? pi_Nx : 0;
    is_ima = false;
  } else {
    pto_TabNaxis[0] = (pi_Ny != 0) ? pi_Ny : 0;
    pto_TabNaxis[1] = (pi_Nx != 0) ? pi_Nx : 0;  
    is_ima = true;
  }
  i_NbAxis = (pi_Nx != 0) ? 1 : 0;
  i_NbAxis = (pi_Ny != 0) ? 2 : i_NbAxis;
  i_NbAxis = (pi_Nz != 0) ? 3 : i_NbAxis;

  memset (po_Buffer, 0, i_NbElem*sizeof(PARAM_TYPE));
  //if (Name != NULL) strcpy(tc_NameArray, Name);
  if (Name != NULL) o_NameArray = Name;

//if (ARRAY_TYPE==true) 
//  cout << "new Ifloat : Nlignes=" << nl() << ", Ncol=" << nc() << endl;
//else
//  cout << "new fltarr : Nlignes=" << nl() << ", Ncol=" << nc() << endl;

}   

//------------------------------------------------------------------------------
// alloc (PARAM_TYPE *BuffData, int Nbr_Line, int Nbr_Col, const char *Name, bool MemManag)
//------------------------------------------------------------------------------
template <class PARAM_TYPE, bool ARRAY_TYPE>
void to_array<PARAM_TYPE,ARRAY_TYPE>::alloc (PARAM_TYPE *BuffData, int Nbr_Line, int Nbr_Col, 
                                  const char *Name, bool MemManag) {
  if (i_NbElem != 0) {
    if (e_UseClassMemAlloc == true) {
#ifdef _USEMEM
       MemMg_free (po_Buffer); 
#endif
    } else if (e_GetBuffer == false) delete [] po_Buffer;  
  }			  
  e_GetBuffer = true;
  e_UseClassMemAlloc = MemManag;
  po_Buffer = BuffData;
  i_NbElem = Nbr_Line * Nbr_Col;
  if (ARRAY_TYPE==true) {
     pto_TabNaxis[1] = Nbr_Col;
     pto_TabNaxis[0] = Nbr_Line;
     is_ima = false;
  } else {
     pto_TabNaxis[0] = Nbr_Col;
     pto_TabNaxis[1] = Nbr_Line; 
     is_ima = true;
 }
  i_NbAxis = 2;				  
}

//------------------------------------------------------------------------------
// alloc (PARAM_TYPE *BuffData, int Nx, int Ny, int Nz, const char *Name, bool MemManag)
//------------------------------------------------------------------------------
template <class PARAM_TYPE, bool ARRAY_TYPE>
void to_array<PARAM_TYPE,ARRAY_TYPE>::alloc (PARAM_TYPE *BuffData, int _Nx, int _Ny, int _Nz, const char *Name, bool MemManag)
{// doesnt work well, treated as an image
	if (i_NbElem != 0)
	{
		if (e_UseClassMemAlloc == true)
		{
			#ifdef _USEMEM
			MemMg_free (po_Buffer); 
			#endif
		}
		else if (e_GetBuffer == false) delete [] po_Buffer;  
	}			  
	e_GetBuffer = true;
	e_UseClassMemAlloc = MemManag;
	po_Buffer = BuffData;
	i_NbElem = _Nx * _Ny * _Nz;
	
	pto_TabNaxis[0] = _Nx;
	pto_TabNaxis[1] = _Ny;
	pto_TabNaxis[2] = _Nz;
	i_NbAxis = 3;
    is_ima = false;
}

//------------------------------------------------------------------------------
// reform (const int pi_Nx, const int pi_Ny, const int pi_Nz)
//------------------------------------------------------------------------------
template <class PARAM_TYPE, bool ARRAY_TYPE>
void to_array<PARAM_TYPE,ARRAY_TYPE>::reform (const int pi_Nx, const int pi_Ny, 
                             const int pi_Nz) {

  if (i_NbElem == 0) alloc(pi_Nx,pi_Ny,pi_Nz,"alloc resize");
  else {

    int ai_Inter;
    i_NbAxis = 1; ai_Inter = pi_Nx;
    pto_TabNaxis[0] = 0; pto_TabNaxis[1] = 0; pto_TabNaxis[2] = 0;
    
    // test Array type    
    if (ARRAY_TYPE==true) { 
      pto_TabNaxis[0] = pi_Nx;
      if (pi_Ny != 0) {pto_TabNaxis[1] = pi_Ny; i_NbAxis=2; 
                       ai_Inter = pi_Nx*pi_Ny;}
     is_ima = false;
   } else {
      pto_TabNaxis[1] = pi_Nx;
      if (pi_Ny != 0) {pto_TabNaxis[0] = pi_Ny; i_NbAxis=2; 
                       ai_Inter = pi_Nx*pi_Ny;}   
      is_ima = true;
   }
    if (pi_Nz != 0) {pto_TabNaxis[2] = pi_Nz; i_NbAxis=3; 
                     ai_Inter = pi_Nx*pi_Ny*pi_Nz;}   
		     
    // increase buffer size
    if (ai_Inter > i_NbElem) {
      
      // deallocate previous bufferr
      if (e_UseClassMemAlloc == true) {
#ifdef _USEMEM
        MemMg_free (po_Buffer);
#endif
      } else if (e_GetBuffer == false && i_NbElem != 0) delete [] po_Buffer;
      
      // allocate new buffer
      if (ai_Inter > TA_MIN_SIZE_FOR_MEM_ALLOC_CALL) {
#ifdef _USEMEM
        e_UseClassMemAlloc = true;
	PARAM_TYPE Dummy=0;
        po_Buffer = MemMg_alloc (ai_Inter, Dummy);
#else
        e_UseClassMemAlloc = false;
        po_Buffer = new PARAM_TYPE [ai_Inter];
        if (po_Buffer == 0) cout << "Not enought memory " << endl;
#endif
      } else {
        e_UseClassMemAlloc = false;
        po_Buffer = new PARAM_TYPE [ai_Inter];
        if (po_Buffer == 0) cout << "Not enought memory " << endl;
      }
      e_GetBuffer = false;
    }
  i_NbElem = ai_Inter;
  }
}

//------------------------------------------------------------------------------
// resize (const int pi_Nx, const int pi_Ny=0, const int pi_Nz=0)
//------------------------------------------------------------------------------
template <class PARAM_TYPE, bool ARRAY_TYPE>
void to_array<PARAM_TYPE,ARRAY_TYPE>::resize (const int pi_Nx, const int pi_Ny, 
                             const int pi_Nz) {
   reform (pi_Nx,pi_Ny, pi_Nz);			     
}

//------------------------------------------------------------------------------
// operator (int x)
//------------------------------------------------------------------------------
// could be used with 2d or 3d tab, on all the element.... => no border test...
template <class PARAM_TYPE, bool ARRAY_TYPE>
inline PARAM_TYPE& to_array<PARAM_TYPE,ARRAY_TYPE>::operator() (int x)  const {
   // if (naxis() != 1) {cout << "One dim array" << endl; exit(-1);}
//!!!!!!!!!   assert (test_indice_i (tc_NameArray, x, nx()));
   return po_Buffer[x];
}
/*
template <class PARAM_TYPE, bool ARRAY_TYPE>
inline PARAM_TYPE& to_array<PARAM_TYPE,ARRAY_TYPE>::operator[] (int x)  const {
   // if (naxis() != 1) {cout << "One dim array" << endl; exit(-1);}
//!!!!!!!!!   assert (test_indice_i (tc_NameArray, x, nx()));
   return po_Buffer[x];
}
*/

//------------------------------------------------------------------------------
// operator (int x, type_border bord)
//------------------------------------------------------------------------------
// !!!!! convert function from int to PARAM_TYPE must exist
template <class PARAM_TYPE, bool ARRAY_TYPE>
inline PARAM_TYPE to_array<PARAM_TYPE,ARRAY_TYPE>::operator() (int x, type_border bord)  const {
 
  if (naxis() != 1) 
  {
      cout << "Error: naxis = " << naxis()  << " and a one dim array is expected ..." << endl; 
      exit(-1);
  }
 
   if ((x<0) || (x>=nx())) {
    PARAM_TYPE Val;
    int indx=x;
    switch (bord) {
    case I_CONT:
      indx = test_index_cont(x,nx());
      Val = po_Buffer[indx]; break;
    case I_MIRROR:
      indx = test_index_mirror(x,nx());
      Val = po_Buffer[indx]; break;     
    case I_PERIOD:
      indx = test_index_period(x,nx());
      Val = po_Buffer[indx]; break;   
    case I_ZERO: Val=0;break;
      break;
    default:exit(-1);break;
    }
    return Val;
  } else return po_Buffer[x];
}


template <class PARAM_TYPE, bool ARRAY_TYPE>
inline PARAM_TYPE & to_array<PARAM_TYPE,ARRAY_TYPE>::setx (int x, type_border bord)    {
 
  if (naxis() != 1) 
  {
      cout << "Error: naxis = " << naxis()  << " and a one dim array is expected ..." << endl; 
      exit(-1);
  }
 
   if ((x<0) || (x>=nx())) {
    int indx=x;
    switch (bord) {
    case I_CONT:
      indx = test_index_cont(x,nx());
      JunkVar = po_Buffer[indx]; 
      return JunkVar;
      break;
    case I_MIRROR:
      indx = test_index_mirror(x,nx());
      return po_Buffer[indx]; break;     
    case I_PERIOD:
      indx = test_index_period(x,nx());
      return po_Buffer[indx]; break;   
    case I_ZERO: JunkVar=0; return JunkVar; break;
      break;
    default:exit(-1);break;
    }
   } else return po_Buffer[x];
} 


//------------------------------------------------------------------------------
// operator (int x, int y)
//------------------------------------------------------------------------------
template <class PARAM_TYPE, bool ARRAY_TYPE>
inline PARAM_TYPE& to_array<PARAM_TYPE,ARRAY_TYPE>::operator() (int x, int y) const {
 
   if (naxis() != 2) 
   {
      cout << "Error: naxis = " << naxis() << " and a two dim array is expected ..." << endl; 
      exit(-1);
  }
 
   if (ARRAY_TYPE==true) {
     return po_Buffer[y*pto_TabNaxis[0]+x];
   } else {
     return po_Buffer[x*pto_TabNaxis[0]+y];
   }
}

//------------------------------------------------------------------------------
// operator (int x, int y, type_border bord)
//------------------------------------------------------------------------------

template <class PARAM_TYPE, bool ARRAY_TYPE>
inline PARAM_TYPE to_array<PARAM_TYPE,ARRAY_TYPE>::operator() (int x, int y, type_border bord) const {
 
  if (naxis() != 2) 
  {
      cout << "Error: naxis = " << naxis() << " and a two dim array is expected ... " << endl; 
      exit(-1);
  }
 
  int indx=x; int indy=y;
  int Nx = (ARRAY_TYPE==true) ? nx(): ny();
  int Ny = (ARRAY_TYPE==true) ? ny(): nx();
  
  if ((x<0) || (x>=Nx) || (y<0) || (y>=Ny)) {
    PARAM_TYPE Val;
    switch (bord) {
    case I_CONT:
      if (ARRAY_TYPE==true) {
        indx = test_index_cont(x,nx());
        indy = test_index_cont(y,ny());
      } else {
        indx = test_index_cont(x,ny());
        indy = test_index_cont(y,nx());     
      } 
      break;     
    case I_MIRROR:
      if (ARRAY_TYPE==true) {    
        indx = test_index_mirror(x,nx());
        indy = test_index_mirror(y,ny());
      } else {
        indx = test_index_mirror(x,ny());
        indy = test_index_mirror(y,nx());      
      }
      break;
    case I_PERIOD:
      if (ARRAY_TYPE==true) {    
        indx = test_index_period(x,nx());
        indy = test_index_period(y,ny());
      } else {
        indx = test_index_period(x,ny());
        indy = test_index_period(y,nx());      
      }
      break;
    case I_ZERO:
      Val=0; return Val; break;
    default:
        exit(-1);break;
    }
  } 
  if (ARRAY_TYPE==true) {
     return po_Buffer[indy*pto_TabNaxis[0]+indx];
   } else {
     return po_Buffer[indx*pto_TabNaxis[0]+indy];
   }
} 


template <class PARAM_TYPE, bool ARRAY_TYPE>
inline PARAM_TYPE & to_array<PARAM_TYPE,ARRAY_TYPE>::setxy (int x, int y, type_border bord)    {
 
 if (naxis() != 2) 
  {
      cout << "Error: naxis = " << naxis() << " and a two dim array is expected ... " << endl; 
      exit(-1);
  }
 
  int indx=x; int indy=y;
  int Nx = (ARRAY_TYPE==true) ? nx(): ny();
  int Ny = (ARRAY_TYPE==true) ? ny(): nx();
  
 if ((x<0) || (x>=Nx) || (y<0) || (y>=Ny)) {
     switch (bord) 
    {
       case I_ZERO:
           JunkVar = 0;
            return JunkVar;
        case I_CONT:
             if (is_ima==false) {
        indx = test_index_cont(x,nx());
        indy = test_index_cont(y,ny());
        JunkVar = po_Buffer[indy*pto_TabNaxis[0]+indx];
       } else {
        indx = test_index_cont(x,ny());
        indy = test_index_cont(y,nx());     
       JunkVar =  po_Buffer[indx*pto_TabNaxis[0]+indy];
        } 
           return JunkVar;
           break;
       case I_MIRROR:
         if (ARRAY_TYPE==true) {    
        indx = test_index_mirror(x,nx());
        indy = test_index_mirror(y,ny());
        } else {
        indx = test_index_mirror(x,ny());
        indy = test_index_mirror(y,nx());      
       }
       break;
    case I_PERIOD:
      if (ARRAY_TYPE==true) {    
        indx = test_index_period(x,nx());
        indy = test_index_period(y,ny());
      } else {
        indx = test_index_period(x,ny());
        indy = test_index_period(y,nx());      
      }
      break;
     default:
        exit(-1);break;
    }
  } 

  if (ARRAY_TYPE==true) {
     return po_Buffer[indy*pto_TabNaxis[0]+indx];
   } else {
     return po_Buffer[indx*pto_TabNaxis[0]+indy];
   }
} 


 //------------------------------------------------------------------------------
// operator (int x, int y, int z)
//------------------------------------------------------------------------------
template <class PARAM_TYPE, bool ARRAY_TYPE>
inline PARAM_TYPE& to_array<PARAM_TYPE,ARRAY_TYPE>::operator() (int x, int y, int z) const {
#if CHECK_DIM
   if (naxis() != 3) 
   {
     cout << "Error: naxis = " << naxis() << " and a three dimension array is expected ... " << endl; 
     exit(-1);
  }
#endif
   return po_Buffer[z*pto_TabNaxis[0]*pto_TabNaxis[1]+y*pto_TabNaxis[0]+x];
}

 
//------------------------------------------------------------------------------
// operator (int x, int y, int z, type_border bord)
//------------------------------------------------------------------------------

template <class PARAM_TYPE, bool ARRAY_TYPE>
inline PARAM_TYPE to_array<PARAM_TYPE,ARRAY_TYPE>::operator() (int x, int y, int z, type_border bord) const {
#if CHECK_DIM
  if (naxis() != 3) 
  {
     cout << "Error: naxis = " << naxis() << " and a three dimension array is expected ... " << endl; 
     exit(-1);
  }
#endif
  int indx=x; int indy=y; int indz=z;
 if ((x<0) || (x>=nx()) || (y<0) || (y>=ny()) || (z<0) || (z>=nz())) {
    PARAM_TYPE Val;
    switch (bord) {
    case I_CONT:
      indx = test_index_cont(x,nx());
      indy = test_index_cont(y,ny());     
      indz = test_index_cont(z,nz()); 
      break;      
    case I_MIRROR:
      indx = test_index_mirror(x,nx());
      indy = test_index_mirror(y,ny());
      indz = test_index_mirror(z,nz());
      break;
    case I_PERIOD:
      indx = test_index_period(x,nx());
      indy = test_index_period(y,ny());
      indz = test_index_period(z,nz());
      break;
    case I_ZERO: Val=0; return Val; break;
      break;
    default:exit(-1);break;
    }
  } 
  return po_Buffer[indz*pto_TabNaxis[0]*pto_TabNaxis[1]+indy*pto_TabNaxis[0]+indx];
}


template <class PARAM_TYPE, bool ARRAY_TYPE>
inline PARAM_TYPE & to_array<PARAM_TYPE,ARRAY_TYPE>::setxyz (int x, int y, int z, type_border bord)    {
#if CHECK_DIM
  if (naxis() != 3) 
  {
     cout << "Error: naxis = " << naxis() << " and a three dimension array is expected ... " << endl; 
     exit(-1);
  }
#endif
  int indx=x; int indy=y; int indz=z;
  if ((x<0) || (x>=nx()) || (y<0) || (y>=ny()) || (z<0) || (z>=nz())) {
    switch (bord) {
    case I_CONT:
      indx = test_index_cont(x,nx());
      indy = test_index_cont(y,ny());     
      indz = test_index_cont(z,nz()); 
      return JunkVar;
      break;      
    case I_MIRROR:
      indx = test_index_mirror(x,nx());
      indy = test_index_mirror(y,ny());
      indz = test_index_mirror(z,nz());
      break;
    case I_PERIOD:
      indx = test_index_period(x,nx());
      indy = test_index_period(y,ny());
      indz = test_index_period(z,nz());
      break;
    case I_ZERO: JunkVar=0; return JunkVar; break;
      break;
    default:exit(-1);break;
    }
  } 
  return po_Buffer[indz*pto_TabNaxis[0]*pto_TabNaxis[1]+indy*pto_TabNaxis[0]+indx];
}






//------------------------------------------------------------------------------
// operator =
//------------------------------------------------------------------------------
template <class PARAM_TYPE, bool ARRAY_TYPE> const 
to_array<PARAM_TYPE,ARRAY_TYPE>& to_array<PARAM_TYPE,ARRAY_TYPE>::operator = (const to_array<PARAM_TYPE,ARRAY_TYPE>& pro_Mat) {
  reform (pro_Mat.n_elem());
  for (int i=0; i<i_NbElem; i++) po_Buffer[i] = pro_Mat.po_Buffer[i]; 
  i_NbAxis = pro_Mat.naxis();
  for (int j=0; j<i_NbAxis; j++) pto_TabNaxis[j] = pro_Mat.axis(j+1);
  return (*this);
}

//------------------------------------------------------------------------------
// operator +=
//------------------------------------------------------------------------------
template <class PARAM_TYPE, bool ARRAY_TYPE> const 
to_array<PARAM_TYPE,ARRAY_TYPE>& to_array<PARAM_TYPE,ARRAY_TYPE>::operator += (const to_array<PARAM_TYPE,ARRAY_TYPE>& pro_Mat) {
  for (int x=0; x<i_NbElem; x++) po_Buffer[x] += pro_Mat.po_Buffer[x];
  return (*this);
}

//------------------------------------------------------------------------------
// operator *=
//------------------------------------------------------------------------------
template <class PARAM_TYPE, bool ARRAY_TYPE> const 
to_array<PARAM_TYPE,ARRAY_TYPE>& to_array<PARAM_TYPE,ARRAY_TYPE>::operator *= (const to_array<PARAM_TYPE,ARRAY_TYPE>& pro_Mat) {
  for (int x=0; x<i_NbElem; x++) po_Buffer[x] *= pro_Mat.po_Buffer[x];
  return (*this);
}

//------------------------------------------------------------------------------
// operator *=
//------------------------------------------------------------------------------
template <class PARAM_TYPE, bool ARRAY_TYPE> const 
to_array<PARAM_TYPE,ARRAY_TYPE>& to_array<PARAM_TYPE,ARRAY_TYPE>::operator *= (const double coef) {
  for (int i=0; i<i_NbElem; i++) 
    po_Buffer[i] = (PARAM_TYPE) ((double)po_Buffer[i] * coef);
  return (*this);
}

//------------------------------------------------------------------------------
// operator -=
//------------------------------------------------------------------------------
template <class PARAM_TYPE, bool ARRAY_TYPE> const 
to_array<PARAM_TYPE,ARRAY_TYPE>& to_array<PARAM_TYPE,ARRAY_TYPE>::operator -= (const to_array<PARAM_TYPE,ARRAY_TYPE>& pro_Mat) {
  for (int x=0; x<i_NbElem; x++) po_Buffer[x] -= pro_Mat.po_Buffer[x];
  return (*this);
}

//------------------------------------------------------------------------------
// operator /=
//------------------------------------------------------------------------------
template <class PARAM_TYPE, bool ARRAY_TYPE> const 
to_array<PARAM_TYPE,ARRAY_TYPE>& to_array<PARAM_TYPE,ARRAY_TYPE>::operator /= (const to_array<PARAM_TYPE,ARRAY_TYPE>& pro_Mat) {
  for (int x=0; x<i_NbElem; x++) 
     if ((pro_Mat.po_Buffer[x] > 1e-07) || (pro_Mat.po_Buffer[x] < -1e-07))
        po_Buffer[x] /= pro_Mat.po_Buffer[x];
     else po_Buffer[x]=0;
  return (*this);
}

//------------------------------------------------------------------------------
// operator ^
//------------------------------------------------------------------------------
template <class PARAM_TYPE, bool ARRAY_TYPE> const 
// !!!!! convert function from PARAM_TYPE to double must exist
to_array<PARAM_TYPE,ARRAY_TYPE>& to_array<PARAM_TYPE,ARRAY_TYPE>::operator ^ (const double pf_coef) {
  for (int i=0; i<i_NbElem; i++) 
    po_Buffer[i] = (PARAM_TYPE) pow ((double)po_Buffer[i], pf_coef);
  return (*this);
}

//------------------------------------------------------------------------------
// info (string Name)
//------------------------------------------------------------------------------
template <class PARAM_TYPE, bool ARRAY_TYPE>
// !!!!! PARAM_TYPE must accept operator << !!!!!!!!!!!
void to_array<PARAM_TYPE,ARRAY_TYPE>::info(string Name) {
  if (Name=="") cout << "   :" ;
  else cout << "  " << Name;
  cout << ", mean = " << mean() << ", sigma = " << sigma();
  cout << ", min = " << min() << ", max = " << max() << endl;
}

//------------------------------------------------------------------------------
// display (int pi_NbElem)
//------------------------------------------------------------------------------
template <class PARAM_TYPE, bool ARRAY_TYPE>
// !!!!! PARAM_TYPE must accept operator << !!!!!!!!!!!
void to_array<PARAM_TYPE,ARRAY_TYPE>::display (int pi_NbElem) {
  if (pi_NbElem == 0) {
       cout <<"  nx="<<pto_TabNaxis[0]<<", ny="<<pto_TabNaxis[1]<<
              ", nz="<<pto_TabNaxis[2]<<", naxis="<<i_NbAxis<<endl;      
  } else {
    if (pi_NbElem > i_NbElem) pi_NbElem=i_NbElem;
    info();
    cout << "  ";
    for (int i=0; i < pi_NbElem; i++)   
      cout << po_Buffer[i] << " " ;
    cout << endl;
  }
}

//------------------------------------------------------------------------------
// rampgen()
//------------------------------------------------------------------------------
template <class PARAM_TYPE, bool ARRAY_TYPE>
// !!!!! convert function from int to PARAM_TYPE must exist
void to_array<PARAM_TYPE,ARRAY_TYPE>::rampgen() {
  for (int i=0;i<i_NbElem;i++) po_Buffer[i]=(PARAM_TYPE)i;
}

//------------------------------------------------------------------------------
// line()
//------------------------------------------------------------------------------
//template <class PARAM_TYPE, bool ARRAY_TYPE>
//to_array<PARAM_TYPE,ARRAY_TYPE> to_array<PARAM_TYPE,ARRAY_TYPE>::line (int i) {
//}

//------------------------------------------------------------------------------
// column()
//------------------------------------------------------------------------------
//template <class PARAM_TYPE, bool ARRAY_TYPE>
//to_array<PARAM_TYPE,ARRAY_TYPE> to_array<PARAM_TYPE,ARRAY_TYPE>::column (int j) {
//}

//------------------------------------------------------------------------------
// sup_threshold (float ThresholLevel)
//------------------------------------------------------------------------------
template <class PARAM_TYPE, bool ARRAY_TYPE>
void to_array<PARAM_TYPE,ARRAY_TYPE>::sup_threshold (float ThresholLevel) {
  for (int x=0;x<i_NbElem;x++) 
    if ((PARAM_TYPE)po_Buffer[x] > ThresholLevel) 
      po_Buffer[x] = (PARAM_TYPE)ThresholLevel;
}


//------------------------------------------------------------------------------
// inf_threshold (float ThresholLevel)
//------------------------------------------------------------------------------
template <class PARAM_TYPE, bool ARRAY_TYPE>
void to_array<PARAM_TYPE,ARRAY_TYPE>::inf_threshold (float ThresholLevel) {
  for (int x=0;x<i_NbElem;x++) 
    if ((PARAM_TYPE)po_Buffer[x] < ThresholLevel) 
      po_Buffer[x] = (PARAM_TYPE)ThresholLevel;
}

//------------------------------------------------------------------------------
// min ()
//------------------------------------------------------------------------------
template <class PARAM_TYPE, bool ARRAY_TYPE>
PARAM_TYPE to_array<PARAM_TYPE,ARRAY_TYPE>::min () {
  int ai_temp=0;
  return (min (ai_temp));
}

//------------------------------------------------------------------------------
// min (int& pri_ind)
//------------------------------------------------------------------------------
template <class PARAM_TYPE, bool ARRAY_TYPE>
PARAM_TYPE to_array<PARAM_TYPE,ARRAY_TYPE>::min (int& pri_ind) {
  PARAM_TYPE ao_prov=po_Buffer[0];
  pri_ind=0;
  for (int i=1; i<i_NbElem; i++) 
    if (ao_prov>po_Buffer[i]) {
      ao_prov=po_Buffer[i];
      pri_ind = i;
    }
  return ao_prov;
}

//------------------------------------------------------------------------------
// max ()
//------------------------------------------------------------------------------ 
template <class PARAM_TYPE, bool ARRAY_TYPE>
PARAM_TYPE to_array<PARAM_TYPE,ARRAY_TYPE>::max () {
  int ai_temp=0;
  return (max (ai_temp));
}

//------------------------------------------------------------------------------
// max (int& pri_ind)
//------------------------------------------------------------------------------ 
template <class PARAM_TYPE, bool ARRAY_TYPE>
//!!!!!!!!!! PARAM_TYPE must define operator <...
PARAM_TYPE to_array<PARAM_TYPE,ARRAY_TYPE>::max (int& pri_ind) {
  PARAM_TYPE ao_prov=po_Buffer[0];
  pri_ind=0;
  for (int i=1; i<i_NbElem; i++) 
    if (ao_prov<po_Buffer[i]) {
      ao_prov=po_Buffer[i];
      pri_ind = i;
    }
  return ao_prov;
}

//------------------------------------------------------------------------------
// maxfabs ()
//------------------------------------------------------------------------------ 
template <class PARAM_TYPE, bool ARRAY_TYPE>
PARAM_TYPE to_array<PARAM_TYPE,ARRAY_TYPE>::maxfabs () {
  int ai_temp=0;
  return (maxfabs (ai_temp));
}

//------------------------------------------------------------------------------
// maxfabs (int& pri_ind)
//------------------------------------------------------------------------------ 
template <class PARAM_TYPE, bool ARRAY_TYPE>
//!!!!!!!!!! PARAM_TYPE must define operator <...
PARAM_TYPE to_array<PARAM_TYPE,ARRAY_TYPE>::maxfabs (int& pri_ind) {
  PARAM_TYPE ao_prov=0;
  pri_ind=0;
  for (int i=0; i<i_NbElem; i++) 
    if (fabs(ao_prov)<fabs(po_Buffer[i])) {
      ao_prov=po_Buffer[i];
      pri_ind = i;
    }
  return ao_prov;
}

//------------------------------------------------------------------------------
// total ()
//------------------------------------------------------------------------------
template <class PARAM_TYPE, bool ARRAY_TYPE>
// !!!!! convert function from PARAM_TYPE to double must exist
double to_array<PARAM_TYPE,ARRAY_TYPE>::total () const {
  double ao_prov= 0.;
  for (int i=0; i<i_NbElem; i++) ao_prov += po_Buffer[i];
  return ao_prov;
}

//------------------------------------------------------------------------------
// energy ()
//------------------------------------------------------------------------------
template <class PARAM_TYPE, bool ARRAY_TYPE>
// !!!!! convert function from PARAM_TYPE to double must exist
double to_array<PARAM_TYPE,ARRAY_TYPE>::energy () const {
  PARAM_TYPE ao_prov=(PARAM_TYPE)0;
  for (int i=0; i<i_NbElem; i++) ao_prov += po_Buffer[i]*po_Buffer[i];
  return ao_prov;
}

//------------------------------------------------------------------------------
// sigma ()
//------------------------------------------------------------------------------
template <class PARAM_TYPE, bool ARRAY_TYPE>
// !!!!! convert function from PARAM_TYPE to double must exist
double to_array<PARAM_TYPE,ARRAY_TYPE>::sigma () const {
  double ao_moy = mean();
  double ad_sigma=0., ad_val=0;
  for (int i=0; i<i_NbElem; i++) {
    ad_val = po_Buffer[i] - ao_moy;
    ad_sigma += ad_val*ad_val;
  }
  if ((ad_sigma /= i_NbElem) > 1e-07) ad_sigma = sqrt (ad_sigma);
  else ad_sigma = 0.;
  return ad_sigma;
}

//------------------------------------------------------------------------------
// mean ()
//------------------------------------------------------------------------------
template <class PARAM_TYPE, bool ARRAY_TYPE>
// !!!!! convert function from PARAM_TYPE to double must exist
double to_array<PARAM_TYPE,ARRAY_TYPE>::mean () const {
  return (double(total())/i_NbElem);
}
 
//------------------------------------------------------------------------------
// sigma_clip (float& pf_Mean, float &pf_Sigma, int pi_Nit)
//------------------------------------------------------------------------------
template <class PARAM_TYPE, bool ARRAY_TYPE>
// !!!!! convert function from PARAM_TYPE to double must exist
void to_array<PARAM_TYPE,ARRAY_TYPE>::sigma_clip (float& pf_Mean, float &pf_Sigma, 
                                       int pi_Nit) const {

  double ad_s0, ad_s1, ad_s2, ad_sm=0., ad_inter;
  PARAM_TYPE ao_val;
  pf_Mean = 0.;
  for (int it=0; it<pi_Nit; it++) {
    ad_s0=ad_s1=ad_s2=0.;
    for (int i=0; i<i_NbElem; i++) {
      ao_val = po_Buffer[i];
      if ((it==0) || (fabs(double(ao_val)-pf_Mean) < ad_sm)) {
	ad_s0++; ad_s1 += double(ao_val); 
	ad_s2 += double(ao_val)*double(ao_val);
      }
    }
    pf_Mean = ad_s1/ad_s0;
    ad_inter = ad_s2/ad_s0 - pf_Mean*pf_Mean;
    if (ad_inter > 1e-7) pf_Sigma = sqrt (ad_inter);
    ad_sm = 3. * pf_Sigma;
  }
} 

//------------------------------------------------------------------------------
// sigma_clip (int pi_Nit)
//------------------------------------------------------------------------------
template <class PARAM_TYPE, bool ARRAY_TYPE>
float to_array<PARAM_TYPE,ARRAY_TYPE>::sigma_clip (int pi_Nit) const {
  float af_Mean=0., af_Sigma=0.;
  sigma_clip (af_Mean, af_Sigma, pi_Nit);
  return (af_Sigma);
}

//------------------------------------------------------------------------------
// to_array (to_array&)
//------------------------------------------------------------------------------
template <class PARAM_TYPE, bool ARRAY_TYPE>
to_array<PARAM_TYPE,ARRAY_TYPE>::to_array 
  (const to_array<PARAM_TYPE,ARRAY_TYPE>& pro_Obj) {
  i_NbElem=0;
  if (ARRAY_TYPE==true) { 
     alloc (pro_Obj.nx(), pro_Obj.ny(), pro_Obj.nz());
  } else {
     alloc (pro_Obj.ny(), pro_Obj.nx(), pro_Obj.nz());
  }
  for (int i=0;i<n_elem();i++) po_Buffer[i]=pro_Obj.po_Buffer[i];
}

//------------------------------------------------------------------------------
// set_attrib (int pi_Nit)
//------------------------------------------------------------------------------
template <class PARAM_TYPE, bool ARRAY_TYPE>
void to_array<PARAM_TYPE,ARRAY_TYPE>::set_attrib () {
  po_Buffer = (PARAM_TYPE*)NULL;
  i_NbElem = 0;
  i_NbAxis = 0;
  for (int i=0;i<MAX_NBR_AXIS;i++)
    pto_TabNaxis[i]=0;
  //tc_NameArray[0]='\0';
  //o_NameArray = "";
  e_UseClassMemAlloc = false;
  e_GetBuffer = false; 
  is_ima=false;
  JunkVar=0;
}


//------------------------------------------------------------------------------
// operator + (to_array, to_array)
//------------------------------------------------------------------------------
template <class PARAM_TYPE, bool ARRAY_TYPE>
to_array<PARAM_TYPE,ARRAY_TYPE> operator + (const to_array<PARAM_TYPE,ARRAY_TYPE>& pro_obj1, 
                           const to_array<PARAM_TYPE,ARRAY_TYPE>& pro_obj2) {
  //to_array<PARAM_TYPE,ARRAY_TYPE>* apo_array = new to_array<PARAM_TYPE,ARRAY_TYPE>;
  to_array<PARAM_TYPE,ARRAY_TYPE> ao_array;
  //apo_array->init(pro_obj1);
  ao_array.init(pro_obj1);
  for (int i=0; i<pro_obj1.n_elem(); i++)
    ao_array(i) = pro_obj1(i) + pro_obj2(i);
    //(*apo_array)(i) = pro_obj1(i) + pro_obj2(i);
  return (to_array<PARAM_TYPE,ARRAY_TYPE>(ao_array));
}

//------------------------------------------------------------------------------
// operator - (to_array, to_array)
//------------------------------------------------------------------------------
template <class PARAM_TYPE, bool ARRAY_TYPE>
to_array<PARAM_TYPE,ARRAY_TYPE> operator - (const to_array<PARAM_TYPE,ARRAY_TYPE>& pro_obj1, 
                           const to_array<PARAM_TYPE,ARRAY_TYPE>& pro_obj2) {
  to_array<PARAM_TYPE,ARRAY_TYPE> ao_array;
  ao_array.init(pro_obj1);
  for (int i=0; i<pro_obj1.n_elem(); i++) 
    ao_array(i)  = pro_obj1(i) - pro_obj2(i);
  return (to_array<PARAM_TYPE,ARRAY_TYPE>(ao_array));
}

//------------------------------------------------------------------------------
// operator * (to_array, to_array)
//------------------------------------------------------------------------------
template <class PARAM_TYPE, bool ARRAY_TYPE>
to_array<PARAM_TYPE,ARRAY_TYPE> operator * (const to_array<PARAM_TYPE,ARRAY_TYPE>& pro_obj1, 
                           const to_array<PARAM_TYPE,ARRAY_TYPE>& pro_obj2) {
  to_array<PARAM_TYPE,ARRAY_TYPE> ao_array;
  ao_array.init(pro_obj1);
  for (int i=0; i<pro_obj1.n_elem(); i++) 
    ao_array(i) = pro_obj1(i) * pro_obj2(i);
  return (to_array<PARAM_TYPE,ARRAY_TYPE>(ao_array));
}

//------------------------------------------------------------------------------
// operator / (to_array, to_array)
//------------------------------------------------------------------------------ 
template <class PARAM_TYPE, bool ARRAY_TYPE>
to_array<PARAM_TYPE,ARRAY_TYPE> operator / (const to_array<PARAM_TYPE,ARRAY_TYPE>& pro_obj1, 
                           const to_array<PARAM_TYPE,ARRAY_TYPE>& pro_obj2) {
  to_array<PARAM_TYPE,ARRAY_TYPE> ao_array;
  ao_array.init(pro_obj1);
  for (int i=0; i<pro_obj1.n_elem(); i++) 
    if ((pro_obj2(i) > 1e-07) || (pro_obj2(i) < -1e-07)) 
      ao_array(i) = pro_obj1(i) / pro_obj2(i);
    else ao_array(i) = 0;
  return (to_array<PARAM_TYPE,ARRAY_TYPE>(ao_array));
}


//------------------------------------------------------------------------------
// operator * (double , to_array)
//------------------------------------------------------------------------------ 
template <class PARAM_TYPE, bool ARRAY_TYPE>
to_array<PARAM_TYPE,ARRAY_TYPE> operator * (const double mult_coeff, 
											const to_array<PARAM_TYPE,ARRAY_TYPE>& pro_obj) {
	to_array<PARAM_TYPE,ARRAY_TYPE> ao_array;
	ao_array.init(pro_obj);
	for (int i=0; i<pro_obj.nx(); i++) 
		ao_array(i)=mult_coeff*pro_obj(i);
	return (to_array<PARAM_TYPE,ARRAY_TYPE>(ao_array));
}
//------------------------------------------------------------------------------
// operator / (to_array,double)
//------------------------------------------------------------------------------ 
template <class PARAM_TYPE, bool ARRAY_TYPE>
to_array<PARAM_TYPE,ARRAY_TYPE> operator / (const to_array<PARAM_TYPE,ARRAY_TYPE>& pro_obj,
											const double div_coeff) 
{
	to_array<PARAM_TYPE,ARRAY_TYPE> ao_array;
	ao_array.init(pro_obj);
	for (int i=0; i<pro_obj.nx(); i++) 
		ao_array(i)=pro_obj(i)/div_coeff;
	return (to_array<PARAM_TYPE,ARRAY_TYPE>(ao_array));
}

//------------------------------------------------------------------------------
// operator > (to_array,double)
//------------------------------------------------------------------------------ 
template <class PARAM_TYPE, bool ARRAY_TYPE>
to_array<PARAM_TYPE,ARRAY_TYPE> operator > (const to_array<PARAM_TYPE,ARRAY_TYPE>& pro_obj,
											const double bound_coeff) {
	to_array<PARAM_TYPE,ARRAY_TYPE> ao_array;
	ao_array.init(pro_obj);
	for (long int i=0; i<pro_obj.nx(); i++) 
		if(pro_obj(i) > bound_coeff) ao_array(i)=1.0;
		else ao_array(i)=0.0;
	return (to_array<PARAM_TYPE,ARRAY_TYPE>(ao_array));
}

//------------------------------------------------------------------------------
// operator < (to_array,double)
//------------------------------------------------------------------------------ 
template <class PARAM_TYPE, bool ARRAY_TYPE>
to_array<PARAM_TYPE,ARRAY_TYPE> operator < (const to_array<PARAM_TYPE,ARRAY_TYPE>& pro_obj,
											const double bound_coeff) {
	to_array<PARAM_TYPE,ARRAY_TYPE> ao_array;
	ao_array.init(pro_obj);
	for (long int i=0; i<pro_obj.nx(); i++) 
		if(pro_obj(i) < bound_coeff) ao_array(i)=1.0;
		else ao_array(i)=0.0;
	return (to_array<PARAM_TYPE,ARRAY_TYPE>(ao_array));
}

//------------------------------------------------------------------------------
// mult (to_array,to_array)
//------------------------------------------------------------------------------
template <class PARAM_TYPE, bool ARRAY_TYPE>
to_array<PARAM_TYPE,ARRAY_TYPE> mult (const to_array<PARAM_TYPE,ARRAY_TYPE>& pro_obj1, 
											const to_array<PARAM_TYPE,ARRAY_TYPE>& pro_obj2) {
	if(pro_obj1.ny() != pro_obj2.nx()) 
	{
		printf("Can't multiply: 1st matrix number of columns different from 2nd matrix number of rows. \n");
		exit(-1);
	}
	to_array<PARAM_TYPE,ARRAY_TYPE> ao_array;
	if(pro_obj2.ny()>0)
	{
		ao_array.alloc(pro_obj1.nx(),pro_obj2.ny());
		for (int i=0; i<pro_obj1.nx(); i++)
		{
			for (int j=0; j<pro_obj2.ny(); j++) 
			{
				for (int k=0; k<pro_obj1.ny(); k++)
				{
					ao_array(i,j) += pro_obj1(i,k) * pro_obj2(k,j);
				}
			}
		}
	}
	else
	{
		ao_array.alloc(pro_obj1.nx());
		for (int i=0; i<pro_obj1.nx(); i++)
		{
			for (int k=0; k<pro_obj1.ny(); k++)
				ao_array(i) += pro_obj1(i,k) * pro_obj2(k);
		}
	}
    return (to_array<PARAM_TYPE,ARRAY_TYPE>(ao_array));
}

//------------------------------------------------------------------------------
// transpose (to_array)
//------------------------------------------------------------------------------
template <class PARAM_TYPE, bool ARRAY_TYPE>
to_array<PARAM_TYPE,ARRAY_TYPE> transpose (const to_array<PARAM_TYPE,ARRAY_TYPE>& pro_obj) {

	to_array<PARAM_TYPE,ARRAY_TYPE> ao_array;
	ao_array.alloc(pro_obj.ny(),pro_obj.nx());
	for(int i=0;i<pro_obj.nx();i++)
	{
		for(int j=0;j<pro_obj.ny();j++)
		{
			ao_array(j,i)=pro_obj(i,j);
		}
	}
	return (to_array<PARAM_TYPE,ARRAY_TYPE>(ao_array));
}


//------------------------------------------------------------------------------
// invert (to_array)
//------------------------------------------------------------------------------
template <class PARAM_TYPE, bool ARRAY_TYPE>
to_array<PARAM_TYPE,ARRAY_TYPE> invert (const to_array<PARAM_TYPE,ARRAY_TYPE>& pro_obj) {
	if(pro_obj.nx() !=2 || pro_obj.ny()!=2) 
	{
		printf("Matrix must be 2x2 to be inverted. \n");
		exit(-1);
	}
	to_array<PARAM_TYPE,ARRAY_TYPE> ao_array;
	ao_array.alloc(2,2);
	double det=pro_obj(0,0)*pro_obj(1,1)-pro_obj(0,1)*pro_obj(1,0);
	ao_array(0,0)=1/det*pro_obj(1,1);
	ao_array(1,1)=1/det*pro_obj(0,0);
	ao_array(0,1)=-1/det*pro_obj(0,1);
	ao_array(1,0)=-1/det*pro_obj(1,0);
    return (to_array<PARAM_TYPE,ARRAY_TYPE>(ao_array));
}

//------------------------------------------------------------------------------
// mult (to_array,to_array, int)
//------------------------------------------------------------------------------
template <class PARAM_TYPE, bool ARRAY_TYPE>
to_array<PARAM_TYPE,ARRAY_TYPE> mult (const to_array<PARAM_TYPE,ARRAY_TYPE>& pro_obj1, 
									  const to_array<PARAM_TYPE,ARRAY_TYPE>& pro_obj2, int index_z) {
	if(index_z >= pro_obj1.nz()) printf("Wrong z index. \n");
	if(pro_obj1.ny() != pro_obj2.nx()) 
	{
		printf("Can't multiply: 1st matrix number of columns different from 2nd matrix number of rows. \n");
		exit(-1);
	}
	to_array<PARAM_TYPE,ARRAY_TYPE> ao_array;
	if(pro_obj2.ny()>0)
	{
		ao_array.alloc(pro_obj1.nx(),pro_obj2.ny());
		for (int i=0; i<pro_obj1.nx(); i++)
		{
			for (int j=0; j<pro_obj2.ny(); j++) 
			{
				for (int k=0; k<pro_obj1.ny(); k++)
				{
					ao_array(i,j) += pro_obj1(i,k,index_z) * pro_obj2(k,j);
				}
			}
		}
	}
	else
	{
		ao_array.alloc(pro_obj1.nx());
		for (int i=0; i<pro_obj1.nx(); i++)
		{
			for (int k=0; k<pro_obj1.ny(); k++)
				ao_array(i) += pro_obj1(i,k,index_z) * pro_obj2(k);
		}
	}
    return (to_array<PARAM_TYPE,ARRAY_TYPE>(ao_array));
}

//------------------------------------------------------------------------------
// log (to_array)
//------------------------------------------------------------------------------
template <class PARAM_TYPE, bool ARRAY_TYPE>
to_array<PARAM_TYPE,ARRAY_TYPE> log (const to_array<PARAM_TYPE,ARRAY_TYPE>& pro_obj) {
	
	to_array<PARAM_TYPE,ARRAY_TYPE> ao_array;
	ao_array.alloc(pro_obj.nx(),pro_obj.ny(),pro_obj.nz());
	for(int i=0;i<pro_obj.nx();i++)
		for(int j=0;j<pro_obj.ny();j++)
			for(int k=0;k<pro_obj.nz();k++)
				ao_array(i,j,k)=log(pro_obj(i,j,k));
	
    return (to_array<PARAM_TYPE,ARRAY_TYPE>(ao_array));
}


#endif


