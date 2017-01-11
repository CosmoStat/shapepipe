
/******************************************************************************
**                   Copyright (C) 2002 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: J.L. Starck
**
**    Date:  31.8.99
**    
**    File: IM_regul.h
**
*******************************************************************************
**
**    DESCRIPTION  Class for regularization
**    -----------  
**                 
******************************************************************************/

#ifndef _HREGUL_H_
#define _HREGUL_H_

/****************************************************************************/

enum type_oper_grad {OPER_LAPLACIAN,OPER_DIR_LAPLACIAN,OPER_MARKOV_4,OPER_MARKOV_8,OPER_ENERGY};

class RegulIma {
  inline int sgn(float Val) {return ( (Val >= 0) ? 1: -1);}

  inline float laplacian_val_dir(Ifloat &Obj, int i, int j)
  {
     int Ndir=4;
     fltarray GradDir(Ndir);
     fltarray GradParam(Ndir);
     float Total, Val=0.;
     GradDir(0) = Obj(i,j) - 0.5*( Obj(i+1,j,BorderType) + Obj(i-1,j,BorderType));
     GradDir(1) = Obj(i,j) - 0.5*( Obj(i,j+1,BorderType) + Obj(i,j-1,BorderType));
     GradDir(2) = Obj(i,j) - 0.5*( Obj(i+1,j+1,BorderType) + Obj(i-1,j-1,BorderType));
     GradDir(3) = Obj(i,j) - 0.5*( Obj(i+1,j-1,BorderType) + Obj(i-1,j+1,BorderType));
     GradParam(0) =  (Obj(i,j) - Obj(i,j-1,BorderType))*(Obj(i,j) - Obj(i,j-1,BorderType));
     GradParam(1) =  (Obj(i,j) - Obj(i,i-1,BorderType))*(Obj(i,j) - Obj(i,i-1,BorderType));
     GradParam(2) =  (Obj(i,j) - Obj(i+1,j-1,BorderType))*(Obj(i,j) - Obj(i+1,j-1,BorderType));
     GradParam(3) =  (Obj(i,j) - Obj(i-1,j-1,BorderType))*(Obj(i,j) - Obj(i-1,j-1,BorderType));
     Total = GradParam.total();
     if (Total > FLOAT_EPSILON)
     { 
        for (i = 0; i < Ndir; i++) GradParam(i) /= Total;
        for (i = 0; i < Ndir; i++) Val = GradDir(i)*GradParam(i);
     }  
     return Val/2.;
  }		    
  inline float markov_val4(Ifloat &Ima, int i, int j)
  {
      double p = MarkovPowerParam;
      double dp=p-1;
      int N = 4;
      dblarray Tab(N);
      Tab(0) = Ima(i,j) - Ima(i-1,j,BorderType);
      Tab(1) = Ima(i,j) - Ima(i,j-1,BorderType);
      Tab(2) = Ima(i,j) - Ima(i+1,j,BorderType);
      Tab(3) = Ima(i,j) - Ima(i,j+1,BorderType);
      double  Val = 0;
      for (int i=0; i < N; i++) Val += p*sgn(Tab(i))*pow(ABS(Tab(i)),dp);
      return (float) Val / (float) N;
  }
  inline float markov_val8(Ifloat &Ima, int i, int j)
  {
      double p = MarkovPowerParam;
      double dp=p-1;
      int N = 8;
      double S2 = 1./sqrt(2.);
      dblarray Tab(N);
      Tab(0) = Ima(i,j) - Ima(i-1,j,BorderType);
      Tab(1) = Ima(i,j) - Ima(i,j-1,BorderType);
      Tab(2) = S2*(Ima(i,j) - Ima(i-1,j-1,BorderType));
      Tab(3) = Ima(i,j) - Ima(i+1,j,BorderType);
      Tab(4) = Ima(i,j) - Ima(i,j+1,BorderType);
      Tab(5) = S2*(Ima(i,j) - Ima(i+1,j+1,BorderType));
      Tab(6) = S2*(Ima(i,j) - Ima(i-1,j+1,BorderType));
      Tab(7) = S2*(Ima(i,j) - Ima(i+1,j-1,BorderType));
      double  Val = 0;
      for (int i=0; i < N; i++) Val += p*sgn(Tab(i))*pow(ABS(Tab(i)),dp);
      return (float) Val / (float) (4.+4*S2);
  }
  public:   
    inline float ernergy_val(Ifloat &Obj, int i, int j)
    {
       float Val = Obj(i,j);
       return Val;
    }
    inline float laplacian_val(Ifloat &Obj, int i, int j)
    {
       float Val = Obj(i,j) - 0.25*( Obj(i+1,j,BorderType)
		    +  Obj(i-1,j,BorderType) +  Obj(i,j+1,BorderType) 
		    +  Obj(i,j-1,BorderType));
       return Val;
    }
    RegulIma () {GradOperType = OPER_MARKOV_8;BorderType=I_CONT;
                 MarkovPowerParam=1.1;}
    ~RegulIma () {}	
    type_oper_grad GradOperType;
    type_border BorderType;
    float MarkovPowerParam;
    inline void ima_regul(Ifloat &Obj, Ifloat &Grad, float Regul)
    {
       int Nl = Obj.nl();
       int Nc = Obj.nc();
       int i,j;
       if ((Grad.nl() != Nl)  || (Grad.nc() != Nc)) Grad.resize(Nl,Nc);

       if (Regul > 0)
       {
         switch(GradOperType)
         {
          case OPER_LAPLACIAN:
	     for (i = 0; i < Nl; i ++)
             for (j = 0; j < Nc; j ++) Grad(i,j) = Regul*laplacian_val(Obj,i,j);
             break;
	  case OPER_MARKOV_4:
             for (i = 0; i < Nl; i ++)
             for (j = 0; j < Nc; j ++) Grad(i,j) =  Regul*markov_val4(Obj,i,j); 
             break;
         case OPER_MARKOV_8:
             for (i = 0; i < Nl; i ++)
             for (j = 0; j < Nc; j ++) Grad(i,j) =  Regul*markov_val8(Obj,i,j); 
             break;
	 case OPER_DIR_LAPLACIAN:
             for (i = 0; i < Nl; i ++)
             for (j = 0; j < Nc; j ++) Grad(i,j) =  Regul*laplacian_val_dir(Obj,i,j); 
             break;
         case OPER_ENERGY:
             for (i = 0; i < Nl; i ++)
             for (j = 0; j < Nc; j ++) Grad(i,j) =  Regul*ernergy_val(Obj,i,j); 
             break;
	 }
      }
      else Grad.init();
   }
   inline void obj_regul(Ifloat &Obj, Ifloat &Grad, float Regul)
   {
       int Nl = Obj.nl();
       int Nc = Obj.nc();
       Ifloat ImaAux(Nl,Nc,"aux");
       ima_regul(Obj,ImaAux,Regul);
       Grad -= ImaAux;
   }
};



#endif
