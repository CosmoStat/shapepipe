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
**    Date:  5/03/2000 
**    
**    File:  Filter.h
**
**    Modification history :
**
*******************************************************************************
**
**    DESCRIPTION : This class defines a filter bank
**    -----------   It can be created by different way:
**                 
**                  1) By using a predefined filter bank 
**                     ex: F_MALLAT_7_9,F_DAUBE_4, ...
**                  2) By using the file name if the file which
**                     contains the filter bank values to use.
**                     The format of the file is wvf format.
**                     see the Bath Wavelet Warehouse for more information:
**                          http://dmsun4.bath.ac.uk/wavelets/warehouse.html
**                    File Format (.wvf)
**                    General format to describe wavelet filter coefficients.
**
**                    [Range low] [Range high]
**                    [Analysis LP filter coefficients]
**                    .
**                    .
**                    .
**                    [Range low] [Range high]
**                    [Synthesis LP filter coefficients]
**                    .
**                    .
**                    .
**           
**                  3) By using the default user file name: mr1.wvf
**                     In this case, the file name must be in the 
**                     current directory
**                  4) By using the environment variable: CEA_FILTER
**                     
******************************************************************************/

#ifndef _DESIGN_FILTER_H_
#define _DESIGN_FILTER_H_

#include "GlobalInc.h"

// Type of filters (only the first 5 are visible for the users).
// F_MALLAT_9_7 ==>  analysis 9, synthesis 7  
#define NBR_SB_FILTER 14
enum type_sb_filter {SB_UNKNOWN, 
                     F_MALLAT_7_9, 
                     F_DAUBE_4, 
                     F_BI2HAAR, 
                     F_BI4HAAR,
                     F_ODEGARD_7_9,
		     F_5_3,
		     F_LEMARIE_1,
		     F_LEMARIE_3,
		     F_LEMARIE_5,
                     F_USER,
                     F_HAAR,
		     F_3_5,
		     F_4_4,
		     F_5_3_DIV,
		     F_MALLAT_9_7};

#define DEF_SB_FILTER F_MALLAT_7_9            // Default user filter 
#define DEF_USER_FILTER_FILE_NAME  "mr1.wvf"  // Default user filter file name
#define USER_FILER_FILE_NAME    "CEA_FILTER"  // Environment variable name
                                              // containing the filter file name

extern char *UserFilterFileName;

const char * StringSBFilter (type_sb_filter  type);
void sb_usage(type_sb_filter Filter);
type_sb_filter get_filter_bank(char *UserArg);

// *********  1D filter bank one step transform

enum sb_type_norm {NORM_L1, NORM_L2};
#define DEF_SB_NORM NORM_L1  // type of normalization

class FilterAnaSynt
{
    sb_type_norm TypeNorm;
    type_sb_filter TypeFilter;
    float *Analysis, *Synthesis;
    int Size_Ana, Size_Synt;
    int Start_Ana, Start_Synt;
    void read_from_file(char *FileName);
            // read the filter H0 and H1 from a file
            // File format is:
            //  [Range low] [Range high]
            //  [Analysis LP filter coefficients]
            //  .
            //  .
            //  .
            //  [Range low] [Range high]
            //  [Synthesis LP filter coefficients]
   void reset_param();
  public:
    FilterAnaSynt() {reset_param();}
    FilterAnaSynt(type_sb_filter T_Filter);
    FilterAnaSynt(char *FileName);
    ~FilterAnaSynt();

    type_sb_filter type_filter()  {return  TypeFilter;}
    sb_type_norm norm()   {return TypeNorm;}
    float *analysis()     {return Analysis;}
    float *synthesis()    {return Synthesis;}
    int size_analysis()   {return Size_Ana;}
    int size_synthesis()  {return Size_Synt;}
    int start_analysis()  {return Start_Ana;}
    int start_synthesis() {return Start_Synt;}

    void alloc(type_sb_filter T_Filter);
    Bool Verbose;
    char *FilterFileName; // User filter bank file name
                          // This field is used when T_Filter=F_USER
};


#endif
