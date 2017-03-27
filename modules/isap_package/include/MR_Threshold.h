

#ifndef _MR_THRESHOLD_H_
#define _MR_THRESHOLD_H_

 
#define NBR_THRESHOLD 5

enum type_threshold {T_KSIGMA, T_FDR, T_UNIVERSAL,
                     T_SURE, T_MRSURE};

#define DEF_THRESHOLD T_KSIGMA

inline const char * StringThreshold (type_threshold type)
{
    switch (type)
    {
        case T_KSIGMA:
              return ("K-SigmaNoise Threshold");break;
        case T_UNIVERSAL:
              return ("Universal Threshold");break;
        case T_SURE:
              return ("SURE Threshold");break;
        case T_MRSURE:
              return ("Multiscale SURE Threshold");break;
        case T_FDR:
              return ("False Discovery Rate (FDR) Theshold");break;
   }
   return ("Error: bad type of filtering");
}

inline void get_threshold_usage(type_threshold type)
{
    fprintf(OUTMAN, "         [-C Coef_Detection_Method]\n");
    for (int i = 0; i < NBR_THRESHOLD; i++)
    fprintf(OUTMAN, "              %d: %s\n",i+1, StringThreshold( (type_threshold) i));
    fprintf(OUTMAN, "              default is %s.\n", StringThreshold(type));
}

#endif

