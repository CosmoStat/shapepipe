/******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  06/08/99   
**    
**    File:  IM3D_IO.cc
**
*******************************************************************************
**
**    DESCRIPTION   input output routines for 3D signals
**    ----------- 
**                 
**    PARAMETRES    
**    ----------    
** 
**    RESULTS      
**    -------  
**
**
**   
******************************************************************************/

#include "IM3D_IO.h"
 
type_3d_format IO_3D_Format =  F3D_UNKNOWN;
type_data TypeInput3DData = UNKNOWN;

#ifdef USELM
extern DemoLic MRDEMO;
#endif

/**************************************************************************/

type_3d_format  io_which_3d_format(const char *Format) 
{
    char Name_Up[MAXCHAR];
    char Name[MAXCHAR];
    const char * Fits = ".fit";
    const char * Gif = ".gif";
    const char *Jpg = ".jpg"; 
    const char *Tif = ".tif";  
    const char *Tiff = ".tiff";  
    unsigned int i;

    strcpy(Name_Up, Format);
    type_3d_format Val_IO_3D_Format = F3D_UNKNOWN;
    strcpy(Name, Name_Up);
    for (i=0; i < strlen(Name); i++) Name[i] = tolower(Name[i]);
    if (strstr(Name,  Gif) != NULL)  Val_IO_3D_Format =  F3D_GIF;
    else if (strstr(Name,  Jpg) != NULL)  Val_IO_3D_Format = F3D_JPEG;
    else if (strstr(Name, Fits) != NULL)  Val_IO_3D_Format = F3D_FITS;
    else if (strstr(Name, Tif) != NULL)  Val_IO_3D_Format  = F3D_TIFF;
    else if (strstr(Name, Tiff) != NULL)  Val_IO_3D_Format = F3D_TIFF;
    return  Val_IO_3D_Format;
}

/**************************************************************************/

type_3d_format  io_detect_3dformat(const char *Format) 
{
    IO_3D_Format =  F3D_UNKNOWN;
    IO_3D_Format  = io_which_3d_format(Format);
    if (IO_3D_Format  ==  F3D_UNKNOWN) IO_3D_Format  =  DEFAULT_3D_FORMAT;
    return  IO_3D_Format;
}

/**************************************************************************/

void io_3d_set_format(const char *Format) 
{
     IO_3D_Format = io_which_3d_format(Format);
}
 
/****************************************************************************/

void io_select_3d_format(char *Name)
{
    type_3d_format AuxFormat = io_which_3d_format(Name);
    // 1) if name defines a format file, we use it
    // 2) else if a format is already define, we use it
    // 3) else we use the default format return by io_detect_3dformat
    if (AuxFormat != F3D_UNKNOWN)  IO_3D_Format = AuxFormat;
    else if (IO_3D_Format == F3D_UNKNOWN) 
                   IO_3D_Format = io_detect_3dformat (Name);
}

/****************************************************************************/

void io_3d_read_data(char *File_Name, fltarray & Dat, fitsstruct *Header)
{
   fitsstruct HD,*PtrHD;
   int Status;
    io_select_3d_format(File_Name);
    extern type_data TypeInputData;

    switch (IO_3D_Format)
    {
        case  F3D_GIF: 
#if IO_GIF_OK
	  io_read3d_gif(File_Name, Dat);
#else
              fprintf (stderr, "Error:GIF is not active\n");
              exit (-1); 
#endif
  	  break;
	case  F3D_JPEG: 
#if IO_JPEG_OK
	  io_read3d_jpeg(File_Name, Dat);
#else
              fprintf (stderr, "Error:JPEG is not active\n");
              exit (-1); 
#endif
  	  break;
 	case  F3D_FITS: 
#if IO_FITS_OK
             if (Header != NULL) PtrHD = Header;
             else PtrHD = &HD;           
             fits_read_fltarr(File_Name, Dat, PtrHD);
             switch(PtrHD->bitpix)
             {
                 case BP_BYTE: TypeInputData = T_BYTE; break;
                 case BP_SHORT: TypeInputData = T_SHORT; break;
                 case BP_INT: TypeInputData = T_INT; break;
                 case BP_FLOAT: TypeInputData = T_FLOAT; break;
                 case BP_DOUBLE: TypeInputData = T_DOUBLE; break;
                 default: break;
             }
#else
              fprintf (stderr, "Error: FITS is not active\n");
              exit (-1); 
#endif
	  break;
	case  F3D_TIFF: 
#if IO_TIFF_OK  
          Status = io_read3d_tiff(File_Name, Dat);
#else
              fprintf (stderr, "Error: TIFF is not active\n");
              exit (-1); 
#endif
	  break;
	default: cerr << "Error: unknown 3D format ... " << endl;
	         break;
   }
   if (Dat.naxis() != 3)
   {
      cerr << "Error: input data must be 3D data ... " << endl;
      exit(-1);
   }
#ifdef USELM
    MRDEMO.test(Dat.nx(), Dat.ny(), Dat.nz());
#endif
}

/****************************************************************************/


void io_3d_write_data(char *File_Name, fltarray & Dat, fitsstruct *Header)
{
    fitsstruct HD;
   int Status;
   io_select_3d_format(File_Name);

   //  cout << "Format = " << String3DFormat(IO_3D_Format) << endl;
   if (Dat.naxis() != 3)
   {
      cerr << "Error in io_3d_write_data: data must be 3D data ... " << endl;
      exit(-1);
   }
    switch (IO_3D_Format)
    {
        case  F3D_GIF: 
#if IO_GIF_OK
 	  io_write3d_gif(File_Name, Dat);
#else
              fprintf (stderr, "Error:GIF is not active\n");
              exit (-1); 
#endif
  	  break;
	case  F3D_JPEG: 
#if IO_JPEG_OK
   	      io_write3d_jpeg(File_Name, Dat);
#else
              fprintf (stderr, "Error:JPEG is not active\n");
              exit (-1); 
#endif
  	  break;
 	case  F3D_FITS: 
#if IO_FITS_OK
             if (Header != NULL) 
             {
                 Header->naxis = 3;
                 Header->width =  Dat.nx();
                 Header->height = Dat.ny();
		 Header->filename = strdup(File_Name);
                 (Header->TabAxis)[0] = Dat.nx();
                 (Header->TabAxis)[1] = Dat.ny();
                 (Header->TabAxis)[2] = Dat.nz();
                 Header->npix = Dat.n_elem();
                 fits_write_fltarr(File_Name, Dat, Header);
             }
             else 
             {
                  initfield(&HD);
                  HD.naxis = 3;
                  HD.bitpix = -32;
                  HD.width =  Dat.nx();
                  HD.height = Dat.ny();
                  HD.filename = strdup(File_Name);
                  (HD.TabAxis)[0] = Dat.nx();
                  (HD.TabAxis)[1] = Dat.ny();
                  (HD.TabAxis)[2] = Dat.nz();
                  HD.npix = Dat.n_elem();
                  fits_write_fltarr(File_Name, Dat, &HD);
              }
 #else
              fprintf (stderr, "Error: FITS is not active\n");
              exit (-1); 
#endif
	  break;
	case  F3D_TIFF: 
#if IO_TIFF_OK  
          Status = io_write3d_tiff(File_Name, Dat);
#else
              fprintf (stderr, "Error: TIFF is not active\n");
              exit (-1); 
#endif
	  break;
  	default: cerr << "Error: unknown 3D format ... " << endl;
	         break;
   }
}


/****************************************************************************/

void IO3DInfoData::init_reading(char *Name, Bool ReadInPseudo)
{   
    int Status;
    io_select_3d_format(Name);
    Format = IO_3D_Format;
    strcpy(File_Name, Name);
    PseudoColRead = ReadInPseudo;

    // cout << "init_reading " << File_Name << " " << String3DFormat(Format) <<  endl;

    switch (Format)
    {
        case F3D_FITS:
#if IO_FITS_OK
              PtrFits = new fitsstruct;
	      initfield(PtrFits);
              fits_read_header(File_Name, PtrFits);
	      if (PtrFits->naxis > 0) Nc = (PtrFits->TabAxis)[0];
	      if (PtrFits->naxis > 1) Nl = (PtrFits->TabAxis)[1];
	      if (PtrFits->naxis > 2) Nima  = (PtrFits->TabAxis)[2];
	      switch (PtrFits->bitpix)
	      {
	         case BP_BYTE: Type = T_BYTE; break;
		 case BP_SHORT: Type = T_SHORT; break;        
                 case BP_INT: Type = T_INT; break;           
                 case BP_FLOAT: Type = T_FLOAT; break;
                 case BP_DOUBLE: Type = T_DOUBLE; break;
		 default:
		    cerr << "Error: unknown format ... " << endl;
		    exit(-1);
		    break;
	      }
#else
              fprintf (stderr, "Error: FITS is not active\n");
              exit (-1); 
#endif
              break;
         case F3D_GIF: 
#if IO_GIF_OK
              readgif(File_Name, Pinfo);
              Nl = Pinfo.h;
              Nc = Pinfo.w;
              PseudoCol = True;
              Type = T_BYTE;
#else
              fprintf (stderr, "Error: GIF is not active\n");
              exit (-1); 
#endif	       
               break;    
	case  F3D_TIFF: 
#if IO_TIFF_OK  
              Status = io_read3d_tiff(File_Name, &Pinfo);
              Nl = Pinfo.h;
              Nc = Pinfo.w;
              FullCol = True;
              Type = T_BYTE;
#else
              fprintf (stderr, "Error: TIFF is not active\n");
              exit (-1); 
#endif
	  break;
 	 case F3D_JPEG:
#if IO_JPEG_OK
             {
               if (ReadInPseudo == False)
                     RJEP.startread(File_Name, JPEG_READ_IN_RGB);
               else  RJEP.startread(File_Name, JPEG_READ_IN_PSEUDO);
               Nl =  RJEP.ny();
               Nc =  RJEP.nx();
               if (ReadInPseudo == False) FullCol = True;
               else PseudoCol = True;
               Type = T_BYTE;
             }
#else
              fprintf (stderr, "Error: JPEG is not active\n");
              exit (-1); 
#endif	       
               break;	 
         default:
              fprintf (stderr, "Error: bad image format. cannot read ...\n");
              exit (-1);
              break;
    }
    TypeInput3DData = Type;
#ifdef USELM
    MRDEMO.test(Nc, Nl);
#endif
}

/**************************************************************************/

void IO3DInfoData::read_pseudo_block(fltarray &Data, int Indi, int Indj)
{
   byte R,G,B;
   int i,j, Pix;
   byte *ptr = Pinfo.pic;
   int First = Indi*Nc+Indj;
   for (i=0; i < Data.ny(); i++)
   for (j=0; j < Data.nx(); j++) 
   {
      Pix = (int) ptr[First+i*Nc+j];
      IO_RGB.pseudo_to_color(Pix, R,G,B);
      
      Data(j,i,0) = R;
      Data(j,i,1) = G;
      Data(j,i,2) = B;
   }
}

/**************************************************************************/

void IO3DInfoData::read_col_block(intarray &Data, int Indi, int Indj)
{
   int i,j,Pix,First = Indi*Nc+Indj;
   byte *Ptr = Pinfo.pic;
   for (i=0; i < Data.ny(); i++)
   for (j=0; j < Data.nx(); j++) 
   {
      Pix = (First+i*Nc+j)*3;
      Data(j,i,0) = Ptr[Pix];
      Data(j,i,1) = Ptr[Pix+1];
      Data(j,i,2) = Ptr[Pix+2];  
  }
}

/**************************************************************************/

void IO3DInfoData::read_col_block(fltarray &Data, int Indi, int Indj)
{
   int i,j, Pix;
   int First = Indi*Nc+Indj;
   byte *Ptr = Pinfo.pic;
   for (i=0; i < Data.ny(); i++)
   for (j=0; j < Data.nx(); j++) 
   {
      Pix =  (First+i*Nc+j)*3;
      Data(j,i,0) = Ptr[Pix];
      Data(j,i,1) = Ptr[Pix+1];
      Data(j,i,2) = Ptr[Pix+2];
  } 
}

/**************************************************************************/

void IO3DInfoData::read_pseudo_block(intarray &Data, int Indi, int Indj)
{
   byte R,G,B;
   int i,j, Pix;
   byte *ptr = Pinfo.pic;
   int First = Indi*Nc+Indj;
   for (i=0; i < Data.ny(); i++)
   for (j=0; j < Data.nx(); j++) 
   {
      Pix = (int) ptr[First+i*Nc+j];
      IO_RGB.pseudo_to_color(Pix, R,G,B);      
      Data(j,i,0) = R;
      Data(j,i,1) = G;
      Data(j,i,2) = B;
  }
}

/**************************************************************************/

void IO3DInfoData::write_pseudo_block(intarray &Data, int Indi, int Indj)
{
  int i,j,First = Indi * Nc + Indj;
  byte *Ptr = Pinfo.pic;
  for (i=0; i < Data.ny(); i++)
  for (j=0; j < Data.nx(); j++)
  {
     byte ValR =  float_to_byte(Data(j,i,0));
     byte ValG =  float_to_byte(Data(j,i,1));
     byte ValB =  float_to_byte(Data(j,i,2));
     Ptr [First+i*Nc+j] =  IO_RGB.color_to_pseudo(ValR, ValG, ValB);
  }
}

/**************************************************************************/

void IO3DInfoData::write_pseudo_block(fltarray &Data, int Indi, int Indj)
{
  int i,j,First = Indi * Nc + Indj;
  byte *Ptr = Pinfo.pic;
  for (i=0; i < Data.ny(); i++)
  for (j=0; j < Data.nx(); j++)
  {
     byte ValR =  float_to_byte(Data(j,i,0));
     byte ValG =  float_to_byte(Data(j,i,1));
     byte ValB =  float_to_byte(Data(j,i,2));
     Ptr [First+i*Nc+j] =  IO_RGB.color_to_pseudo(ValR, ValG, ValB);
  }
}

/**************************************************************************/

void IO3DInfoData::write_col_block(intarray &Data, int Indi, int Indj)
{
   int i,j, Pix;
   int First = Indi*Nc+Indj;
   byte *Ptr = Pinfo.pic;
   for (i=0; i < Data.ny(); i++)
   for (j=0; j < Data.nx(); j++) 
   {
      Pix = (First+i*Nc+j)*3;
      Ptr[Pix] = float_to_byte(Data(j,i,0));
      Ptr[Pix+1] = float_to_byte(Data(j,i,1));
      Ptr[Pix+2] = float_to_byte(Data(j,i,2));
  }
}

/**************************************************************************/

void IO3DInfoData::write_col_block(fltarray &Data, int Indi, int Indj)
{
   int i,j, Pix;
   int First = Indi*Nc+Indj;
   byte *Ptr = Pinfo.pic;
   for (i=0; i < Data.ny(); i++)
   for (j=0; j < Data.nx(); j++) 
   {
      Pix = (First+i*Nc+j)*3;
      Ptr[Pix] = float_to_byte(Data(j,i,0));
      Ptr[Pix+1] = float_to_byte(Data(j,i,1));
      Ptr[Pix+2] = float_to_byte(Data(j,i,2));
   }
}

/**************************************************************************/

void IO3DInfoData::init_writing(char *fname, int Nlima, int Ncima, int NbrIma)
{
    int i;
    strcpy(File_Name, fname);
    io_select_3d_format(File_Name);
    Format = IO_3D_Format;
    Nl = Pinfo.h = Nlima;
    Nc = Pinfo.w = Ncima; 
    Nima = NbrIma;
    Pinfo.comment = new char [1];
    Pinfo.comment = '\0';
    Pinfo.colType = F_FULLCOLOR;

    switch (Format)
    {
        case F3D_FITS:
#if IO_FITS_OK
              if (PtrFits == NULL)
              {
                  cout << "Error: fits structure is not allocated ... " << endl;
                  exit(-1);
              }
              fits_write_header(File_Name, PtrFits);
#else
              fprintf (stderr, "Error: FITS is not active\n");
              exit (-1); 
#endif
              break;
         case F3D_GIF:
#if IO_GIF_OK
               PseudoCol = True;
               Pinfo.pic = new byte[Nl*Nc];
               Nima = 1;
#else
            fprintf (stderr, "Error:GIF is not active\n");
            exit (-1); 
#endif            
               break; 
 	 case F3D_TIFF:
#if  IO_TIFF_OK
               FullCol = True;
               Pinfo.pic = new byte[Nl*Nc*3];
               Nima = 3;
#else
            fprintf (stderr, "Error:TIFF is not active\n");
            exit (-1); 
#endif     
            break;
	 case F3D_JPEG:
#if IO_JPEG_OK
               FullCol = True;
               Nima = 3;
               WJEP.Color = True;
#else
            fprintf (stderr, "Error: JPEG is not active\n");
            exit (-1); 
#endif               
               break;
         default:
              fprintf (stderr, "Error: bad image format. cannot read ...\n");
              exit (-1);
              break;
    }

    // allocated LUT for pseudo color images
    if ((PseudoCol == True) || (FullCol == True))
    {
      IO_RGB.alloc_color();
      for (i=0; i< IO_RGB.lut_size(); i++)
      {
        // cout << IO_RGB.red(i) << " " << IO_RGB.green(i) << " " << IO_RGB.blue(i) << endl;
        Pinfo.r[i] =  IO_RGB.red(i);
        Pinfo.g[i] =  IO_RGB.green(i);
        Pinfo.b[i] =  IO_RGB.blue(i);
      }
    }
#if IO_JPEG_OK
   if (Format == F3D_JPEG) WJEP.startwrite(File_Name, Nl, Nc);
#endif               

}

/**************************************************************************/

void IO3DInfoData::end_reading()
{
       FullCol = PseudoCol = False;
       PseudoColRead=False;
       Nl = Nc = 0;
       Nima = 1;
       Format = F3D_UNKNOWN;
       Type = UNKNOWN;
       File_Name[0]='\0';
       Pinfo.type = 0;
       Pinfo.comment = new char [1];
       Pinfo.comment = '\0';
       Pinfo.colType = F_FULLCOLOR;
       if (Pinfo.pic  != NULL) delete [] Pinfo.pic;
       PtrFits = NULL;
       Pinfo.pic = NULL;
       AllocBufferSize=0;
       NbrReadLine = LastBlockLineNbr = 0;
       switch (Format)
       {
        case F3D_FITS:
        case F3D_GIF: 
                break;  
	 case F3D_JPEG:
#if IO_JPEG_OK
              RJEP.end_read();
#else
              fprintf (stderr, "Error: JPEG is not active\n");
              exit (-1); 
#endif
              break;	
         case F3D_TIFF: 
	  break;
         default:
              fprintf (stderr, "Error: bad image format. cannot read ...\n");
              exit (-1);
              break;
    }
}

/**************************************************************************/

void IO3DInfoData::end_writing()
{
    switch (Format)
    {
        case F3D_FITS:
                 break;  
        case F3D_GIF: 
#if IO_GIF_OK
              writegif(File_Name, Pinfo);
#else
              fprintf (stderr, "Error: GIF is not active\n");
              exit (-1); 
#endif
                break;  
	 case F3D_JPEG:
#if IO_JPEG_OK
              WJEP.end_write();
#else
              fprintf (stderr, "Error: JPEG is not active\n");
              exit (-1); 
#endif
              break;	
         case F3D_TIFF: 
#if IO_TIFF_OK  
              io_write3d_tiff(File_Name, &Pinfo);
#else
              fprintf (stderr, "Error: TIFF is not active\n");
              exit (-1); 
#endif
	  break;
         default:
              fprintf (stderr, "Error: bad image format. cannot read ...\n");
              exit (-1);
              break;
    }
}

/**************************************************************************/

void io_3d_read_block_ima(char *File_Name, fltarray &Image, 
                       int Indi, int Indj, IO3DInfoData &InfoDat, Bool NoBscale)
{
    if ((Image.ny() + Indi > InfoDat.Nl) ||
        (Image.nx() + Indj > InfoDat.Nc))
    {
       cerr << "Error: this block cannot be extracted from file: "  <<  File_Name << endl;
       cerr << "       Xs  = " <<        Indj << " Ys  = " <<  Indi       << endl;
       cerr << "       Ncb = " <<  Image.nx() << " Nlb = " <<  Image.ny() << endl;
       cerr << "       Nc  = " <<  InfoDat.Nc << " Nl  = " <<  InfoDat.Nl << endl;
    }
    
    switch (InfoDat.Format)
    {
        case F3D_FITS:
#if IO_FITS_OK
              fits_read_block(File_Name, Image, Indi, Indj, NoBscale);
#else
              fprintf (stderr, "Error: FITS is not active\n");
              exit (-1); 
#endif
              break;
         case F3D_GIF: 
#if IO_GIF_OK
              InfoDat.read_pseudo_block(Image, Indi, Indj);
#else
              fprintf (stderr, "Error: GIF is not active\n");
              exit (-1); 
#endif
                break;  
	case F3D_JPEG:   
#if IO_JPEG_OK
              // Test if the block is in memory
              if (Indi < InfoDat.NbrReadLine-InfoDat.LastBlockLineNbr)
              {
                 cout << "Error: block alredy read  " <<  Indi << endl;                                  
                 exit(-1);
              }
              else if (Indi >= InfoDat.NbrReadLine) // Block not in memory
              {
                 int BlockNl = Image.ny();
                 if (Indi >= InfoDat.NbrReadLine+BlockNl) 
                 {          
                     // next block cannot contains the needed lines
                     cout << "Error: next block to small to read line " << Indi << endl;                                  
                     exit(-1);
                 }
                 Bool RGBIma = (InfoDat.PseudoColRead == False) ? True: False;

                 // Allocate InfoDat.Pinfo if it is necessary
                 unsigned int BuffSize=BlockNl*InfoDat.Nc;
                 if (RGBIma==True) BuffSize *= 3;
                 // cout << "BlockNl  = " << BlockNl  <<   " Nc = " << InfoDat.Nc  <<  "BuffSize = " << BuffSize << endl;
                 if (InfoDat.AllocBufferSize < BuffSize)
                 {
                     if (InfoDat.AllocBufferSize > 0) delete [] InfoDat.Pinfo.pic;
                     InfoDat.Pinfo.pic = new byte [BuffSize];
                     InfoDat.AllocBufferSize = BuffSize;
                 }
                 
                 // Read next buffer
                 InfoDat.RJEP.readbuff(BlockNl, &(InfoDat.Pinfo), RGBIma);
                 InfoDat.NbrReadLine += BlockNl;
                 InfoDat.LastBlockLineNbr = BlockNl;
              }
	      if (InfoDat.PseudoColRead == False)
                   InfoDat.read_col_block(Image, 0, Indj);
              else InfoDat.read_pseudo_block(Image, 0, Indj);
 	      // if (InfoDat.PseudoColRead == True) cout << "Pseudo " << endl;
#else
              fprintf (stderr, "Error: JPEG is not active\n");
              exit (-1); 
#endif	       
               break;
         case F3D_TIFF: 
#if IO_TIFF_OK  
              InfoDat.read_col_block(Image, Indi, Indj);
#else
              fprintf (stderr, "Error: TIFF is not active\n");
              exit (-1); 
#endif
	  break;
        default:
              fprintf (stderr, "Error: bad image format. cannot read ...\n");
              exit (-1);
              break;
    }
}

/**************************************************************************/

void io_3d_write_block_ima(char *File_Name, fltarray &Image, 
                       int Indi, int Indj, IO3DInfoData &InfoDat, Bool NoBscale)
{
    Bool RGBIma = (InfoDat.PseudoColRead == False) ? True: False;
    int BlockNl = Image.ny();
    int BlockNc = Image.nx();
    if ((Image.ny() + Indi > InfoDat.Nl) ||
        (Image.nx() + Indj > InfoDat.Nc))
    {
       cerr << "Error: this block cannot be inserted in file: "  <<  File_Name << endl;
       cerr << "       Xs  = " <<        Indj << " Ys  = " <<  Indi       << endl;
       cerr << "       Ncb = " <<  Image.nx() << " Nlb = " <<  Image.ny() << endl;
       cerr << "       Nc  = " <<  InfoDat.Nc << " Nl  = " <<  InfoDat.Nl << endl;
    }
    
    switch (InfoDat.Format)
    { 
        case F3D_FITS:
#if IO_FITS_OK
              fits_write_block(File_Name, Image, Indi, Indj, NoBscale);
#else
              fprintf (stderr, "Error: FITS is not active\n");
              exit (-1); 
#endif
              break;
         case F3D_GIF: 
#if IO_GIF_OK
              InfoDat.write_pseudo_block(Image, Indi, Indj);
#else
              fprintf (stderr, "Error: GIF is not active\n");
              exit (-1); 
#endif
                break;  
	case F3D_JPEG:   
#if IO_JPEG_OK
             // Test if the block is in memory
              if (Indi < InfoDat.NbrWriteLine)
              {
                 cout << "Error: block alredy written  " <<  Indi << endl;                                  
                 exit(-1);
              }
              else if (Indi >= InfoDat.NbrWriteLine+BlockNc) 
              {
                 cout << "Error: cannot skip blocks  " <<  Indi << endl;                                  
                 exit(-1);
              }
              else
              {
                InfoDat.PseudoColRead = False;
                if (Indj == 0)  // New bloc line
                {
                    // Allocate InfoDat.Pinfo if it is necessary
                    unsigned int BuffSize=BlockNl*InfoDat.Nc;
                    if (RGBIma==True) BuffSize *= 3;
                 // cout << "BlockNl  = " << BlockNl  <<   " Nc = " << InfoDat.Nc  <<  "BuffSize = " << BuffSize << endl;
                    if (InfoDat.AllocBufferSize < BuffSize)
                    {
                       if (InfoDat.AllocBufferSize > 0) delete [] InfoDat.Pinfo.pic;
                       InfoDat.Pinfo.pic = new byte [BuffSize];
                       InfoDat.AllocBufferSize = BuffSize;
                    }
                 }
              }
	      if (InfoDat.PseudoColRead == False)
                   InfoDat.write_col_block(Image, 0, Indj);
              else InfoDat.write_pseudo_block(Image, 0, Indj);

 	     // if (InfoDat.PseudoColRead == False) cout << "FCOL " << endl;

             // Write next buffer
             if (Indj + BlockNc >= InfoDat.Nc)
             {
                 InfoDat.WJEP.writebuff(BlockNl, &(InfoDat.Pinfo), RGBIma);
                 InfoDat.NbrWriteLine += BlockNl;
                 InfoDat.LastBlockLineNbr = BlockNl;
             }
#else
              fprintf (stderr, "Error: JPEG is not active\n");
              exit (-1); 
#endif	       
               break;
         case F3D_TIFF: 
#if IO_TIFF_OK  
              InfoDat.write_col_block(Image, Indi, Indj);
#else
              fprintf (stderr, "Error: TIFF is not active\n");
              exit (-1); 
#endif
	  break;

        default:
              fprintf (stderr, "Error: bad image format. cannot read ...\n");
              exit (-1);
              break;
    }
}

/**************************************************************************/
