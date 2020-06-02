[Home](../shapepipe.md) | [Modules](../module_docs.md)

# SExtractor package documentation

## General

This package is a wrapper for the [Astromatic](https://www.astromatic.net/) [SExtractor](https://www.astromatic.net/software/sextractor) software.  
All the possible input file are handle. It use only one config file for the whole parametrisation.  
It is also possible to use values from the image header to set SExtractor input parameters.

## Run SExtractor package

**Note :**
- All the default files needed to run SExtractor are in `SExtractor_package/config/SExtractor_default`. If you are using other files you can add them here except for : default.sex and default.param.
- The default.sex as been arrange in a way to make the usage easier. It's better to not change this file.
- After a run of SExtractor_package all the files used by SExtractor can be find in :
run*****/logs/sextractor_config_inputs/.

### Config file

All the parametrisation is done in the package_config_file.cfg in the section : `[CODE]`, `[SEXTRACTOR_INPUT]` and `[SEXTRACTOR_OUTPUT]`.
```text
# ----------------------------- Code configuration ---------------------------
[CODE]

# Path to executable
EXEC_PATH = sextractor          -> SExtractor executable can be `sextractor' or
                                   `sex' depending on your personal build
# Name of executable configuration file
DEFAULT_FILENAME = NONE          -> No config file are provide
# Input file name formats
INPUT_FILENAME_FORMATS = ['CFIS.fits', 'weight.fits', 'flag.fits']    -> See section
                                                                         input files
# Output file name ending
OUTPUT_CATALOG_FILE_ENDING = .cat

WEIGHT = True                    -> If "True" weight file will be used
FLAG = True
PSF = False
ASSOC = False

# ---------------------------- *.sex configuration ---------------------------
[SEXTRACTOR_INPUT]               -> You can set here any input parameter. They will
                                    overwrite the value set in default.sex
DETECT_THRESH    10              -> You can set directly the value
DETECT_MINAREA    @3.141592*IQFINAL/0.186/2.*IQFINAL/0.186/2. -> You can make operation
                                                                 with header value by
                                                                 by starting the line
                                                                 with "@"
PIXEL_SCALE    0.186
SEEING_FWHM    @IQFINAL

GAIN    @GAIN
SATUR_LEVEL    @SATURATE

WEIGHT_TYPE    MAP_WEIGHT
WEIGHT_GAIN    Y

# --------------------------- *.param configuration --------------------------
[SEXTRACTOR_OUTPUT]              -> This section is the *.param file. You can copy/paste
                                    here the content of your *.param file
NUMBER
X_IMAGE
X_WORLD
Y_IMAGE
Y_WORLD

IMAFLAGS_ISO(1)

VIGNET(51,51)
```

All other parameters are set the usual way.

#### Input files

All the input files are provide in the section `[CODE]` to the key `INPUT_FILENAME_FORMATS` following this specific order :

 ```
INPUT_FILENAME_FORMATS = ['image.fits', 'weight.fits', 'flag.fits', 'psf.psf', assoc.list']
```
The corresponding type of files need to be set to `True` in the same section.
All input file provide in the section `[SEXTRACTOR_INPUT]` are ignored.

#### Value from image header

You can set value in the section `[SEXTRACTOR_INPUT]` using value from the header. It is also possible to make linear operation between them. To do so, the parameter value has to start with "@".
```text
GAIN @GAIN*12

DETECT_THRESH @ 3*4.56 + 78
```

### Run the package

Once the parametrisation is done the package is run the usual way.
