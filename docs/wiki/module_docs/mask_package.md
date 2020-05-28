[Home](../shapepipe.md) | [Modules](../module_docs.md)

# Mask package documentation

## General

Most of the code is a Python adaptation of a script from the THELI pipeline.   
The code use the GSC 2.2 as star catalog. Data are obtain with the CDSclient.

This package use *.reg template mask for spike and halo. They are then fitted on selected stars depending of their magnitude. It's also possible to mask a border all around the image.

## Create masks

All the mask properties are set in the config file *.mask.   
An example follow :
```text
# ----------------------------- Path ----------------------------------
[PROGRAM_PATH]

WW_PATH = ww                           -> Full path to WeightWatcher executable
                                          (WeightWatcher executable can be `weightwatcher' or
                                          `ww' depending on your personal build)
WW_CONFIG_FILE = default.ww            -> Path to WW config file *.ww
CDSCLIENT_PATH = findgsc2.2            -> Full path to CDSclient/findgsc2.2 executable


# ----------------------------- Border parameters ----------------------------------
[BORDER_PARAMETERS]

BORDER_MAKE = True                     -> If 'True' make this mask

BORDER_WIDTH = 300                     -> Width in pixels of the border to mask
BORDER_FLAG_VALUE = 4                  -> Flag value for this masking


# ----------------------------- Halo parameters ----------------------------------
[HALO_PARAMETERS]

HALO_MAKE = True

HALO_MASKMODEL_PATH = halo_mask.reg    -> Path to the DS9-region_file format template
                                          for halo
HALO_MAG_LIM = 13.                     -> Higher stars' magnitude to mask
HALO_SCALE_FACTOR = 0.05               -> see following section
HALO_FLAG_VALUE = 1
HALO_REG_FILE = halo.reg               -> Name of the region file created


# ----------------------------- Spike parameters ----------------------------------
[SPIKE_PARAMETERS]

SPIKE_MAKE = True

SPIKE_MASKMODEL_PATH = MEGAPRIME_star_i_13.8.reg   -> Path to the DS9-region_file format    
                                                      template for
SPIKE_MAG_LIM = 18.
SPIKE_SCALE_FACTOR = 0.3
SPIKE_FLAG_VALUE = 2
SPIKE_REG_FILE = spike.reg

# ---------------------------- Messier parameters ---------------------------------
[MESSIER_PARAMETERS]

MESSIER_MAKE = True

MESSIER_CAT_PATH = $HOME/ShapePipe/modules/mask_package/config/mask_default/Messier_catalog.npy
MESSIER_PIXEL_SCALE = 0.187
MESSIER_SIZE_PLUS = 1.
MESSIER_FLAG_VALUE = 8


# -------------------------------- External flag ----------------------------------
[EXTERNAL_FLAG]

EF_MAKE = True


# -------------------------- missing data parameters ------------------------------
[MD_PARAMETERS]

MD_MAKE = True

MD_THRESH_FLAG = 0.3
MD_THRESH_REMOVE = 0.75
MD_REMOVE = False




# ----------------------------- Other parameters ----------------------------------
[OTHER]

KEEP_REG_FILE = False                  -> If 'True' will keep temporary region files
KEEP_INDIVIDUAL_MASK = False           -> If 'True' will keep individual mask in fits
                                          format
TEMP_DIRECTORY = .temp/                -> Directory where temporary file are store
                                          It need to be provided even if you don't keep
                                          temporary files and must already exist
```

Once this config file is complete, we need to update some variable in the package_config_file.cfg  :
```text
[CODE]
DEFAULT_FILENAME = config.mask                            -> The mask config file to use
INPUTE_FILENAME_FORMATS = ['image.fits', 'weight.fits']   -> Follow this specific order
                                                             to provide files
```
(The rest of the parametrisation is done as usual)

Now the package can be run.

## How the template mask fitting is done ?

At the moment an empirical formula is used to fit the mask template to every stars depending of the magnitude. Here is the scaling function :
```text
scaling = 1 - [TYPE]_SCALE_FACTOR * (star_mag - [TYPE]_MAG_PIVOT)
```

## Tiles versus single exposures

For tiles, external flags do not exist. Make sure to set the keyword in the `[EXTERNAL_FLAG]` section to False:

```
EF_MAKE = True
```

Exposures are dealt with on a CCD-basis. The much reduced image size requires a smaller border mask, e.g in the `[BORDER_PARAMETERS]` section, set

```
BORDER_WIDTH = 50
```






## Future plan

Some idea to improve/re-built this package :
- Star catalog. Replace GSC by Gaia.
- Use another wrapper than CDSclient which provide too complex output really hard to handle.
- Avoid the .reg syntax (from DS9).
- Avoid the usage of WeightWatcher (WW) which create temporary file possibly problematic (size ~ 100Mo for CFIS images).
- Need weight image and do nothing with it only needed for WW.
- Fit automaticly mask on object without external parametrization.
- No possibility to use an external mask.
