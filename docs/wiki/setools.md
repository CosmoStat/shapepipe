[Home](../shapepipe.md) | [Modules](../module_docs.md)

# SETools package documentation

## General

This package allow you to easily handle files in SExtractor FITS_LDAC catalog format.
At the moment only masking is possible (thresholding on output parameters, thereby
allowing to create selected sub-samples of the data).

## Mask config file

The processing is controlled via a .setools config file.
Create a new text file with the extension : <name>.setools.

### Mask processing

Create a mask section:
```text
[MASK]
```
The default name in the code is "mask[1,2,...]".  
To change the name, proceed this way :
```text
[MASK:name_of_your_mask]
```
Now you can create your thresholding using those comparisons : `[<, >, <=, >=, ==, !=]`  
Basic operations between parameters/numbers are possible: `[+, -, *, /]`  
There is also some functions implemented: `[mean, median, mode]`  
To use a function proceed in this way:
```text
mode(FWHM_IMAGE * 4.2) < 12
```
**Warning :**  
At the moment it is NOT possible to perform operation with a function:
```text
mean(X_IMAGE + 2) + 3 < Y_IMAGE
```
You have to proceed this way instead:
```text
mean(X_IMAGE + 2) < Y_IMAGE - 3
```

### Post-processing.

Optional sections for statistics and plotting can be created.

#### Statistics

For each mask defined in the .setools file, a statistics section can be added:

[STAT:mask1]
...
[STAT:mask2]
...
...

The possible entries in each statistics section are:
- NUMBER_OBJ
  Number of objects in selected sub-sample

For each section, an output ascii file BASE_OUTPUT_DIR/run_<date.time>/st


### Update package config file

Once this config file is complete, we need to update some variables in the package_config_file.cfg, e.g.:
```text
DEFAULT_FILENAME = star_galaxy_selection.setools      -> name of the <name>.setools file
```
(The rest of the parametrisation is done as usual)

Now the package can be run.


## Example

*.setools file to make a star/galaxy selection:

```text
[MASK:star_selection]
MAG_AUTO > 18.
MAG_AUTO < 23.
FLAGS == 0
IMAFLAGS_ISO == 0
FWHM_IMAGE - 0.05 <= mode(FWHM_IMAGE)
FWHM_IMAGE + 0.05 >= mode(FWHM_IMAGE)
CLASS_STAR != 0

[MASK:galaxy_selection]
MAG_AUTO > 20.
FLAGS == 0
IMAFLAGS_ISO == 0
```

## Future plan

- Add more functions : `[sqrt, pow, log, e, ...]`
- Add more functionalities :
  - Plot
  - Histogram
  - Statistic
  - ...
- Add support of other catalog format

- Make it compatible with other Astromatic software
