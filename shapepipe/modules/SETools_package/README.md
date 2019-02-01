# SETools package documentation

## General

This package allow you to easly handle SExtractor FITS_LDAC catalog.  
At the moment only the masking is possible (thresholding on output parameters).

All the process is done using one *.setools config file.

## Create a mask

First create a new text file with the extension : *.setools.

Then you need to create a mask section :
```text
[MASK]
```
The default name is "mask_[1,2,...]".  
To change the name, proceed this way :
```text
[MASK:name_of_your_mask]
```
Now you can create your thresholding using those comparisons : `[<, >, <=, >=, ==, !=]`  
Basic operation between parameters/numbers are possible : `[+, -, *, /]`  
There is also some functions : `[sqrt, pow, mean, median, mode]`  
To use a function proceed this way :
```text
mode(FWHM_IMAGE * 4.2) < 12
```
**Warning :**  
At the moment it is NOT possible to make operation with a function :
```text
mean(X_IMAGE + 2) + 3 < Y_IMAGE
```
You have to proceed this way instead :
```text
mean(X_IMAGE + 2) < Y_IMAGE - 3
```



Once this config file is complete, we need to update some variable in the package_config_file.cfg  :
```text
DEFAULT_FILENAME = star_galaxy_selection.setools      -> name of your .setools file
```
(The rest of the parametrisation is done as usual)

Now the package can be run.


## Example

*.setools file to make a star/galaxy selection :

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

- Add more functions : `[log, e, ...]`
- Add more functionalities :
  - Plot
  - Histogram
  - Statistic
  - ...
- Add support of other catalog format

- Make it compatible with other Astromatic software
