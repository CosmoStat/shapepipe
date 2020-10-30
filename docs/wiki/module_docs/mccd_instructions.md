[Home](../shapepipe.md) | [Modules](../module_docs.md)

# MCCD instructions

One of the main differences in the handling of files in the pipeline context of
the MCCD and PSFEx modelling methods is that the first one has one model per
exposure while the latter has one model per exposure per CCD. 

More details on the MCCD package can be found in its 
[repository](https://github.com/CosmoStat/mccd) and in its
[documentation](https://cosmostat.github.io/mccd/).


#### MCCD configuration file

The MCCD package has its own configuration file. The default parameters can
be found on the MCCD repository
[here](https://github.com/CosmoStat/mccd/blob/master/config_MCCD.ini) or in
the example folder. When using the MCCD package integrated into the shapepipe
there is no need to change the sections ``[INPUTS], [VALIDATION]`` as the
pipeline will only use the other configuration sections.

Information about MCCD-specific parameters can be found directly in the 
MCCD config file. More information can be found on the documentation page
and on the MCCD article.

## MCCD: preprocessing

---
> Input runner: ``setools_runner``  
> Runner: ``mccd_preprocessing_runner.py``  
> Parameter section: ``[MCCD]``  
---

**Input:** it can take as input the output from ``setools_runner``
and runs in serial.

**Comment:** This runner is in charge of merging and do some preprocessing on
the stars becasuse the training of the MCCD model needs to use the whole
exposure.

Depending on the parameter ``MODE`` on the ``[MCCD]`` parameters, there could
be several options.

- ``MODE: FIT`` There should be only one type of file on the inputs.
Meaning only one pattern and one extension on the ``[FILE]`` parameters, that
should correspond to the training stars.

- ``MODE: FIT_VALIDATION`` There should be two type of files on the inputs on
``[FILE]``, the first one for the training stars and the second one for the
validation stars. If only one type of file is given at the input, it will be
taken as it was for training.

- ``MODE: VALIDATION`` In this case only one type of file is expected in 
``[FILE]``, and it should be the validation stars.

**Output:** the files are saved in ``.fits`` format with two name
patterns. ``train_star_selection`` for the training catalog and
``test_star_selection`` for the validation catalog. The patterns are followed
by the separator, generally ``-``, and the exposure id number.


## MCCD: fit only

---
> Input runner: ``mccd_preprocessing_runner``  
> Runner: ``mccd_fit_runner.py``  
> Parameter section: ``[MCCD]``  
---

**Input:** takes the output from ``mccd_preprocessing_runner`` and trains the 
MCCD model in parallel. Needs as input the ``train_star_selection``.

**Comment:** The parameter ``MODE`` on the ``[MCCD]`` must be set to ``FIT``.

**Output:** the trained models that are saved in ``.npy`` format and by default
by the pattern ``fitted_model`` followed by a separator, typically ``-``, and
the exposure id number.


## MCCD: validation only

---
> Runner: ``mccd_val_runner.py``  
> Parameter section: ``[MCCD]``  
---

**Input:** takes as first input the ``fitted_model`` from ``mccd_fit_runner``
and as second input the ``test_star_selection`` from the preprocessing runner.

**Comment:** Naturally, this model validates the trained MCCD model with the
corresponding test stars.
The parameter ``MODE`` on the ``[MCCD]`` must be set to ``VALIDATION``.

**Output:** the validation star catalog which is saved under the pattern 
``validation_psf`` followed by the separator, typically ``-``, and
the exposure id number on a ``.fits`` file.
This file contains the recovered PSF at the position of the test stars, and 
also the test stars. The momements are computed and saved for both of them,
the PSFs and the stars.


## MCCD: fit and validation

---
> Input runner: ``mccd_preprocessing_runner``  
> Runner: ``mccd_fit_val_runner.py``  
> Parameter section: ``[MCCD]``
---

**Input:** takes the output from ``mccd_preprocessing_runner``. Trains the 
MCCD model in parallel with  the ``train_star_selection`` and then validates
it with the ``test_star_selection``.

**Comment:** The parameter ``MODE`` on the ``[MCCD]`` must be set to
``FIT_VALIDATION``.
Performs the task of fitting the models and then validating it with the test
stars.

**Output:** same as ``mccd_val_runner``.


## MCCD: interpolation

---
> Runner: ``mccd_interp_runner.py``  
> Parameter section: ``[MCCD_INTERP_RUNNER]``  
---

**Comment:** This runner mimics the psfex interpolation runner so that the
existent pipeline structure can be maintained without major changes.
The interpolation files that contain the interpolation position can be the
same ones as for PSFEx, meaning that there is no need to do a preprocessing
over them.

The MCCD models need to be already trained. So the
``mccd_fit_runner`` or the ``mccd_fit_val_runner`` should be launched before.

The parameter ``PSF_MODEL_DIR`` should be the directory containing all the 
trained PSF models that follow the ``PSF_MODEL_PATTERN`` with the separator
``PSF_MODEL_SEPARATOR``.


### Classic

**Comment:** The ``MODE`` parameter should be set up in ``CLASSIC``.

**Output:** ``.fits`` file containing the interpolated PSF at the requested
positions. If ``GET_SHAPES`` is set to ``True`` the measured shapes of the PSFs
are also stored in the output. The pattern followed is ``galaxy_psf``.

### Multi-epoch

**Comment:** The ``MODE`` parameter should be set up in ``MULTI-EPOCH``.

**Output:** ``.sqlite`` file containing the same information as in the 
``CLASSIC`` mode. The pattern followed is ``galaxy_psf``.


## MCCD: merge validation catalog

---
> Input runner: ``mccd_fit_val_runner``  
> Runner: ``mccd_merge_starcat_runner.py``  
> Parameter section: No parameter inputs  
---

**Input:** this module collects all the ``validation_psf`` files from 
``mccd_fit_val_runner`` or ``mccd_validation_runner``.

**Comment:** he module merges the input validation files into a single catalog
containing the shape information about the stars and the PSFs.
It also computes some statistics over the validation data and prints it on the
log file.

**Output:** A single ``.fits`` file containing the merged catalog under the
name ``full_starcat-0000000.fits``.

## MCCD: plot utilities

---
> Input runner: ``mccd_merge_starcat_runner``  
> Runner: ``mccd_plots_runner.py``  
> Parameter section: ``[MCCD_PLOTS]``  
---

**WARNING!** For using this module  and compute the rho statistics 
you need to install more dependencies.
The python packages [treecorr](https://github.com/rmjarvis/TreeCorr),
[Stile](https://github.com/msimet/Stile), and their dependencies
should be installed. 

**Stile bug:** There is a known bug on the Stile package.
Fist install the treecorr package by ``pip install treecorr ``.
Then there are two options to fix the problem. You can a fork of the Stile
package that includes the fix as follows:

```bash
git clone https://github.com/tobias-liaudat/Stile.git
cd Stile
python setup.py install
```

The other option is to download and fix the bug yourself and then install the
package. 
To fix the bug what one need to do is to clone the package
``git clone https://github.com/msimet/Stile.git``. Then, go to the file
``treecorr_utils.py`` and add one line:

```python
102: DEPRECIATED_fields = f.readline().split() # [TL] BUG FIX
103: fields = f.readline().split() 
```

Then just do:
```bash
cd Stile
python setup.py install
```

Everything should be fine.

---

**Input:** The merged catalog ``full_starcat-0000000.fits`` from the merge
runner.

**Comment:** Draw different plots like, meanshapes over the focal plane,
histograms of the shape errors or compute the rho statistics.

**Output:** Several plots in ``.pdf`` and/or ``.png``.

