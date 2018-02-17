# Python Scripts

These scripts are used to process pipeline products and do not require job
handling (*i.e.* they can be run in serial).

## Programs

* `cfis_field_select.py`
  Handle and select CFIS fields, pointings, tiles, weight images, for book keeping
* `scp_CFIS_cc.py`
  Copy files to and from cc@in2p3
* `cfis_get_coord_exposures.py`
  Print coordinates extracted from FITS exposure CFIS images 


## Modules/libraries


### cfis

* cfis.py
  CFIS-specific classes and routines


### generic

* stuff.py

## jupyer notebooks

* `cfis_field_select.ipybn`
  Examples and tests for `cfis_field_select.py`


