[Home](../shapepipe.md) | [Modules](../module_docs.md)

# File IO package documentation

This documentation is only for the class FITSCatalog.

## Open a file

The first step is to create a FITSCatalog object (call 'example' here):

```python
example = scatalog.FITSCatalog(fullpath, hdu_no=None,
          open_mode=BaseCatalog.OpenMode.ReadOnly,
          memmap=False, SEx_catalog=False)
```

**NOTE :** the file doesn't need to exist to do that. This is always the first step.

  - The path :
    - The path to an already existing fits file or where you want to create a new one.
  - HDU :
    - By default : `hdu_no = 1`
  - Possible open monde :
    - `BaseCatalog.OpenMode.ReadOnly` : allow you to read the contents of your fits file but none modification can be done.
    - `BaseCatalog.OpenMode.ReadWrite` : allow you to modify the fits file. For example, you can add HDUs or create a new fits file.
    - `BaseCatalog.OpenMode.Append`
  - SEx_catalog :
    - This option allow you to set some default parameters in the case you're working on a SExtractor output.
    - By default `SEx_catalog = False`
    - Set `hdu_no = 2` (where the data are in a SExtractor FITS_LDAC file).
    - Needed to create a new SEXtractor-like fits file.

The second step is to actually open the file (if the file already exist):

`example.open()`

**Examples :**

- Open a fits table with a PrimaryHDU :
```python
example = scatalog.FITSCatalog('path_to_the_table')
example.open()
```
- Open a fits image :
```python
example = scatalog.FITSCatalog('path_to_the_image', hdu_no=0)
example.open()
```
- Open a SExtractor LDAC catalog :
```python
example = scatalog.FITSCatalog('path_to_SExCat', SEx_catalog=True)
example.open()
```

## Getting the data and the header :

Functions to use :
```python
scatalog.FITSCatalog.get_data(hdu_no)
scatalog.FITSCatalog.get_header(hdu_no)
```

- hdu_no :
  - By default, use the value set when the file was opened.

If the HDU is not specified it's set to the value used when the object was created.
The data are returned as an astropy.io.fits.fitsrec.FITS_rec format, which is similar to a numpy.recarray (every column as a label).
The header is returned as a python dictionary.

**Example :**

```python
data = example.get_data()
data['X_POSITION']
out :
  array([ 1.23, 4.56, ..., 7.89], dtype=float32)

header = example.get_header()
header['BITPIX']
out :
  8
```

## Create a new file

**NOTE :** When you create a fits table file it always have a PrimaryHDU. To create an image fits file with the image in the first HDU see section 'save data to a fits file'.

The first step is to create a FITSCatalog object (see 'Open a file' section).

The second step is to create the file using :

```python
example.create(ext_name=None, overwrite=True, s_hdu=True, sex_cat_path=None)
```

- ext_name :
  - The name of the "main" extansion/HDU.
  - The option is necessary only if you create a secondary HDU.
- overwrite :
  - If `overwrite = True` and the file already exist, will overwrite it.
- s_hdu :
  - If `s_hdu = True` will create a secondary HDU.
- sex_cat_path :
  - If you create a file opened with `SEx_catalog = True` you need to provide a path to a SExtractor catalog to mimic. It will copy the secondary HDU from the path you provide.

**Example :**

- Create a fits table with a PrimaryHDU :
```python
example = scatalog.FITSCatalog('path_of_the_file')
example.create()
```

- Create a SExtractor LDAC like catalog :
```python
example = scatalog.FITSCatalog('path_of_the_file', SEx_catalog=True)
example.create(sex_cat_path='path_to_SExtractor_catalog_to_mimic')
```

## Save data to a fits file

The following function allow one to save data in an already existing file or in a new one. It will depend weither the FITSCatalog object point to an already existing file or not.

List of the compatible formats for the data :
- dict
- list*
- numpy.ndarray*
- numpy.recarray
- astropy.io.fits.fitsrec.FITS_rec (data format in an astropy fits file)

\* : a list with the names of the different columns is required.

The first step is to create a FITSCatalog object (see 'Open a file' section).

**WARNING !  :** The file need to be open with the option `open_mode = scatalog.BaseCatalog.OpenMode.ReadWrite`. Otherwise, no changement will be applied.

Function to use :
```python
example.save_as_fits(self, data=None, names=None, ext_name=None, sex_cat_path=None, image=False, overwrite=False)
```

- data :
  - The data to save in one of the format mentioned previously.
- names :
  - Name to give to the different columns.
- ext_name :
  - Name of the new HDU. Default : 'NEW'
- sex_cat_path :
  - If you create a file opened with `SEx_catalog = True` you need to provide a path to a SExtractor catalog to mimic. It will copy the secondary HDU from the path you provide.
- image :
  - if `image = True` it will create a new file with the image as PrimaryHDU.
- overwrite :
  - Only use when you create an image.
  - If `overwrite = True` and the file already exist, will overwrite it.

**Example :**

- Save a list into a new file :

 ```python
  example = scatalog.FITSCatalog('path_to_the_new_file', open_mode=scatalog.BaseCatalog.OpenMode.ReadWrite)
  data = [[1,2,3],[4,5,6],[7,8,9]]
  names = ['a','b','c']
  example.save_as_fits(data,names)

  example.get_info()
  out :
    Filename: example.fits
    No.    Name      Ver    Type      Cards   Dimensions   Format
      0  PRIMARY       1 PrimaryHDU       4   ()      
      1  NEW           1 BinTableHDU     15   3R x 3C   [K, K, K]
  ```

- Save a dictionary in an already existing file :

  ```python
  example = scatalog.FITSCatalog('path_of_the_file', open_mode=scatalog.BaseCatalog.OpenMode.ReadWrite)
  data = {'a': [1,2,3], 'b': [4,5,6], 'c': [7,8,9]}
  example.save_as_fits(data, ext_name='new_hdu')

  example.get_info()
  out :
    Filename: example.fits
    No.    Name      Ver    Type      Cards   Dimensions   Format
      0  PRIMARY       1 PrimaryHDU       6   ()      
      1  LDAC_IMHEAD    1 BinTableHDU     12   1R x 1C   [8640A]   
      2  LDAC_OBJECTS    1 BinTableHDU     66   1573R x 28C   [J, E, E, E, E, E, E, E, E, D, D, D, D, D, I, J, E, D, D, E, E, E, E, E, 2500E, E, E, E]   
      3  NEW_HDU       1 BinTableHDU     15   3R x 3C   [K, K, K]
  ```

- Save an array as image :

  ```python
  example = example = scatalog.FITSCatalog('path_to_the_new_file', open_mode=scatalog.BaseCatalog.OpenMode.ReadWrite)
  data = numpy.array([[1,2,3],[4,5,6],[7,8,9]])
  example.save_as_fits(data, image=True)

  example.get_info()
  out :
    Filename: test.fits
    No.    Name      Ver    Type      Cards   Dimensions   Format
      0  PRIMARY       1 PrimaryHDU       6   (3, 3)   int64
  ```
