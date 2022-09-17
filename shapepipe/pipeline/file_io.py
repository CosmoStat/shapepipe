"""FILE IO.

This file contains methods for file I/O handling.

:Author: Marc Gentile and Axel Guinot

"""

import os

import numpy as np
from astropy.io import fits
from astropy import units
from astropy.table import Table


class BaseCatalogue(object):
    """Base catalogue.

    Base catalogue management class.

    Parameters
    ----------
    fullpath : str
      Full path to catalogue

    Notes
    -----
    This class is not meant to be used directly

    """

    def __init__(self, fullpath):

        # catalogue file path
        self._directory, self._filename = os.path.split(fullpath)
        # input/output format
        self._format = BaseCatalogue.InputFormat.Undefined
        # catalogue internal data (e.g. AsciiData class)
        self._cat_data = None

    @property
    def fullpath(self):
        """Set Full Path.

        Get the full path of the catalogue file

        """
        return os.path.join(self._directory, self._filename)

    @property
    def directory(self):
        """Set Directory.

        Get the directory of the catalogue file

        """
        return self._directory

    @property
    def filename(self):
        """Set File Name.

        Get the name of the catalogue file

        """
        return self._filename

    @property
    def format(self):
        """Format.

        Get the default input/output format of the catalogue
        (e.g. Text, SExtractor, FITS)

        """
        return self._format

    def get_nb_rows(self):
        """Get Number of Rows.

        Get the number of rows in the catalogue

        Returns
        -------
        int
            Number of rows

        """
        raise BaseCatalogue.FeatureNotImplemented('get_nb_rows()')

    def get_nb_cols(self):
        """Get Number of Columns.

        Get the number of columns in the catalogue

        Returns
        -------
        int
            Number of columns
        """
        raise BaseCatalogue.FeatureNotImplemented("get_nb_cols()")

    def get_col_names(self):
        """Get Column Names.

        Get the list of column names in the catalogue

        Returns
        -------
        list
            list of column names
        """
        raise BaseCatalogue.FeatureNotImplemented('get_col_names()')

    def get_col_formats(self):
        """Get Column Formats.

        Get the list of column formats in the order of columns

        """
        raise BaseCatalogue.FeatureNotImplemented('get_col_names()')

    def add_col(
        self,
        col_name,
        col_format=None,
        col_comment=None,
        col_data=None,
    ):
        """Add Column.

        Add a Column to the catalogue

        Parameters
        ----------
        col_name : str
            Column name
        col_format : str
            Column python format
        col_comment : str
            Column comment (ignored in for class)
        col_data : numpy.ndarray
            Column data as a numpy array

        """
        raise BaseCatalogue.FeatureNotImplemented('add_col()')

    def _file_exists(self, filepath):
        """File Exists.

        Check if input file path is a valid file.

        Parameters
        ----------
        filepath : str
            Path to file

        """
        if not os.path.isfile(filepath):
            return False
        else:
            return True

    class InputFormat:
        """Supported input catalogue formats."""

        Undefined = 0
        TabulatedText = 1
        SExtractor = 2
        FITS = 4
        FITS_LDAC = 5

    class OpenMode:
        """Supported input catalogue open modes."""

        ReadOnly = 'readonly'
        ReadWrite = 'update'
        Append = 'append'

    class Column(object):
        """Column.

        Represents a column in the catalogue.

        """

        def __init__(self):
            self._cat_col = None

        @property
        def name(self):
            """Set Name.

            Get the name of the column
            """
            raise BaseCatalogue.FeatureNotImplemented('column.name')

        @property
        def format(self):
            """Format.

            Get the format of the column

            """
            raise BaseCatalogue.FeatureNotImplemented('column.format')

        @property
        def data(self):
            """Set Data.

            Get the data associated with the column

            """
            raise BaseCatalogue.FeatureNotImplemented('column.data')

        def get_nb_rows(self):
            """Get Number of Rows.

            Retrieve the number of rows of the column.

            """
            raise BaseCatalogue.FeatureNotImplemented('get_nb_rows()')

        def get_info(self):
            """Get Information.

            Retrieve information about the column.

            """
            raise BaseCatalogue.FeatureNotImplemented('get_info()')

        def get_type(self):
            """Get Type.

            Get the data type of the column

            """
            raise BaseCatalogue.FeatureNotImplemented('get_type()')

    class FeatureNotImplemented(NotImplementedError):
        """Feature Not Implemented.

        Parameters
        ----------
        msg : str
            Message

        """

        def __init__(self, msg):
            self._msg = msg

        def __str__(self):
            return (
                f'File IO *** ERROR ***: Feature: {self._msg} is not '
                + 'implemented in this class'
            )

    class catalogueNotOpen(Exception):
        """Catalogue Not Open.

        Catalogue has not been open yet.

        Parameters
        ----------
        filepath : str
            Path to file

        """

        def __init__(self, filepath):
            self._filepath = filepath

        def __str__(self):
            return (
                f'File IO *** ERROR ***: catalogue: {self._filepath} '
                + 'is not open'
            )

    class DataNotFound(Exception):
        """Data Not Found.

        No data found (at given hdu).

        Parameters
        ----------
        filepath : str
            Path to file
        hdu : int
            HDU number

        """

        def __init__(self, filepath, hdu):
            self._filepath = filepath
            self._hdu = hdu

        def __str__(self):
            return (
                f'File IO *** ERROR ***: File \'{self._filepath}\', '
                + f'hdu={self._hdu}: data not found'
            )

    class catalogueFileNotFound(Exception):
        """Catalogue File Not Found.

        Exception thrown when a catalogue file is not found on disk.

        Parameters
        ----------
        filepath : str
            Path to file

        """

        def __init__(self, filepath):
            self._filepath = filepath

        def __str__(self):
            """Set string representation of the exception object."""
            return f'File IO *** ERROR ***: file {self._filepath} no found'

    class ColumnNotFound(Exception):
        """Column Not Found.

        Exception thrown when a named catalogue column could not be found

        Parameters
        ----------
        col_name : str
            Column name

        """

        def __init__(self, col_name):
            self._col_name = col_name

        def __str__(self):
            """Set string representation of the exception object."""
            return f'File IO *** ERROR ***: column {self._col_name} no found'

    class catalogueNotCreated(Exception):
        """Catalogue Not Created.

        Catalogue could not be created.

        Parameters
        ----------
        filepath : str
            Path to file

        """

        def __init__(self, filepath):
            self._filepath = filepath

        def __str__(self):
            return (
                f'File IO *** ERROR ***: catalogue: {self._filepath} could '
                + 'not be created'
            )

    class OpenModeNotSupported(Exception):
        """Open Mode Not Supported.

        Opening mode not supported by this version.

        Parameters
        ----------
        filepath : str
            Path to file
        open_mode : OpenMode
            File opening mode

        """

        def __init__(self, filepath, open_mode):
            self._filepath = filepath
            self._open_mode = open_mode

        def __str__(self):
            return (
                f'File IO *** ERROR ***: catalogue: {self._filepath} '
                + 'Open Mode {self._open_mode} not supported'
            )

    class OpenModeConflict(Exception):
        """Open Mode Conflict.

        Opening mode is in conflict with an action.

        Parameters
        ----------
        open_mode : OpenMode
            File opening mode
        open_mode_needed : OpenMode
            File opening mode needed

        """

        def __init__(self, open_mode, open_mode_needed):
            self._open_mode = open_mode
            self._open_mode_needed = open_mode_needed

        def __str__(self):
            return (
                'File IO *** ERROR ***: catalogue has to be open as : '
                + f'{self._open_mode_needed} , Mode used : {self._open_mode}'
            )


class FITSCatalogue(BaseCatalogue):
    """FITS Catalogue.

    Catalogues management in .FITS format.

    Parameters
    ----------
    fullpath : str
        Full path to file
    hdu_no : int
        HDU number
    open_mode : OpenMode
        File opening mode
    memmap : Bool
        Option to use memory mapping
    SEx_catalogue : bool
        Option to specify if the input is a SExtractor catalogue

    """

    def __init__(
        self,
        fullpath,
        hdu_no=None,
        open_mode=BaseCatalogue.OpenMode.ReadOnly,
        memmap=False,
        SEx_catalogue=False,
    ):
        BaseCatalogue.__init__(self, fullpath)

        # default input/output format
        self._format = BaseCatalogue.InputFormat.FITS

        # opening mode (see FITSCatalogue.OpenMode)
        self._open_mode = open_mode
        # Work with SExtractor fits format or not
        self._SEx_catalogue = SEx_catalogue
        # HDU number of the underlying .FITS table
        if hdu_no is None:
            # Default is 1 (or 2 if you are using )
            if SEx_catalogue:
                self._hdu_no = 2
            else:
                self._hdu_no = 1
        else:
            self._hdu_no = hdu_no
        # use memory mapping or not
        self._use_memmap = memmap

    def __str__(self):
        if self._cat_data is not None:
            info = f'{self.get_info()}'
        else:
            info = 'No information'
        return info

    @property
    def hdu_no(self):
        """Set HDU Number.

        Get the HDU index of the table.

        """
        return self._hdu_no

    @property
    def open_mode(self):
        """Set Open Mode.

        Return the catalogue opening mode.

        See Also
        --------
        FITSCatalogue.OpenMode

        """
        return self._open_mode

    @property
    def use_memmap(self):
        """Use Memory Map.

        If True, use memory mapping.

        """
        return self._use_memmap

    @hdu_no.setter
    def hdu_no(self, hdu_no):
        self._hdu_no = hdu_no

    @open_mode.setter
    def open_mode(self, open_mode):
        self._open_mode = open_mode

    @use_memmap.setter
    def use_memmap(self, use_memmap):
        self._use_memmap = use_memmap

    def open(self):
        """Open.

        Open an existing catalogue.

        """
        if self._file_exists(self.fullpath):

            # Open catalogue file
            self._cat_data = fits.open(
                self.fullpath,
                mode=self.open_mode,
                memmap=self.use_memmap,
                ignore_missing_end=True,
            )
        else:
            raise BaseCatalogue.catalogueFileNotFound(self.fullpath)

    def create(self, ext_name=None, s_hdu=True, sex_cat_path=None):
        """Create.

        Create an empty catalogue in FITS format.

        Parameters
        ----------
        ext_name : str
            Extension name or number
        s_hdu : bool
            If true add a secondary HDU
        sex_cat_path : str
            Path to SEXtractor catalogue

        """
        primary_hdu = fits.PrimaryHDU()
        if self._SEx_catalogue:
            if sex_cat_path is not None:
                if self._file_exists(sex_cat_path):
                    sex_cat = FITSCatalogue(sex_cat_path, hdu_no=1)
                    sex_cat.open()
                    secondary_hdu = sex_cat._cat_data[1]
                    self._cat_data = fits.HDUList([primary_hdu, secondary_hdu])
                    self._cat_data.writeto(self.fullpath, overwrite=True)
                    sex_cat.close()
                    del sex_cat
                else:
                    raise BaseCatalogue.catalogueFileNotFound(sex_cat_path)
            else:
                raise ValueError(
                    'sex_cat_path needs to be provided to create a '
                    + 'SEXtractor catalogue'
                )
        elif s_hdu:
            secondary_hdu = fits.BinTableHDU(
                data=None,
                header=None,
                name=ext_name,
            )
            self._cat_data = fits.HDUList([primary_hdu, secondary_hdu])
            self._cat_data.writeto(self.fullpath, overwrite=True)
        else:
            self._cat_data = fits.HDUList([primary_hdu])
            self._cat_data.writeto(self.fullpath, overwrite=True)

    def copy_hdu(self, fits_file=None, hdu_no=None, hdu_name=None):
        """Copy HDU.

        Copy an HDU from a FITS file.

        Parameters
        ----------
        fits_file : str
            Catalogue containing the HDU
        hdu_no : int
            HDU number
        hdu_name : str
            HDU name

        """
        if fits_file is None:
            fits_file = self

        if self._cat_data is None:
            raise BaseCatalogue.catalogueNotOpen(self.fullpath)
        if fits_file._cat_data is None:
            raise BaseCatalogue.catalogueNotOpen(fits_file.fullpath)

        if hdu_no is None:
            hdu_no = fits_file.hdu_no

        if hdu_name is None:
            hdu_name = fits_file._cat_data[hdu_no].name

        self._cat_data.append(
            fits.BinTableHDU(fits_file.get_data(hdu_no=hdu_no), name=hdu_name)
        )

    def apply_mask(
        self,
        fits_file=None,
        hdu_no=None,
        mask=None,
        hdu_name=None,
    ):
        """Apply Mask.

        Apply a mask to data for a specified HDU.

        Parameters
        ----------
        fits_file : str
            FITS file containing the data to be masked
        hdu_no : int
            HDU number
        mask : numpy.ndarray
            Array of booleans or an array of indices
        hdu_name : str
            HDU name

        Notes
        -----
        This will create a new HDU containing data[mask]

        """
        if fits_file is None:
            fits_file = self
        if self._cat_data is None:
            raise BaseCatalogue.catalogueNotOpen(self.fullpath)
        if fits_file._cat_data is None:
            raise BaseCatalogue.catalogueNotOpen(fits_file.fullpath)
        if mask is None:
            raise ValueError('Mask not provided')
        if type(mask) is not np.ndarray:
            raise TypeError('Mask need to be a numpy.ndarray')
        if hdu_no is None:
            hdu_no = fits_file.hdu_no

        if hdu_name is None:
            hdu_name = fits_file._cat_data[hdu_no].name

        if mask.dtype == bool:
            mask = np.where(mask is True)
            self._cat_data.append(
                fits.BinTableHDU(
                    fits_file.get_data(hdu_no)[:][mask],
                    name=hdu_name,
                )
            )
        elif mask.dtype == int:
            self._cat_data.append(
                fits.BinTableHDU(
                    fits_file.get_data(hdu_no)[:][mask],
                    name=hdu_name,
                )
            )
        else:
            raise TypeError('Mask type must be of type int or bool')

    def save_as_fits(
        self,
        data=None,
        names=None,
        ext_name=None,
        sex_cat_path=None,
        image=False,
        image_header=None,
        overwrite=False,
    ):
        """Save as FITS.

        Save data into an existing FITS file or into a new one.

        Save data from dict, list, numpy.ndarray, numpy.recarray or
        astropy.io.fits.fitsrec.FITS_rec (data format in an astropy fits file)
        When creating a new FITS to store BinTable data it create a PrimaryHDU.
        When creating a new FITS to store Image there is no PrimaryHDU.
        You can create a SExtractor format FITS by specifying a SExtractor
        catalogue from where data come from.

        Parameters
        ----------
        data : numpy.ndarray
            Data to be stored
        names : list
            List of column names
        ext_name : str
            Name of the HDU where data are stored
        sex_cat_path : str
            Path of the existing SExtractor catalogue to mimic
        image : bool
            If true create a fits image
        image_header : astropy.io.fits.header
            Header to use when saving an image
        overwrite : bool
            Option to overwrite an existing image, only used when creating a
            FITS image

        Notes
        -----
        To create a SExtractor-like FITS file you need to specify
        ``SEx_catalogue=True`` when declaring the FITSCatalogue object.

        """
        if self.open_mode != FITSCatalogue.OpenMode.ReadWrite:
            raise BaseCatalogue.OpenModeConflict(
                open_mode=self.open_mode,
                open_mode_needed=FITSCatalogue.OpenMode.ReadWrite,
            )

        if data is None:
            raise ValueError('Data not provided')

        if not image:
            if type(data) is dict:
                names = list(data.keys())
                it = list(range(len(names)))
                if len(names) == 1:
                    data = np.array(data[names[0]])
                else:
                    data = [np.array(data[i]) for i in names]
                self._save_to_fits(
                    data,
                    names,
                    it,
                    ext_name,
                    sex_cat_path,
                    overwrite=overwrite,
                )

            elif type(data) is np.recarray:
                names = list(data.dtype.names)
                it = names
                self._save_to_fits(
                    data,
                    names,
                    it,
                    ext_name,
                    sex_cat_path,
                    overwrite=overwrite,
                )

            elif type(data) is fits.fitsrec.FITS_rec:
                self._save_from_recarray(data, ext_name, sex_cat_path)

            elif type(data) is np.ndarray:
                if names is None:
                    if data.dtype.names is not None:
                        names = data.dtype.names
                        it = names
                    else:
                        raise ValueError('Names not provided')
                else:
                    it = range(len(names))
                self._save_to_fits(
                    data,
                    names,
                    it,
                    ext_name,
                    sex_cat_path,
                    overwrite=overwrite,
                )

            elif type(data) is list:
                if names is None:
                    raise ValueError('Names not provided')
                it = range(len(names))
                data = np.asarray(data)
                self._save_to_fits(
                    data,
                    names,
                    it,
                    ext_name,
                    sex_cat_path,
                    overwrite=overwrite,
                )

            elif type(data) is Table:
                if names is None:
                    raise ValueError('Names not provided')
                it = names
                self._save_to_fits(
                    data,
                    names,
                    it,
                    ext_name,
                    sex_cat_path,
                    overwrite=overwrite,
                )

        else:
            if type(data) is np.ndarray:
                self._save_image(
                    data=data,
                    header=image_header,
                    overwrite=overwrite,
                )
            else:
                raise TypeError('Data need to be a numpy.ndarray')

    def create_from_numpy(
        self,
        matrix,
        col_names,
        ext_name=None,
        ext_ver=None,
        header=None,
    ):
        """Create from Numpy.

        Create a new catalogue from a two-dimensional numpy array.

        Parameters
        ----------
        matrix : numpy.ndarray
            Two-dimensional numpy array
        col_names : list
            List of column names to use as the header
        ext_name : str
            Extension name or number
        ext_ver : str
            Extension version
        header : list
            List of dictionaries with keys: 'card', name', 'value',
            'value_orig', 'comment'

        """
        col_list = []
        for col_name in col_names:
            icol = col_names.index(col_name)
            col_type = self._get_fits_col_type(matrix[:, icol])
            col_data = fits.Column(
                name=col_name,
                format=col_type,
                array=np.ravel(matrix[:, icol]),
            )
            col_list.append(col_data)

        fits_header = None
        if header is not None:
            fits_header = fits.Header()
            for (k, v) in header.items():
                fits_header[k] = v

        primary_hdu = fits.PrimaryHDU()
        secondary_hdu = fits.BinTableHDU.from_columns(
            col_list,
            header=fits_header,
        )
        if ext_name is not None:
            secondary_hdu.name = ext_name

        self._cat_data = fits.HDUList(hdus=[primary_hdu, secondary_hdu])
        self._cat_data.writeto(self.fullpath, overwrite=True)

    def close(self):
        """Close."""
        if self._cat_data is not None:
            if self.open_mode == FITSCatalogue.OpenMode.ReadWrite:
                self.save()
            self._cat_data.close()
            del self._cat_data
            self._cat_data = None
        else:
            raise BaseCatalogue.catalogueNotOpen(self.fullpath)

    def save(self):
        """Save."""
        if self.open_mode == FITSCatalogue.OpenMode.ReadWrite:
            self._cat_data.flush()
        else:
            raise BaseCatalogue.OpenModeConflict(
                open_mode=self.open_mode,
                open_mode_needed=FITSCatalogue.OpenMode.ReadWrite,
            )

    def get_nb_rows(self, hdu_no=None):
        """Get Number of Rows.

        Get the number of rows in the catalogue.

        Parameters
        ----------
        hdu_no : int
            HDU index

        Returns
        -------
        int
            Number of rows or 0 if None

        """
        nb_rows = 0
        if self._cat_data is None:
            raise BaseCatalogue.catalogueNotOpen(self.fullpath)

        if len(self._cat_data) > 0:
            if hdu_no is None:
                hdu_no = self.hdu_no

            if self._cat_data[hdu_no].size > 0:
                nb_rows = len(self._cat_data[hdu_no].data)

        return nb_rows

    def get_nb_cols(self, hdu_no=None):
        """Get Number of Columns.

        Get the number of columns in the catalogue.

        Parameters
        ----------
        hdu_no : int
            HDU index

        Returns
        -------
        int
            Number of columns or 0 if None

        """
        if self._cat_data is None:
            raise BaseCatalogue.catalogueNotOpen(self.fullpath)
        if hdu_no is None:
            hdu_no = self.hdu_no
        return len(self.get_col_names(hdu_no))

    def get_col_names(self, hdu_no=None):
        """Get Column Names.

        Get the list of column names in the catalogue.

        Parameters
        ----------
        hdu_no : int
            HDU index

        Returns
        -------
        list
            list of column names

        """
        if self._cat_data is None:
            raise BaseCatalogue.catalogueNotOpen(self.fullpath)
        else:
            if hdu_no is None:
                hdu_no = self.hdu_no
            return self._cat_data[hdu_no].columns.names

    def get_info(self):
        """Get Information.

        Retrieve some information about the catalogue.

        Returns
        -------
        dict
            A dictionary with detailed information

        Notes
        -----
        See the fitsio documentation of the info() function for the details

        """
        if self._cat_data is not None:
            return self._cat_data.info()
        else:
            return fits.info(self.fullpath)

    def get_ext_name(self, hdu_no=None):
        """Get Extension Name.

        Return the name of an extansion or all of them.

        Parameters
        ----------
        hdu_no : int
            index of the hdu to return, if None return all of them

        Returns
        -------
        str or list
            The hdu name or a list of strings

        """
        if hdu_no is None:
            n = [self._cat_data[i].name for i in range(len(self._cat_data))]
        else:
            n = self._cat_data[int(hdu_no)].name

        return n

    def col_exists(self, col_name, hdu_no=None):
        """Column Exists.

        Check whether a named column exists.

        Parameters
        ----------
        col_name : str
            Column name
        hdu_no : int
            HDU index

        Returns
        -------
        bool
            True if column exists, False otherwise
        """
        return col_name in self.get_col_names()

    def get_col_index(self, col_name, hdu_no=None):
        """Get Column Index.

        Get the column index from a column name.

        Parameters
        ----------
        col_name : str
            Column name
        hdu_no : int
            HDU index

        Returns
        -------
        int
            Column index (zero-based)

        """
        col_names = self.get_col_names()
        if col_name not in col_names:
            raise BaseCatalogue.ColumnNotFound(col_name)
        return col_names.index(col_name)

    def get_col_data(self, col_index, hdu_no=None):
        """Get Column Data.

        Return the data of a column from its index.

        Parameters
        ----------
        col_name : str
            Column name
        hdu_no : int
            HDU index

        Returns
        -------
        numpy.ndarray
            Data associated with the column

        """
        if self._cat_data is not None:
            if hdu_no is None:
                hdu_no = self.hdu_no

            col_data = self._cat_data[hdu_no].columns[col_index].array
            if isinstance(col_data, fits.column.Delayed):
                # trick to force data to be loaded into memory
                _ = self._cat_data[hdu_no].data
            return self._cat_data[hdu_no].columns[col_index].array
        else:
            raise BaseCatalogue.catalogueNotOpen(self.fullpath)

    def get_named_col_data(self, col_name, hdu_no=None):
        """Get Named Column Data.

        Return the data of a column from its index (zero-based).

        Parameters
        ----------
        col_name : str
            Column name
        hdu_no : int
            HDU index

        Returns
        -------
        numpy.ndarray
            Data associated with the column

        """
        if self._cat_data is not None:
            if hdu_no is None:
                hdu_no = self.hdu_no

            col_data = self._cat_data[hdu_no].columns[col_name].array
            if isinstance(col_data, fits.column.Delayed):
                # trick to force data to be loaded into memory
                _ = self._cat_data[hdu_no].data
            return self._cat_data[hdu_no].columns[col_name].array
        else:
            raise BaseCatalogue.catalogueNotOpen(self.fullpath)

    def get_data(self, hdu_no=None):
        """Get Data.

        Return data of the specified hdu.

        Parameters
        ----------
        hdu_no : int
            HDU index

        Returns
        -------
        numpy.ndarray
            Data associated with the HDU

        """
        if self._cat_data is not None:
            if hdu_no is None:
                hdu_no = self.hdu_no

            dat = self._cat_data[hdu_no].data
            if dat is None:
                raise BaseCatalogue.DataNotFound(self.fullpath, self.hdu_no)

            return dat

        else:
            raise BaseCatalogue.catalogueNotOpen(self.fullpath)

    def get_header(self, hdu_no=None):
        """Get Header.

        Return the catalogue header as a list of dictionaries.

        Parameters
        ----------
        hdu_no : int
            HDU index

        Returns
        -------
        astropy.io.fits.header
            FITS header

        Notes
        -----
        See astropy documentation

        """
        if self._cat_data is not None:
            if hdu_no is None:
                hdu_no = self.hdu_no
            return dict(self._cat_data[hdu_no].header.items())
        else:
            raise BaseCatalogue.catalogueNotOpen(self.fullpath)

    def get_header_value(self, request, hdu_no=None):
        """Get Header Value.

        Return the value of a parameters or a linear combination of parameters
        and/or numbers.

        Parameters
        ----------
        request : str
            Request parameter or a linear combination of parameters
        hdu_no : int
            HDU index

        Returns
        -------
        float
            Result of the request

        """
        if request is None:
            raise ValueError('request not provided')
        if type(request) is not str:
            raise TypeError('request has to be a string')

        if hdu_no is None:
            hdu_no = self._hdu_no

        header = self.get_header(hdu_no=hdu_no)
        if header is None:
            raise ValueError(f'Empty header in the hdu : {hdu_no}')

        return interpreter(string=request, catalogue=header).result

    def add_header_card(self, key, value=None, comment=None, hdu_no=None):
        """Add Header Card.

        Add a card in the header of the specified HDU.

        Parameters
        ----------
        key : str
            The key to add
        value : any
            The value of the key
        comment : str
            Comment for the key
        hdu_no : int
            HDU index

        """
        if self.open_mode != FITSCatalogue.OpenMode.ReadWrite:
            raise BaseCatalogue.OpenModeConflict(
                open_mode=self.open_mode,
                open_mode_needed=FITSCatalogue.OpenMode.ReadWrite,
            )

        if self._cat_data is None:
            raise BaseCatalogue.catalogueNotOpen(self.fullpath)

        if hdu_no is None:
            hdu_no = self._hdu_no

        card = []
        if key is None:
            raise ValueError('key not provided')
        else:
            card.append(key)

        if value is not None:
            card.append(value)
        else:
            if comment is not None:
                card.append('')

        if comment is not None:
            card.append(comment)

        card = tuple(card)

        self._cat_data[hdu_no].header.append(card, end=True)

    def get_headers(self):
        """Get Headers.

        Return the catalogue header as a list of dictionaries.

        Returns
        -------
        list
            list of headers

        """
        headers = []
        try:
            for hdu in self._cat_data:
                headers.append(dict(hdu.header.items()))
        except Exception:
            pass

        return headers

    def get_comments(self, hdu_no=None):
        """Get Comments.

        Return the list catalogue comments.

        Parameters
        ----------
        hdu_no : int
            HDU index

        Returns
        -------
        list
            List of catalogue comments

        """
        if self._cat_data is not None:
            if hdu_no is None:
                hdu_no = self.hdu_no
            return [
                self._cat_data[hdu_no].header.comments[c]
                for c in self._cat_data[hdu_no].header.keys()
            ]
        else:
            raise BaseCatalogue.catalogueNotOpen(self.fullpath)

    def get_col_comments(self, hdu_no=None):
        """Get Column Comments.

        Return the list of column comments.

        Parameters
        ----------
        hdu_no : int
            HDU index

        Returns
        -------
        list
            List of column comments

        """
        if self._cat_data is not None:
            if hdu_no is None:
                hdu_no = self.hdu_no
            hdr_col_types = [
                tt for tt in self._cat_data[hdu_no].header.keys()
                if 'TTYPE' in tt
            ]
            return [
                self._cat_data[hdu_no].header.comments[c]
                for c in hdr_col_types
            ]
        else:
            raise BaseCatalogue.catalogueNotOpen(self.fullpath)

    def get_col_formats(self, hdu_no=None):
        """Get Column Formats.

        Get the list of python column formats in the order of columns.

        Parameters
        ----------
        hdu_no : int
            HDU index

        Returns
        -------
        list
            List of column formats

        """
        if self._cat_data is not None:
            if hdu_no is None:
                hdu_no = self.hdu_no
            return self._cat_data[hdu_no].columns.formats
        else:
            raise BaseCatalogue.catalogueNotOpen(self.fullpath)

    def add_col(
        self,
        col_name,
        col_data,
        hdu_no=None,
        ext_name=None,
        new_cat=False,
        new_cat_inst=None,
    ):
        """Add Column.

        Add a Column to the catalogue.

        Parameters
        ----------
        col_name : str
            Column name
        col_data : numpy.ndarray
            Column data
        hdu_no : int
            HDU index
        ext_name : str, optional
            Change the name of the extansion
        new_cat : bool, optional
            If true will save the changes into a new catalogue
        new_cat_inst : io.FITSCatalogue
            New catalogue object

        """
        if new_cat:
            open_mode = new_cat_inst.open_mode
            output_path = new_cat_inst.fullpath
        else:
            open_mode = self.open_mode
            output_path = self.fullpath

        if self._cat_data is None:
            raise BaseCatalogue.catalogueNotOpen(self.fullpath)

        if open_mode != FITSCatalogue.OpenMode.ReadWrite:
            raise BaseCatalogue.OpenModeConflict(
                open_mode=open_mode,
                open_mode_needed=FITSCatalogue.OpenMode.ReadWrite,
            )

        if type(col_data) != np.ndarray:
            TypeError('col_data must be a numpy.ndarray')

        if hdu_no is None:
            hdu_no = self.hdu_no
        if ext_name is None:
            ext_name = self._cat_data[hdu_no].name

        n_of_hdu = len(self._cat_data)
        old_hdu_prev = []
        for i in range(0, hdu_no):
            old_hdu_prev.append(self._cat_data[i])
        old_hdu_next = []
        for i in range(hdu_no + 1, n_of_hdu):
            old_hdu_next.append(self._cat_data[i])

        new_fits = fits.HDUList(old_hdu_prev)

        col_list = self._cat_data[hdu_no].data.columns

        data_type = self._get_fits_col_type(col_data)
        data_shape = col_data.shape[1:]
        dim = str(tuple(data_shape))
        mem_size = 1
        if len(data_shape) != 0:
            for k in data_shape:
                mem_size *= k
            data_format = f'{mem_size}{data_type}'
            new_col = fits.ColDefs([fits.Column(
                name=col_name,
                format=data_format,
                array=col_data,
                dim=dim,
            )])
            col_list += new_col
        elif data_type == 'A':
            mem_size *= len(max(col_data, key=len))
            data_format = f'{mem_size}{data_type}'
            new_col = fits.ColDefs([fits.Column(
                name=col_name,
                format=data_format,
                array=col_data,
                dim=str((mem_size,)),
            )])
            col_list += new_col
        else:
            data_format = f'{mem_size}{data_type}'
            new_col = fits.ColDefs([fits.Column(
                name=col_name,
                format=data_format,
                array=col_data,
            )])
            col_list += new_col

        new_fits.append(fits.BinTableHDU.from_columns(col_list, name=ext_name))

        new_fits += fits.HDUList(old_hdu_next)

        new_fits.writeto(output_path, overwrite=True)

        if not new_cat:
            self._cat_data.close()
            del self._cat_data
            self._cat_data = fits.open(
                self.fullpath,
                mode=self.open_mode,
                memmap=self.use_memmap,
            )

    def remove_col(self, col_index):
        """Remove Column.

        Delete a column from its index.

        Parameters
        ----------
        col_index : int
            Index of the column to delete

        """
        raise BaseCatalogue.FeatureNotImplemented('remove_col()')

    def remove_named_col(self, col_name):
        """Remove Named Column.

        Delete a column from its index.

        Parameters
        ----------
        col_name : str
            Name of the column to delete

        """
        if self._cat_data is not None:
            col_index = self.get_col_index(col_name)
            self.remove_col(col_index)
        else:
            raise BaseCatalogue.catalogueNotOpen(self.fullpath)

    def _append_col(self, column, hdu_no=None):
        """Append Column.

        Append a Column object.

        Parameters
        ----------
        column : ?
            An object derived from BaseCatalogue.Column
        hdu_no : int
            HDU index

        """
        if self._cat_data is not None:
            new_col = fits.Column(
                name=column.name,
                format=column.format,
                array=column.data,
            )

            if hdu_no is None:
                hdu_no = self.hdu_no

            orig_table = fits.open(self.fullpath)[hdu_no].data
            orig_cols = orig_table.columns

            new_col = fits.ColDefs([fits.Column(
                name=column.name,
                format=column.format,
                array=np.zeros(len(orig_table)),
            )])
            col_list = orig_cols + new_col
            hdu = fits.BinTableHDU.from_columns(col_list)
            hdu.data[column.name] = column.data
            hdu.writeto(self.fullpath, overwrite=True)

        else:
            raise BaseCatalogue.catalogueNotOpen(self.fullpath)

    def _get_fits_col_type(self, col_data):
        """Get FITS Column Type.

        Get the FITS data type of a given column.

        Parameters
        ----------
        col_data : any
            Column data

        Returns
        -------
        str
            Column FITS data type

        """
        if col_data is None or len(col_data) == 0:
            col_type = 'D'
        elif type(col_data[0]) in [np.int16]:
            col_type = 'I'
        elif type(col_data[0]) in [np.int32]:
            col_type = 'J'
        elif type(col_data[0]) in [int, np.int64]:
            col_type = 'K'
        elif type(col_data[0]) in [float, np.float16, np.float32, np.float64]:
            col_type = 'D'
        elif type(col_data[0]) is bool:
            col_type = 'L'
        elif type(col_data[0]) in [str, np.str, np.str_, np.str0]:
            col_type = 'A'
        else:
            col_type = 'D'

        return col_type

    def _get_python_col_type(self, col_type):
        """Get Python Column Type.

        Get the Python data type of a given column.

        Parameters
        ----------
        col_data : any
            Column data

        Returns
        -------
        str
            Column Python data type

        """
        if col_type in ['B', 'I', 'J', 'K']:
            pcol_type = '%d'
        elif col_type in ['D', 'E']:
            pcol_type = '%f'
        elif col_type in ['A', 'C', 'M']:
            pcol_type = '%s'
        elif col_type == 'L':
            pcol_type = '%s'
        else:
            pcol_type = '%f'

        return pcol_type

    def _save_to_fits(
            self,
            data,
            names,
            it,
            ext_name=None,
            sex_cat_path=None,
            overwrite=False,
    ):
        """Save to FITS.

        Save array of data as fits with their associated column names.

        Parameters
        ----------
        data : numpy.ndarray
            Array with the data
        names : list
            List of the column names
        it : iterator
            ?
        ext_name : str
            Name of the HDU where data are stored
        sex_cat_path : str
            Path of the existing SExtractor catalogue to mimic
        overwrite : bool
            Option to overwrite an existing catalogue

        """
        if data is None:
            raise ValueError('Data not provided')

        if self._file_exists(self.fullpath) and not overwrite:
            if self._cat_data is None:
                self.open()
            if ext_name is None:
                ext_name = 'new'
        else:
            if self._SEx_catalogue:
                self.create(s_hdu=False, sex_cat_path=sex_cat_path)
                self.open()
                if ext_name is None:
                    ext_name = 'LDAC_OBJECTS'
            else:
                self.create(s_hdu=False)
                self.open()
                if ext_name is None:
                    ext_name = 'new'

        if len(names) == 1:
            data = np.array([data])
        col_list = []
        for idx in it:
            data_shape = data[idx].shape[1:]
            dim = str(tuple(data_shape))
            name = names[it.index(idx)]
            data_type = self._get_fits_col_type(data[idx])
            mem_size = 1
            if len(data_shape) != 0:
                for shape in data_shape:
                    mem_size *= shape
                data_format = f'{mem_size}{data_type}'
                col_list.append(fits.Column(
                    name=name,
                    format=data_format,
                    array=data[idx],
                    dim=dim,
                ))
            elif data_type == 'A':
                mem_size *= len(max(data[idx], key=len))
                data_format = f'{mem_size}{data_type}'
                col_list.append(fits.Column(
                    name=name,
                    format=data_format,
                    array=data[idx],
                    dim=str((mem_size,)),
                ))
            else:
                data_format = f'{mem_size}{data_type}'
                col_list.append(fits.Column(
                    name=name,
                    format=data_format,
                    array=data[idx],
                ))

        self._cat_data.append(
            fits.BinTableHDU.from_columns(col_list, name=ext_name)
        )
        self.close()

    def _save_from_recarray(
        self,
        data=None,
        ext_name=None,
        sex_cat_path=None,
        overwrite=False,
    ):
        """Save From Record Array.

        Save a numpy.recarray or astropy.io.fits.fitsrec.FITS_rec into a FITS
        file.

        Parameters
        ----------
        data : numpy.ndarray
            Array with the data
        ext_name : str
            Name of the HDU where data are stored
        sex_cat_path : str
            Path of the existing SExtractor catalogue to mimic
        overwrite : bool
            Option to overwrite an existing catalogue

        """
        if data is None:
            raise ValueError('Data not provided')

        if self._file_exists(self.fullpath) and not overwrite:
            if self._cat_data is None:
                self.open()
            if ext_name is None:
                ext_name = 'new'
            self._cat_data.append(fits.BinTableHDU(data, name=ext_name))
            self.close()
        else:
            if self._SEx_catalogue:
                self.create(s_hdu=False, sex_cat_path=sex_cat_path)
                self.open()
                if ext_name is None:
                    ext_name = 'LDAC_OBJECTS'
                self._cat_data.append(fits.BinTableHDU(data, name=ext_name))
                self.close()
            else:
                self.create(s_hdu=False)
                self.open()
                if ext_name is None:
                    ext_name = 'new'
                self._cat_data.append(fits.BinTableHDU(data, name=ext_name))
                self.close()

    def _save_image(self, data=None, header=None, overwrite=False):
        """Save Image.

        Save an image into a FITS file.

        Parameters
        ----------
        data : numpy.ndarray
            Array with the data
        header : astropy.io.fits.header
            FITS header
        overwrite : bool
            Option to overwrite an existing catalogue

        """
        if (data is not None):
            fits.PrimaryHDU(data, header).writeto(
                self.fullpath,
                overwrite=overwrite,
            )
        else:
            raise ValueError('Data or names not provided')

    class Column(BaseCatalogue.Column):
        """Column.

        Represents a column in the catalogue.

        Parameters
        ----------
        name : str
            Name name of the column
        format : str
            Python format of the column
        comment : str
            Comment strin
        data : numpy.ndarray
            Associated column data

        """

        def __init__(self, name, format=None, comment=None, data=None):
            BaseCatalogue.Column.__init__(self)

            self._name = name

            if format is None:
                format = 'D'
            self._format = format

            if comment is None:
                comment = name
            self._comment = comment

            if data is None:
                data = np.asarray([])
            else:
                if isinstance(data, list):
                    self._data = np.asarray(data)
                else:
                    self._data = data

        def __str__(self):
            info = f'{self._cat_col}'
            return info

        @property
        def name(self):
            """Set Name.

            Get the name of the column.

            """
            return self._name

        @property
        def format(self):
            """Format.

            Get the format of the column.

            """
            return self._format

        @property
        def comment(self):
            """Comment.

            Get the comment string associated with the column.

            """
            return self._comment

        @property
        def data(self):
            """Set Data.

            Get the data associated with the column.

            """
            return self._data

        @name.setter
        def name(self, name):
            self._name = name

        @format.setter
        def format(self, format):
            self._format = format

        @comment.setter
        def comment(self, comment):
            self._comment = comment

        @data.setter
        def data(self, data):
            self._data = data


def get_unit_from_fits_header(header, key):
    """Get Unit From FITS Header.

    Return coordinate unit corresponding to column with name ``key``.

    Parameters
    ----------
    header : FITS header
        Header information
    key : str
        Column name

    Raises
    ------
    IndexError
        If key not in header

    Returns
    -------
    astropy.units.Unit
        Unit object

    """
    # Loop over column names to find key
    idx = 1
    idx_found = -1
    while True:
        ttype_idx = f'TTYPE{idx}'
        if ttype_idx not in header:
            # Reached beyond last column
            break
        else:
            if header[ttype_idx] == key:
                # Found correct column
                idx_found = idx
                break
        idx += 1

    if idx_found == -1:
        raise IndexError(f'Column \'{key}\' not found in FITS header')

    # Extract coordinate unit string from header
    tcunit_idx = f'TCUNI{idx}'
    if tcunit_idx not in header:
        raise IndexError(
            f'No coordinate unit found for column \'{key}\''
            ' in FITS header'
        )
    unit_str = header[tcunit_idx]

    # Transform to Unit object
    u = units.Unit(unit_str)

    return u
