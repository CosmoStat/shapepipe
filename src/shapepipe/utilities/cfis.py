"""CFIS TOOLS.

CFIS module.

:Authors: Martin Kilbinger

"""

import glob
import os
import re
import sys

import astropy.coordinates as coords
import numpy as np
import pylab as plt
from astropy.io import ascii

from shapepipe.utilities.file_system import mkdir

unitdef = 'degree'

# Maybe define class for these constants?
size = {}
size['tile'] = 0.5
size['weight'] = 0.5
size['exposure'] = 1.0

# Cut criteria for exposures
exp_time_min = 95
flag_valid = 'V'


class param:
    """Param Class.

    General class to store (default) variables.

    """

    def __init__(self, **kwds):
        self.__dict__.update(kwds)

    def print(self, **kwds):
        """Print."""
        print(self.__dict__)

    def var_list(self, **kwds):
        """Get Variable List."""
        return vars(self)


class CfisError(Exception):
    """Cfis Error.

    Generic error that is raised by this script.

    """

    pass


class image():
    """Image Class.

    Class to store and create image information.

    Parameters
    ----------
    name : str
        file name
    ra : Angle
        right ascension
    dec : Angle
        declination
    exp_time : int, optional, default=-1
        exposure time
    valid : str, optional, default='Unknown'
        validation flag

    """

    def __init__(self, name, ra, dec, exp_time=-1, valid='Unknown'):
        self.name = name
        self.ra = ra
        self.dec = dec
        if exp_time is None:
            self.exp_time = -1
        else:
            self.exp_time = exp_time
        if valid is None:
            self.valid = 'Unknown'
        else:
            self.valid = valid

    def cut(self, no_cuts=False):
        """Cut.

        Return True (False) if image does (not) need to be cut from selection.

        Parameters
        ----------
        no_cuts : bool, optiona, default=False
            do not cut if True

        Returns
        -------
        bool
            True (False) if image is (not) cut

        """
        # Do not cut if no_cuts flag is set
        if no_cuts:
            return False

        # Cut if exposure time smaller than minimum (and not flagged as
        # unknown or n/a)
        if self.exp_time < exp_time_min and self.exp_time != -1:
            return True

        # Cut if validation flag is not valid (and not unknown)
        if self.valid != flag_valid and self.valid != 'Unknown':
            return True

        return False

    def get_ID(self):
        """Get ID.

        Return image ID.

        Returns
        -------
        str
            image iD

        Raises
        ------
        ValueError
            if name does not match to ID pattern

        """
        m = re.search(r'(\d{3}).{1}(\d{3})', self.name)
        if m is None:
            raise ValueError(f'No ID match in file name {name}')
        else:
            return f'{m[1]}.{m[2]}'

    def print(
        self,
        file=sys.stdout,
        base_name=False,
        name_only=True,
        ID_only=False
    ):
        """Print.

        Print image information as ascii Table column.

        Parameters
        ----------
        file : file, optional, default=sys.stdout
            output file handle
        base_name : bool, optional, default=False
            if True (False), print image base name (full path)
        name_only : bool, optional, default=False
            if True, do not print metainfo
        ID_only : bool, optional, default=False
            if True, only print file ID instead of entire name

        Raises
        ------
        ValueError
            if name does not match to ID pattern

        """
        if base_name:
            name = os.path.basename(self.name)
        else:
            name = self.name

        if ID_only:
            m = re.search(r'\d{3}.\d{3}', name)
            if m is None:
                raise ValueError(f'No ID match in file name {name}')
            else:
                name = m[0]
        print(name, end='', file=file)

        if not name_only:
            if self.ra is not None:
                ra_unit = getattr(self.ra, unitdef)
                print(f' {ra_unit:10.2f}', end='', file=file)
            if self.dec is not None:
                dec_unit = getattr(self.dec, unitdef)
                print(f' {dec_unit:10.2f}', end='', file=file)
            print(f' {self.exp_time:5d} {self.valid:8s}', end='', file=file)
        print(file=file)

    def print_header(self, file=sys.stdout):
        """Print Header.

        Print header for ascii Table output

        Parameters
        ----------
        file : file, optional, default=sys.stdout
            output file handle

        """
        print(
            f'#Name ra[{unitdef}] dec[{unitdef}] exp_time[s] validation',
            file=file
        )


def log_command(argv, name=None, close_no_return=True):
    """Log Command.

    Write command with arguments to a file or stdout.
    Choose name = 'sys.stdout' or 'sys.stderr' for output on sceen.

    Parameters
    ----------
    argv : numpy.ndarray of str
        Command line arguments
    name : str
        Output file name (default: 'log_<command>')
    close_no_return : bool
        If True (default), close log file. If False, keep log file open
        and return file handler

    Returns
    -------
    filehandler
        log file handler (if close_no_return is False)

    """
    if name is None:
        name = 'log_' + os.path.basename(argv[0])

    if name == 'sys.stdout':
        f = sys.stdout
    elif name == 'sys.stderr':
        f = sys.stderr
    else:
        f = open(name, 'w')

    for a in argv:

        # Quote argument if special characters
        if ']' in a or ']' in a:
            a = f'\"{a}\"'

        print(a, end='', file=f)
        print(' ', end='', file=f)

    print('', file=f)

    if not close_no_return:
        return f

    if name != 'sys.stdout' and name != 'sys.stderr':
        f.close()


def my_string_split(string, num=-1, verbose=False, stop=False, sep=None):
    """My String Split.

    Split a *string* into a list of strings. Choose as separator
    the first in the list [space, underscore] that occurs in the string.
    (Thus, if both occur, use space.)

    Parameters
    ----------
    string : str
        Input string
    num : int
        Required length of output list of strings, -1 if no requirement.
    verbose : bool
        Verbose output
    stop : bool
        Stop programs with error if True, return None and continues otherwise
    sep : bool
        Separator, try ' ', '_', and '.' if None (default)

    Raises
    ------
    CfisError
        If number of elements in string and num are different, for stop=True
    ValueError
        If no separator found in string

    Returns
    -------
    list
        List of string on success, and None if failed

    """
    if string is None:
        return None

    if sep is None:
        has_space = string.find(' ')
        has_underscore = string.find('_')
        has_dot = string.find('.')

        if has_space != -1:
            my_sep = ' '
        elif has_underscore != -1:
            my_sep = '_'
        elif has_dot != -1:
            my_sep = '.'
        else:
            # no separator found, does string consist of only one element?
            if num == -1 or num == 1:
                my_sep = None
            else:
                raise Valueerror(
                    'No separator (\' \', \'_\', or \'.\') found in string'
                    + f' \'{string}\', cannot split'
                )
    else:
        if not string.find(sep):
            raise ValueError(
                f'No separator \'{sep}\' found in string \'{string}\', '
                + 'cannot split'
            )
        my_sep = sep

    res = string.split(my_sep)

    if num != -1 and num != len(res) and stop:
        raise CfisError(
            f'String \'{len(res)}\' has length {num}, required is {num}'
        )

    return res


def get_file_pattern(pattern, band, image_type, want_re=True, ext=True):
    """Get File Pattern.

    Return file pattern of CFIS image file.

    Parameters
    ----------
    pattern : str
        input pattern, can be ''
    band : str
        band, one of 'r', 'u'
    image_type : str
        image type, one of 'exposure', 'exposure_flag', 'exposure_flag.fz',
        'exposure_weight', 'exposure_weight.fz', 'tile', 'cat',
        'weight', 'weight.fz'
    want_re : bool, optional, default=True
        return regular expression if True
    ext : bool, optional, default=True
        if True add file extention to pattern

    Returns
    -------
    str
        output pattern

    """
    if pattern == '':
        if image_type in (
            'exposure',
            'exposure_flag',
            'exposure_flag.fz',
            'exposure_weight',
            'exposure_weight.fz'
        ):
            pattern_base = r'\d{7}p'
        else:
            pattern_base = rf'CFIS.*\.{band}'
    else:
        pattern_base = pattern

    if ext:
        if image_type == 'exposure':
            pattern_out = rf'{pattern_base}\.fits\.fz'
        elif image_type == 'exposure_flag':
            pattern_out = rf'{pattern_base}\.flag\.fits'
        elif image_type == 'exposure_flag.fz':
            pattern_out = rf'{pattern_base}\.flag\.fits\.fz'
        elif image_type == 'exposure_weight':
            pattern_out = rf'{pattern_base}\.weight\.fits'
        elif image_type == 'exposure_weight.fz':
            pattern_out = rf'{pattern_base}\.weight\.fits\.fz'
        elif image_type == 'tile':
            pattern_out = rf'{pattern_base}\.fits'
        elif image_type == 'cat':
            pattern_out = rf'{pattern_base}\.cat'
        elif image_type == 'weight':
            pattern_out = rf'{pattern_base}\.weight\.fits'
        elif image_type == 'weight.fz':
            pattern_out = rf'{pattern_base}\.weight\.fits\.fz'
        else:
            raise CfisError(f'Invalid type \'{image_type}\'')
    else:
        pattern_out = pattern_base

    if not want_re:
        pattern_out = pattern_out.replace('\\', '')

    return pattern_out


def get_tile_number_from_coord(ra, dec, return_type=str):
    """Get Tile Number From Coord.

    Return CFIS stacked image tile number covering input coordinates.
    This is the inverse to get_tile_coord_from_nixy.

    Parameters
    ----------
    ra : Angle
        right ascension
    dec : Angle
        declination
    return type : <type 'type'>
        return type, int or str

    Raises
    ------
    CfisError
        for invalid return type

    Returns
    -------
    tuple
        tile number for x and tile number for y

    """
    y = (dec.degree + 90) * 2.0
    yi = int(np.rint(y))

    x = ra.degree * np.cos(dec.radian) * 2.0
    xi = int(np.rint(x))
    if xi == 720:
        xi = 0

    if return_type == str:
        nix = f'{xi:03d}'
        niy = f'{yi:03d}'
    elif return_type == int:
        nix = xi
        niy = yi
    else:
        raise CfisError(f'Invalid return type {return_type}')

    return nix, niy


def get_tile_coord_from_nixy(nix, niy):
    """Get Tile Coord From Nixy.

    Return coordinates corresponding to tile with number (nix,niy).
    This is the inverse to get_tile_number_from_coord.

    Parameters
    ----------
    nix : str or int
        tile number for x, can be list
    niy : str or int
        tile number for y, can be list

    Returns
    -------
    tuple
        right ascension and declination

    """
    if not np.isscalar(nix):
        # Transform to int, necessary if input is string
        xi = np.array(nix).astype(int)
        yi = np.array(niy).astype(int)
    else:
        xi = int(nix)
        yi = int(niy)

    d = yi / 2 - 90
    dec = coords.Angle(d, unit='deg')
    r = xi / 2 / np.cos(dec.radian)
    ra = coords.Angle(r, unit='deg')

    return ra, dec


def get_tile_name(nix, niy, band, image_type='tile', input_format='full'):
    """Get Tile Name.

    Return tile name for given tile numbers.

    Parameters
    ----------
    nix : str
        tile number for x
    niy : str
        tile number for y
    band : str
        band, one in 'r' or 'u'
    image_type : str optional, default='tile'
        image type
    input_format : str optional, default='full'
        'full' (name) or 'ID_only' input format for image names

    Raises
    ------
    CfisError
        for invalid type

    Returns
    -------
    str
        tile name

    """
    if type(nix) is int and type(niy) is int:
        if input_format == 'ID_only':
            tile_base = f'{nix:03d}.{niy:03d}'
        else:
            tile_base = f'CFIS.{nix:03d}.{niy:03d}.{band}'
    elif type(nix) is str and type(niy) is str:
        if input_format == 'ID_only':
            tile_base = f'{nix}.{niy}'
        else:
            tile_base = f'CFIS.{nix}.{niy}.{band}'
    else:
        raise CfisError(f'Invalid type for input tile numbers {nix}, {niy}')

    if input_format == 'ID_only':
        tile_name = tile_base
    else:
        if image_type == 'tile':
            tile_name = f'{tile_base}.fits'
        elif image_type == 'weight':
            tile_name = f'{tile_base}.weight.fits'
        elif image_type == 'weight.fz':
            tile_name = f'{tile_base}.weight.fits.fz'
        else:
            raise CfisError(f'Invalid image type {image_type}')

    return tile_name


def get_tile_number(tile_name):
    """Get Tile Number.

    Return tile number of given image tile name.

    Parameters
    ----------
    tile_name : str
        tile name

    Raises
    ------
    CfisError
        if tile name does not match expected pipeline numbering scheme

    Returns
    -------
    tuple
        tile number for x and tile number for y

    """
    m = re.search(r'(\d{3})[\.-](\d{3})', tile_name)
    if m is None or len(m.groups()) != 2:
        raise CfisError(
            f'Image name \'{tile_name}\' does not match tile name syntax'
        )

    nix = m.groups()[0]
    niy = m.groups()[1]

    return nix, niy


def get_tile_number_list(tile_name_list):
    """Get Tile Number List.

    Return tile numbers of given image tiles.

    Parameters
    ----------
    tile_name_list : list of str
        tile names

    Returns
    -------
    tuple
        tile numbers for x and tile numbers for y

    """
    nix_list = []
    niy_list = []
    for tile_name in tile_name_list:
        nix, niy = get_tile_number(tile_name)
        nix_list.append(nix)
        niy_list.append(niy)

    return nix_list, niy_list


def get_log_file(path, verbose=False):
    """Get Log File.

    Return log file content.

    Parameters
    ----------
    path : str
        log file path
    verbose : bool, optional, default=False
        verbose output if True

    Raises
    ------
    CfisError
        if input path does not exist

    Returns
    -------
    list of str
        log file lines

    """
    if not os.path.isfile(path):
        raise CfisError(f'Log file \'{path}\' not found')

    f_log = open(path, 'r')
    log = f_log.readlines()
    if verbose:
        print(f'Reading log file, {len(log)} lines found')
    f_log.close()

    return log


def check_ra(ra):
    """Check RA.

    Range check of right ascension.

    Parameters
    ----------
    ra : Angle
        right ascension

    Raises
    ------
    CfisError
        for invalid ra range on input

    Returns
    -------
    bool
        result of check (True if pass, False if fail)

    """
    print(ra.deg)
    if ra.deg < 0 or ra.deg > 360:
        raise CfisError('Invalid ra, valid range is 0 < ra < 360 deg')
        return 1

    return 0


def check_dec(dec):
    """Check Dec.

    Range check of declination.

    Parameters
    ----------
    dec : Angle
        declination

    Raises
    ------
    CfisError
        for invalid dec range on input

    Returns
    -------
    bool
        result of check (True if pass, False if fail)

    """
    if dec.deg < -90 or dec.deg > 90:
        raise CfisError('Invalid dec, valid range is -90 < dec < 90 deg')
        return 1

    return 0


def get_Angle(str_coord):
    """Get Angle.

    Return Angles ra, dec from coordinate string

    Parameters
    ----------
    str_coord : string
        string of input coordinates

    Returns
    -------
    tuple
        right ascension and declination

    """
    ra, dec = my_string_split(str_coord, num=2, stop=True)

    a_ra = coords.Angle(ra)
    a_dec = coords.Angle(dec)

    return a_ra, a_dec


def get_Angle_arr(str_coord, num=-1, wrap=True, verbose=False):
    """Get Angle Arr.

    Return array of Angles from coordinate string

    Parameters
    ----------
    str_coord : str
        string of input coordinates
    num : int, optional, default=-1
        expected number of coordinates (even number)
    wrap : bool, optional, default=True
        if True, wrap ra to [0; 360)
    verbose : bool, optional, default=False
        verbose output

    Returns
    -------
    numpy.ndarray of SkyCoord
        array of sky coordinates (pairs ra, dec)

    """
    angles_mixed = my_string_split(
        str_coord,
        num=num,
        verbose=verbose,
        stop=True
    )
    n = len(angles_mixed)
    n = int(n / 2)

    angles = []
    for idx in range(n):
        ra = angles_mixed[2 * idx]
        dec = angles_mixed[2 * idx + 1]
        if wrap:
            c = coords.SkyCoord(ra, dec)
        else:
            c = param(ra=coords.Angle(ra), dec=coords.Angle(dec))
        angles.append(c)

    return angles


def read_list(fname, col=None):
    """Read List.

    Read list from ascii file.

    Parameters
    ----------
    fname : str
        ascii file name
    col : str, optional, default=None
        column name

    Returns
    -------
    list of str
        list of file names

    """
    if col is None:
        f = open(fname, 'rU', encoding='latin1')
        file_list = [x.strip() for x in f.readlines()]
        f.close()
    else:
        import pandas as pd
        dat = pd.read_csv(fname, sep=r'\s+', dtype='string', header=None)
        if col not in dat:
            col = int(col)
        file_list = dat[col]

    file_list.sort()

    return file_list


def create_image_list(fname, ra, dec, exp_time=[], valid=[]):
    """Create Image List.

    Return list of image information.

    Parameters
    ----------
    fname : list of str
        file names
    ra : list of str
        right ascension
    dec : list of str
        declination
    exp_time : list of int, optional, default=[]
        exposure time
    valid : list of str, optional, default=[]
        QSO exposure validation flag

    Returns
    -------
    list of images
        list of image information

    """
    nf = len(fname)
    nr = len(ra)
    nd = len(dec)
    if nf == 0:
        raise CfisError('No entries in file name list')
    if (nf != nr or nf != nd) and nr != 0 and nd != 0:
        raise CfisError(
            f'Lists fname, ra, dec have not same length ({nf}, {nr}, {nd})'
        )

    images = []
    for i in range(nf):
        if nr > 0 and nd > 0:
            r = Angle(f'{ra[i]} {unitdef}')
            d = Angle(f'{dec[i]} {unitdef}')
        else:
            r = None
            d = None
        if len(exp_time) > 0:
            e = exp_time[i]
        else:
            e = -1
        if len(valid) > 0:
            v = valid[i]
        else:
            v = None
        im = image(str(fname[i]), r, d, exp_time=e, valid=v)
        images.append(im)

    return images


def get_image_list(
    inp,
    band,
    image_type,
    col=None,
    input_format='full',
    verbose=False
):
    """Get Image List.

    Return list of images.

    Parameters
    ----------
    inp : str
        file name or direcory path
    band : str
        optical band
    image_type : str
        image type ('tile', 'exposure', 'cat', 'weight', 'weight.fz')
    col : str, optionalm default=None
        column name for file list input file
    input_format : str, optional, default='full'
        'full' (name) or 'ID_only' input format for image names
    verbose : bool, optional, default=False
        verbose output if True

    Returns
    -------
    list of class images
        image list

    """
    file_list = []
    ra_list = []
    dec_list = []
    exp_time_list = []
    valid_list = []

    if os.path.isdir(inp):
        if col is not None:
            raise CfisError(
                'Column name (-c option) only valid if input is file'
            )

        # Read file names from directory listing
        inp_type = 'dir'
        file_list = glob.glob(f'{os.path.abspath(inp)}/*')

    elif os.path.isfile(inp):
        if image_type in ('tile', 'weight', 'weight.fz'):
            # File names in single-column ascii file
            inp_type = 'file'
            file_list = read_list(inp, col=col)
        elif image_type == 'exposure':
            # File names and coordinates in ascii file
            inp_type = 'file'
            dat = ascii.read(inp)

            if len(dat.keys()) == 3:
                # File is exposure + coord list
                # (obtained from get_coord_CFIS_pointings.py)
                file_list = dat['Pointing']
                ra_list = dat['R.A.[degree]']
                dec_list = dat['Declination[degree]']
            elif len(dat.keys()) == 12:
                # File is log file, e.g. from
                # http://www.cfht.hawaii.edu/Science/CFIS-DATA
                #  /logs/MCLOG-CFIS.r.qso-elx.log
                # Default file separator is '|'
                for d in dat:
                    file_list.append(f'd{d["col1"]}p.fits.fz')
                    ra = re.split(r'\s*', d['col4'])[0]
                    dec = re.split(r'\s*', d['col4'])[1]
                    ang = coords.Angle('{ra} hours')
                    ra_list.append(ang.degree)
                    dec_list.append(dec)
                    exp_time = int(d['col5'])
                    exp_time_list.append(exp_time)
                    valid = re.split(r'\s*', d['col11'])[2]
                    valid_list.append(valid)
            else:
                raise CfisError(
                    f'Wrong file format, #columns={len(dat.keys())},'
                    + f' has to be 3 or 12'
                )
        else:
            raise CfisError(f'Image type \'{image_type}\' not supported')

    # Create list of objects, coordinate lists can be empty
    image_list = create_image_list(
        file_list,
        ra_list,
        dec_list,
        exp_time=exp_time_list,
        valid=valid_list
    )

    # Filter file list to match CFIS image pattern
    img_list = []
    if input_format == 'ID_only':
        pattern = get_file_pattern(r'\d{3}.\d{3}', band, image_type, ext=False)
    else:
        pattern = get_file_pattern(
            rf'CFIS.\d{{3}}.\d{{3}}\.{band}',
            band,
            image_type
        )

    for img in image_list:

        # Get link source name if symbolic link
        try:
            link_src = os.readlink(img.name)
            name = link_src
        except Exception:
            # No link, continue
            name = img.name

        m = re.findall(pattern, name)
        if len(m) != 0:
            img_list.append(img)

    if verbose and len(img_list) > 0:
        print(
            f'{len(img_list)} image files found in input {inp_type} \'{inp}\''
        )

    return img_list


def exclude(f, exclude_list):
    """Exclude.

    Return True if f is on ``exclude_list``.

    Parameters
    ----------
    f : str
        file name
    exclude_list : list of str
        list of files

    Returns
    -------
    bool
        True (False) if f is in list

    """
    return f in exclude_list


def find_image_at_coord(
    images,
    coord,
    band,
    image_type,
    no_cuts=False,
    input_format='full',
    verbose=False
):
    """Find Image At Coordinates.

    Return image covering given coordinate.

    Parameters
    ----------
    images : list of class image
        list of images
    coord : str
        coordinate ra and dec with units
    band : str
        optical band
    image_type : str
        image type ('tile', 'weight', 'weight.fz', 'exposure',
        'exposure_weight', 'exposure_weight.fz', 'exposure_flag',
        'exposure_flag.fz', 'cat')
    no_cuts : bool, optional, default=False
        no cuts (of short exposure, validation flag) if True
    input_format : str, optional, default='full'
        one of 'full', 'ID_only'
    verbose : bool, optional
        verbose output if True, default=False

    Raises
    ------
    CfisError
        if image type is not 'tile';
        if more than one tile matches
        if input image coordinates are already set

    Returns
    -------
    list of images
        Found image(s), None if none found.

    """
    ra, dec = get_Angle(coord)

    if verbose:
        print(f'Looking for image at coordinates {ra}, {dec}')

    if image_type in ('tile', 'weight', 'weight.fz'):
        nix, niy = get_tile_number_from_coord(ra, dec, return_type=int)
        tile_name = get_tile_name(
            nix,
            niy,
            band,
            image_type,
            input_format=input_format
        )

        img_found = []
        for img in images:
            if os.path.basename(img.name) == tile_name:
                # Update coordinate in image for tiles with central coordinates
                ra_c, dec_c = get_tile_coord_from_nixy(nix, niy)
                if img.ra is not None or img.dec is not None:
                    raise CfisError(
                        'Coordinates in image are already '
                        + f'set to {img.ra}, {img.rec}, '
                        + f'cannot update to {ra_c}, {ra_dec}'
                    )
                img.ra = ra_c
                img.dec = dec_c
                img_found.append(img)

        if len(img_found) != 0:
            pass
        else:
            if verbose:
                print(f'Tile with numbers ({nix}, {niy}) not found')

        if len(img_found) > 1:
            raise CfisError(f'More than one tile ({img_found}) found')

    elif image_type == 'exposure':
        sc_input = coords.SkyCoord(ra, dec)

        img_found = []
        for img in images:
            # Check distance along ra and dec from image center
            sc_img_same_ra = coords.SkyCoord(ra, img.dec)
            sc_img_same_dec = coords.SkyCoord(img.ra, dec)
            distance_ra = sc_input.separation(sc_img_same_dec)
            distance_dec = sc_input.separation(sc_img_same_ra)
            if (
                distance_ra.degree < size[image_type] / 2
                and distance_dec.degree < size[image_type] / 2
            ):
                if not img.cut(no_cuts=no_cuts):
                    img_found.append(img)

        if len(img_found) != 0:
            pass
        else:
            if verbose:
                print('No exposure image found')

    else:
        raise CfisError('Only implemented for image_type=tile')

    return img_found


def find_images_in_area(
    images,
    angles,
    band,
    image_type,
    no_cuts=False,
    verbose=False
):
    """Fine Images In Area.

    Return image list within coordinate area (rectangle)

    Parameters
    ----------
    images : list of class image
        list of images
    angles : array(2) of SkyCoord
        coordinates [(ra0, dec0), (ra1, dec1)]
    band : string
        optical band
    image_type : str
        image type ('tile', 'exposure', 'cat', 'weight', 'weight.fz')
    no_cuts : bool, optional, default=False
        no cuts (of short exposure, validation flag) if True
    verbose : bool, optional, default=False
        verbose output if True

    Returns
    -------
    list of images
        found images

    """
    found = []

    # if coordinates extend over 360:
    #  check ranges [ra_min, 360] and [0, ra_max-360]
    # if not:
    #  check range [ra_min, ra_max]
    ra_bounds = []
    threesixty = coords.Angle(360, unitdef)
    if angles[1].ra > threesixty:
        ra_bounds = [
            [angles[0].ra, threesixty],
            [coords.Angle(0, unitdef), angles[1].ra - threesixty]
        ]
    else:
        ra_bounds = [[angles[0].ra, angles[1].ra]]

    if verbose:
        print(
            'Looking for all images within rectangle, '
            + f'dec=({angles[0].dec}, {angles[1].dec}), ',
            end=''
        )
        for ra_min_max in ra_bounds:
            print(f'ra=[({ra_min_max[0]}, {ra_min_max[1]}) ', end='')
        print()

    if image_type in ('tile', 'weight', 'weight.fz'):
        for img in images:
            nix, niy = get_tile_number(img.name)
            ra, dec = get_tile_coord_from_nixy(nix, niy)

            if dec.is_within_bounds(angles[0].dec, angles[1].dec):
                within = False

                # Check whether image is in any of the ra bound pairs
                for (ra_min, ra_max) in ra_bounds:
                    if ra.is_within_bounds(ra_min, ra_max):
                        within = True
                        break
                if within:
                    if img.ra is None or img.dec is None:
                        img.ra = ra
                        img.dec = dec

                    found.append(img)

    elif image_type == 'exposure':
        for img in images:
            if img.dec.is_within_bounds(angles[0].dec, angles[1].dec):
                within = False

                # Check whether image is in any of the ra bound pairs
                for (ra_min, ra_max) in ra_bounds:
                    if ra.is_within_bounds(ra_min, ra_max):
                        within = True
                        break
                if within:
                    if not img.cut(no_cuts=no_cuts):
                        found.append(img)

    else:
        raise CfisError(f'Image type \'{image_type}\' not implemented yet')

    if verbose:
        print(f'{len(found)} images found in area')

    return found


def plot_init():
    """Plot Init.

    Initialize a plot

    """
    fs = 12
    fig, ax = plt.subplots()

    ax = plt.gca()
    ax.yaxis.label.set_size(fs)
    ax.xaxis.label.set_size(fs)

    plt.tick_params(axis='both', which='major', labelsize=fs)

    plt.rcParams.update({'figure.autolayout': True})

    return ax


def plot_area(
    images,
    angles,
    image_type,
    outbase,
    interactive,
    col=None,
    show_numbers=False,
    show_circle=True,
    show_area_border=False,
    ax=None,
    lw=None,
    save=True,
    dxy=0
):
    """Plot Area.

    Plot images within area.

    Parameters
    ----------
    images : numpy.ndarray of image
        images
    angles : numpy.ndarray (SkyCoord, 2)
        Corner coordinates of area rectangle
    image_type : str
        image type ('tile', 'exposure', 'cat', weight')
    outbase : str
        output file name base
    interactive : bool
        show plot if True
    col : str, optional, default=None
        color
    show_numbers : bool, optional, default=False
        show tile numbers if True
    show_circle : bool, optional, default=True
        plot circle center and circumference around area if True
    show_area_border : bool, optional, default=False
        plot rectangle around area
    ax : axes, optional, default None
        Init axes if None
    lw : float, optional, default None
        line width
    save : bool, optional, default=True
        save plot to pdf file if True
    dxy : float, optional, default=0
        shift

    """
    if outbase is None:
        outname = 'plot.pdf'
    else:
        outname = f'{outbase}.pdf'

    if not lw:
        my_lw = 0.1
    else:
        my_lw = lw
    color = {'tile': 'b', 'exposure': 'g', 'weight': 'r'}

    if not ax:
        ax = plot_init()

    # Field center
    n_ima = len(images)
    if n_ima > 0:
        ra_c = sum([img.ra for img in images]) / float(n_ima)
        dec_c = sum([img.dec for img in images]) / float(n_ima)
        cos_dec_c = np.cos(dec_c)
    else:
        ra_c = 0
        dec_c = 0
        cos_dec_c = 1

    if show_circle:
        # Circle around field
        dx = abs(angles[0].ra - angles[1].ra)
        dy = abs(angles[0].dec - angles[1].dec)
        dx = getattr(dx, unitdef)
        dy = getattr(dy, unitdef)
        radius = (
            max(dx, dy) / 2 + (size['exposure'] + size['tile']) * np.sqrt(2)
        )
        circle = plt.Circle(
            (ra_c.deg, dec_c.deg),
            radius,
            color='r',
            fill=False,
        )
        plt.plot(ra_c, dec_c, 'or', mfc='none', ms=3)
        ax.add_artist(circle)
    else:
        radius = 0

    if col:
        c = col
    else:
        c = color[image_type]

    for img in images:
        x = img.ra.degree
        y = img.dec.degree
        nix, niy = get_tile_number(img.name)
        if show_numbers:
            plt.text(
                x,
                y,
                f'{nix}.{niy}',
                fontsize=3,
                horizontalalignment='center',
                verticalalignment='center'
            )

        # Image boundary
        dx = size[image_type] / 2 / cos_dec_c
        dy = size[image_type] / 2
        cx, cy = square_from_centre(x, y, dx, dy, dxy=dxy)
        ax.plot(cx, cy, '-', color=c, linewidth=my_lw)

    if show_area_border:
        cx, cy = square_from_corners(angles[0], angles[1])
        ax.plot(cx, cy, 'r-.', linewidth=my_lw)

    plt.xlabel('R.A. [degree]')
    plt.ylabel('Declination [degree]')
    if outbase is not None:
        plt.title(outbase)

    # Limits
    border = 0.25
    if angles[1].ra.degree > angles[0].ra.degree:
        xm = (angles[1].ra.degree + angles[0].ra.degree) / 2
        dx = angles[1].ra.degree - angles[0].ra.degree
    else:
        dx = max(360 - angles[0].ra.deg, angles[1].ra.deg) * 2
        xm = ((angles[0].ra.deg - 360) + angles[1].ra.deg) / 2
    ym = (angles[1].dec.degree + angles[0].dec.degree) / 2
    dy = angles[1].dec.degree - angles[0].dec.degree
    lim = max(dx, dy)
    plt.xlim(xm - lim / 2 - border, xm + lim / 2 + border)
    plt.ylim(ym - lim / 2 - border, ym + lim / 2 + border)

    if save:
        print(f'Saving plot to {outname}')
        plt.savefig(outname)

    if interactive:
        plt.show()

    return ra_c, dec_c, radius


def square_from_centre(x, y, dx, dy, dxy=0):
    """Square From Centre.

    Return coordinate vectors of corners that define a closed
    square for plotting.

    Parameters
    ----------
    x : float
        x-coordinate centre
    y : float
        y-coordinate centre
    dx : float
        size in x
    dy : float
        size in y
    dxy : float, optional, default=0
        constant offset

    Returns
    -------
    tuple
        square coordinates in x and y

    """
    a = dxy
    cx = [x - dx + a, x + dx + a, x + dx + a, x - dx + a, x - dx + a]
    cy = [y - dy + a, y - dy + a, y + dy + a, y + dy + a, y - dy + a]

    return cx, cy


def square_from_corners(ang0, ang1):
    """Square From Corners.

    Return coordinate vectors of corners that define a closed
    square for plotting.

    Parameters
    ----------
    ang0 : Angle
        lower-left square coordinates
    ang1 : Angle
        upper-right square coordinates

    Returns
    -------
    tuple
        square coordinates in x and y, in unit 'unitdef'

    """
    cx = [ang0.ra, ang1.ra, ang1.ra, ang0.ra, ang0.ra]
    cy = [ang0.dec, ang0.dec, ang1.dec, ang1.dec, ang0.dec]

    cxd = [getattr(x, unitdef) for x in cx]
    cyd = [getattr(y, unitdef) for y in cy]

    return cxd, cyd
