#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Module stuff.py

Various routines

:Authors: Martin Kilbinger

:Date: 19/01/2018
"""


# Compability with python2.x for x>6
from __future__ import print_function


import sys
import os
import re
import subprocess
import shlex
import errno
import glob

from optparse import IndentedHelpFormatter
import textwrap



class param:
    """General class to store (default) variables
    """

    def __init__(self, **kwds):
        self.__dict__.update(kwds)

    def print(self, **kwds):
        print(self.__dict__)

    def var_list(self, **kwds):
        return vars(self)



class MyError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)



def error(str, val=1, stop=True, verbose=True):
    """Print message str and exits program with code val.
       See [ABC:]covest.py for 2.7 version.

    Parameters
    ----------
    str: string
        message
    val: integer
        exit value, default=1
    stop: boolean
        stops program if True (default), continues if False
    verbose: boolean
        verbose output if True (default)

    Returns
    -------
    None
    """

    if verbose is True:
        print_color('red', str, file=sys.stderr, end='')

    if stop is False:
        if verbose is True:
            print_color('red', ', continuing', file=sys.stderr)
    else:
        if verbose is True:
            print_color('red', '', file=sys.stderr)
        sys.exit(val)



def mkdir_p(path, verbose=False):
    """Create diretorcy by calling os.makedirs. Emulate shell function 'mkdir -p'.

    Parameters
    ----------
    path: string
        Directory name
    verbose: boolean
        Verbose mode, default False

    Returns
    -------
    None
    
    """

    if verbose is True:
        print('Creating directory \'{}\''.format('{}'.format(path)))

    try:
        os.makedirs(str(path))
    except OSError as exc:
        if exc.errno == errno.EEXIST:
            pass
        else: raise



def print_color(color, txt, file=sys.stdout, end='\n'):
    """Print text with color. If not supported, print standard text.

    Parameters
    ----------
    color: string
        color name
    txt: string
        message
    file: file handler
        output file handler, default=sys.stdout
    end: string
        end string, default='\n'

    Returns
    -------
    None
    """

    try:
        import colorama
        colors = {'red'    : colorama.Fore.RED,
                  'green'  : colorama.Fore.GREEN,
                  'blue'   : colorama.Fore.BLUE,
                  'yellow' : colorama.Fore.YELLOW,
                  'black'  : colorama.Fore.BLACK,
                 }

        if colors[color] is None:
            col = colorama.Fore.BLACK
        else:
            col = colors[color]

        print(col + txt + colors['black'] + '', file=file, end=end)

    except ImportError:
        print(txt, file=file, end=end)



def my_string_split(string, num=-1, verbose=False, stop=False):
    """Split a *string* into a list of strings. Choose as separator
        the first in the list [space, underscore] that occurs in the string.
        (Thus, if both occur, use space.)

    Parameters
    ----------
    string: string
        Input string
    num: int
        Required length of output list of strings, -1 if no requirement.
    verbose: bool
        Verbose output
    stop: bool
        Stop programs with error if True, return None and continues otherwise

    Raises
    ------
    MyError
        If number of elements in string and num are different, for stop=True

    Returns
    -------
    list_str: string, array()
        List of string on success, and None if failed.
    """

    if string is None:
        return None

    has_space      = string.find(' ')
    has_underscore = string.find('_')
    has_dot        = string.find('.')

    if has_space != -1:
        sep = ' '
    elif has_underscore != -1:
        sep = '_'
    elif has_dot != -1:
        sep = '.'
    else:
        # no separator found, does string consist of only one element?
        if num == -1 or num == 1:
            sep = None
        else:
            error('No separator (\' \', \'_\', or \'.\') found in string \'{}\', cannot split'.format(string))

    #res = string.split(sep=sep) # python v>=3?
    res = string.split(sep)

    if num != -1 and num != len(res) and stop==True:
        raise MyError('String \'{}\' has length {}, required is {}'.format(string, len(res), num))

    return res



def substitute(dat, key, val_old, val_new, sep='='):
    """Performs a substitution val_new for val_old as value corresponding to key.
       See run_csfisher_cut_bins.py

    Parameters
    ----------
    dat: string
        file content
    key: string
        key
    val_old: n/a
        old value
    val_new: n/a
        new value
    sep: character
        separater between key and value, optional, default '='

    Returns
    -------
    dat: string
        file content after substitution
    """

    str_old = '{}\s*{}\s*{}'.format(key, sep, val_old)
    str_new = '{}\t\t{} {}'.format(key, sep, val_new)

    #print('Replacing \'{}\' -> \'{}\''.format(str_old, str_new))

    dat, n  = re.subn(str_old, str_new, dat)

    if n != 1:
        msg = 'Substitution {} -> {} failed, {} entries replaced'.format(str_old, str_new, n)
        error(msg, val=1)

    return dat



def substitute_arr(dat, key, val_old, val_new):
    """Performs a substitution val_new for val_old as values in an array corresponding to key.

    Parameters
    ----------
    dat: string
        file content
    key: string
        key
    val_old: n/a
        old value
    val_new: n/a
        new value

    Returns
    -------
    dat: string
        file content after substitution
    """

    n = 0
    str_old = '({}\s*=\s*[\[,].*){}(.*[,\]])'.format(key, val_old)
    str_new = r'\1{}\2'.format(val_new)

    print('Replacing \'{}\' -> \'{}\''.format(str_old, str_new))

    n_tries = 0
    while True:
        dat, n  = re.subn(str_old, str_new, dat)

        # If first time there is no substitution -> error
        if n_tries == 0 and n != 1:
            msg = 'Substitution {} -> {} failed, {} entries replaced'.format(str_old, str_new, n)
            error(msg, val=1)

        n_tries = n_tries + 1

        # After first time, break when no more subsitutions
        if n == 0:
            break

    return dat



def add_to_arr(dat, key, val, empty=False):
    """Adds a value to an existing value array for a given key.

    Parameters
    ----------
    dat: string
        file content
    key: string
        key
    val: n/a
        value to add
    empty: bool, optional, default=False
        if True, original array assumed empty

    Returns
    -------
    dat: string
        file content after substitution
    """

    n = 0

    str_old = '({}\s*=.*)\]'.format(key)
    if empty:
        str_new = r'\1 {}]'.format(val)
    else:
        str_new = r'\1, {}]'.format(val)
    dat, n  = re.subn(str_old, str_new, dat)

    if n != 1:
        msg = 'Substitution {} -> {} failed, {} entries replaced'.format(str_old, str_new, n)
        error(msg, val=1)

    return dat



class IndentedHelpFormatterWithNL(IndentedHelpFormatter):
  """Allows newline to have effect in option help.
     From https://groups.google.com/forum/#!msg/comp.lang.python/bfbmtUGhW8I/sZkGryaO8gkJ
     Usage: parser = OptionParser(usage=usage, formatter=stuff.IndentedHelpFormatterWithNL())
  """
  def format_description(self, description):
    if not description: return ""
    desc_width = self.width - self.current_indent
    indent = " "*self.current_indent
# the above is still the same
    bits = description.split('\n')
    formatted_bits = [
      textwrap.fill(bit,
        desc_width,
        initial_indent=indent,
        subsequent_indent=indent)
      for bit in bits]
    result = "\n".join(formatted_bits) + "\n"
    return result

  def format_option(self, option):
    # The help for each option consists of two parts:
    #   * the opt strings and metavars
    #   eg. ("-x", or "-fFILENAME, --file=FILENAME")
    #   * the user-supplied help string
    #   eg. ("turn on expert mode", "read data from FILENAME")
    #
    # If possible, we write both of these on the same line:
    #   -x    turn on expert mode
    #
    # But if the opt string list is too long, we put the help
    # string on a second line, indented to the same column it would
    # start in if it fit on the first line.
    #   -fFILENAME, --file=FILENAME
    #       read data from FILENAME
    result = []
    opts = self.option_strings[option]
    opt_width = self.help_position - self.current_indent - 2
    if len(opts) > opt_width:
      opts = "%*s%s\n" % (self.current_indent, "", opts)
      indent_first = self.help_position
    else: # start help on same line as opts
      opts = "%*s%-*s  " % (self.current_indent, "", opt_width, opts)
      indent_first = 0
    result.append(opts)
    if option.help:
      help_text = self.expand_default(option)
# Everything is the same up through here
      help_lines = []
      for para in help_text.split("\n"):
        help_lines.extend(textwrap.wrap(para, self.help_width))
# Everything is the same after here
      result.append("%*s%s\n" % (
        indent_first, "", help_lines[0]))
      result.extend(["%*s%s\n" % (self.help_position, "", line)
        for line in help_lines[1:]])
    elif opts[-1] != "\n":
      result.append("\n")
    return "".join(result)



def log_command(argv, name=None, close_no_return=True):
    """Write command with arguments to a file or stdout.
       Choose name = 'sys.stdout' or 'sys.stderr' for output on sceen.

    Parameters
    ----------
    argv: array of strings
        Command line arguments
    name: string
        Output file name (default: 'log_<command>')
    close_no_return: bool
        If True (default), close log file. If False, keep log file open
        and return file handler

    Returns
    -------
    log: filehandler
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
            a = '\"{}\"'.format(a)

        print(a, end='', file=f)
        print(' ', end='', file=f)

    print('', file=f)

    if close_no_return == False:
        return f

    if name != 'sys.stdout' and name != 'sys.stderr':
        f.close()



def run_cmd(cmd_list, run=True, verbose=True, stop=False, parallel=True, file_list=None, devnull=False, env=None):
    """Run shell command or a list of commands using subprocess.Popen().

    Parameters
    ----------

    cmd_list: string, or array of strings
        list of commands
    run: bool
        If True (default), run commands. run=False is for testing and debugging purpose
    verbose: bool
        If True (default), verbose output
    stop: bool
        If False (default), do not stop after command exits with error.
    parallel: bool
        If True (default), run commands in parallel, i.e. call subsequent comands via
        subprocess.Popen() without waiting for the previous job to finish.
    file_list: array of strings
        If file_list[i] exists, cmd_list[i] is not run. Default value is None
    devnull: boolean
        If True, all output is suppressed. Default is False.
    env: bool, optional, default=None
        Modified environment in which shell command will runi

    Returns
    -------
    sum_ex: int
        Sum of exit codes of all commands
    """

    if type(cmd_list) is not list:
        cmd_list = [cmd_list]

    if verbose is True and len(cmd_list) > 1:
        print('Running {} commands, parallel = {}'.format(len(cmd_list), parallel))


    ex_list   = []
    pipe_list = []
    out_list  = []
    err_list  = []

    if env is None:
        env = os.environ.copy()

    for i, cmd in enumerate(cmd_list):

        ex = 0
        out = ''
        err = ''

        if run is True:
            # Check for existing file
            if file_list is not None and os.path.isfile(file_list[i]):
                if verbose is True:
                    print_color('blue', 'Skipping command \'{}\', file \'{}\' exists'.format(cmd, file_list[i]))
            else:
                if verbose is True:
                        print_color('green', 'Running command \'{0}\''.format(cmd))

                # Run command
                try:
                    cmds = shlex.split(cmd)
                    if devnull is True:
                        pipe = subprocess.Popen(cmds, stdout=subprocess.DEVNULL, env=env)
                    else:
                        pipe = subprocess.Popen(cmds, stdout=subprocess.PIPE, env=env)
                        #out, err = pipe.communicate()
                        #print(out)
                        # See https://www.endpoint.com/blog/2015/01/28/getting-realtime-output-using-python
                        while True:
                            output = pipe.stdout.readline()
                            if output == '' and pipe.poll() is not None:
                                break
                            if output:
                                print(output.strip())
                            ex = pipe.poll()

                    if parallel is False:
                        # Wait for process to terminate
                        pipe.wait()

                    pipe_list.append(pipe)

                    # If process has not terminated, ex will be None
                    #ex = pipe.returncode
                except OSError as e:
                    print_color('red', 'Error: {0}'.format(e.strerror))
                    ex = e.errno

                    check_error_stop([ex], verbose=verbose, stop=stop)

        else:
            if verbose is True:
                print_color('yellow', 'Not running command \'{0}\''.format(cmd))

        ex_list.append(ex)
        out_list.append(out)
        err_list.append(err)


    if parallel is True:
        for i, pipe in enumerate(pipe_list):
            pipe.wait()

            # Update exit code list
            ex_list[i] = pipe.returncode
    s = check_error_stop(ex_list, verbose=verbose, stop=stop)

    #return s
    return s, out_list, err_list



def check_error_stop(ex_list, verbose=True, stop=False):
    """Check error list and stop if one or more are != 0 and stop=True

    Parameters
    ----------
    ex_list: list of integers
        List of exit codes
    verbose: boolean
        Verbose output, default=True
    stop: boolean
        If False (default), does not stop program

    Returns
    -------
    s: integer
        sum of absolute values of exit codes
    """

    if ex_list is None:
        s = 0
    else:
        s = sum([abs(i) for i in ex_list])


    # Evaluate exit codes
    if s > 0:
        n_ex = sum([1 for i in ex_list if i != 0])
        if verbose is True:
            if len(ex_list) == 1:
                print_color('red', 'The last command returned sum|exit codes|={}'.format(s), end='')
            else:
                print_color('red', '{} of the last {} commands returned sum|exit codes|={}'.format(n_ex, len(ex_list), s), end='')
        if stop is True:
            print_color('red', ', stopping')
        else:
            print_color('red', ', continuing')

        if stop is True:
            sys.exit(s)


    return s



def ln_s(orig, new, orig_to_check=False, verbose=False, force=False):
    """Create symbolic link.

    Parameters:
    -----------
    orig: string
        Name of original file
    new: string
        Name of new, to be created, link
    orig_to_check: string, optional, default=False
        Original file to check for existance (can be different than orig
        is link is not created in same directory where this function is called)
    verbose: bool
        Verbose output
    force: bool
        If True, link creation is forced even if file exists

    Returns:
    --------
    None
    """

    if orig_to_check is None or os.path.isfile(orig_to_check) or os.path.isdir(orig_to_check):
        if os.path.isfile(new) or os.path.islink(new):
            if force == False:
                if verbose:
                    print('File \'{}\' exists, skipping link creation...'.format(new))
            else:
                if verbose:
                    print('File \'{0}\' exists, deleting file and creating new link {1} <- {0}'.format(new, orig))
                os.remove(new)
                os.symlink(orig, new)
        else:
            if verbose:
                print('Creating link \'{}\' <- \'{}\''.format(orig, new))
            os.symlink(orig, new)
    else:
        if verbose:
            print('Original file \'{}\' does not exist, skipping...'.format(orig_to_check))




def get_file_list(input_dir, pattern_base, ext='.fits', verbose=False):
    """Return list of all files in directory whose names follow base pattern and pipeline numbering scheme

    Parameters
    ----------
    input_dir: string
        input directory
    pattern_base: string
        base of file name pattern
    ext: string, optional, default='.fits'
        file extension
    verbose: bool, optional, default=False
        verbose output if True

    Returns
    -------
    dst_list: list of string
        file name list
    """

    files = glob.glob('{}/*{}'.format(input_dir, ext))

    dst_list = []
    for f in files:

        # Test if file matches pattern_base
        m = re.findall(pattern_base, f)
        if len(m) != 0:
            dst_list.append(f)

    if len(dst_list) == 0:
        stuff.error('No files found in \'{}\' that matches pattern_base \'{}\''.format(input_dir, pattern_base))
    if verbose == True:
        print('Found {} tiles'.format(len(dst_list)))

    return dst_list



def get_pipe_file_number(pattern_base, file_name):
    """Returns pipeline file number.

    Parameters
    ----------
    pattern_base: string
        base file name to match
    file_name: string
        file name from where to extract number

    Returns
    -------
    num: int
        pipeline file number
    """

    pattern = re.compile('{}(\d*)-0'.format(pattern_base))
    m       = re.search(pattern, file_name)
    if m is not None:
        num = m.group(1)
    else:
        stuff.error('Could not extract number from tile file \'{}\''.format(file_name))

    return num


