"""! 
   @package mks.mks_helper helper class
   @author Marc Gentile
   @file mks_help.py
   Helper class
""" 

# -- Python imports
try:
    from thread import get_ident as _get_ident
except ImportError:
    from dummy_thread import get_ident as _get_ident
from _abcoll import KeysView, ValuesView, ItemsView

# -- External import
from scatalog import *                   # catalog management
from mpfg.mp_helper import *             # base Helper
from mks_version import __version__      # version information


# -------------------------------------------------------------------------------------------------
class MksHelper(Helper):

   """! Convenient utility functions """

   # -----------------------------------------------------------------------------------------------
   def __init(self):
      Helper.__init__(self)

   # -----------------------------------------------------------------------------------------------
   def get_version(self):
      """! 
          Get the version number of this mks code as text 
          @return version number of this mks code as text
      """
      return __version__ 


   # -----------------------------------------------------------------------------------------------
   def show_config_summary(self, master):
      
      if master.logging_enabled():

         # --- gfit version
         master.logger.log_info_p("\n*** mks {0} ***\n".format(self.get_version()))

         # --- Python modules
         master.logger.log_info_p("Standard Python modules:")
         try:
            np = __import__('numpy', globals(), locals(), [], -1)
            master.logger.log_info_p("- numpy {0}\t\t{1}".format(np.__version__, np.__file__))        
            sp = __import__('scipy', globals(), locals(), [], -1)
            master.logger.log_info_p("- scipy {0}\t\t{1}".format(sp.__version__, sp.__file__))           
            pyfits = __import__('pyfits', globals(), locals(), [], -1)
            master.logger.log_info_p("- pyfits {0}\t\t{1}".format(pyfits.__version__, pyfits.__file__))     
            galsim = __import__('galsim', globals(), locals(), [], -1)
            master.logger.log_info_p("- galsim {0}\t\t{1}".format(galsim.version, galsim.__file__))  
         except Exception as detail:
            master.logger.log_error_p("- some modules could not be found: {0}\n".format(detail)) 

         master.logger.log_info_p("\n")
         master.logger.log_info_p("MPF Python modules:")
         try:
            mpfg = __import__('mpfg', globals(), locals(), [], -1)
            master.logger.log_info_p("- mpfg {0}\t\t{1}".format(mpfg.__version__, mpfg.__file__))  
            mpfx = __import__('mpfx', globals(), locals(), [], -1)
            master.logger.log_info_p("- mpfx {0}\t\t{1}".format(mpfx.__version__, mpfx.__file__))  
            slog = __import__('slogger', globals(), locals(), [], -1)
            master.logger.log_info_p("- slogger {0}\t\t{1}".format(slog.__version__, slog.__file__))  
            sconf = __import__('sconfig', globals(), locals(), [], -1)
            master.logger.log_info_p("- sconfig {0}\t\t{1}".format(sconf.__version__, sconf.__file__))  
            scat = __import__('scatalog', globals(), locals(), [], -1)
            master.logger.log_info_p("- scatalog {0}\t{1}".format(scat.__version__, scat.__file__)) 
         except Exception as detail:
            master.logger.log_error_p("- some modules could not be imported: {0}\n".format(detail)) 

         master.logger.log_info_p("\n") 


   # -----------------------------------------------------------------------------------------------
   def create_from_list_dico(self, data_dico, output_directory, output_filename, \
                                   job, worker, col_list=[], 
                                   key_index_map={}, key_fmt_map={}, default_fmt="%.6e",
                                   is_ascii=True, is_sextractor=False, hdu_no=1):
      """! 
         Save a dictionary to disk as a catalog file. It is assumed a list of values is attached to
         each first-level key in the dictionary (a "list" dictionary)
         @param data_dico dictionary with the data
         @param output_directory directory where to create the file
         @param output_filename name of the file to create
         @param job an object of class MksJob to process
         @param worker instance of the worker process  
         @param col_list list of column names. If empty, take all the keys of the dictionary 
         @param key_index_map if not empty, contains a map with the preferred order for some keys
         @param key_fmt_map if not empty, contains the preferred output format of some key values
         @param default_fmt default format if not explicitly specified
         @param is_ascii if True tells the dictionary to create as ASCII format
         @param is_sextractor if set to True and catalog has text format, create SExtractor file
                              with ASCII_HEAD format
         @param hdu_no HDU index                     
      """

      try:

         # --- Create output file
         output_filepath = os.path.join(output_directory, output_filename)

         # --- Build the list of columns in the required order
         cat_col_list = []
         cat_col_fmt = []
         cat_col_comments = ""

         if len(col_list) == 0:
            col_names = data_dico.keys()
         else:
            col_names = col_list

         # --- Sort column names according to their index
         if len(key_index_map) > 0:
            sorted_index_tuples = sorted(key_index_map.items(), key=itemgetter(1, 0))
            sorted_col_names = [ sorted_index_tuples[i][0] for (k, i) in sorted_index_tuples ]
            left_col_names = [c for c in col_names if not c in sorted_col_names]
            col_names = sorted_col_names
            col_names.extend(left_col_names)

         # --- Check the column names are indeed in the dictionary
         col_names_to_remove = []
         for col_name in col_names:
            if not col_name in data_dico:
               self.print_warning( "column: {0} not found in the dictionary".format(col_name) )
               col_names_to_remove.append(col_name)

         for col_name in col_names_to_remove:
            self.print_warning( "Removing column: {0}".format(col_name) )
            col_names.remove(col_name)

         for col_name in col_names:
            col_no = col_names.index(col_name)
            ### print col_name, len(data_dico[col_name])

            if col_name in key_index_map.keys():
               cat_col_list.insert(key_index_map[col_name], col_name)

               if col_name in key_fmt_map:
                  cat_col_fmt.insert(key_index_map[col_name], key_fmt_map[col_name])
               else:
                  cat_col_fmt.insert(key_index_map[col_name], default_fmt)
            else:
               cat_col_list.append(col_name)
               cat_col_fmt.append(default_fmt)

         #print "col_list:", col_list

         # --- Insert the columns in the catalog
         col_data_list = []   
         for col_name in cat_col_list:
            #print col_name
            cat_col_comments += col_name + " " 
            col_data_list.append(numpy.asarray([]))

         for col_name in cat_col_list:
            #print col_name
            col_no = cat_col_list.index(col_name)
            col_data_list[col_no] = numpy.concatenate(  
                                                    (col_data_list[col_no], data_dico[col_name] ) ) 
         data_matrix = numpy.asmatrix(col_data_list).transpose().squeeze()

         #print "matrix:", data_matrix

         # --- Create catalog file of the correct format
         catalog = None
         if not is_ascii:
            # --- Assume FITS format
            catalog = FITSCatalog(output_filepath, hdu_no=hdu_no)
            catalog.create_from_numpy(data_matrix, cat_col_list, ext_name=None)
         else:   
            if is_sextractor:
               # --- ASCII_HEAD SExtractor format
               catalog = SExCatalog(output_filepath)
               catalog.create_from_numpy(data_matrix, cat_col_list, col_formats=cat_col_fmt)
            else:         
               # --- ASCII tabulated format
#                catalog = TextCatalog(output_filepath)
#                catalog.create_from_numpy(data_matrix, cat_col_list)
               numpy.savetxt(output_filepath, data_matrix, 
                                              fmt=cat_col_fmt, header=cat_col_comments)


      except:
         self.print_error("could not create catalog from dictionary ({0})".format(sys.exc_info()))

      finally:
         if not catalog is None:
            catalog.close()   

#    # -----------------------------------------------------------------------------------------------
#    def _generate_shear_value(self, mean, sigma):
# 
#       base_data = numpy.arange(0, 10, 1)
#       base_shear = (base_data - numpy.mean(base_data)) / 100.0
#       random_shear = numpy.random.normal(mean, sigma, 10)
#       shear_values = (base_shear + random_shear)
#       numpy.random.shuffle(shear_values)
#       return shear_values[0]
# 
#    # -----------------------------------------------------------------------------------------------
#    def _generate_shear_values(self, nb_shear, mean, sigma):
# 
#       return [self._generate_shear_value(mean, sigma) for i in xrange(0, nb_shear)]


   # ~~~~~~~~~~~~~~~
   # Private methods 
   # ~~~~~~~~~~~~~~~


# --------------------------------------------------------------------------------------------------
class OrderedDict(dict):
    'Dictionary that remembers insertion order'
    # An inherited dict maps keys to values.
    # The inherited dict provides __getitem__, __len__, __contains__, and get.
    # The remaining methods are order-aware.
    # Big-O running times for all methods are the same as for regular dictionaries.

    # The internal self.__map dictionary maps keys to links in a doubly linked list.
    # The circular doubly linked list starts and ends with a sentinel element.
    # The sentinel element never gets deleted (this simplifies the algorithm).
    # Each link is stored as a list of length three:  [PREV, NEXT, KEY].

    def __init__(self, *args, **kwds):
        '''Initialize an ordered dictionary.  Signature is the same as for
        regular dictionaries, but keyword arguments are not recommended
        because their insertion order is arbitrary.

        '''
        if len(args) > 1:
            raise TypeError('expected at most 1 arguments, got %d' % len(args))
        try:
            self.__root
        except AttributeError:
            self.__root = root = []                     # sentinel node
            root[:] = [root, root, None]
            self.__map = {}
        self.__update(*args, **kwds)

    def __setitem__(self, key, value, dict_setitem=dict.__setitem__):
        'od.__setitem__(i, y) <==> od[i]=y'
        # Setting a new item creates a new link which goes at the end of the linked
        # list, and the inherited dictionary is updated with the new key/value pair.
        if key not in self:
            root = self.__root
            last = root[0]
            last[1] = root[0] = self.__map[key] = [last, root, key]
        dict_setitem(self, key, value)

    def __delitem__(self, key, dict_delitem=dict.__delitem__):
        'od.__delitem__(y) <==> del od[y]'
        # Deleting an existing item uses self.__map to find the link which is
        # then removed by updating the links in the predecessor and successor nodes.
        dict_delitem(self, key)
        link_prev, link_next, key = self.__map.pop(key)
        link_prev[1] = link_next
        link_next[0] = link_prev

    def __iter__(self):
        'od.__iter__() <==> iter(od)'
        root = self.__root
        curr = root[1]
        while curr is not root:
            yield curr[2]
            curr = curr[1]

    def __reversed__(self):
        'od.__reversed__() <==> reversed(od)'
        root = self.__root
        curr = root[0]
        while curr is not root:
            yield curr[2]
            curr = curr[0]

    def clear(self):
        'od.clear() -> None.  Remove all items from od.'
        try:
            for node in self.__map.itervalues():
                del node[:]
            root = self.__root
            root[:] = [root, root, None]
            self.__map.clear()
        except AttributeError:
            pass
        dict.clear(self)

    def popitem(self, last=True):
        '''od.popitem() -> (k, v), return and remove a (key, value) pair.
        Pairs are returned in LIFO order if last is true or FIFO order if false.

        '''
        if not self:
            raise KeyError('dictionary is empty')
        root = self.__root
        if last:
            link = root[0]
            link_prev = link[0]
            link_prev[1] = root
            root[0] = link_prev
        else:
            link = root[1]
            link_next = link[1]
            root[1] = link_next
            link_next[0] = root
        key = link[2]
        del self.__map[key]
        value = dict.pop(self, key)
        return key, value

    # -- the following methods do not depend on the internal structure --

    def keys(self):
        'od.keys() -> list of keys in od'
        return list(self)

    def values(self):
        'od.values() -> list of values in od'
        return [self[key] for key in self]

    def items(self):
        'od.items() -> list of (key, value) pairs in od'
        return [(key, self[key]) for key in self]

    def iterkeys(self):
        'od.iterkeys() -> an iterator over the keys in od'
        return iter(self)

    def itervalues(self):
        'od.itervalues -> an iterator over the values in od'
        for k in self:
            yield self[k]

    def iteritems(self):
        'od.iteritems -> an iterator over the (key, value) items in od'
        for k in self:
            yield (k, self[k])

    def update(*args, **kwds):
        '''od.update(E, **F) -> None.  Update od from dict/iterable E and F.

        If E is a dict instance, does:           for k in E: od[k] = E[k]
        If E has a .keys() method, does:         for k in E.keys(): od[k] = E[k]
        Or if E is an iterable of items, does:   for k, v in E: od[k] = v
        In either case, this is followed by:     for k, v in F.items(): od[k] = v

        '''
        if len(args) > 2:
            raise TypeError('update() takes at most 2 positional '
                            'arguments (%d given)' % (len(args),))
        elif not args:
            raise TypeError('update() takes at least 1 argument (0 given)')
        self = args[0]
        # Make progressively weaker assumptions about "other"
        other = ()
        if len(args) == 2:
            other = args[1]
        if isinstance(other, dict):
            for key in other:
                self[key] = other[key]
        elif hasattr(other, 'keys'):
            for key in other.keys():
                self[key] = other[key]
        else:
            for key, value in other:
                self[key] = value
        for key, value in kwds.items():
            self[key] = value

    __update = update  # let subclasses override update without breaking __init__

    __marker = object()

    def pop(self, key, default=__marker):
        '''od.pop(k[,d]) -> v, remove specified key and return the corresponding value.
        If key is not found, d is returned if given, otherwise KeyError is raised.

        '''
        if key in self:
            result = self[key]
            del self[key]
            return result
        if default is self.__marker:
            raise KeyError(key)
        return default

    def setdefault(self, key, default=None):
        'od.setdefault(k[,d]) -> od.get(k,d), also set od[k]=d if k not in od'
        if key in self:
            return self[key]
        self[key] = default
        return default

    def __repr__(self, _repr_running={}):
        'od.__repr__() <==> repr(od)'
        call_key = id(self), _get_ident()
        if call_key in _repr_running:
            return '...'
        _repr_running[call_key] = 1
        try:
            if not self:
                return '%s()' % (self.__class__.__name__,)
            return '%s(%r)' % (self.__class__.__name__, self.items())
        finally:
            del _repr_running[call_key]

    def __reduce__(self):
        'Return state information for pickling'
        items = [[k, self[k]] for k in self]
        inst_dict = vars(self).copy()
        for k in vars(OrderedDict()):
            inst_dict.pop(k, None)
        if inst_dict:
            return (self.__class__, (items,), inst_dict)
        return self.__class__, (items,)

    def copy(self):
        'od.copy() -> a shallow copy of od'
        return self.__class__(self)

    @classmethod
    def fromkeys(cls, iterable, value=None):
        '''OD.fromkeys(S[, v]) -> New ordered dictionary with keys from S
        and values equal to v (which defaults to None).

        '''
        d = cls()
        for key in iterable:
            d[key] = value
        return d

    def __eq__(self, other):
        '''od.__eq__(y) <==> od==y.  Comparison to another OD is order-sensitive
        while comparison to a regular mapping is order-insensitive.

        '''
        if isinstance(other, OrderedDict):
            return len(self)==len(other) and self.items() == other.items()
        return dict.__eq__(self, other)

    def __ne__(self, other):
        return not self == other

    # -- the following methods are only used in Python 2.7 --

    def viewkeys(self):
        "od.viewkeys() -> a set-like object providing a view on od's keys"
        return KeysView(self)

    def viewvalues(self):
        "od.viewvalues() -> an object providing a view on od's values"
        return ValuesView(self)

    def viewitems(self):
        "od.viewitems() -> a set-like object providing a view on od's items"
        return ItemsView(self)


# -- EOF mks_helper.py
