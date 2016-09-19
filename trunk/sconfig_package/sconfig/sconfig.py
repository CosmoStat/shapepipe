"""! 
   @package sconfig.sconfig Simple section-based configuration file parser
   @author Marc Gentile
   @file sconfig.py
   Simple section-based configuration file parser
   Configuration file parser
"""

# -- Python imports
import os, sys
import string
import re
import imp
from operator import itemgetter, attrgetter
import numpy as np

# --- Module-specific imports
import sconfig_helper as sch

# ------------------------------------------------------------------------------
class SConfig(object):

   """!
      A simple configuration file parser 
   """

   ROOT_SECTION_NAME = "ROOT"             # name of the 'ROOT' section

   def __init__(self, filepath, from_config=None):
      """! 
         Create a new configuration file object
         @param filepath filepath of the configuration file
         @param from_config [optional] SConfig object to merge with
         @note the @c from_config parameter is not used in this version 
      """

      # Input parameters
      self._input_config_path = filepath
      self._ext_config_obj = from_config

      self._next_section_rank = 0 # Section ranking

      self._root_section = SConfig.Section(SConfig.ROOT_SECTION_NAME, None, self._next_section_rank)
      self._config_file_fd = None   # file descriptor (if file open)

      # Internal section name map: section name -> section object
      self._section_index = {SConfig.ROOT_SECTION_NAME:self._root_section}

      # Input parameter checking
      self._valid_parameters()

      # Read & Parse input configuration file
      self._read_config()
      if not from_config is None:
         # TODO
         self._update_config()

#         # Internal section dictionary (section tree)
#         self._section_dico = {SConfig.ROOT_SECTION_NAME:None}

   def __str__(self):
      """! String displayed when printing a SConfig object """
      return "Configuration file path: {0}".format(self._input_config_path)  

   # ~~~~~~~~~~~
   # Properties 
   # ~~~~~~~~~~~

   @property
   def config_path(self):
      """! @return the path of the configuration file """
      return self._input_config_path

   @property
   def config_external_config(self):
      """! @return the SConfig object to merge with
           @note this is not implemented in this version
      """
      return self._ext_config_obj

   # ~~~~~~~~~~~~~~~
   # Public methods 
   # ~~~~~~~~~~~~~~~

   def has_section(self, section_name):
      """ 
         !tells if @c section_name has been defined 
         @param section_name section name to check for existence
         @retval True if the section exists
         @retval False section not found      
      """
      return self._find_section(section_name) is not None

   def has_key(self, key, section=ROOT_SECTION_NAME):
      """! 
         Tells if @c key belongs to @c section  
         @param key key to check
         @param section section to which the key belongs to 
         @retval True if the section has the specifield key
         @retval False key not found in the specified section
         @throws SectionNotFound
      """
      section_obj = self._find_section(section)
      if not section_obj is None:
         return section_obj.has_key(key)
      else:
         raise SConfig.SectionNotFound(section)

   def get_as_string(self, key, section=ROOT_SECTION_NAME, evaluate=False):
      """! 
         Retrieve, as a Python string, the value associated with a key 
         @param key key for which the value is to be retrieved
         @param section section to which the key belongs to
         @param evaluate if @c True, invoke @c eval() on the value before returning
         @throws KeyNotFound
         @throws SectionNotFound
      """
      section_obj = self._find_section(section)
      if not section_obj is None:
         if key in section_obj.config_dico:
            value = section_obj.config_dico[key]
            if evaluate:
               return str(eval(value)) 
            else:   
               return value
         else:   
            raise SConfig.KeyNotFound(key, section_obj.name)
      else:
         raise SConfig.SectionNotFound(section)

   def get_as_boolean(self, key, section=ROOT_SECTION_NAME, evaluate=False):
      """! 
         Retrieve, as a Python boolean, the value associated with a key 
         @param key key for which the value is to be retrieved
         @param section section to which the key belongs to
         @param evaluate if @c True, invoke @c eval() on the value before returning
         @throws KeyNotFound
         @throws SectionNotFound
      """
      return self.get_as_string(key, section, evaluate).lower() in ("yes", "true", "t", "1")

   def get_as_int(self, key, section=ROOT_SECTION_NAME, evaluate=False):
      """! 
         Retrieve, as an Python integer, the value associated with a key 
         @param key key for which the value is to be retrieved
         @param section section to which the key belongs to
         @param evaluate if @c True, invoke @c eval() on the value before returning
         @throws KeyNotFound
         @throws SectionNotFound
      """
      return int(self.get_as_string(key, section, evaluate))

   def get_as_float(self, key, section=ROOT_SECTION_NAME, evaluate=False):
      """! 
         Retrieve, as a Python floating-point number, the value associated with a key 
         @param key key for which the value is to be retrieved
         @param section section to which the key belongs to
         @param evaluate if @c True, invoke @c eval() on the value before returning
         @throws KeyNotFound
         @throws SectionNotFound
      """
      return float(self.get_as_string(key, section, evaluate))

   def get_as_dict(self, key, section=ROOT_SECTION_NAME, evaluate=False):
      """! 
         Retrieve, as a Python dictionary, the value associated with a key 
         @param key key for which the value is to be retrieved
         @param section section to which the key belongs to
         @param evaluate if @c True, invoke @c eval() on the value before returning
         @throws KeyNotFound
         @throws SectionNotFound
         @note if @c evaluate is specified, @c eval() is invoked on each value in the dictionary
      """
      tokens = self.get_as_string(key, section, False)
      if evaluate:
         token_dico = eval(tokens)
         for k in token_dico.keys():
            token_dico[k] = eval(token_dico[k])
         return token_dico
      else:
         return eval(tokens)

   def get_as_list(self, key, section=ROOT_SECTION_NAME, evaluate=False):
      """! 
         Retrieve, as a Python list, the value associated with a key 
         @param key key for which the value is to be retrieved
         @param section section to which the key belongs to
         @param evaluate if @c True, invoke @c eval() on the value before returning
         @throws KeyNotFound
         @throws SectionNotFound
         @note if @c evaluate is specified, @c eval() is invoked on each element of the list
      """
      tokens = self.get_as_string(key, section, False)
      elem_list = eval(tokens)
      if evaluate:
         return [eval(e) for e in elem_list]
      else:
         return elem_list             

   def get_as_array(self, key, section=ROOT_SECTION_NAME, evaluate=False):
      """! 
         Retrieve, as a Python numpy array, the value associated with a key 
         @param key key for which the value is to be retrieved
         @param section section to which the key belongs to
         @param evaluate if @c True, invoke @c eval() on the value before returning
         @throws KeyNotFound
         @throws SectionNotFound
         @note if @c evaluate is specified, @c eval() is invoked on each element of the array
      """
      return np.asarray(self.get_as_list(key, section, evaluate))

   def get_as_tuple(self, key, section=ROOT_SECTION_NAME, evaluate=False):
      """! 
         Retrieve, as a Python tuple, the value associated with a key 
         @param key key for which the value is to be retrieved
         @param section section to which the key belongs to
         @param evaluate if @c True, invoke @c eval() on the value before returning
         @throws KeyNotFound
         @throws SectionNotFound
         @note if @c evaluate is specified, @c eval() is invoked on each element of the tuple
      """      
      return tuple(self.get_as_list(key, section, evaluate))
      
   def get_as_module(self, key, section=ROOT_SECTION_NAME):
      """! 
         Retrieve, as a Python module, the value associated with a key 
         @param key key for which the value is to be retrieved
         @param section section to which the key belongs to
         @throws KeyNotFound
         @throws SectionNotFound
      """   

      [module_name, module_path] = self.get_as_list(key, section, evaluate=False)   
      file_obj, filename, data = imp.find_module(module_name, [module_path])
      imp.load_module(module_name, file_obj, filename, data)
      return imp.load_module(module_name, file_obj, filename, data)

   def get_section_data(self, section=ROOT_SECTION_NAME):
      """! 
         Return a Python dictionary with the key and values defined in the specified section
         @param section section whose data is to be retrieved     
         @throws SectionNotFound  
      """
      section_obj = self._find_section(section)
      if not section_obj is None:
         #print("### rank_key_tuples of {0}: {1}".format(section, section_obj.rank_key_tuples))
         return section_obj.config_dico
      else:
         raise SConfig.SectionNotFound(section)

   def get_section_keys(self, section=ROOT_SECTION_NAME):
      """!
         Return the keys defined in the specified section (sorted by natural order) 
         @param section section whose data is to be retrieved     
         @throws SectionNotFound  
      """      
      section_obj = self._find_section(section)
      if not section_obj is None:
         return section_obj.keys
      else:
         raise SConfig.SectionNotFound(section)

   def get_section_values(self, section=ROOT_SECTION_NAME):
      """!
         Return the values of the keys defined in a section (sorted by natural order of keys) 
         @param section section whose data is to be retrieved     
         @throws SectionNotFound 
      """
      section_obj = self._find_section(section)
      if not section_obj is None:
         return section_obj.values
      else:
         raise SConfig.SectionNotFound(section)

   def get_key_value(self, key, section=ROOT_SECTION_NAME):
      """!
         Return, a an Python object, the value bound to a key in a section 
         @param key key for which the value is to be retrieved
         @param section section to which the key belongs to
         @throws SectionNotFound 
      """
      section_obj = self._find_section(section)
      if not section_obj is None:
         return section_obj.config_dico[key]
      else:
         raise SConfig.SectionNotFound(section)        

   def set_key_value(self, key, value, section=ROOT_SECTION_NAME, create_section=True):
      """! 
         Record in memory a value of arbitray type 
         @param key key for which the value is to be set
         @param value value to set
         @param section section to which the key belongs to
         @param create_section if @c True, also create the section to hold the (key, value) pair
         @throws InvalidSection 
         @throws SectionNotFound 
      """   
      section_obj = self._find_section(section)
      if not section_obj is None:
         section_obj.add_key_value(key, value)
      else:
         if create_section:
            section_obj = self._add_section(section, rank=self._next_section_rank)
            if not section_obj is None:            
               section_obj.add_key_value(key, value)
            else:
               raise SConfig.InvalidSection(section)
         else:
            raise SConfig.SectionNotFound(section)      

   def open(self, input_filepath):
      """!
         Open the cnfiguration file
         @param input_filepath file path of the configuration file
         @throws InputOutputError 
      """
      try:
         self._config_file_fd = open(input_filepath)
         return self._config_file_fd
      except(InputOutputError):
         raise(SConfig.InputOutputError(sys.exc_info()[1]))
      finally:
         self.close()

   def close(self):
      """! Close open configuration file """   
      if not self._config_file_fd is None:
         self._config_file_fd.close()

   def dump_to_string(self, section=ROOT_SECTION_NAME, recurse=True):
      """!
         Dump configuration file content to a string starting from a specified section
         @param section section whose content to dump
         @param recurse if @c True, walk down all the sub-sections
         @note if @c section is @c ROOT_SECTION_NAME and @c recurse is @c True, dump the entire 
               configuration file content
         @see dump_to_file
         @throws SectionNotFound
      """

      str = ""    
      subsections = []
      section_obj = self._find_section(section)
      if not section_obj is None:
         #if not section_obj.is_root(): 
            # Only section <section> (if not root)
            subsections = [section_obj]
      else:
         raise SConfig.SectionNotFound(section)

      if recurse:
         # Include subsections
         subsections.extend(self._get_subsections(section, recurse=recurse))

      for subsection in subsections:

         str += ("[{0}]\n".format(subsection.name))
         for (rank, key) in subsection.rank_key_tuples:
            str += ("{0} = {1}\r\n".format(key, subsection.config_dico[key]))   
         str += "\n"
          
      return str


   def dump_to_file(self, fd=sys.stdout, section=ROOT_SECTION_NAME, recurse=True):
      """!
         Dump configuration file content to a string starting from a specified section
         @param fd file descriptor of an open file
         @param section section whose content to dump
         @param recurse if @c True, walk down all the sub-sections
         @note if @c section is @c ROOT_SECTION_NAME and @c recurse is @c True, dump the entire 
               configuration file content
         @note if @c fd corresponds to the standard output, dump to the console instead of a file on 
               disk
         @see dump_to_string
         @throws SectionNotFound
      """
      fd.write(self.dump_to_string(section, recurse))

#   def write_section(self, section=ROOT_SECTION_NAME, flush=False):
#      """ Write content to open configuration file """   
#      pass

#   def write_to_file(output_filepath, ROOT_SECTION_NAME):
#      """ Write content to non-open configuration file """   
#      pass

#   def write_to_file(output_dir, output_filepath, ROOT_SECTION_NAME):
#      write(os.path.join(output_dir, output_filepath, ROOT_SECTION_NAME))
#         print("SConfig: file {0} no found".format(self._input_config_path))

   def get_subsections(self, section=ROOT_SECTION_NAME, recurse=False):
      """!
         Return the list of subsection names of a specified section (any level)
         @param section section whose subsection names have to be returned
         @param recurse if @c True, walk down all the sub-sections
         @note if @c section is @c ROOT_SECTION_NAME, list all the names of the top-level sections 
               of the configuration file. 
         @note If @c recurse is also @c True, return the names all the sections and subsections 
               of the configuration file at all levels 
         @throws SectionNotFound
      """
      return [s.name for s in self._get_subsections(section, recurse=recurse)]

   def is_subsection(self, section_name, parent_name):
      """! 
         Tells whether a section is a subsection of a given section 
         @param section_name name of the supposedly-child section
         @param parent_name name of the parent section 
         @retval True if @c section is a subsection of @c section
         @throws SectionNotFound
      """
      section_obj = self._find_section(section_name)
      if not section_obj is None:
         parent_section_obj = self._find_section(parent_name)
         if not parent_section_obj is None:
            return section_obj.is_subsection_of(parent_section_obj)
         else:
            raise SConfig.SectionNotFound(parent_name) 
      else:
         raise SConfig.SectionNotFound(section_name)  

   def has_subsection(self, section_name, subsection_name):
      """!
         Tells whether a section has a given subsection 
         @param section_name name of the owner section
         @param subsection_name of the supossedly-owned subsection
         @retval True if section with name @c section_name contains a subsection with name
                 @c subsection_name
         @throws SectionNotFound
      """
      section_obj = self._find_section(section_name)
      if not section_obj is None:
         subsection_obj = self._find_section(subsection_name)
         if not subsection_obj is None:
            return section_obj.has_subsection(subsection_obj)
         else:
            raise SConfig.SectionNotFound(subsection_name) 
      else:
         raise SConfig.SectionNotFound(section_name)  

   def section_inherits_from(self, section_name, mother_name):
      """!
         Tells whether a section inherits from another section (any level) 
         @param section_name name of a section
         @param mother_name name of the mother subsection
         @throws SectionNotFound
      """
      section_obj = self._find_section(section_name)
      if not section_obj is None:
         mother_obj = self._find_section(mother_name)
         if not mother_obj is None:
            return section_obj.inherits_from(mother_obj)
         else:
            raise SConfig.SectionNotFound(mother_name) 
      else:
         raise SConfig.SectionNotFound(section_name)   


   # ~~~~~~~~~~~~~~~~
   # Private methods 
   # ~~~~~~~~~~~~~~~~

   def _valid_parameters(self):
      """! 
         Validate user-specified parameters 
         @throws ConfigFileNotFound
      """
      if not sch.filepath_exists(self._input_config_path):
         raise SConfig.ConfigFileNotFound(self._input_config_path)

   def _update_config(self):
      """!
         @todo not currently implemented
      """
      pass

   def _read_config(self):
      """! Read configuration file and create data structures """
      # Parse configuration file
      if self._parse_config():
         # Interpret content
         self._interpret_config()   

   def _parse_config(self):
      """! 
         Parse each configuration line and create appropriate data structures 
      """

      with open(self._input_config_path) as fd:
   
         section_stack = []   # to keep track of the current section

         bracket_pattern = re.compile("\[(.+)\]")  # section pattern matching

         current_section = self._root_section
         current_section_invalid = False
         
         # --- Parse each line...
         token_stack = []

         for line in fd.readlines():
            line = line.strip()

            # --- Ignore blank lines and comments
            if len(line) != 0 and not line.startswith('#'):

               # --- Split line into (name,value) pairs or sections

               is_continued = '\\' in line   # continuation character found

               token_list = re.split('[#]', line, maxsplit=1)  # strip comment
               token_list = re.split('[=]', token_list[0], maxsplit=1) # parse '='

               # --- Locate section declaration if any in the line
               section_match = bracket_pattern.search(token_list[0])
               if not section_match is None and len(token_stack) == 0:

                  # A definition for a section was found
                  section_name = self._format_section_name(section_match.group(1).strip())
                     
                  # Adding section
                  self._next_section_rank += 1
                  section = self._add_section(section_name, self._next_section_rank)
                  current_section_invalid = (section is None)
                  if not section is None:
                     current_section = section
               else:
                  if is_continued:
                     if len(token_list) == 2:
                        # Remove continuation character if any
                        if '\\' in token_list[1]:
                           token_list[1] = token_list[1].strip().replace('\\', "", 1)
                        token_stack.append(token_list)
                        continue
                     elif len(token_list) == 1:
                        # partial line, without the '=' separator
                        token_list = self._join_line_tokens(token_list, token_stack)
                        token_stack.append(token_list)
                        continue
                  else:
                     if len(token_stack) > 0:
                        if len(token_list) == 1:
                           # partial line, without the '=' separator
                           token_list = self._join_line_tokens(token_list, token_stack)
                           token_stack = []
                     
                  # --- Adding (key,value) pairs
                  if not current_section_invalid:
                     # A valid [key,pair,comment] definition
                     if len(token_list) == 2:
                        current_section.add_key_value(token_list[0].strip(), 
                                                      sch.expand_env_vars(token_list[1]).strip())
                        #print 'adding (key,value) pair ({0}, {1})'.format(token_list[0].strip(), sch.expand_env_vars(token_list[1]).strip())
                     elif len(token_list) == 1:
                        raise SConfig.ValueNotFound(token_list[0], current_section.name)
                  else:
                     raise SConfig.InvalidSection(current_section.name)
                     #print("SConfig *** ERROR ***: section [{0}] is invalid => corresponding data not recorded".format(current_section.name))

                  #print(current_section)                     

         # --- end for

   def _join_line_tokens(self, token_list, token_stack):
      token_list[0] = token_list[0].replace('\\', '', 1)
      last_token_list = token_stack.pop()
      return [last_token_list[0], last_token_list[1] + token_list[0]]

   def _interpret_config(self):
      """!
         @todo not currently implemented
      """
      pass  

   def _section_exists(self, section_name):
      """! Check if a section has already been defined """
      return section_name in self._section_index

   def _get_subsections(self, parent_name, recurse=False):
      """! Return the tree of subsections of a parent section as a list """
      parent_obj = self._find_section(parent_name)
      if not parent_obj is None:
         if recurse:
            children = []
            self._get_section_tree(parent_obj, children)
            return sorted(children, key=attrgetter('rank'))  # sort / rank
         else:
            return parent_obj.subsections
      else:
         raise SConfig.SectionNotFound(parent_name)  
     
   def _get_section_tree(self, section_obj, children):
      """! Enumerate all children sections (subsections) of a section """
      for obj in section_obj.subsections:
         children.extend([self._get_section_tree(obj, children)])
      return section_obj

#      if section_obj.subsections != []:
#         for obj in section_obj.subsections:
#            #print "### section_obj:", obj.name, "children:", obj.subsections
#            children.extend([self._get_section_tree(obj, children)])   
#         return section_obj
#      else:
#         #print "section obj: {0} has NO children".format(section_obj.name)
#         return section_obj

   def _format_section_name(self, section_name):
      """! Make sure a section name is properly formatted"""
      if '.' in section_name:
         section_tokens = [token.strip() for token in re.split("[.]", section_name)]
         new_section_name = "" 
         for token in section_tokens:
            new_section_name += "{0}.".format(token)     
         return new_section_name[:-1]
      else:
         return section_name.strip()         

   def _add_section(self, section_expr, rank):
      """! Create and register a new section """
      name, parent_name, section_name, base_section_names = self._parse_section_name(section_expr)
#      print("@@@ PARSE: section: [{0}]  name: [{1}] parent name: [{2}] section_name: [{3}] base_section_names: {4}".format(section_expr, name, parent_name, section_name, base_section_names))
      section = None

      # Check parent and create section  
      parent_section = self._find_section(parent_name)
      if not parent_section is None:   
         if parent_section.name == parent_name:
#            print("### Creating section: [{0}]...".format(section_name))   
            section = SConfig.Section(section_name, parent_section, rank)   # create section object
            parent_section.subsections.append(section)
            self._section_index[section_name] = section        # record section in index
         else:
            print("SConfig *** ERROR ***: parents for section [{0}] do not match ([{1}] versus [{2}]) => section creation aborted".format(section.name, parent_section.name, parent_name))
            raise SConfig.InvalidSection(parent_name) 
      else:
         #print("SConfig *** ERROR ***: subsection: [{0}] has unknown parent [{1}] => section creation aborted".format(name, parent_name)) 
         raise SConfig.SectionNotFound(parent_name)

      # Copy (key,value) pairs from any specified base section(s) to inherit from
      if len(base_section_names) > 0:
         #print("### copy_section(): {0} to [{1}]".format(base_section_names, section_name) )
         self._copy_sections(base_section_names, section_name)

#      print("created: {0}".format(section))

      return section

   def _copy_sections(self, src_section_names, dst_section_name):
      """!
         Copy the missing (key, value) pairs in destination section from a list of source sections
      """
      dst_section = self._find_section(dst_section_name)
      if not dst_section is None:
#         dst_section = self._add_section(dst_section_name)
         for src_section_name in src_section_names:
            # Copy missing (key,value) pairs from source section(s)
#            print("### Copying {0}...".format(src_section_name))
            src_section = self._find_section(src_section_name)
            if not src_section is None:
               dst_section._copy_data_from(src_section)  # copy mother's dico
               dst_section.mothers.append(src_section)   # mother section obj
#               print('*** Mothers of : {0}: {1}'.format(dst_section.name, [s.name for s in dst_section.mothers] )) 
            else:
               raise SConfig.SectionNotFound(src_section_name)
      else:
         raise SConfig.SectionNotFound(dst_section_name)

   def _parse_section_name(self, section_def):
      """! Parse section expression and extract (section, parent section) """
      base_section_names = []      
      if ':' in section_def:
         # Section declared as inheriting from a base section
         section_tokens = [token.strip() for token in section_def.rsplit(':', 1)]
         if len(section_tokens) == 2:
            [section_expr, base_section_expr] = section_tokens
            base_section_names = [token.strip() for token in re.split("[,]", base_section_expr.strip())]
      else:
         section_expr = section_def.strip()

      section_tokens = [SConfig.ROOT_SECTION_NAME] + section_expr.rsplit('.', 1)
      return section_tokens[-1], section_tokens[-2], section_expr, base_section_names
         
   def _find_section(self, section_name):
      """! Find and return an already-defined section. Return None if not found """
      return self._section_index.get(self._format_section_name(section_name), None)

   def _find_section_safe(self, section_name):   
      """! Find and return an already-defined section. Raise SectionNotFound exception if not found """
      section_obj = self._find_section_safe(section_name)
      if not section_obj is None:
         raise SConfig.SectionNotFound(section_name)
      else:
         return section_obj

   def _get_root_section(self):
      """! Returns the 'ROOT' section object """   
      return self._root_section    


   # ------------------------------------------------------------------------------
   class ConfigFileNotFound(Exception):
      """!
         Exception thrown when a configuration file is not found on disk 
      """
      def __init__(self, filepath):
         """!
            Exception constructor
            @param filepath file path of the configuration file 
         """
         self._filepath = filepath

      def __str__(self):
         """! String representation of the exception object """
         return "SConfig *** ERROR ***: configuration file {0} no found".format(self._filepath)

   # ------------------------------------------------------------------------------
   class KeyNotFound(Exception):
      """!
         Exception thrown when a key could not be found in a given section
      """
      def __init__(self, key, section_name):
         """!
            Exception constructor
            @param key that could not be found
            @param section_name name of the section
         """
         self._key = key
         self._section_name = section_name

      def __str__(self):
         """! String representation of the exception object """
         return "SConfig *** ERROR ***: key {0} no found in section [{1}]".format(
                                                                    self._key, self._section_name)

   # ------------------------------------------------------------------------------
   class ValueNotFound(Exception):
      """!
         Exception thrown when a value associated with a key could not be found in a given section
      """
      def __init__(self, key, section_name):
         """!
            Exception constructor
            @param key whose value could not be found
            @param section_name name of the section
         """
         self._key = key
         self._section_name = section_name

      def __str__(self):
         """! String representation of the exception object """
         return "SConfig *** ERROR ***: no value found for key {0} in section [{1}]".format(
                                                                    self._key, self._section_name)

   # ------------------------------------------------------------------------------
   class SectionNotFound(Exception):
      """!
         Exception thrown when a section could not be found in the list of registered sections 
      """
      def __init__(self, section_name):
         """!
             @param section_name name of the section that could not be found
         """
         self._section_name = section_name

      def __str__(self):
         """! String representation of the exception object """
         return "SConfig *** ERROR ***: section [{0}] no found".format(self._section_name)

   # ------------------------------------------------------------------------------
   class SectionInvalidFormat(Exception):
      """!
         Exception thrown when a section has an invalid format
      """
      def __init__(self, section_name):
         """!
             @param section_name name of the section with an invalid format
         """
         self._section_name = section_name

      def __str__(self):
         """! String representation of the exception object """
         return "SConfig *** ERROR ***: invalid format for section [{0}]".format(self._section_name)

   # ------------------------------------------------------------------------------
   class InvalidSection(Exception):
      """!
         Exception thrown in case a section is invalid
      """
      def __init__(self, expr):
         """!
             @param expr expression for the error
         """
         self._expr = expr

      def __str__(self):
         """! String representation of the exception object """
         return "SConfig *** ERROR ***: invalid section [{0}]".format(self._expr)

   # ------------------------------------------------------------------------------
   class InputOutputError(Exception):
      """!
         Exception thrown in case of input/output error
      """
      def __init__(self, error):
         """!
             @param error text of the error
         """
         self._error = error

      def __str__(self):
         """! String representation of the exception object """
         return "SConfig *** ERROR ***: IO error {0}".format(self._error)

   # ------------------------------------------------------------------------------
   class Section(object):
      
      """! 
         Inner class that represents a section in the configuration file

      """

      def __init__(self, name="ROOT", parent=None, rank=0):
         """!
            Construct a Section object
            @param name name of the section to create, @c ROOT_SECTION_NAME if the root section
            @param parent name of the parent section, @c ROOT_SECTION_NAME if the root section    
            @param rank order of creation of the section
            @note this constructor is internal and not supposed to be invoked from user programs
         """
         
         self._name = name             # section name
         self._parent = parent         # parent section object
         self._config_dico = {}        # holds the key-value pairs for the section
         self._subsections = []        # defined subsections
         self._mothers = []            # mother sections this section is inheriting from
         self._rank = rank             # section order of creation
         self._next_key_rank = 0       # ranking counter
         self._rank_key_tuples = []    # to track key natural order         
         self._referenced_keys = []    # keys referenced with the '&' REF operator

         self._prefixes = [SConfig.Section.KeyPrefix.REF,\
                           SConfig.Section.KeyPrefix.DEREF]  # special key prefixes

      # ~~~~~~~~~~~
      # Properties 
      # ~~~~~~~~~~~

      @property
      def name(self):
         """! @return the name of the section """
         return self._name

      @property
      def config_dico(self):
         """! @return the section internal dictionary """
         return self._config_dico 

      @property
      def parent(self):
         """! @return the parent Section object """
         return self._parent 

      @property
      def subsections(self):
         """! @return the list of subsection Section objects """
         return self._subsections 

      @property
      def mothers(self):
         """! @return the list of parent Section objects of this Section """
         return self._mothers 

      @property
      def rank(self):
         """! @return the rank of this Section """
         return self._rank
       
      @rank.setter 
      def rank(self, rank):
         """! 
            Set the rank of this Section
            @param rank the rank of this Section 
         """
         self._rank = rank     

      @property
      def rank_key_tuples(self):
         """! @return the sorted list used to track key ordering """
         return sorted(self._rank_key_tuples)

      @property
      def keys(self):
         """! @return the list of keys of the Section sorted by ranks """
         return [key for (rank,key) in self.rank_key_tuples] 
   
      @property  
      def values(self):
         """! @return the list of the keys of the Section sorted by ranks """
         return [self._config_dico[key] for (rank, key) in self.rank_key_tuples]

      # ~~~~~~~~~~~~~~~ 
      # Public methods 
      # ~~~~~~~~~~~~~~~ 

      def is_subsection_of(self, parent):
         """! 
            Tells if this section belongs to a parent section
            @param parent parent Section object
            @retval True if this section belongs to @c parent 
            @retval False otherwise
         """
         return self in parent.subsections

      def has_subsection(self, child):
         """! 
            Tells if this section has @c child as a subsection
            @param child child Section object
            @retval True if this section has @c child as a subsection
            @retval False otherwise
         """
         return child.parent == self

      def has_key(self, key):
         """! 
            Tells whether a specified key belongs to this section
            @param key key to check
            @retval True if this section contains @c key
            @retval False otherwise
         """
         return self._strip_key(key) in self._config_dico

      def is_root(self):
         """!
            Tells if this section is the @c ROOT_SECTION_NAME Section 
            @retval True if the top-level @c ROOT_SECTION_NAME Section
            @retval False otherwise
         """   
         return self._name == "ROOT" 

      def inherits_from(self, mother):
         """! 
            Tells if this section is an instance of a mother Section 
            @param mother mother Section 
            @retval True if this section is a child of @c mother
            @retval False otherwise
         """
#         print("mothers of {0}: {1}".format(self.name, mother.name))
         if mother.is_root():
            return True
         else:
            return mother in self.mothers

      def add_key_value(self, key, value):
         """! 
            Add a (key,value) pair 
            @param key key to add
            @param value value to associate with @c key
         """
         self._config_dico[key] = value
         self._rank_key_tuples.append((self._next_key_rank, key))
         self._next_key_rank += 1

      def write(self, output_file):
         """!
            @todo not implemented in this section 
         """
         pass

      # ~~~~~~~~~~~~~~~~ 
      # Private methods 
      # ~~~~~~~~~~~~~~~~ 

      def _copy_data_from(self, section):
         # Update (key,value) pairs
         d = dict(self._config_dico.items())
         d.update(section.config_dico)
         self._config_dico = d

         # Update ranking info
         mother_rank_key_tuples = section.rank_key_tuples
         self_ranK_key_tuples = zip(np.arange(len(mother_rank_key_tuples), len(self._config_dico)), self.keys)
         #print("=> self_ranK_key_tuples:", self.name, self.rank_key_tuples, self_ranK_key_tuples)
         self._rank_key_tuples = mother_rank_key_tuples
         self._rank_key_tuples.extend(self_ranK_key_tuples) 
         self._next_key_rank = len(self._rank_key_tuples)
         #print("*** COPY: rank_key_tuples of {0}: {1} {2} <-> {3}".format(self.name, self._rank_key_tuples, len(self._config_dico), len(self._rank_key_tuples)))

      def _strip_key(self, key):
         """ Return a key without prefix if any """
         if key[0] in self._prefixes:
            return key[1:]
         else:
            return key 

      # ~~~~~~~~~~ 
      # Operators
      # ~~~~~~~~~~ 

      def __str__(self):
         """! Format of a Section object for printing """
         output_str = ""   
         if self.is_root():
            output_str += ("[{0}] #keys: {1}".format(self.name, len(self._config_dico))) 
         else:   
            output_str += ("rank: {0} [{1}] #keys: {2} subsection of: [{3}] mother: {4}".format(self.rank, self.name, len(self._config_dico), self.parent.name, [s.name for s in self._mothers]))
         
         for (s,t) in self._config_dico.items():
            output_str += ("\n\t{0}={1}".format(s,t))      

   #      [output_str += ("\t{0}={1}".format(s,t)) for (s,t) in self._config_dico.items()]      

         return output_str   

      def __eq__(self, other):
         """! Equality comparator between Section objects """
         return self._name == other.name.strip() and self._parent.name == other.parent.name.strip()      

      def __ne__(self, other):
         """! Non-Equality comparator between Section objects """
         return not self.__eq__(other)

      # ~~~~~~~~~~~~~ 
      # Inner Classes
      # ~~~~~~~~~~~~~ 

      # -----------------------------------------------------------------------------------------------
      class KeyPrefix:
         """!
            @todo For a future implementation
         """
         (REF, DEREF) = ('&','*')


def main():
   
   config = SConfig("test.cfg")


if __name__ == "__main__":
    main()
