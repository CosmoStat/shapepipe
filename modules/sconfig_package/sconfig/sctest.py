#!/usr/bin/env python

"""! 
   @package sconfig.sctest Simple section-based configuration file management
   @author Marc Gentile
   @file sctest.py
   Testing program for SConfig
"""

import sys
from sconfig import SConfig

def test():

   # --- Reading configuration file
   config = SConfig("./test/config.cfg")

   # --- Getting & setting values
   
   print (config.get_as_string("GLOBAL_1"))    

   # --- Testing continuation character
   print('Continuation line 1: {0}'.format(config.get_as_string("LONG_LINE")) )    

   # --- Testing module loading 
   print("Loading module numpy.fft...")
   print("numpy.linalg module: {0}".format(config.get_as_module("FFT_MODULE")) )

   # --- Testing Section management
   print(config.get_as_string("KEY_HOME", "section 1")) 
   print(config.get_as_dict("DICT_1", "section 1")) 
   print(config.get_as_dict("DICT_2", "section 1")) 
      
   print("Section {0} has key {1}: {2}".format("KEY_HOME", "section 1", 
                                               config.has_key("KEY_HOME", "section 1"))) 

   print(config.get_as_string("KEY_1_SECTION_3_SUBSECTION_3", "section 2. subsection 3"))
   print(config.get_as_string("KEY_1_SECTION_3_SUBSECTION_3_SUBSECTION_4",
                              "section 2.subsection 3.subsection 4")) 
   l = config.get_as_list("LIST_1", "section 1")
   print(l, type(l))
   a = config.get_as_array("ARRAY_1", "section 1", False)
   print(a, type(a))
   print(config.get_as_boolean("BOOL_1", "section 1"))

   print(config.get_section_data(" section 2.   subsection 3 "))
   print(config.get_section_data())
   print("Values of [{0}]: {1}".format(" section 2.   subsection 3 ", 
                                       config.get_section_values(" section 2.   subsection 3 ")))
   print("Keys of [{0}]: {1}".format("section 4", config.get_section_keys("section 4")))
   print("Values of [{0}]: {1}".format("section 4", config.get_section_values("section 4")))

   config.set_key_value("KEY_1", 1000, "section 1")
   print(config.get_key_value("KEY_1", "section 1"))   
   print(config.get_key_value("KEY_1_SECTION_1", "section 1"))   

   # --- Testing section inheritance

   l = config.set_key_value("LIST_1", [5,6,7,8], "section 1")
   print(l, type(l))

   print("Data for section: {0}: {1} ".format("section 1", config.get_section_data("section 1")))
   print("Data for section: {0}: {1} ".format("section 4", config.get_section_data("section 4")))
   print("Data for section: {0}: {1} ".format("section 5", config.get_section_data("section 4"))) 

   print("Section [{0}] inherits from [{1}]: {2}".format("section 2", "ROOT", 
                               config.section_inherits_from("section 2", "ROOT"))) 
   print("Section [{0}] inherits from [{1}]: {2}".format("section 4", "section 1", 
                               config.section_inherits_from("section 4", "section 1"))) 
   print("Section [{0}] inherits from [{1}]: {2}".format("section 2", "section 2.subsection 3", 
                               config.section_inherits_from("section 2", "section 2.subsection 3"))) 

   # --- Opening/closing configuration file
   fd = config.open("./test/config.cfg")
   if fd is not None:
      fd.close()
   else:
      print("not opened")

   # --- Sub-Section hierarchy

   print("TOP subsection names of [{0}]: {1}".format("ROOT", config.get_subsections("ROOT", recurse=True)))

   # Finding children
   print("TOP Children of [{0}]: {1}".format("ROOT", 
                         [child.name for child in config._get_subsections("ROOT", recurse=False)]))      
   for child in config._get_subsections("ROOT"):
      print("Children of [{0}]: {1}".format(child.name, 
                         [child.name for child in config._get_subsections(child.name)]))      

   print("TOP Children of [{0}]: {1}".format("section 1", 
                         [child.name for child in config._get_subsections("section 1")] ))      
   print("TOP Children of [{0}]: {1}".format("section 4", 
                         [child.name for child in config._get_subsections("section 4")] ))      

   print("ALL Children of [{0}]: {1}".format("ROOT", [child.name for child in config._get_subsections("ROOT", True)] ))        

   # Child relationship
   print("Section [{0}] is subsection of [{1}]: {2}".format("section 2.subsection 3", "section 2", 
                                  config.is_subsection("section 2.  subsection 3", "section 2")) )
   print("Section [{0}] is subsection of [{1}]: {2}".format("section 2.subsection 3", "section 1", 
                                  config.is_subsection("section 2.  subsection 3", "section 1")) )
   print("Section [{0}] is subsection of [{1}]: {2}".format("section 2.subsection 3.subsection 4",\
                                                            "section 2.subsection 3", 
            config.is_subsection("section 2.subsection 3.subsection 4", "section 2.subsection 3")) )

   # Parent relationship
   print("Section [{0}] has subsection [{1}]: {2}".format("section 2", "section 2.subsection 3", 
                                    config.has_subsection("section 2", "section 2.subsection 3")) )   
   print("Section [{0}] has subsection [{1}]: {2}".format("section 1", "section 2.subsection 3", 
                                    config.has_subsection("section 1", "section 2.subsection 3")) )

   # --- Dump configularation file
   print("\nDump of [{0}]: {1}".format("ROOT", 
                                       config.dump_to_string(section="ROOT", recurse=True)))
   print("\nDump of [{0}]: {1}".format("section 4", 
                                       config.dump_to_string(section="section 4", recurse=False)))

   print("Full recursive dump:")
   config.dump_to_file(recurse=True)

if __name__ == "__main__":
   test()
