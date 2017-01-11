"""! 
   sp_spatial.py - Computation of distances - Search of Nearest-neighbors over a 2D space
"""

# -- Python imports
import os, sys
import string
import time
import imp
import numpy
import scipy
from scipy.spatial import KDTree, cKDTree

# --- Module-specific imports

# --- External imports


# -------------------------------------------------------------------------------------------------
class SpatialCalculator(object):
   
   """! 
      Computation of distances - Search of Nearest-neighbors over a 2D space
   """

   def __init__(self, coords):
      """! 
         SpatialCalculator constructor 
         @param coords coordinate tuple (x_coords, y_coords)
         @param master master process
      """

      # --- Public attributes
      self._coords = coords    # coordinate area to make calculations from

      # --- Private attributes
      self.__kdtree = KDTree(coords)  # private KDTree object

   # ~~~~~~~~~~
   # Properties 
   # ~~~~~~~~~~

   @property
   def coords(self):
      """! @return the underlying coordinates """
      return self._coords



   # ~~~~~~~~~~~~~~~
   # Public methods 
   # ~~~~~~~~~~~~~~~


   # -----------------------------------------------------------------------------------------------
   def locate_neighbors_within_range(self, location, min_nb_neighbors, max_nb_neighbors, 
                                           max_dist, dist_tol, dist_incr, dist_type=2):

      
      neighbors = []
      neighbor_indice = []

      # --- Find the minimum and maximum neighbor counts within <max_distance>
      guess_dists = numpy.arange(0.0, max_dist + dist_tol, dist_incr)

      #print "min_nb_neighbors", min_nb_neighbors, "max_nb_neighbors", max_nb_neighbors

      neighbor_counts = self.__kdtree.count_neighbors(KDTree([numpy.array(location)]), guess_dists)

      #print "neighbor_counts:", neighbor_counts

      matching_indice = numpy.where(neighbor_counts >= min_nb_neighbors)[0]

      #print "matching_indice:", matching_indice

      if len(matching_indice) > 0:
         matching_dists = guess_dists[matching_indice]

         #print "matching_dists:", matching_dists

         neighbor_indice = self.__kdtree.query_ball_point(location, matching_dists.min(), 
                                                          p=dist_type, eps=0)         

         #print "neighbor_indice:", neighbor_indice

         x_obj, y_obj = numpy.hsplit(self.__kdtree.data[neighbor_indice], 2)
         x_obj, y_obj = x_obj.squeeze(), y_obj.squeeze()

         return [zip(x_obj, y_obj), neighbor_indice]      
      else:

         # --- No neighbors found within <max_dist>
         return [[], []]         

      

   # -----------------------------------------------------------------------------------------------
   def build_distance_matrix(self, object_list):
      """! 
         Construct a distance matrix 
         @param object_list array of coordinates [x,y]
      """   

      return scipy.spatial.distance.squareform(scipy.spatial.distance.pdist(
                                                                      numpy.asarray(object_list)))


   # -----------------------------------------------------------------------------------------------
   def build_distance_vector_from_location(self, location_tuple, object_list):
      """! 
         Construct an array containing distances from @c location to objects whose coordinates
         are found in @c object_list    
         @param coordinates [x0,y0] to compute distances from
         @param object_list array of coordinates [x,y]
      """   

      if not location_tuple in object_list:
         object_list.insert(0, list(location_tuple))
         dist_vector = self.build_distance_matrix(object_list)[1:,0]
      else:
         dist_vector = self.build_distance_matrix(object_list)[object_list.index(location_tuple)]

      return dist_vector


# -- EOF sp_spatial.py
