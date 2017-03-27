#!/usr/bin/env python

"""! 
   multifit - Unit Testing program
""" 

import numpy
import numpy.random

from multifit import *



# --------------------------------------------------------------------------------------------------
def profile_func(x_array, y_array, *params):

   [p1, p2, p3] = params
   r = numpy.hypot(x_array, y_array)   

   profile = func_exp(r, p1, p2, p3)
   
   #profile = p1**2 * numpy.exp(-r**((p2)*2)) + p3 
   #profile = p1**3 + p2**2 + p3 

   return profile      

# --------------------------------------------------------------------------------------------------
def create_profile(grid, **kwargs):

   [grid_size] = kwargs['args']
   params = kwargs['params']

   xc, yc = grid_size/2, grid_size/2
   x, y = grid
   x -= xc   
   y -= yc 

   stamp = profile_func(x, y, *params)

   return stamp   
   
# # --------------------------------------------------------------------------------------------------
# def objective_func_1(grid, **kwargs):
# 
#    guess_stamp = make_stamp(create_profile, **kwargs)      
# 
#    residuals = ref_stamp**2 - guess_stamp**2
#    diff = numpy.sum(residuals)
# 
#    #print 'In objective_func_1:', params
# 
#    return diff

# --------------------------------------------------------------------------------------------------
def make_stamp(func, **kwargs):
   """ Create a postage stamp with a profile given by func(), a set of parameters and args """

   #print 'kwargs:', kwargs['args']

   [stamp_size] = kwargs['args']
   [p1, p2, p3] = kwargs['params'] 

   lin_space = numpy.linspace(0, stamp_size -1, stamp_size)
   grid = numpy.meshgrid(lin_space, lin_space)

   return func(grid, **kwargs)

# --------------------------------------------------------------------------------------------------
def func_exp(r, p1, p2, p3):
#   return p1*r**2 + p2*r + p3
   return p1*r * numpy.exp(-r**(p2)) + p3

# --------------------------------------------------------------------------------------------------
def make_profile(x, y, func, params):
      
   r = numpy.hypot(x, y)   
   return func(r, *params)

# --------------------------------------------------------------------------------------------------
def create_profile2d(model_func, true_params, grid_size):

   lin_space = numpy.linspace(0, grid_size -1, grid_size)
   x, y = numpy.meshgrid(lin_space, lin_space)
   xc, yc = grid_size/2, grid_size/2
   x -= xc   
   y -= yc 
   true_profile = make_profile(x, y, model_func, true_params)   

   # --- add some noise
   noise = numpy.random.random((grid_size, grid_size)) / 10.0
   true_profile += noise    

   return true_profile

# --------------------------------------------------------------------------------------------------
def profile2d_objective_func(params, model_func, true_profile, profile_func, grid_size):
   """ Example of objective function for a profile2d model """
   
   return (true_profile.flatten() - model_func(params, profile_func, grid_size).flatten())**2

# --------------------------------------------------------------------------------------------------
def profile2d_objective_func_2(params, model_func, true_profile, profile_func, grid_size):
   """ Example of objective function for a profile2d model """
   
   return numpy.sum((true_profile.flatten() - model_func(params, profile_func, grid_size).flatten())**2)

# --------------------------------------------------------------------------------------------------
def profile2d_objective_func_3(params, data):
   """ Example of objective function for a profile2d model """
   
   model_func, true_profile, profile_func, grid_size = data

   #print "params:", params

   return (true_profile.flatten() - model_func(params, profile_func, grid_size).flatten())

# --------------------------------------------------------------------------------------------------
def profile2d_fitting_func(params, model, profile_func, grid_size):

   return model.func(params, profile_func, grid_size)


# --------------------------------------------------------------------------------------------------
def test():

   # ---- Model1

   # Make Reference stamp
   stamp_size = 31
   back_mean  = 0.0
   back_sigma = 1e-10
   params = [1, 1 , 1]
   #ref_stamp = make_stamp(create_profile, args=[stamp_size], params=params)
   #print ref_stamp

   # --- Two-dimensional least square fitting, arbitrary model
   grid_size = 31
   true_params = [10,2,5]
   true_profile = create_profile2d(func_exp, true_params, grid_size)
   guess_params = [ModelParam('p'+str(n), 1.0) for n in [1,2,3]]

#   # --- Scipy least square fitting

#   fitter = ModelFitter(model_name="profile2d", method_name="scileastsq")                       

#   fitter.model.params = guess_params   
#   fitter.method.objective_func = profile2d_objective_func


#   args = [fitter.model.func, true_profile, func_exp, grid_size]
#   status = fitter.fit(args)

#   print "*** scipy leastsq ***"
#   print "status:", status
#   print "Fitting info:", fitter.fitting_info
#   print "Fitting options:", fitter.method.fitting_options
#   print "Diag options:", fitter.method.diag_options

#   # --- Scipy least square fitting - Version 2

#   fitter = ModelFitter(model_name="profile2d", method_name="scileastsq")                       

#   fitter.model.params = guess_params   
#   fitter.method.objective_func = fitter.method.get_chi2_objective_func()

#   args = [profile2d_fitting_func, fitter.model, true_profile, 1.0, func_exp, grid_size]
#   #args = [fitter.model.func, true_profile, func_exp, grid_size]

#   #args = [fitter.model.func, true_profile, func_exp, grid_size]
#   status = fitter.fit(args)

#   print "*** scipy leastsq ***"
#   print "status:", status
#   print "Fitting info:", fitter.fitting_info
#   print "Fitting options:", fitter.method.fitting_options
#   print "Diag options:", fitter.method.diag_options

   # --- KMPFIT fitting

   # Model: profile2d

   fitter = ModelFitter(model_name="profile2d", method_name="kmpfit", 
                      model_config_dir="./config/models", 
                      method_config_dir="./config/methods", 
                      fitting_config_dir = "./config/methods")             

   #fitter.model.params = guess_params   
   fitter.method.objective_func = profile2d_objective_func_3

   args = [fitter.model.func, true_profile, func_exp, grid_size]
   status = fitter.fit(args)

   print "*** KMPFIT Profile2D ***"
   print "status:", status
   print "Fitting info:", fitter.fitting_info
   print "Fitting options:", fitter.method.fitting_options
   print "Diag options:", fitter.method.diag_options


   # --- KMPFIT fitting (default chi2 objective functiuon)

   # Model: profile2d

   fitter = ModelFitter(model_name="profile2d", method_name="kmpfit", 
                      model_config_dir="./config/models", 
                      method_config_dir="./config/methods", 
                      fitting_config_dir = "./config/methods")             

   fitter.model.params = guess_params   
   #fitter.method.objective_func = fitter.method.get_chi2_objective_func()

   sigmas = numpy.ones(grid_size)

   #args = [fitter.model.func, true_profile, func_exp, grid_size]
   args = [profile2d_fitting_func, fitter.model, true_profile, sigmas, [func_exp, grid_size]]

   status = fitter.fit(args)

   print "*** KMPFIT Profile2D ***"
   print "status:", status
   print "Fitting info:", fitter.fitting_info
#   print "Fitting options:", fitter.method.fitting_options
#   print "Diag options:", fitter.method.diag_options

##   # --- SCDM fitting

#   # Model: profile2d

#   fitter = ModelFitter(model_name="profile2d", method_name="scdmin", 
#                      model_config_dir="./config/models", 
#                      method_config_dir="./config/methods", 
#                      fitting_config_dir = "./config/methods")             

#   fitter.model.params = guess_params   
#   fitter.method.objective_func = profile2d_objective_func_2

#   args = [fitter.model.func, true_profile, func_exp, grid_size]
#   status = fitter.fit(args)

#   print "*** SCDM Profile2D ***"
#   print "status:", status
#   print "Fitting info:", fitter.fitting_info
#   print "Fitting options:", fitter.method.fitting_options
#   print "Diag options:", fitter.method.diag_options

#   # Model: profile2d - version 2

#   fitter = ModelFitter(model_name="profile2d", method_name="scdmin",
#                        model_config_dir="./config/models", 
#                        method_config_dir="./config/methods", 
#                        fitting_config_dir = "./config/methods")                

#   fitter.model.params = guess_params   
#   fitter.method.objective_func = fitter.method.get_chi2_objective_func()

#   sigmas = numpy.ones(grid_size)
#   args = [profile2d_fitting_func, fitter.model, true_profile, sigmas, [func_exp, grid_size]]
#   status = fitter.fit(args)

#   print "*** SCDM Profile2D bis ***"
#   print "status:", status
#   print "Fitting info:", fitter.fitting_info
#   print "Fitting options:", fitter.method.fitting_options
#   print "Diag options:", fitter.method.diag_options


if __name__ == "__main__":
   test()
