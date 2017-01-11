#!/usr/bin/env python

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
# test program for scdm.py 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	

import numpy.random
import time
import pylab
import matplotlib
import matplotlib.figure
from matplotlib.font_manager import FontProperties
import mpl_toolkits.mplot3d.axes3d as pylab3

from scdm import *


#
# make_stamp
#
def make_stamp(stamp_size, func, args):
   lin_space = numpy.linspace(0, stamp_size -1, stamp_size)
   grid = numpy.meshgrid(lin_space, lin_space)
   return func(grid, *args)



def make_stamp_func1(grid, grid_size, a, b, c):
   
   #stamp_size = data.shape[0]
   xc, yc = grid_size/2, grid_size/2
   x, y = grid
   x -= xc   
   y -= yc   

   #stamp = numpy.random.normal(0, d, (grid_size, grid_size))

   
   #stamp =  numpy.exp(-(x**2 + y**2) / (2 * d**2)) / (d * math.sqrt(2*math.pi))

   r = numpy.hypot(x, y)   
   stamp = a * numpy.exp(-x**2) + b * y**2 - c
   #stamp = (a * x**2 + b * y**2) * c
   #stamp = (a * x**2 + b * y**2) * c
#   stamp =  c * numpy.exp(-r**((a)*2)) + b 
   #stamp = a * y**3 + b * x**2 + c * x * y
   #stamp = -(a * x**2 + b *y**2 + c)**0.5 
   #stamp = x**2 / a - b *y**2 + c 
   #stamp = (x**2 / a + y**2/ b) + c
 
   #stamp =  c * numpy.exp( -(2.0 * a) * (r/b)**(1/a))
   #stamp =  c  - (2.0 * a - 0.331) * (r/b)**(1/a)
   #stamp =  numpy.exp(-(r/b)**(1/c) + a) 
   #r = numpy.hypot(a*x, b*y)   
   #stamp = numpy.exp(-r**2) * c

#   back_mean  = 0.0
#   back_sigma = 1e-10
#   stamp += numpy.random.normal(back_mean, back_sigma, stamp.shape)

   return stamp   


def test_func_1(param_values, ref_stamp, stamp_size):   

   #print 'Evaluating for parameter values', param_values

   #[guess_a, guess_b, guess_c] = param_values
   args = [stamp_size] + param_values
   guess_stamp = make_stamp(stamp_size, make_stamp_func1, args)

   #residuals = numpy.absolute(stamp - guess_stamp)
   #residuals = (ref_stamp - guess_stamp)**2
   residuals = numpy.absolute((ref_stamp - guess_stamp))
   diff = numpy.sum(residuals)
   #print numpy.sum(stamp), numpy.sum(guess_stamp), diff

   #plot_contour_stamp(guess_stamp, plot_title='Guess Contour Stamp: ' + str(diff), has_grid=False, output_dir='../plots', output_file='guess_stamp_contour_plot_'+ str(diff)+'.png')
   #plot_contour_stamp(residuals, plot_title='Residuals Contour Stamp: ' + str(diff), has_grid=False, output_dir='../plots', output_file='residuals_stamp_contour_plot_'+ str(diff)+'.png')

   param_values = [('%9.3f' %p) for p in param_values]

   #plot_stamp_3D(residuals, plot_title='Residuals 3D', has_grid=False, output_dir='../plots', output_file=('residuals_plot_3D_%s_%9.3e.png' %(param_values, diff) ), show=False)
   #plot_stamp_3D(guess_stamp, plot_title='Estimated Stamp 3D', has_grid=False, output_dir='../plots', output_file=('fitted_stamp_plot_3D_%s_%9.3e.png' %(param_values, diff) ), show=False)

   return [diff, residuals]


def test_out_of_bound_handler(iiter, param, paramValue, boundIndex):
      
   print '*** test_out_of_bound_handler() called for:', param.name, paramValue, boundIndex
   status = 0
   if param.name == 'param1':
      status = -1
   return [status, paramValue, param.param_range[boundIndex]]                              

def test_numerical_error_handler(iiter, param, value, valueName):
      
   print '*** test_numerical_error_handler() called for:', param.name, value, valueName
   status = -1
   if param.name == 'param1':
      status = -1
   return status                         


#array1 = numpy.arange(0,100)
#tuple1 = [-1, 1]
#args = [array1, tuple1]
#test_func_1(params, *args)

# Make Reference stamp
stamp_size = 99
back_mean  = 0.0
back_sigma = 1e-10
stamp_args = [stamp_size] + [1, 2,  3]
ref_stamp = make_stamp(stamp_size, make_stamp_func1, stamp_args)
ref_stamp += numpy.random.normal(back_mean, back_sigma, ref_stamp.shape)

args = [ref_stamp, stamp_size]
param1 = Param(name='param1', param_range=[None, None], guess_value= 0,  step_size=1e-1, step_range=[None, 1e-1], step_factor=2, max_tol=1e-8)
param2 = Param(name='param2', param_range=[None, None], guess_value= 0,  step_size=1e-1, step_range=[None, 1e-1], step_factor=2, max_tol=1e-8)
param3 = Param(name='param3', param_range=[None, None], guess_value= 0,  step_size=1e-1, step_range=[None, 1e-1], step_factor=2, max_tol=1e-8)

#param4 = Param(name='param4', param_range=[1, None],     guess_value= 500,  step_size=100,  step_range=[1e-8, None], step_factor=2, maxFuncTol=1e-8, scaling=ParamScaling.Scale)

param1.out_of_bound_handler = test_out_of_bound_handler
param1.numerical_error_handler = test_numerical_error_handler

params = Params([param1, param2, param3])
#target_params = Params([param1, param2])
target_params = params


#print 'param names:', params.getParamNames()   
#print 'param ranks:', params.getParamRanks() 
#print params.getParamRank('param3')

#print params.getParamByName('param1').name
#print params.getParamByName('param2').name
#print params.getParamByName('param3').name

#print 'param at rank 1:', params.getParamByRank(1).name 
#print 'param at rank 2:', params.getParamByRank(2).name 
#print 'param at rank 3:', params.getParamByRank(3).name 
#print 'param tuples: ', params.getParamTuples()

#minimization_type = MinimizationType.minimize_param_sigma
minimization_type = MinimizationType.minimize_func_value | MinimizationType.minimize_param_delta | MinimizationType.minimize_param_sigma
fitting_options = FittingOptions(max_iter=20000, max_fev=10000, max_func_tol=1e-10, all_minima=False, mini_type=minimization_type)

base_output_directory = '.'
diag_options = DiagOptions(diag_flags=0, base_output_directory=base_output_directory)
#diag_options = DiagOptions(diag_flags=DiagFlag.logs | DiagFlag.plots | DiagFlag.stats, base_output_directory=base_output_directory)
#diag_options = DiagOptions(diag_flags=DiagFlag.logs | DiagFlag.plots, base_output_directory=base_output_directory) 
#diag_options = DiagOptions(diag_flags=DiagFlag.logs, base_output_directory=base_output_directory)

minimizer = Minimizer(params=params, target_params=target_params, func=test_func_1, args=args, fitting_options=fitting_options, diag_options=diag_options)

plot_helper = minimizer.plot_helper

start_time = time.clock()

[fitted_param_values, fitting_results] = minimizer.minimize()

elapsed_time = time.clock() - start_time

# Fitted stamp
stamp_args = [stamp_size] + fitted_param_values
fitted_stamp = make_stamp(stamp_size, make_stamp_func1, stamp_args)

print("Fitted Parameters: {0} - Time: {1}".format(fitted_param_values, elapsed_time))

#if minimizer.isPlottingEnabled():
#   #plot_stamp(stamp, plot_title='Stamp', has_grid=False, output_dir='.', output_file='stamp_plot')
#   #plot_contour_stamp(ref_stamp, plot_title='Ref Contour Stamp', has_grid=False, output_dir='../plots', output_file='ref_stamp_contour_plot')
#   plot_helper.plotContourStamp(ref_stamp, plot_title='Reference Contour Stamp', has_grid=False, output_dir='../plots', output_file='ref_stamp_contour_plot')
#   plot_helper.plotContourStamp(fitted_stamp, plot_title='Fitted Contour Stamp', has_grid=False, output_dir='../plots', output_file='fitted_stamp_contour_plot')
#   plot_helper.plotStamp3D(ref_stamp, plot_title='Reference Surface to fit', x_label='parameter 1', y_label='parameter 2', z_label='parameter 3', 
#                              has_grid=False, output_dir='../plots', output_file='ref_surface_3D', show=False)
#   plot_helper.plotStamp3D(fitted_stamp, plot_title='Fitted Surface to fit', x_label='parameter 1', y_label='parameter 2', z_label='parameter 3', 
#                          has_grid=False, output_dir='../plots', output_file='fitted_surface_3D', show=False)

#   plot_helper.plotContourStamp(fitting_results.getCustomData()[0],  plot_title='First residuals', has_grid=False, output_dir='../plots', output_file='first_residuals_contour_plot')
#   plot_helper.plotContourStamp(fitting_results.getCustomData()[-1], plot_title='Best-Fit residuals', has_grid=False, output_dir='../plots', output_file='best-fit_residuals_contour_plot')
#   plot_helper.plotStamp3D(fitting_results.getCustomData()[0],  plot_title='First residuals',    x_label='parameter 1', y_label='parameter 2', z_label='parameter 3',
#                          has_grid=False, output_dir='../plots', output_file='first_residuals_3D', show=False)
#   plot_helper.plotStamp3D(fitting_results.getCustomData()[-1], plot_title='Best-fit residuals', x_label='parameter 1', y_label='parameter 2', z_label='parameter 3',
#                          has_grid=False, output_dir='../plots', output_file='best-fit_residuals_3D', show=False)


#   # Parameter Surface plots
#   fittingData = fitting_results.getFittingData()
#   fittingDataDico = fittingData.getFittingDataDico()   
#   #print fittingData.getParamFittingDataFor('param1').getFittingDataDico()

#   prettyParamValues = [string.strip('%9.3f' %v) for v in fitted_param_values]

#   plot_helper.plotParameterSet3D(['param1', 'param2', 'param3'], 'feval',
#                                  'feval_3D', '../plots', plot_title='Function evaluation evolution Surface', extra_title=prettyParamValues, pfilter=None, first_index=-1, last_index=-1,
#                                  x_label='parameter 1', y_label='parameter 2', z_label='parameter 3', has_grid=True, pfigure=None, is_normed=False, has_equalaxes=False, show=False)

#   plot_helper.plotParameterSet3D(['param1', 'param2', 'param3'], 'steps',
#                                  'steps_3D', '../plots', 'Step evolution Surface for %s' %prettyParamValues, pfilter=None, first_index=-1, last_index=-1,
#                                  x_label='parameter 1', y_label='parameter 2', z_label='parameter 3', has_grid=True, pfigure=None, is_normed=False, has_equalaxes=False, show=False)

#   plot_helper.plotParameterSet3D(['param1', 'param2', 'param3'], 'deval',
#                                  'deval_3D', '../plots', 'Derivative evolution Surface for %s' %prettyParamValues, pfilter=None, first_index=-1, last_index=-1,
#                                  x_label='parameter 1', y_label='parameter 2', z_label='parameter 3', has_grid=True, pfigure=None, is_normed=False, has_equalaxes=False, show=False)

#   plot_helper.plotParameterSet3D(['param1', 'param2', 'param3'], 'params',
#                                  'peval_3D', '../plots', 'Parameter evolution Surface for %s' %prettyParamValues, pfilter=None, first_index=-1, last_index=-1,
#                                  x_label='parameter 1', y_label='parameter 2', z_label='parameter 3', has_grid=True, pfigure=None, is_normed=False, has_equalaxes=False, show=False)

#   plot_helper.plotParameterSet3D(['param1', 'param2', 'param3'], 'dcount',
#                                  'dcount_3D', '../plots', 'Direction count evolution Surface for %s' %prettyParamValues, pfilter=None, first_index=-1, last_index=-1,
#                                  x_label='parameter 1', y_label='parameter 2', z_label='parameter 3', has_grid=True, pfigure=None, is_normed=False, has_equalaxes=False, show=False)

