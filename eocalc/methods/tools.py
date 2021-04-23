# -*- coding: utf-8 -*-
"""Basic toolbox for space emission calculator used in emission calcs."""
import numpy as np
import pandas as pd
import scipy.special
import binas
# formulas gaussian plumes
# based on (Fioletov et al., 2017; Dammers et al., 2021)
# https://acp.copernicus.org/articles/17/12597/2017/
# any location can be described as VCD(x,y) = sum(flow_contribution_i(x,y) * emis_i) + bias
# for this we need VCD's, added windfields to calculate flow, and potential source locations (grid)
# x,y are in km x km space, not degrees, so rotation and adjustment needed compared to lat lon.

# winddir/speed back to u,v 
def winddir_speed_to_u_v(wind_spd, wind_dir):
    dir_rp1 = wind_dir * 2 * np.pi / 360.
    dir_rp = np.array(dir_rp1) + np.pi / 2
    u = wind_spd * np.cos(dir_rp)
    v = -wind_spd * np.sin(dir_rp)
    return u, v
        
# wind speed
def calc_wind_speed(u,v):
    return np.sqrt(np.array(u)**2.+np.array(v)**2.)

# winddirection
def calc_wind_direction(u, v):
    # calc wind direction from u,v
    # ensure array form
    u = np.array(u)
    v = np.array(v)
    wind_dir = np.zeros(len(u),float)
    
    sel_lzero = (v >= 0)
    sel_neg_both = ((u < 0) & ( v < 0))
    sel_negv_posu = ((u >= 0) &( v < 0))
    
    if len(wind_dir[sel_lzero])>0:
        wind_dir[sel_lzero] = ((180. / np.pi) * np.arctan(u[sel_lzero] / v[sel_lzero]) + 180.)
    if len(wind_dir[sel_neg_both])>0:
        wind_dir[sel_neg_both] = ((180. / np.pi) * np.arctan(u[sel_neg_both] / v[sel_neg_both]) + 0.)
    if len(wind_dir[sel_negv_posu])>0:
        wind_dir[sel_negv_posu] = ((180. / np.pi) * np.arctan(u[sel_negv_posu] / v[sel_negv_posu]) + 360.)

    return wind_dir
    
# rotation
def rotate_plume_around_point(reflon, reflat, lon, lat, wind_direction):
    # rotate lat/lon grid plume in relation to a point
    dtr = np.pi / 180. # conversion degrees to rad
    x_globe = binas.earth_radius * dtr * (lon - reflon) * np.cos(reflat * dtr)
    y_globe = binas.earth_radius * dtr * (lat - reflat)
    cos_wd = np.cos(-wind_direction * dtr)
    sin_wd = np.sin(-wind_direction * dtr)
    x_grid = x_globe * coswd + y_globe * sin_wd
    y_grid = -x_globe * sinwd + y_globe * cos_wd
    return x_grid, y_grid
    
# footprint mapper --> for if we want to oversample individual observations to create more contrast/sharpness


# plumewidth adjustment, shown to fit real life plumes better than basic function.
# plumewidth can also be calculated in more detail using radiation etc (Nassar, et al., 201?, CO2 powerplant paper), not used here
def function_adjust_plumewidth(y, plumewidth):
    if np.size(y) > 1:
        plumewidth_adj = np.full(len(y), plumewidth)
        sel_y = np.where(y <= 0)
        plumewidth_adj[sel_y] = np.sqrt((plumewidth ** 2) - 1.5 * y[sel_y])
    else:
        if y <= 0:
            plumewidth_adj = np.sqrt((plumewidth ** 2) - 1.5 * y)
        else:
            plumewidth_adj = plumewidth
    return plumewidth_adj

# crosswind diffusion
def flow_function_f(x, y, plumewidth):
    val = 1. / (function_adjust_plumewidth(y, plumewidth) * np.sqrt(2. * np.pi))
    val2 = val * np.exp(-(x ** 2.) / (2. * function_adjust_plumewidth(y, plumewidth) ** 2.))
    return val2

# downwind advection / diffusion
def flow_function_g(y, s, plumewidth, decay):
    decay_adj = decay / s
    val = (decay_adj / 2.) * np.exp((decay_adj * (decay_adj * plumewidth ** 2. + 2. * y)) / 2.)
    var1 = (decay_adj * (plumewidth ** 2.) + y) / (np.sqrt(2) * plumewidth)
    val2 = val * scipy.special.erfc(var1)
    return val2

# i.e. flowfun to calculate the factors in arrays etc
# convolution of f/g essentially.
def constant_flowfun(x1, x2, s1, decay, plumewidth): 
    fout = flow_function_f(x1, x2, plumewidth) * flow_function_g(x2, s1, plumewidth, decay)
    return fout




# legendre bias function or other

# setup_lin_matrix --> create the matrix for AX=B

# fitting operator?

# multi source helper --> if multisourcing operations


# flatten lists --> little helper to correct obsession with lists
def flatten_list(list_of_lists):
    return [item for sublist in list_of_lists for item in sublist]

# make nc file for list of variables, dimensions, attributes

# create basic emission fields from a inventory for comparison?

