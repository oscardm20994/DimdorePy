#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Jan 11 17:13:01 2020

@author: Oscar Dimdore-Miles

script that contains taylor made tools for processing Iris data cubes

iris reference: https://scitools-iris.readthedocs.io/en/stable/
"""

import iris
import numpy as np
import iris.coord_categorisation as coord_cat


def add_times(self, cube, time):
		"""Adds time coordinate metadata to an input iris cube. 

		Args
		----------
		cube : iris.cube.Cube object
	    An iris cube with a time dimension.

		time : str
		String indicating the standard name of the time dimension of the
		input cube. E.G 'time' or 't' etc

		Returns
		----------
		cube: iris.cube.Cube object
		same iris cube as was input with added time coordinate metadata

		"""

		# add dimensions for month, clim_season (indicating season e.g. 'DJF'),
		# year, day_number (day in the year) and season_year.
		coord_cat.add_month(cube, time, name='month')
		coord_cat.add_season(cube, time, name='clim_season')
		coord_cat.add_year(cube, time, name='year')
		coord_cat.add_day_of_year(cube, time, name='day_number')
		coord_cat.add_season_year(cube, time, name='season_year')
		
		return cube



def seasonally_detrend(self, cube, tname, year_len = 360):
		"""Function that seasonally detrends the data in an input iris cube. 
		This is done by estimating the seasonal component to data in the cube
		(the saesonal cycle) and subtracting this component from the data.

		Args
		----------
		cube : iris.cube.Cube object
    		An iris cube with a time dimension.
    
		year_len : int
    		integer indicating the number of timepoints in each year of the
    		cube's data. Default value is 360 day calendar in pi-control model data.
    		Also works with 366 for daily observation data.

		tname : str
    		String indicating the standard name of the time dimension of the
    		input cube. E.G 'time' or 't' etc.

		Returns
		----------
		cube: iris.cube.Cube object
    		same iris cube as was input with added time coordinate metadata

		"""
        
		# get the number of time points in the cube   
		ntime = len(cube.coord(tname).points)

		# find seasonal cycle by aggregating
		cycle = cube.aggregated_by('day_number', iris.analysis.MEAN).data,

		# repeat this seasonal cycle over the time dimension for the purpose
		# of subtracting
		cycle_repeated = np.tile(cycle, [int(math.ceil(ntime/year_len)), 1])[0,:]

		# subtract off repeated saeasonal cycle from the cube's data
		detrended = iris.analysis.maths.subtract(cube, cycle_repeated)

		# return the detrended cube
		return detrended