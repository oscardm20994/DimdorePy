#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sep 20 13:22:27 2019

@author: Oscar Dimdore-Miles

script to define the strat class. Instances of this class
contain attributes of key climate indices which track variability in
atmospheric circulation of the stratosphere. This class is wrapped around the
Metoffice analysis package iris (reference link below) and many methods of the
class require pre-loaded iris data cubes to run. 

iris reference: https://scitools-iris.readthedocs.io/en/stable/
"""

# import necessary packages
import cubetools.add_times as add_times
import numpy as np
import iris
import cubetools.seasonally_detrend as seasonally_detrend
from iris.analysis.cartography import cosine_latitude_weights
import iris.analysis.stats
from eofs.iris import Eof


## define strat class
class strat:
	"""
    A class used to represent a set of stratospheric variability indices 
    derived from a given dataset (e.g. a model simulation or reanalysis).

    ...

    Attributes
    ----------
    dataset_name : str
        a string representing the dataset on which the instance is based. 
        
    U_QBO : iris.cube.Cube object
        Iris cube containing the monthyl QBO zonal mean zonal winds.
        
    QBO_spectrum : numpy array
        Numpy array of the Fourier Power spectrum of QBO winds at
        different pressure levels.
        
    QBO_freqs: numpy array
        1D numpy array with the sample frequencies corresponding to the 
        frequency dimension of QBO_spectrum. 
        
    NAM_pc : dict 
       dictionary containing the timeseries of the Northern Annular Mode, the 
       1st principle component of NH daily geopotential height. Keys of this 
       dict represent different pressure levels on which the NAM is calculated
    
    NAM_eof : dict
        dictionary containing the 1st empirical orthogonal function of
        daily geopotential height. Keys of this dict represent different
        pressure levels on which the EOF is calculated.
    
    
    Methods
    -------    
    
    QBO_indices(self, U)
        Calculates QBO wind timeseries and associated Fourier power spectrum.
        Assigns these as attributes to the instance
    
    NAM(self, GPH)
        Calculates the northern annular mode index on a given pressure level.
        Returns the 1st Empitical Orthogonal Function (EOF) and Principle 
        component timeseries of daily NH geopotential height.
    
    calculate_NAMS(self, GPH, pname, tname)
        Evaluates the PC and EOFs of GPH data on a range of pressure levels
        using the NAM method. Adds the PC and EOF cubes to NAM_eof and NAM_pc
        attributes of the instance.
        
        """
    
    
    ## init method for the strat class asks user to assign dataset name 
    ## and the years over which the data are evaluated
	def __init__(self, dataset_name, years):
        
        #raise  type error if dataset_name is not a string
		if not isinstance(dataset_name, str):
			raise TypeError("dataset_name must be a string")    
        
        # initialise the attributes of the instance
		self.dataset_name = dataset_name
		self.years = years
		self.QBO_spectrum = None
		self.QBO_freqs = None
		self.NAM_pc = {}
		self.NAM_eof = {}
        

	def QBO_indices(self, U):
		"""Calculates the Quasi Biennial Oscillation timeseries and 
		corresponding fourier power spectrum at a range of pressure levels. 
		Assigns this data to attributes of the instance. The QBO qinds are
		defined in Baldwin et al. 2001 (link below) as the equatorial zonal 
		mean zonal wind in the stratosphere.
		(https://agupubs.onlinelibrary.wiley.com/doi/10.1029/1999RG000073)

		Args
		----------
		U : iris.cube.Cube object
			A cube provided by the user containing the monthly zonal mean zonal 
			wind (ZMZW) from the dataset. Dimensions of this cube should be time, 
			longitude and pressure level. 

		Raises
		------
		TypeError if argument U is not an instance of iris.cube.Cube
		"""
 
		# raise type error if U is not an iris cube
		if not isinstance(U, iris.cube.Cube):
			raise TypeError("U must be an iris cube")  

		# define latitude constraint between 5N and 5S
		lat_constraint=iris.Constraint(latitude=lambda cell: -5.1 < cell < 5.1)
		
		# restrcit latitudes to 5S-5N 
		U_QBO = U.extract(lat_constraint)

		# average over latitudes using the cosine of the latitude of each
		# gridbox as weights to take into account gridbox area differences
		# across the range.
		U.coord('latitude').guess_bounds()
		cos_weights = cosine_latitude_weights(U_QBO)
		U_QBO = U_QBO.collapsed('latitude', iris.analysis.MEAN, 
                  weights = cos_weights) 
	    
		# assign QBO winds to instance
		self.U_QBO = U_QBO

		# dt in months should always be 1 for monthly data. 
		dt = 1

		# compute fast fourier transform
		fft_U = np.fft.fft(U_QBO.data, axis = 0)
 		
		# fetch associated frequencies
		freq = np.fft.fftfreq(U_QBO.data[:,0].size, d=dt)
	
		# throw away negative frequencies
		keep = freq>=0
		fft_U = fft_U[keep,:]
		freq = freq[keep]

		# calculate the spectrum (absolute value of FFT)
		# and normalise by the standard deviation of the QBO timeseries
		# at each pressure level.
		fourier_spectrum = np.array([np.abs(fft_U[:,i])/np.std(U_QBO.data, 
                        axis = 0)[i] for i in range(len(fft_U[0,:]))])

		# assign numpy arrays for the fourier power spectrum 
		# and frequency array
		self.QBO_spectrum = fourier_spectrum
		self.QBO_freqs = freq
		
		return
        
	def NAM(self, GPH, tname):
		"""Calculates the Northern Annular Mode timeseries 

		Args
		----------
		GPH : iris.cube.Cube object
    		A cube provided by the user containing the daily zonal mean
    		geopotential height data on a single pressure level from the 
    		dataset.
    
		tname: str
    		String indicating the standard name of the time dimension of the
    		input cube. E.G 'time' or 't' etc.
    

		Raises
		------
		TypeError if argument GPH is not an instance of iris.cube.Cube

		TypeError if argument tname is not a string

		"""
        
		# raise  type error if tname is not a string
		if not isinstance(tname, str):
			raise TypeError("tname must be a string")  

		# raise  type error if GPH is not an iris cube
		if not isinstance(GPH, iris.cube.Cube):
			raise TypeError("GPH must be an iris cube")  

		# add extra time coordinates to the GPH cube
		GPH = add_times(GPH, tname)
        
		# define latitude constraint object between 20N and 90N
		lat_constraint=iris.Constraint(latitude=lambda cell: 19.1 < cell < 91)
		
		# restrcit latitudes
		GPH = GPH.extract(lat_constraint)

		# seasonally detrend the GPH cube
		GPH_anom = seasonally_detrend(GPH, 'time')

		# extract winter months from the GPH data
		winter_anom = GPH_anom.extract(iris.Constraint(month = ['Dec',
                                                'Jan', 'Feb','Mar']))
        
		# calculate the 1st empirical orthogonal function of the GPH anomaly
		# field. This operation collapses the time dimension and provides the 
		# latitude pattern that maximises the GPH variance. The 1st EOF in GPH
		# gives an indication of the strenght of atmospheric circulation which
		# exhibits significant NH winter variability.
		# These lines uses the iris eof analysis package.      
		solver = Eof(winter_anom)

		# extract just the 1st EOF pattern. This is a 1D numpy array in the
		# latitude dimension.
		eof1 = solver.eofs(neofs = 1)[0,:]

		# The EOF pattern produced either indicates a strong or weak
		# circulation pattern. For ease of use, check which type of pattern
		# is calculated by evaluating the NS gradient in the pattern.
		# If the pattern corresponds to negative circulation, switch the sign.
		if eof1[-1].data > 0:
			eof1 = -1*eof1

		# wrap EOF pattern around time dimension.
		ntime = len(GPH_anom.coord('time').points)
		eof_wrap = np.tile(eof1.data[np.newaxis,...], (ntime,1))

		# calculate the principle component timeseries by projecting the
		# wrapped eof pattern onto anomaly field.
		PC_non_collapsed = iris.analysis.maths.multiply(GPH_anom, eof_wrap)
		PC = PC_non_collapsed.collapsed(['latitude'], iris.analysis.SUM)

		#normalise the PC timeseries by the mean and STD from winter months
		PC_DJFM = PC.extract(iris.Constraint(month = ['Dec', 'Jan', 
                                                'Feb', 'Mar']))
		mean = np.mean(PC_DJFM.data)
		std = np.std(PC_DJFM.data)
		PC = (PC - mean)/std

		#return PC timeseries and EOF pattern in iris cubes
		return PC, eof1

	def calculate_NAMS(self, GPH, pname, tname):
		"""calls the NAM method to calculate the Northern Annular Mode 
		timeseries on a range of pressure levels. Fill the NAM_eof and NAM_pc
		dictionaries of the instance. 

		Args
		----------
		GPH : iris.cube.Cube object
			A cube provided by the user containing the daily zonal mean
    		geopotential height data on a range of pressure levels from the 
    		dataset.
    
		tname: str
    		String indicating the standard name of the time dimension of the
    		GPH input cube. E.G 'time' or 't' etc.

		pname: str
    		String indicating the standard name of the pressure dimension of 
    		the input cube. E.G 'pressure' or 'pressure_level' etc.
    
		"""

		# get the pressure coordinate values out of the GPH cube
		Ps = GPH.coord(pname).points

		# loop over the pressure levels
		for current_p in Ps:
    
    		#constrain the cube to a single level 
			GPH_single_level = GPH.extract(iris.Constraint(
        	air_pressure = current_p))
    
    		# call the NAM function for GPH on the current level. 
			PC, EOF = self.NAM(GPH_single_level, tname)
    
			# assign the PC and EOF for the current pressure level to the
			# instance's dicts. assign the key as a string of the
			# current pressure level.
			self.NAM_eof[str(current_p) + 'hPa'] = EOF
			self.NAM_pc[str(current_p) + 'hPa'] = PC


		return  


