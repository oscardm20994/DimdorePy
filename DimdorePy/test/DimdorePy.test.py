# testing script for DimdorePy 
# conducts unit testing with pytest on the package

import pytest
from DimdorePy.cubetools import add_times, seasonally_detrend
import iris

#load a sample test iris cube to conduct tests on
testing_cube = iris.load('../data/test_cube.nc')[0]

def test_time_metadata_added():
    
    #call the add_times function
    new_cube = add_times(testing_cube, 't')
    
    #assert expected metadata has been added to the cube
    assert new_cube.coord('month') != None
    assert new_cube.coord('clim_season') != None
    assert new_cube.coord('year') != None
    assert new_cube.coord('day_number') != None
    assert new_cube.coord('season_year') != None
