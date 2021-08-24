from setuptools import setup

setup(
   name='DimdorePy',
   version='0.0.1',
   author='Oscar Dimdore-Miles',
   author_email='oscar.dimdoremiles@gmail.com',
   packages=['DimdorePy', 'DimdorePy.test'],
   description='Analysis package for climate data',
   long_description=open('README.md').read(),
   install_requires=[
       "numpy",
       "pytest",
       "iris",
       "math"
       
   ],
)