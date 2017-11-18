# GeoMat
It is a collection of Matlab functions for exploration geophysics. This project also includes Chen Qi's research codes related to his findings in stratigraphic filtering.

Installation
============
To use all functionalities provided by GeoMat, you only need to add {ROOT}/src folder with subfolders into Matlab search path. {ROOT} represents the path where you copied your GeoMat project.

File Structure
===============
src
  - plotting: all figure related functions are in this folder, for example, wiggle-variable-area plot function.
  - signal_processing_utility: all common signal processing functions are in this folder, for example, bandpass filter and convolution.
  - stratigraphic_filtering: source code for Chen and Hilterman paper: Well ties for severe stratigraphic filtering. Please see examples.m for more details.
