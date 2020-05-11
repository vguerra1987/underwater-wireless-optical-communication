function [T, S] = getParamsFromDepth(D)
% This function returns bathymetry-related information about temperature
% and salinity respect to depth. In this version, information from NOAA was
% downloaded and stored in a separate file. If there were an API, this
% could be dynamically retrieved according to a geographical position.
% Current location corresponds to the Gran Canaria's area.
% Salinity is in parts per million, and the other functions receive it in
% parts per thousand. Hence, proper scaling is applied.

data = readtable('NOAA.csv');

T = interp1(data.Depth, data.Temperature, D);
S = interp1(data.Depth, data.Salinity, D)/1000;