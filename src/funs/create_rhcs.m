function [i,j,k,origin]= create_rhcs(m1,m2,m3)
% This function creates an arbitrarily oriented Right Hand Coordinate
% System using the data from 3 markers. 
% The dimension of output unit vectors (i,j and k) is the same as the
% dimension of the input data points. That is if m1 m2 and m3 are given as
% 1x3 vectors then i j and k will be each 1x3 unit vectors

v1= m1 - m2;
v2= m3 - m2; 
v3= cross(v1,v2); 
v4= cross(v3,v1);
i= v1/norm(v1);
j= v4/norm(v4);
k= v3/norm(v3); 
origin= m2;