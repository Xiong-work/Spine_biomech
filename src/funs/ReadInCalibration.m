%% ASSIGNMENT 4 - 3d kinematics
% Read in data and plot to see how it works
% Start the process of developping coordinate systems etc.
%% Anatomical Markers Labels
% Right Heel = 1:3
% Right Metatarsal 1 = 4:6
% Right Metatarsal 5 = 7:9
% Right Lateral Ankle = 10:12
% Right Medial Ankle = 13:15 
% Right Lateral Knee = 16:18
% Right Medial Knee = 19:21
% Right Greater Trochanter = 34:36
% Left Iliac Crest = 49:51
% Right Iliac Crest = 52:54
% Right PSIS = 55:57
% Left PSIS = 58:60
% Left ASIS = 61:63
% Right PSIS = 64:66; 
%% Cluster Markers Labels
% RSHPA = 22:24
% RSHPP = 25:27
% RSHDP = 28:30
% RSHDA = 31:33
% RTHPP = 37:39
% RTHPA = 40:42
% RTHDA = 43:45
% RTHDP= 46:48
%% Read in data
data= csvread('CalibrationGood.csv',5,0);
markers= data(:,3:end);
%% Split into markers
%% Anatomical
frame= 1; 
RH = markers(frame,1:3);
RMT1 = markers(frame,4:6);
RMT5 = markers(frame,7:9);
RLA = markers(frame,10:12);
RMA = markers(frame,13:15);
RLK = markers(frame,16:18);
RMK = markers(frame,19:21);
RGT = markers(frame,34:36);
LIC= markers(frame,49:51);
RIC = markers(frame,52:54);
RPSIS= markers(frame,55:57);
LPSIS = markers(frame,58:60);
LASIS = markers(frame,61:63);
RASIS = markers(frame,64:66);
%% Clusters
RSHPA = markers(frame,22:24); 
RSHPP = markers(frame,25:27);
RSHDP = markers(frame,28:30);
RSHDA = markers(frame,31:33);
RTHPP = markers(frame,37:39);
RTHPA = markers(frame,40:42);
RTHDA = markers(frame,43:45);
RTHDP= markers(frame,46:48);
%% Plot the markers together on one graph;
plot3Dpoints(markers(1,:));
axis equal; 
