function [TMTAthigh, TMTAshank,...
          alphakneebias,betakneebias,gammakneebias]= CoordSystemsRobertson2
% This function determines the anatomical coordinate systems for the
% bodies, as well as the relationships between the marker and anatomical
% coordinate systems!

% This function also plots the coordinate systems and markers to check if
% they look right, and determines the Euler Angle Offsets from a standing
% calibration trial
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
%% Construct technical coordinate systems for the right thigh and shank
[ithigh,jthigh,kthigh, originthigh]= create_rhcs(RTHDP,RTHDA,RTHPP);
% [ishank,jshank,kshank, originshank]= create_rhcs(RSHDP,RSHDA,RSHPP);
[ishank,jshank,kshank, originshank]= create_rhcs(RSHPA,RSHPP,RSHDP);

thigh1= [ithigh' jthigh' kthigh'];
det(thigh1)
shank1= [ishank' jshank' kshank'];
det(shank1)
%% Create [T 4x4] matrices for each technical coordinate system
%Thigh
TMthigh = eye(4);
TMthigh(1:3,1:3)=thigh1;
TMthigh(1:3,4) = originthigh';
%Shank
TMshank = eye(4);
TMshank(1:3,1:3)=shank1;
TMshank(1:3,4) = originshank';
%Plot
showcs(TMthigh,50); 
showcs(TMshank,50);
%% Determine Joint Centres
%% Do pelvis calculations as per Robertson Equations 2.11 - 2.16
originpelvis= 0.5*(RASIS+LASIS); 
i= RASIS - originpelvis;
v= originpelvis - 0.5 * (RPSIS+LPSIS); 
k= cross(i,v); 
j= cross(k,i); 
ipelv= i/norm(i); 
jpelv= j/norm(j);
kpelv= k/norm(k); 
pelvis1= [ipelv' jpelv' kpelv']; % equat 2.16
det(pelvis1)
%Pelvis
TApelvis = eye(4);
TApelvis(1:3,1:3)=pelvis1;
TApelvis(1:3,4) = originpelvis';
showcs(TApelvis,50); 
% Do thigh segment calculations as per Robertson Equations 2.17-2.23
HIPnorm = norm(RASIS-LASIS);
PHIP = [0.36*HIPnorm -0.19*HIPnorm -0.30*HIPnorm];
HipJC= pelvis1*PHIP'+originpelvis';
%%
hold on;
plot3(HipJC(1),HipJC(2),HipJC(3),'ro','MarkerSize',10,'MarkerFaceColor','r'); % plot a big red dot at Hip JC
KneeJC=RMK+((RLK-RMK)/2);
AnkleJC=RMA +((RLA-RMA)/2);
plot3(KneeJC(1),KneeJC(2),KneeJC(3),'mo','MarkerSize',10,'MarkerFaceColor','m'); % plot a big magenta dot at Knee JC
plot3(AnkleJC(1),AnkleJC(2),AnkleJC(3),'co','MarkerSize',10,'MarkerFaceColor','c'); % plot a big cyan dot at Ankle JC
%% Create anatomical CS, where X= superiorly, Y = anteriorly, Z= medially
%Thigh
k = (HipJC' - KneeJC)/norm(HipJC' - KneeJC);
v = (RLK - RMK)/ norm(RLK - RMK);
j = cross(k,v);
i = cross(j,k);
Rthigh = [i' j' k'];

det(Rthigh); % ans = 1 = right handed

TAthigh = eye(4);
TAthigh(1:3,1:3) = Rthigh;
TAthigh(1:3,4) = HipJC;
%Shank
k = (KneeJC - AnkleJC)/norm(KneeJC - AnkleJC);
v = (RLK - RMK)/norm(RLK - RMK);
j = cross(k,v);
i = cross(j,k);

Rshank = [i' j' k'];

det(Rshank) %ans = 1 = right handed

TAshank = eye(4);
TAshank(1:3,1:3) = Rshank;
TAshank(1:3,4) = KneeJC;
showcs(TAthigh,50);
showcs(TAshank,50);
%%
%% Find transformation matrices (T) between technical and anatomical coordinate systems
[TMTAthigh]=(inv(TMthigh))*TAthigh; 
[TMTAshank]=(inv(TMshank))*TAshank; 
%% Find the tranformation matrices between the fixed and moving ACS for
%% each segment (with the origin fixed at the COM)
[KNEE]= (inv(TAthigh(1:3,1:3))* TAshank(1:3,1:3)); 
%% Calculate the Angle Bias
[alphakneebias, betakneebias, gammakneebias] = euleranglesRobertson(KNEE);
