function [alphaknee2,betaknee2,gammaknee2]...
          =FindInstantDataRobertson(data,TMTAthigh,TMTAshank)

%This function determines the anatomical coordinate systems at every point
%in time using the markers for each triad and the transformation between
%marker and anatomical coordinate systems determined previously in
%'CoordSystems.m'

%This function also determines the euler angles at each point in time:
%Theta= Adduction(+'ve)/ Abduction
%Phi= Flexion(+'ve)/ Extension
%Psi= Internal(+'ve)/ External Rotation
%% Split the data (using the different format for Squat.csv
RSHPA = data(86:88); 
RSHPP = data(89:91); 
RSHDP = data(83:85); 
RSHDA = data(80:82); 
RTHPP = data(101:103); 
RTHPA = data(98:100); 
RTHDA = data(92:94); 
RTHDP= data(95:97); 
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
%% Determine [T 4x4] matrices for each anatomical coordinate system
%Thigh
TAthigh=TMthigh*TMTAthigh;
%Shank
TAshank=TMshank*TMTAshank;
%% Find the tranformation matrices between the fixed and moving ACS for
%% each segment
[KNEE]= (inv(TAthigh(1:3,1:3))* TAshank(1:3,1:3)); 
%% Calculate the Angles 
[alphaknee2, betaknee2, gammaknee2] = euleranglesRobertson(KNEE);
%%
% figure(); 
% plot3Dpoints(data(20:end)); 
% showcs(TAthigh,50); 
% showcs(TAshank,50); 