%THIS PROGRAM CALCULATES 3-D ANKLE ANGLES
clear all;
close all;
%% Addpath
here= pwd;
% whereto= 'C:\Users\rgraham\Desktop\APA 6903 - Grad Biomech\2017\3D Kinematics Code';
whereto= 'E:\Rgraham\rgraham\Desktop\APA 6903 - Grad Biomech\2017\3D Kinematics Code';

cd(whereto);
%% Input Constants
kinematicrate= 120; %hz
%% Create the anatomical coordinate systems (+Plot) and determine the Euler
%% angle biases
         [TMTAthigh, TMTAshank,...
          alphakneebias,betakneebias,gammakneebias]= CoordSystemsRobertsonGood;
%% Load in the squat trial
filename= 'Squat.csv';
data= csvread(filename,5,0);
time= data(:,1)/kinematicrate; 
%% Split out V3D knee angles
kneeangles= data(:,14:16); 
%% Determine the instantaneous tranformation matrices for each segment
%Initialize the variables
x= size(data,1); 
alphaknee2= zeros(1,x); 
betaknee2= zeros(1,x); 
gammaknee2= zeros(1,x); 

%Run them through a loop
hw = waitbar(0,'Calculating the Euler Angles');

for i= 1:x
         [alphaknee2(i),betaknee2(i),gammaknee2(i)]...
          =FindInstantDataRobertson(data(i,:)*1000,TMTAthigh,TMTAshank);
  waitbar(i/x,hw) ;
end
close(hw); 
%% Find the Euler angles at each point in time 
% Knee
alphaknee= rad2deg(unwrap(alphaknee2-alphakneebias));
betaknee= rad2deg(unwrap(betaknee2-betakneebias));
gammaknee= rad2deg(unwrap(gammaknee2-gammakneebias));
% Knee
alphakneebias= rad2deg((alphaknee2));
betakneebias= rad2deg((betaknee2));
gammakneebias= rad2deg((gammaknee2));

%% Plot the angles
figure(2); hold on; 
plot(alphaknee,'b');
plot(betaknee,'r');
plot(gammaknee,'g'); 
legend({'FE','Ab/Adduction','Rotation'}); 
title('Standing Bias Removed'); 
% 
%% Plot the angles
figure(3); hold on; 
plot(alphakneebias,'b');
plot(betakneebias,'r');
plot(gammakneebias,'g'); 
legend({'FE','Ab/Adduction','Rotation'}); 
% 
%% Plot the angles
figure(4); hold on; 
plot(time,alphaknee-alphaknee(1),'b');
plot(time,betaknee-betaknee(1),'r');
plot(time,gammaknee-gammaknee(1),'g'); 
plot(time,kneeangles(:,1)-kneeangles(1,1),'b:');
plot(time,kneeangles(:,2)-kneeangles(1,2),'r:');
plot(time,kneeangles(:,3)-kneeangles(1,3),'g:'); 
legend({'FE','Ab/Adduction','Rotation'}); 
xlabel('Time(s)'); 
ylabel('Degrees'); 
title('Matlab (solid) vs. Visual 3D (dashed) - no biases'); 
