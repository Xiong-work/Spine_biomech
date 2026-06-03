function [h]=plot3Dpoints(dat);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Points
% Plots a group of points and returns a handle
% [h]=plotpoints(dat);
% dat: is a 1xnpoints matrix of points described by x-y-z coordinates
% [h]=plotpoints(dat) 

npoints=size(dat,2)/3;
xdat=dat(1:3:(npoints-1)*3+1);
ydat=dat(2:3:(npoints-1)*3+2);
zdat=dat(3:3:(npoints-1)*3+3);

plot3(xdat,ydat,zdat,'.'); 
xlabel('1');
ylabel('2');
zlabel('3');

for i=1:npoints;
   text(xdat(i),ydat(i),zdat(i),[' m' num2str(i)])
end;


