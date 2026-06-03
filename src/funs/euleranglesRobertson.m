function [alpha,beta,gamma] = euleranglesRobertson(RFAMA) 

% This function determines the euler angles (radians) from a matrix using a
% ZYX rotation sequence

% alpha= atan2(-RFAMA(3,2),RFAMA(3,3)); 
% beta= atan2(RFAMA(3,1),sqrt((RFAMA(1,1).^2+RFAMA(2,1).^2))); 
% gamma= atan2(-RFAMA(2,1),RFAMA(1,1)); 

alpha= atan2(-RFAMA(3,2),RFAMA(3,3)); 
beta= atan2(RFAMA(3,1),(sqrt(RFAMA(1,1).^2+RFAMA(2,1).^2))); 
gamma= atan2(-RFAMA(2,1),RFAMA(1,1)); 


