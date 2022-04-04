
% This function computes the asymptotic Algorithmic Complexity per estimation for the Kalman filter.

function[f]=AC_Kalman_estim()
    
Mw=1; % writing in memory
Mr=1; % reading from memory
Mm=1; % maintaining in memory
S=1; % additions+subtractions
P=1; % product
D=1; % division
E=1; % exponential
F=1; % floor 
R=1; % square root
OPP=1; % opposite 

f = (5*S+3*P + D) + 12*Mr + 2*Mm + 4*Mw;

end


