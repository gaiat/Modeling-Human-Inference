
% This function computes the asymptotic Algorithmic Complexity per estimation for the Hierarchical Gaussian filter.
% Nlevels --> Number of levels in the model.

function[f]=AC_HGF_estim(Nlevels)
    
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
POW=1; % power

f = Mm+E+P+S+3*Mr+3*Mm+Mw+Mm+Mr+D+S+Mr+Mm+Mw+S+Mr+S+2*Mr+Mw+Mm+S+P+D+3*Mr+Mw+(Nlevels-1)*(4*S+POW+D+5*Mr+D+2*Mr+Mw+Mw+Mm+D+S+Mr+D+2*Mr+Mw+Mm+Mr+...
    3*S+6*P+POW+Mm+Mw+10*Mr+S+4*P+D+7*Mr+Mw) + (Nlevels-2)*(E+S+P+3*Mr+3*Mm+Mw) + Mm+Mr ;

end


