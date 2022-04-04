
% This function computes the mean squared estimation error for the Estimate, 
% over a Ttot=5000 time-long instance (iterat) of a Gaussian change-point process 
% with volatility h and noise sigma (=std of the observations around the source).

function f=MSE_Estim_Evidence(h,sigma,iterat) 

warning('off','all')
sigma0=1;
Ttot=5*10^3;

load(['samples/process_s',num2str(sigma*100),'_s0',num2str(sigma0*100),'_h',num2str(h*1000),'_it',num2str(iterat)],'process') 

f = (mean((process([1:Ttot],1)-process([1:Ttot],2)).^2))/(sigma0^2);


end

