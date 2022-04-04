
% This function computes the mean squared estimation error for a single Sliding Window with learning rate x, 
% over a Ttot=5000 time-long instance (iterat) of a Gaussian change-point process 
% with volatility h and noise sigma (=std of the observations around the
% source).

function f=MSE_Estim_SA_1(x,h,sigma,iterat) 

warning('off','all')
sigma0=1;
mu0=0;
Ttot=5*10^3;

load(['samples/process_s',num2str(sigma*100),'_s0',num2str(sigma0*100),'_h',num2str(h*1000),'_it',num2str(iterat)],'process') 
nup=(sigma/sigma0)^2;

runl_integer=floor((1/x)-nup);
d=(1/x)-nup-runl_integer;
if(runl_integer+d>=Ttot/2)    
    f = 1;    
else   
    mul=zeros(Ttot,1);
    mul(runl_integer+1,1)=x*(mu0*nup+sum(process([2:runl_integer+1],1))+d*(process(1,1)));
    for t=runl_integer+2:Ttot
        mul(t,1)=mul(t-1,1)+ x*(process(t,1)-(1-d)*process(t-runl_integer,1)-d*process(t-runl_integer-1,1));  
    end
    f = (mean(( mul([runl_integer+1:Ttot],1)-process([runl_integer+1:Ttot],2) ).^2))/(sigma0^2); 
end


end






