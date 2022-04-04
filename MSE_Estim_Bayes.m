
% This function computes the mean squared estimation error for the Bayesian
% model, over a Ttot=5000 time-long instance (iterat) of a Gaussian change-point process 
% with volatility h and noise sigma (=std of the observations around the source).

function [f]=MSE_Estim_Bayes(h,sigma,iterat)

warning('off','all')
sigma0=1;
mu0=0;
Ttot=5*10^3;

load(['samples/process_s',num2str(sigma*100),'_s0',num2str(sigma0*100),'_h',num2str(h*1000),'_it',num2str(iterat)],'process') 
nup=(sigma/sigma0)^2;
chip=mu0*nup;
nu_t=nup; 
chi_t=chip; 
mu_Bayes=zeros(Ttot,1);

for t=1:Ttot
    mu_LH=chi_t(1,:)./nu_t(1,:); 
    varl_t=(sigma^2).*(1+(1./nu_t(1,:)));     
    chi_t(1,t+1)=chi_t(1,t)+process(t,1);
    chi_t(1,[2:t])=chi_t(1,[2:t])+process(t,1)-process([t-1:-1:1],1)';
    nu_t(1,t+1)=nu_t(1,t)+1;
    mul_t=chi_t(1,[2:t+1])./nu_t(1,[2:t+1]); 
    if(t==1)
        pl=1; 
    else
        num_pl(1,1)=h*normpdf(process(t,1),mu_LH(1,1),sqrt(varl_t(1,1))); 
        num_pl(1,[2:t])=(1-h)*normpdf(process(t,1),mu_LH(1,[2:t]),sqrt(varl_t(1,[2:t]))).*pl; 
        pl=num_pl./sum(num_pl); 
    end
    mu_Bayes(t,1)=mul_t*pl';   
end

f=mean((mu_Bayes-process([1:Ttot],2)).^2)/(sigma0^2); 


end









