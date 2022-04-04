
% This function computes the mean squared estimation error for a Hierarchical Gaussian Filter with 2 levels (Nlevels=2), 
% over a Ttot=5000 time-long instance (iterat) of a Gaussian change-point process 
% with volatility h and noise sigma (=std of the observations around the source); 
% x is a column vector that specifies (in this order): prior mean (first level), prior
% variance (first level), prior mean (second level), prior variance (second level), 
% k (first level), omega (first level), theta (as defined in Mathys et al. Front Hum Neurosci, 2014).

function f=MSE_Estim_HGF(x,h,sigma,iterat)  

warning('off','all')
Nlevels=2; 
sigma0=1;
Ttot=5*10^3;

load(['samples/process_s',num2str(sigma*100),'_s0',num2str(sigma0*100),'_h',num2str(h*1000),'_it',num2str(iterat)],'process') 

mui_hat=zeros(Ttot,Nlevels);
pii_hat=zeros(Ttot,Nlevels);
mui=zeros(Ttot,Nlevels);
pii=zeros(Ttot,Nlevels);
nui=zeros(Ttot,Nlevels);
deltai=zeros(Ttot,Nlevels-1);
sigmai=zeros(Ttot,Nlevels);
deltau=zeros(Ttot,1);
ki=zeros(1,Nlevels-1);
omegai=zeros(1,Nlevels-1);

piu_hat=1/(sigma^2);
for ii=1:Nlevels
    mui(1,ii)=x(2*(ii-1)+1,1);  
    sigmai(1,ii)=x(2*(ii-1)+2,1);
    pii(1,ii)=1/sigmai(1,ii);    
end
for ii=1:Nlevels-1
    ki(1,ii)=x(2*Nlevels+2*(ii-1)+1,1);
    omegai(1,ii)=x(2*Nlevels+2*(ii-1)+2,1);
end
theta=x(end,1);

for tp=2:Ttot+1
    mui_hat(tp,1)=mui(tp-1,1);
    nui(tp,1)=exp(ki(1,1)*mui(tp-1,2)+omegai(1,1));
    pii_hat(tp,1)=1/(nui(tp,1)+sigmai(tp-1,1));
    deltau(tp,1)=process(tp-1,1)-mui_hat(tp,1);
    pii(tp,1)=pii_hat(tp,1)+piu_hat;
    mui(tp,1)=mui_hat(tp,1)+(piu_hat/pii(tp,1))*deltau(tp,1);
    for ii=1:Nlevels-1    
        sigmai(tp,ii)=(pii(tp,ii))^(-1);
        deltai(tp,ii)=((sigmai(tp,ii)+(mui(tp,ii)-mui_hat(tp,ii))^2)/(sigmai(tp-1,ii)+nui(tp,ii)))-1;
        mui_hat(tp,ii+1)=mui(tp-1,ii+1);      
        if(ii==Nlevels-1)
            nui(tp,ii+1)=theta;
        else
            nui(tp,ii+1)=exp(ki(1,ii+1)*mui(tp-1,ii+2)+omegai(1,ii+1));
        end
        pii_hat(tp,ii+1)=1/(sigmai(tp-1,ii+1)+nui(tp,ii+1));
        pii(tp,ii+1)=pii_hat(tp,ii+1)+(1/2)*((ki(1,ii)*nui(tp,ii)*pii_hat(tp,ii))^2)*(1+(1-1/(nui(tp,ii)*pii(tp-1,ii)))*deltai(tp,ii));
        mui(tp,ii+1)=mui_hat(tp,ii+1)+(1/2)*ki(1,ii)*nui(tp,ii)*((pii_hat(tp,ii))/(pii(tp,ii+1)))*deltai(tp,ii);
    end
end
estimate_model=mui([2:Ttot+1],1);

f = (mean((estimate_model-process([1:Ttot],2)).^2))/(sigma0^2); 


end


