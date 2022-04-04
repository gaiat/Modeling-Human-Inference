
% This function computes the mean squared estimation error for a Mixture of Delta Rules with N units, 
% over a Ttot=5000 time-long instance (iterat) of a Gaussian change-point process 
% with volatility h and noise sigma (=std of the observations around the source); 
% x is a column vector that specifies the model learning rates.

function f=MSE_Estim_DR_N(x,h,sigma,iterat,N) 

warning('off','all')
sigma0=1;
mu0=0;
Ttot=5*10^3;
load(['samples/process_s',num2str(sigma*100),'_s0',num2str(sigma0*100),'_h',num2str(h*1000),'_it',num2str(iterat)],'process') 
nup=(sigma/sigma0)^2;

mul=mu0*ones(N,1);
varl=zeros(N,1);
runlength=zeros(N,1);
logLH=zeros(N,1);
mu_pred=zeros(Ttot,1);
log_num=zeros(N,1);
logA=zeros(N,1);
x=sort(x,'descend');
for indl=1:N
    runlength(indl,1)=(1/x(indl,1))-nup;
    varl(indl,1)=(sigma)^2*(1+x(indl,1));
end
pc=zeros(N,N);
pc(1,[1:N])=ones(1,N);
pnc=zeros(N,N);
pnc(N,N)=1;

for indl=2:N
    if(runlength(indl,1)>runlength(indl-1,1)+1)
        pnc(indl,indl-1)=1/(runlength(indl,1)-runlength(indl-1,1));
    else
        pnc(indl,indl-1)=1;
    end
end
for indl=1:N-1
    if(runlength(indl+1,1)>runlength(indl,1)+1)
        pnc(indl,indl)=(runlength(indl+1,1)-runlength(indl,1)-1)/(runlength(indl+1,1)-runlength(indl,1));
    else
        pnc(indl,indl)=0;
    end
end

for tp=1:Ttot
    for indl=1:N
        mul(indl,1)=mul(indl,1)+x(indl,1)*(process(tp,1)-mul(indl,1));
        logLH(indl,1)= - log(sqrt(varl(indl,1)*2*pi)) - (((process(tp,1)-mul(indl,1))^2)/(2*varl(indl,1))); 
    end
    if(tp==1)
        pl=zeros(N,1); 
        pl(1,1)=1;
    else
        for indl=N:-1:1 
            log_num(indl,1)=logLH(indl,1)+log((h*pc(indl,:)+(1-h)*pnc(indl,:))*pl);
            logA(indl,1)=log_num(indl,1)-log_num(N,1); 
        end
        logANan=find(isnan(logA));
        if(isempty(logANan))
            [logA_max,Ind_max]=max(logA); 
            if((1/exp(logA_max))==0)
                pl=zeros(N,1);
                pl(Ind_max,1)=1;
            else
                for indl=1:N
                    pl(indl,1)= (exp(logA(indl,1)))/(sum(exp(logA)));
                end
            end
        else
            pl=zeros(N,1);
            pl(logANan(1,1),1)=1;
        end
    end
    mu_pred(tp,1)=pl'*mul;
end

f = (mean((mu_pred-process([1:Ttot],2)).^2))/(sigma0^2); 


end

