
% This function computes the mean squared estimation error for a Mixture of Sliding Windows with N units, 
% over a Ttot=5000 time-long instance (iterat) of a Gaussian change-point process 
% with volatility h and noise sigma (=std of the observations around the
% source); x is a column vector that specifies the model learning rates.


function f=MSE_Estim_SA_N(x,h,sigma,iterat,N) 

warning('off','all')
sigma0=1;
mu0=0;
Ttot=5*10^3;

load(['samples/process_s',num2str(sigma*100),'_s0',num2str(sigma0*100),'_h',num2str(h*1000),'_it',num2str(iterat)],'process') 
nup=(sigma/sigma0)^2;
mul=mu0*ones(N,1);
varl=zeros(N,1);
runlength=zeros(N,1);
runl_integer=zeros(N,1);
d=zeros(N,1);
logLH=zeros(N,1);
pl=(1/N)*ones(N,1); 
estim=zeros(Ttot,1);
log_num=zeros(N,1);
logA=zeros(N,1);

x=sort(x,'descend');
for indl=1:N
    runlength(indl,1)=(1/x(indl,1))-nup;
    varl(indl,1)=(sigma^2)*(1+x(indl,1));
    runl_integer(indl,1)=floor(runlength(indl,1));
    d(indl,1)=runlength(indl,1)-runl_integer(indl,1);
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

rows_small=find(runlength(:,1)<Ttot/2);
[~,max_rows_small]=max(runlength(rows_small,1));
for ww=1:length(rows_small)
    mul(rows_small(ww,1),1)=x(rows_small(ww,1),1)*(mu0*nup+sum(process([runl_integer(max_rows_small,1)-runl_integer(rows_small(ww,1),1)+2:runl_integer(max_rows_small,1)+1],1))+...
        d(rows_small(ww,1),1)*(process(runl_integer(max_rows_small,1)+1-runl_integer(rows_small(ww,1),1),1))); 
end

for tp=runl_integer(max_rows_small,1)+2:Ttot
    for indl=1:N
        if(any(rows_small==indl)) 
            mul(indl,1)=mul(indl,1)+ x(indl,1)*(process(tp,1)-(1-d(indl,1))*process(tp-runl_integer(indl,1),1)-d(indl,1)*process(tp-runl_integer(indl,1)-1,1)); 
        end
        logLH(indl,1)= - log(sqrt(2*pi*varl(indl,1))) - (((process(tp,1)-mul(indl,1))^2)/(2*varl(indl,1)));        
    end
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
    estim(tp,1)=pl'*mul;
end

if(isempty(max_rows_small))
    f = 1;
else
    f = (mean(( estim([runl_integer(max_rows_small,1)+2:Ttot],1)-process([runl_integer(max_rows_small,1)+2:Ttot],2) ).^2))/(sigma0^2); 
end



end





