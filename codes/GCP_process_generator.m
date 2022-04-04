
% This code generates samples from Gaussian change-point processes with
% different values of noise R and volatility h (examples in Fig.2).

nb_iterat=10;
sigma0=1;
mu0=0;
Ttot=5*10^3;
h_all_sim=[0.02:0.02:1]; 
R_all_sim=[0.01:0.01:0.09,0.1:0.1:6]; 


for q=1:length(R_all_sim)
    R=R_all_sim(1,q);
    sigma=R*sigma0;
    for k=1:length(h_all_sim)
        h=h_all_sim(1,k);
        for iterat=1:nb_iterat
            process=zeros(Ttot,2);
            process(1,2)=mu0+sigma0*randn;
            process(1,1)=process(1,2)+sigma*randn;
            for t=2:Ttot
                rand_nb=rand;
                if(rand_nb<h)
                    process(t,2)=mu0+sigma0*randn;
                else
                    process(t,2)=process(t-1,2);
                end
                process(t,1)=process(t,2)+sigma*randn;
            end
            save(['samples/process_s',num2str(sigma*100),'_s0',num2str(sigma0*100),'_h',num2str(h*1000),'_it',num2str(iterat)],'process')        
        end        
    end    
end

