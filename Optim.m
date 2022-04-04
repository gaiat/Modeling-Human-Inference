
% This code computes the optimal parameter value(s) of the Delta-Rule model (with N units), 
% the associated minimum mean squared error and Hessian at the optimum, for each of the 
% volatility/noise values specified in h_all_sim and R_all_sim, in the estimation task. 
% The code can be adapted to find the optimal configurations of the other parametric models 
% by replacing the "funfmincon" function with that of the corresponding model and 
% setting the appropriate lower (lb) and upper (ub) bounds on parameter values.
% It can also be adapted to find the optimal model configurations in the
% prediction task, by a simple modification of the model mean squared error (MSE*) functions: 
% Prediction(t)=(1-h)*estimation(t) + h*mu0, MSE=mean((prediction(t)-process(t+1))^2)/(sigma0^2).
% Note that the optimal learning-rate values are the same in the estimation and prediction tasks (except for h=1), 
% and are typically estimated with higher accuracy using the estimation error functions.

N=1;
nb_restarts=10; 
nb_iterat=10;
h_all_sim=[0.02:0.02:1]; 
R_all_sim=[0,0.1:0.1:6]; 
sigma0=1;

A=[];
b=[];
Aeq=[];
beq=[];
lb=0*ones(N,1);
ub=1*ones(N,1);
nonlcon=[];
options = optimset('Algorithm','interior-point','Display', 'off');
warning('off','all')

for Rind=1:length(R_all_sim)    
    R=R_all_sim(1,Rind);
    sigma=R*sigma0;
    for hind=1:length(h_all_sim)
        h=h_all_sim(1,hind);
        for restarts=1:nb_restarts
            x_start= ub.*rand(N,1);
            for iterat=1:nb_iterat
                if(N==1)
                    funfmincon = @(x)MSE_Estim_DR_1(x,h,sigma,iterat);
                else
                    funfmincon = @(x)MSE_Estim_DR_N(x,h,sigma,iterat,N);
                end
                [xmin_f,f_min,exitflag,output,lambda,grad,hessian] = fmincon(funfmincon,x_start,A,b,Aeq,beq,lb,ub,nonlcon,options);
                MSErr{hind,Rind}(iterat,restarts)=f_min;
                [xmin,Ind]=sort(xmin_f,1,'ascend');
                Hes{hind,Rind}{iterat,restarts}=hessian(Ind,Ind);
                for nb_node=1:N
                    LearningRates{nb_node,1}{hind,Rind}(iterat,restarts)=xmin(nb_node,1);
                end
            end
        end
    end
end
save(['out_optim_estimation/LearningRates_DR_N',num2str(N)],'LearningRates')
save(['out_optim_estimation/MSErr_DR_N',num2str(N)],'MSErr')
save(['out_optim_estimation/Hes_DR_N',num2str(N)],'Hes')





