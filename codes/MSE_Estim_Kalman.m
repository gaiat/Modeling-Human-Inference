
% This function computes the mean squared estimation error for the Kalman filter with process noise x, 
% over a Ttot=5000 time-long instance (iterat) of a Gaussian change-point process 
% with volatility h and noise sigma (=std of the observations around the source).

function f=MSE_Estim_Kalman(x,h,sigma,iterat) 

warning('off','all')
sigma0=1;
mu0=0;
Ttot=5*10^3;

load(['samples/process_s',num2str(sigma*100),'_s0',num2str(sigma0*100),'_h',num2str(h*1000),'_it',num2str(iterat)],'process') 
A=1;
H=1;
B=[];

kalmanFilter = vision.KalmanFilter(A,H,B,'MeasurementNoise',sigma^2,'ProcessNoise',x,'State',mu0,'StateCovariance',sigma0^2);
mu_pred=zeros(Ttot,1);
for tp=1:Ttot
    predict(kalmanFilter);
    [~,mu_pred(tp,1),~]= correct(kalmanFilter,process(tp,1));
end 

f = (mean((mu_pred-process([1:Ttot],2)).^2))/(sigma0^2); 


end


