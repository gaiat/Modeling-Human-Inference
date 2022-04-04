
% This code computes the algorithmic complexity and inaccuracy of the
% models for different volatility/noise conditions, and the parameters of the
% linear fits of logComplexity versus logInaccuracy. 
% With task='prediction', relative_inaccuracy=1, NSP=0, KalmanHGF=0, it plots Fig4 and Supplementary Fig3A-C; it also computes complexity and inaccuracy used in Fig5A and Supplementary Fig4A-C, and fit parameters used in Fig6A.
% With task='estimation', relative_inaccuracy=1, NSP=0, KalmanHGF=0, it computes complexity and inaccuracy used in Fig5B and Supplementary Fig4G-I, and fit parameters used in Fig6B.
% With task='prediction', relative_inaccuracy=1, NSP=1, KalmanHGF=0, it plots Supplementary Figs2,3D-F.
% With task='prediction', relative_inaccuracy=0, NSP=0, KalmanHGF=0, it computes complexity and non-normalized inaccuracy used in Supplementary Fig4D-F.
% With task='estimation', relative_inaccuracy=0, NSP=0, KalmanHGF=0, it computes complexity and non-normalized inaccuracy used in Supplementary Fig4L-N.
% With KalmanHGF=1, it computes complexity and inaccuracy used in Supplementary Fig9A-C.

task='prediction';
relative_inaccuracy=1; % 1 --> Inaccuracy=(MSE(model)-MSE(Bayesian))/MSE(Bayesian); 0 --> Inaccuracy=(MSE(model)-MSE(Bayesian))/sigma0^2
NSP=0; % NSP=1 -> Neural-Spiking-P system implementation of the algorithms; NSP=0 -> Each algorithmic operation counts 1.
KalmanHGF=0; % 1 --> add Kalman filter and Hierarchical Gaussian filter to the model list; 0 -> include only models of Fig1.

if(KalmanHGF)
    task='estimation';
    relative_inaccuracy=1;
    NSP=0;
    Nlevels_HGF=2;
end
nb_restarts=10; 
nb_iterat=10;
sigma0=1;
nb_models=8; 
min_nb_fit_datapoints=2; 
Nmixt=2;
AC_action_pred=0.15; 
AC_action_estim=AC_action_pred; 
R_all_sim=[0,0.1:0.1:6];
h_all_sim=[0.02:0.02:1];
load(['out_optim_',task,'/Hes_DR_N1'],'Hes')
HesDR1=Hes;
load(['out_optim_',task,'/Hes_DR_N2'],'Hes')
HesDR2=Hes;
load(['out_optim_',task,'/Hes_SA_N1'],'Hes')
HesSA1=Hes;
load(['out_optim_',task,'/Hes_SA_N2'],'Hes')
HesSA2=Hes;
load(['out_optim_',task,'/LearningRates_SA_N1'],'LearningRates') 
LearningRateSA1=LearningRates;
load(['out_optim_',task,'/LearningRates_SA_N2'],'LearningRates')  
LearningRatesSA2=LearningRates;
load(['out_optim_',task,'/Hes_ML'],'Hes')
HesML=Hes;
load(['out_optim_',task,'/MSErr_DR_N1'],'MSErr')
MSErrDR1=MSErr;
load(['out_optim_',task,'/MSErr_DR_N2'],'MSErr')
MSErrDR2=MSErr;
load(['out_optim_',task,'/MSErr_SA_N1'],'MSErr')
MSErrSA1=MSErr;
load(['out_optim_',task,'/MSErr_SA_N2'],'MSErr')
MSErrSA2=MSErr;
load(['out_optim_',task,'/MSErr_ML'],'MSErr')
MSErrML=MSErr;
load(['out_optim_',task,'/MSErr_Prior'],'MSErr')
MSErrPrior=MSErr;
load(['out_optim_',task,'/MSErr_Evidence'],'MSErr')
MSErrEvidence=MSErr;
load(['out_optim_',task,'/MSErr_Bayes'],'MSErr')
MSErrBayes=MSErr;
if(KalmanHGF)
    load(['out_optim_',task,'/Hes_Kalman'],'Hes')
    HesKalman=Hes;
    load(['out_optim_',task,'/Hes_HGF'],'Hes')
    HesHGF=Hes;
    load(['out_optim_',task,'/MSErr_Kalman'],'MSErr')
    MSErrKalman=MSErr;
    load(['out_optim_',task,'/MSErr_HGF'],'MSErr')
    MSErrHGF=MSErr;
end
    

% Compute Inaccuracy, Accuracy and Complexity for all models

AC_Evidence_pred=0;
AC_Evidence_estim=AC_Evidence_pred;
Complexity_all_models=cell(length(h_all_sim),length(R_all_sim));
AccuracyLoss_all_models=cell(length(h_all_sim),length(R_all_sim));
Accuracy_all_models=cell(length(h_all_sim),length(R_all_sim));
for hind=1:length(h_all_sim)
    h=h_all_sim(1,hind);
    for Rind=1:length(R_all_sim)
        R=R_all_sim(1,Rind);
        if(R==0)
            if(strcmp(task,'prediction')) 
                % The model minimizing MSE is the Memoryless model: 
                MSE_Evidence=zeros(nb_iterat,1);
                MSE_Prior=zeros(nb_iterat,1);
                MSE_ML=zeros(nb_iterat,1);
                for iterat=1:nb_iterat
                    MSE_Evidence(iterat,1)=MSErrEvidence{hind,Rind}(iterat,1);
                    MSE_Prior(iterat,1)=MSErrPrior{hind,Rind}(iterat,1);
                    MSE_ML(iterat,1)=MSErrML{hind,Rind}(iterat,1);
                end
                Complexity_all_models{hind,Rind}=AC_action_pred+[AC_Evidence_pred;AC_Prior_pred();AC_ML_pred(NSP)];  
                if(relative_inaccuracy)
                    denom=mean(MSE_ML);
                else
                    denom=sigma0^2;
                end
                Loss_Evidence=(mean(MSE_Evidence)-mean(MSE_ML))/denom;
                Loss_Prior=(mean(MSE_Prior)-mean(MSE_ML))/denom;
                Loss_ML=0;
                Accuracy_Evidence=mean(MSE_ML)/mean(MSE_Evidence);
                Accuracy_Prior=mean(MSE_ML)/mean(MSE_Prior);
                Accuracy_ML=1;
                AccuracyLoss_all_models{hind,Rind}=[Loss_Evidence;Loss_Prior;Loss_ML];  
                Accuracy_all_models{hind,Rind}=[Accuracy_Evidence;Accuracy_Prior;Accuracy_ML];  
            elseif(strcmp(task,'estimation'))
                % The model minimizing MSE is the Evidence: 
                Complexity_all_models{hind,Rind}=AC_action_estim+AC_Evidence_estim;
                Loss_Evidence=0;
                Accuracy_Evidence=1;                
                AccuracyLoss_all_models{hind,Rind}=Loss_Evidence; 
                Accuracy_all_models{hind,Rind}=Accuracy_Evidence;
            end            
        else
            [~,indminMSEDR1]=min(MSErrDR1{hind,Rind},[],2);
            [~,indminMSEDR2]=min(MSErrDR2{hind,Rind},[],2);
            [~,indminMSESA1]=min(MSErrSA1{hind,Rind},[],2);
            [~,indminMSESA2]=min(MSErrSA2{hind,Rind},[],2); 
            [~,indminMSEML]=min(MSErrML{hind,Rind},[],2);
            LR_SA1=[];
            LR_SA2=[];
            MSE_Evidence=[];
            MSE_Prior=[];
            MSE_ML=[];
            MSE_DR1=[];
            MSE_DR2=[];
            MSE_SA1=[];
            MSE_SA2=[];
            MSE_Bayes=[];
            q=0;
            if(KalmanHGF)
                [~,indminMSEKalman]=min(MSErrKalman{hind,Rind},[],2);
                [~,indminMSEHGF]=min(MSErrHGF{hind,Rind},[],2);
                MSE_Kalman=[];
                MSE_HGF=[];
            end
            for iterat=1:nb_iterat
                hsDR1=HesDR1{hind,Rind}{iterat,indminMSEDR1(iterat)};
                hsDR2=HesDR2{hind,Rind}{iterat,indminMSEDR2(iterat)};
                hsSA1=HesSA1{hind,Rind}{iterat,indminMSESA1(iterat)};
                hsSA2=HesSA2{hind,Rind}{iterat,indminMSESA2(iterat)};
                hsML=HesML{hind,Rind}{iterat,indminMSEML(iterat)};
                [~,D_DR2] = eig(hsDR2); 
                [~,D_SA2] = eig(hsSA2);
                if(KalmanHGF)
                    hsKalman=HesKalman{hind,Rind}{iterat,indminMSEKalman(iterat)};
                    hsHGF=HesHGF{hind,Rind}{iterat,indminMSEHGF(iterat)};
                    [~,D_HGF] = eig(hsHGF);
                    min_reached=(~any(diag(D_DR2)<0))&&(~any(diag(D_SA2)<0))&&(hsDR1>0)&&(hsSA1>0)&&(hsML>0)&&(~any(diag(D_HGF)<0))&&(hsKalman>0); 
                else
                    min_reached=(~any(diag(D_DR2)<0))&&(~any(diag(D_SA2)<0))&&(hsDR1>0)&&(hsSA1>0)&&(hsML>0); 
                end                
                if(min_reached)
                    q=q+1;
                    LR_SA1(q,1)=LearningRateSA1{1,1}{hind,Rind}(iterat,indminMSESA1(iterat));
                    LR_SA2(q,1)=LearningRatesSA2{1,1}{hind,Rind}(iterat,indminMSESA2(iterat));
                    LR_SA2(q,2)=LearningRatesSA2{2,1}{hind,Rind}(iterat,indminMSESA2(iterat));
                    MSE_Evidence(q,1)=MSErrEvidence{hind,Rind}(iterat,1);
                    MSE_Prior(q,1)=MSErrPrior{hind,Rind}(iterat,1);
                    MSE_Bayes(q,1)=MSErrBayes{hind,Rind}(iterat,1);
                    MSE_ML(q,1)=MSErrML{hind,Rind}(iterat,indminMSEML(iterat));
                    MSE_DR1(q,1)=MSErrDR1{hind,Rind}(iterat,indminMSEDR1(iterat));
                    MSE_DR2(q,1)=MSErrDR2{hind,Rind}(iterat,indminMSEDR2(iterat));
                    MSE_SA1(q,1)=MSErrSA1{hind,Rind}(iterat,indminMSESA1(iterat));
                    MSE_SA2(q,1)=MSErrSA2{hind,Rind}(iterat,indminMSESA2(iterat));
                    if(KalmanHGF)
                        MSE_Kalman(q,1)=MSErrKalman{hind,Rind}(iterat,indminMSEKalman(iterat));
                        MSE_HGF(q,1)=MSErrHGF{hind,Rind}(iterat,indminMSEHGF(iterat));
                    end
                end                
            end            
            RunLength_SA1=(1/mean(LR_SA1))-R^2;
            RunLengths_SA2(1,1)=(1/mean(LR_SA2(:,1)))-R^2;
            RunLengths_SA2(2,1)=(1/mean(LR_SA2(:,2)))-R^2;
            if(strcmp(task,'prediction'))
                Complexity_all_models{hind,Rind}=AC_action_pred+[AC_Evidence_pred;AC_Prior_pred();AC_ML_pred(NSP);...
                    AC_DR_pred(NSP);AC_SA_pred(RunLength_SA1,NSP);AC_NDRs_pred(Nmixt,NSP);AC_NSAs_pred(RunLengths_SA2,Nmixt,NSP)];
            elseif(strcmp(task,'estimation'))
                if(KalmanHGF)
                    Complexity_all_models{hind,Rind}=AC_action_estim+[AC_Evidence_estim;AC_Prior_estim();AC_ML_estim(NSP);...
                        AC_DR_estim(NSP);AC_SA_estim(RunLength_SA1,NSP);AC_Kalman_estim();AC_HGF_estim(Nlevels_HGF);AC_NDRs_estim(Nmixt,NSP);AC_NSAs_estim(RunLengths_SA2,Nmixt,NSP)]; 
                else
                    Complexity_all_models{hind,Rind}=AC_action_estim+[AC_Evidence_estim;AC_Prior_estim();AC_ML_estim(NSP);...
                        AC_DR_estim(NSP);AC_SA_estim(RunLength_SA1,NSP);AC_NDRs_estim(Nmixt,NSP);AC_NSAs_estim(RunLengths_SA2,Nmixt,NSP)]; 
                end
            end
            if(relative_inaccuracy)
                denom=mean(MSE_Bayes);
            else
                denom=sigma0^2;
            end
           Loss_Evidence=(mean(MSE_Evidence)-mean(MSE_Bayes))/denom;
           Loss_Prior=(mean(MSE_Prior)-mean(MSE_Bayes))/denom;
           Loss_ML=(mean(MSE_ML)-mean(MSE_Bayes))/denom;
           Loss_DR1=(mean(MSE_DR1)-mean(MSE_Bayes))/denom;
           Loss_DR2=(mean(MSE_DR2)-mean(MSE_Bayes))/denom;
           Loss_SA1=(mean(MSE_SA1)-mean(MSE_Bayes))/denom;      
           Loss_SA2=(mean(MSE_SA2)-mean(MSE_Bayes))/denom;
           Accuracy_Evidence=mean(MSE_Bayes)/mean(MSE_Evidence);
           Accuracy_Prior=mean(MSE_Bayes)/mean(MSE_Prior);
           Accuracy_ML=mean(MSE_Bayes)/mean(MSE_ML);
           Accuracy_DR1=mean(MSE_Bayes)/mean(MSE_DR1);
           Accuracy_DR2=mean(MSE_Bayes)/mean(MSE_DR2);
           Accuracy_SA1=mean(MSE_Bayes)/mean(MSE_SA1);       
           Accuracy_SA2=mean(MSE_Bayes)/mean(MSE_SA2);               
           if(KalmanHGF)
               Loss_Kalman=(mean(MSE_Kalman)-mean(MSE_Bayes))/denom;
               Loss_HGF=(mean(MSE_HGF)-mean(MSE_Bayes))/denom;
               Accuracy_Kalman=mean(MSE_Bayes)/mean(MSE_Kalman);
               Accuracy_HGF=mean(MSE_Bayes)/mean(MSE_HGF);
               AccuracyLoss_all_models{hind,Rind}=[Loss_Evidence;Loss_Prior;Loss_ML;Loss_DR1;Loss_SA1;Loss_Kalman;Loss_HGF;Loss_DR2;Loss_SA2];      
               Accuracy_all_models{hind,Rind}=[Accuracy_Evidence;Accuracy_Prior;Accuracy_ML;Accuracy_DR1;Accuracy_SA1;Accuracy_Kalman;Accuracy_HGF;Accuracy_DR2;Accuracy_SA2]; 
           else
               AccuracyLoss_all_models{hind,Rind}=[Loss_Evidence;Loss_Prior;Loss_ML;Loss_DR1;Loss_SA1;Loss_DR2;Loss_SA2]; 
               Accuracy_all_models{hind,Rind}=[Accuracy_Evidence;Accuracy_Prior;Accuracy_ML;Accuracy_DR1;Accuracy_SA1;Accuracy_DR2;Accuracy_SA2]; 
           end
        end
    end  
end
                
% Linear fits in log-log scale and goodness-of-fit statistics
R_acceptableFits=[0.1:0.1:6]; 
h_acceptableFits=[0.02:0.02:0.86]; 
[~,h_sel,~]=intersect(round(100*h_all_sim),round(100*h_acceptableFits));
[~,R_sel,~]=intersect(round(10*R_all_sim),round(10*R_acceptableFits));
first_zero=zeros(length(h_all_sim),length(R_all_sim));
b_all_models=NaN(length(h_all_sim),length(R_all_sim));
loga_all_models=NaN(length(h_all_sim),length(R_all_sim));
pValueModelvsConst=NaN(length(h_all_sim),length(R_all_sim));
Resid_all=[];
Loss_fit_all=[];
for hind=1:length(h_all_sim)
    for Rind=1:length(R_all_sim)
        models0=find(AccuracyLoss_all_models{hind,Rind}<=0);
        if(isempty(models0))
            first_zero(hind,Rind)=nb_models;
        else
            first_zero(hind,Rind)=min(models0);
        end
        if(first_zero(hind,Rind)>min_nb_fit_datapoints) 
            if(ismember(hind,h_sel)&&ismember(Rind,R_sel))
                x=Complexity_all_models{hind,Rind}([1:first_zero(hind,Rind)-1],1);
                lm = fitlm(log10(x),log10(AccuracyLoss_all_models{hind,Rind}([1:first_zero(hind,Rind)-1],1)),'linear');
                b_all_models(hind,Rind)=-lm.Coefficients.Estimate(2,1);
                loga_all_models(hind,Rind)=lm.Coefficients.Estimate(1,1);
                % Goodness-of-fit statistics:
                ANOV=anova(lm,'summary');
                pValueModelvsConst(hind,Rind)=table2array(ANOV(2,5)); 
                Resid_all=[Resid_all;lm.Residuals.Raw];
                Loss_fit_all=[Loss_fit_all;10.^(loga_all_models(hind,Rind)-b_all_models(hind,Rind)*log10(x))];
            end
        end
    end       
end

% Conditions to be excluded:
nofit=[];
posFitSlope=[];
for hind=1:length(h_all_sim)
    for Rind=1:length(R_all_sim)
        if(first_zero(hind,Rind)<=min_nb_fit_datapoints)
            nofit=[nofit;hind,Rind,h_all_sim(1,hind),R_all_sim(1,Rind)];
        else
            if(b_all_models(hind,Rind)<0)
                posFitSlope=[posFitSlope;hind,Rind,h_all_sim(1,hind),R_all_sim(1,Rind)];
            end
        end
    end
end
if((isempty(nofit))&&(isempty(posFitSlope)))
    disp('At least 2 models have Loss>0 and the fit slope is negative in all tested conditions')
else
    if(~isempty(nofit))
        disp('No fit possible: volatility, noise')
        nofit(:,[3,4])
    end
    if(~isempty(posFitSlope))
        disp('Positive fit slope: volatility, noise')
        posFitSlope(:,[3,4])
    end
end

if(~NSP)
    if(relative_inaccuracy)
        if(KalmanHGF)
            save(['DiminReturns_',task,'/Complexity_all_models_KalmanHGF'],'Complexity_all_models')
            save(['DiminReturns_',task,'/AccuracyLoss_all_models_KalmanHGF'],'AccuracyLoss_all_models') 
        else
            save(['DiminReturns_',task,'/Complexity_all_models'],'Complexity_all_models')
            save(['DiminReturns_',task,'/AccuracyLoss_all_models'],'AccuracyLoss_all_models') 
            save(['DiminReturns_',task,'/b_all_models'],'b_all_models')
            save(['DiminReturns_',task,'/loga_all_models'],'loga_all_models')
        end
    else
        save(['DiminReturns_',task,'/Complexity_all_models_abs'],'Complexity_all_models')
        save(['DiminReturns_',task,'/AccuracyLoss_all_models_abs'],'AccuracyLoss_all_models') 
    end
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if((relative_inaccuracy)&&strcmp(task,'prediction'))
    
% Fig4:

h_examples=0.1;
R_examples=1;
h_examples_multR=[0.1;0.1;0.1];
R_examples_multR=[0.2;1;3];
h_examples_multh=[0.06;0.3;0.8];
R_examples_multh=[1;1;1];
tol_I=0.1;
tol_A=1/(tol_I+1);
minCplot=0.1;
if(NSP)
    maxCplot=350;
else
    maxCplot=200;
end
myMap=[0.68,0.92,1;0,0,0; 51/255,51/255,1; 0,204/255,0; 1,0,0; 1,160/255,1; 1,0.97,0; 0.8,0.8,0.8]; 
Colors_fits_multR={[183/255,232/255,255/255];[0/255,153/255,255/255];[0/255,0/255,153/255]};
Colors_fits_multh={[255/255,153/255,0];[255/255,51/255,0];[130/255,51/255,0/255]};
Markers_models={'^';'v';'s';'d';'p';'h';'o';'*'};
pos_title=[0.805 186.061 0];
figure1=figure;
hinde=find(round(h_all_sim*100)==round(100*h_examples));
Rinde=find(round(R_all_sim*100)==round(100*R_examples));
hinde_multR=zeros(length(h_examples_multR),1);
Rinde_multR=zeros(length(h_examples_multR),1);
for k=1:length(h_examples_multR)
    hinde_multR(k,1)=find(round(h_all_sim*100)==round(100*h_examples_multR(k,1)));
    Rinde_multR(k,1)=find(round(R_all_sim*100)==round(100*R_examples_multR(k,1)));
end
hinde_multh=zeros(length(h_examples_multh),1);
Rinde_multh=zeros(length(h_examples_multh),1);
for k=1:length(h_examples_multh)
    hinde_multh(k,1)=find(round(h_all_sim*100)==round(100*h_examples_multh(k,1)));
    Rinde_multh(k,1)=find(round(R_all_sim*100)==round(100*R_examples_multh(k,1)));
end

% Fig4A

axes1 = axes('Position',[0.12007874015748 0.768282351025766 0.290682414698163 0.189053148568057]);
hold(axes1,'on');
hold on
size=150;
pp=zeros(length(Complexity_all_models{hinde,Rinde}),1);
for i=1:length(Complexity_all_models{hinde,Rinde})
     pp(i)=scatter(i,Complexity_all_models{hinde,Rinde}(i,1),size,Markers_models{i,1},'MarkerEdgeColor','k');
end
plot([1:length(Complexity_all_models{hinde,Rinde})],Complexity_all_models{hinde,Rinde},'-','Color','k','LineWidth',2)
hold off
xlabel('Model ranking','FontSize',22)
ylabel('Complexity','FontSize',22)
xlim(axes1,[0.5 7.5]);
ylim(axes1,[0 maxCplot]);
title('A','FontSize',22)
legend([pp(1),pp(2),pp(3),pp(4),pp(5),pp(6),pp(7)],{'Evidence','Prior Model','Memoryless Model','Delta Rule','Sliding Window','Mixture of Delta Rules','Mixture of Sliding Windows'},...
     'Position',[0.12012951191162 0.826514338950981 0.181933837784004 0.130787973525948],'FontSize',13);
box(axes1,'on');
set(axes1,'FontSize',18,'XTick',[1 2 3 4 5 6 7]);

% Fig4B

axes2 = axes('Position', [0.484438810018009 0.766433793663686 0.278028381583041 0.18922851377866]);
hold(axes2,'on');
hold on
x=[minCplot:(maxCplot-minCplot)/1000:maxCplot];
for i=1:first_zero(hinde,Rinde)-1
     plot(Complexity_all_models{hinde,Rinde}(i,1),AccuracyLoss_all_models{hinde,Rinde}(i,1),'Marker',Markers_models{i,1},'MarkerEdgeColor','k','MarkerSize',12)
end
plot(x,10.^(loga_all_models(hinde,Rinde)-b_all_models(hinde,Rinde)*log10(x)),'-','Color','k','LineWidth',2)
hold off
xlim(axes2,[-3 maxCplot]);
ylim(axes2,[0 2.5]);
xlabel('Complexity','FontSize',22)
ylabel('Inaccuracy','FontSize',22)
title('B','FontSize',22)
box(axes2,'on');
if(NSP)
    set(axes2,'FontSize',18,'XTick',...
    [0 100 200 300],'XTickLabel',{'0','100','200','300'},...
    'YTick',[0 1 2],'YTickLabel',{'0','1','2'});
else
    set(axes2,'FontSize',18,'XTick',...
    [0 50 100 150 200],'XTickLabel',{'0','50','100','150','200'},...
    'YTick',[0 1 2],'YTickLabel',{'0','1','2'});
end
% Inset
axes3 = axes('Position', [0.595800524934383 0.845150284321687 0.166666666666667 0.10978521608637]);
hold(axes3,'on');
hold on
for i=1:first_zero(hinde,Rinde)-1
     loglog(Complexity_all_models{hinde,Rinde}(i,1),AccuracyLoss_all_models{hinde,Rinde}(i,1),'Marker',Markers_models{i,1},'MarkerEdgeColor','k','MarkerSize',9)
end
loglog(x,10.^(loga_all_models(hinde,Rinde)-b_all_models(hinde,Rinde)*log10(x)),'-','Color','k','LineWidth',1)
hold off
xlim(axes3,[0.05 1000]);
ylim(axes3,[0.01 15]);
box(axes3,'on');
set(axes3,'FontSize',10,'XMinorTick','on','XScale','log','XTick',...
    [0.1 1 10 100],'XTickLabel',{'10^{-1}','10^{0}','10^{1}','10^{2}'},...
    'YMinorTick','on','YScale','log');

% Fig4C

axes4 = axes('Position',[0.118110236220472 0.472786352558894 0.291994750656168 0.198212835093421]);
hold(axes4,'on');
pp=zeros(length(h_examples_multR),1);
x=[minCplot:(maxCplot-minCplot)/1000:maxCplot];
hold on
for k=1:length(h_examples_multR)
    for i=1:first_zero(hinde_multR(k,1),Rinde_multR(k,1))-1
        loglog(Complexity_all_models{hinde_multR(k,1),Rinde_multR(k,1)}(i,1),AccuracyLoss_all_models{hinde_multR(k,1),Rinde_multR(k,1)}(i,1),'Marker',Markers_models{i,1},'MarkerEdgeColor',Colors_fits_multR{k,1},'MarkerFaceColor',Colors_fits_multR{k,1},'MarkerSize',7)
    end
    pp(k)=loglog(x,10.^(loga_all_models(hinde_multR(k,1),Rinde_multR(k,1))-b_all_models(hinde_multR(k,1),Rinde_multR(k,1))*log10(x)),'-','Color',Colors_fits_multR{k,1},'LineWidth',2);
end
hold off
xlabel('Complexity','FontSize',22)
ylabel('Inaccuracy','FontSize',22)
title('C','FontSize',22)
xlim(axes4,[0.05 maxCplot+300]);
ylim(axes4,[10^-3 10^2]);
box(axes4,'on');
set(axes4,'FontSize',18,'XMinorTick','on','XScale','log','XTick',[0.1 1 10 100],'XTickLabel',{'10^{-1}','10^{0}','10^{1}','10^{2}'},'YMinorTick','on','YScale','log','YTick',[0.01 1 100],'YTickLabel',{'10^{-2}','10^{0}','10^{2}'});

% Fig4D

axes5 = axes('Position', [0.482283464566929 0.472786352558895 0.284776902887139 0.198387978068495]);
hold(axes5,'on');
hold on
for k=1:length(h_examples_multh)
    for i=1:first_zero(hinde_multh(k,1),Rinde_multh(k,1))-1
        loglog(Complexity_all_models{hinde_multh(k,1),Rinde_multh(k,1)}(i,1),AccuracyLoss_all_models{hinde_multh(k,1),Rinde_multh(k,1)}(i,1),'Marker',Markers_models{i,1},'MarkerEdgeColor',Colors_fits_multh{k,1},'MarkerFaceColor',Colors_fits_multh{k,1},'MarkerSize',7)
    end
    loglog(x,10.^(loga_all_models(hinde_multh(k,1),Rinde_multh(k,1))-b_all_models(hinde_multh(k,1),Rinde_multh(k,1))*log10(x)),'-','Color',Colors_fits_multh{k,1},'LineWidth',2)
end
hold off
xlabel('Complexity','FontSize',22)
ylabel('Inaccuracy','FontSize',22)
title('D','FontSize',22)
xlim(axes5,[0.05 maxCplot+300]);
ylim(axes5,[3*10^-6 5*10^2]);
box(axes5,'on');
set(axes5,'FontSize',18,'XMinorTick','on','XScale','log','XTick',[0.1 1 10 100],'XTickLabel',{'10^{-1}','10^{0}','10^{1}','10^{2}'},'YMinorTick','on','YScale','log','YTick',[0.0001 1],'YTickLabel',{'10^{-4}','10^{0}'});

% Fig4E, top panel

axes6 = axes('Position',[0.118110236220472 0.259285134037363 0.293963254593176 0.114541023558083]);
hold(axes6,'on');
pp=zeros(length(h_examples_multR),1);
hold on
for k=1:length(h_examples_multR)
    pp(k)=plot(x,10.^(loga_all_models(hinde_multR(k,1),Rinde_multR(k,1))-b_all_models(hinde_multR(k,1),Rinde_multR(k,1))*log10(x)),'-','Color',Colors_fits_multR{k,1},'LineWidth',2);
    a=10.^(loga_all_models(hinde_multR(k,1),Rinde_multR(k,1)));
    b=b_all_models(hinde_multR(k,1),Rinde_multR(k,1));
    C_intercept=(a/tol_I)^(1/b);
    plot(C_intercept,a*C_intercept^(-b),'*','MarkerEdgeColor',Colors_fits_multR{k,1},'MarkerFaceColor',Colors_fits_multR{k,1},'MarkerSize',10,'LineWidth',3)
    plot([min(x),max(x)],[tol_I,tol_I],'-k','LineWidth',0.5)
end
hold off
ylabel('Inaccuracy','FontSize',22)
xlim(axes6,[0 maxCplot]);
ylim(axes6,[0 0.25]);
title('E','FontSize',22)
legend([pp(1) pp(2) pp(3)],{['h=',num2str(h_examples_multR(1,1)),'; R=',num2str(R_examples_multR(1,1))],['h=',num2str(h_examples_multR(2,1)),'; R=',num2str(R_examples_multR(2,1))],['h=',num2str(h_examples_multR(3,1)),'; R=',num2str(R_examples_multR(3,1))]},...
    'Position',[0.302416195247556 0.613504020529618 0.108142491382649 0.057676684023596],'FontSize',13);
box(axes6,'on');
if(NSP)
    set(axes6,'FontSize',18,'XTick',[0 100 200 300],'XTickLabel',{'' '' '' ''},'YTick',[0 tol_I],'YTickLabel',{'0',num2str(tol_I)});
else
    set(axes6,'FontSize',18,'XTick',[0 50 100 150 200],'XTickLabel',{'' '' '' ''},'YTick',[0 tol_I],'YTickLabel',{'0',num2str(tol_I)});
end

% Fig4F, top panel

axes7 = axes('Position', [0.47244094488189 0.258326563769293 0.293963254593176 0.114541023558083]);
hold(axes7,'on');
pp=zeros(length(h_examples_multh),1);
hold on
for k=1:length(h_examples_multh)
    pp(k)=plot(x,10.^(loga_all_models(hinde_multh(k,1),Rinde_multh(k,1))-b_all_models(hinde_multh(k,1),Rinde_multh(k,1))*log10(x)),'-','Color',Colors_fits_multh{k,1},'LineWidth',2);
    a=10.^(loga_all_models(hinde_multh(k,1),Rinde_multh(k,1)));
    b=b_all_models(hinde_multh(k,1),Rinde_multh(k,1));
    C_intercept=(a/tol_I)^(1/b);
    plot(C_intercept,a*C_intercept^(-b),'*','MarkerEdgeColor',Colors_fits_multh{k,1},'MarkerFaceColor',Colors_fits_multh{k,1},'MarkerSize',10,'LineWidth',3)
    plot([min(x),max(x)],[tol_I,tol_I],'-k','LineWidth',0.5)
end
hold off
ylabel('Inaccuracy','FontSize',22)
xlim(axes7,[0 maxCplot]);
ylim(axes7,[0 0.25]);
title('F','FontSize',22)
legend([pp(1) pp(2) pp(3)],{['h=',num2str(h_examples_multh(1,1)),'; R=',num2str(R_examples_multh(1,1))],['h=',num2str(h_examples_multh(2,1)),'; R=',num2str(R_examples_multh(2,1))],['h=',num2str(h_examples_multh(3,1)),'; R=',num2str(R_examples_multh(3,1))]},...
  'Position',[0.663079875077088 0.613002724392741 0.104325697603268 0.057676684023596],'FontSize',13);
box(axes7,'on');
if(NSP)
    set(axes7,'FontSize',18,'XTick',[0 100 200 300],'XTickLabel',{'' '' '' ''},'YTick',[0 tol_I],'YTickLabel',{'0',num2str(tol_I)});
else
    set(axes7,'FontSize',18,'XTick',[0 50 100 150 200],'XTickLabel',{'' '' '' ''},'YTick',[0 tol_I],'YTickLabel',{'0',num2str(tol_I)});
end

% Fig4E, bottom panel

axes8 = axes('Position',[0.11745406824147 0.115353371242892 0.293963254593176 0.114541023558083]);
hold(axes8,'on');
hold on
step=10^(-4);
pp=zeros(length(h_examples_multR),1);
for k=1:length(h_examples_multR)
    x=[step:step:maxCplot];
    AA=1./(1+10.^(loga_all_models(hinde_multR(k,1),Rinde_multR(k,1))-b_all_models(hinde_multR(k,1),Rinde_multR(k,1))*log10(x)));
    x=[0,x];
    AA=[0,AA];
    pp(k)=plot(x,AA,'-','Color',Colors_fits_multR{k,1},'LineWidth',2);
    a=10.^(loga_all_models(hinde_multR(k,1),Rinde_multR(k,1)));
    b=b_all_models(hinde_multR(k,1),Rinde_multR(k,1));
    C_intercept=((a*tol_A)/(1-tol_A))^(1/b);
    plot(C_intercept,1/(1+a*C_intercept^(-b)),'*','MarkerEdgeColor',Colors_fits_multR{k,1},'MarkerFaceColor',Colors_fits_multR{k,1},'MarkerSize',10,'LineWidth',3)
    plot([0,max(x)],[tol_A,tol_A],'-k','LineWidth',0.5)
end
hold off
xlabel('Complexity','FontSize',22)
ylabel('Accuracy','FontSize',22)
xlim(axes8,[0 maxCplot]);
ylim(axes8,[0.7 1]);
box(axes8,'on');
if(NSP)
    set(axes8,'FontSize',18,'XTick',[0 100 200 300],'XTickLabel',{'0' '100' '200' '300'},'YTick',[0.7 tol_A 1],'YTickLabel',{'0.7',num2str(tol_A),'1'});
else
    set(axes8,'FontSize',18,'XTick',[0 50 100 150 200],'XTickLabel',{'0' '50' '100' '150' '200'},'YTick',[0.7 tol_A 1],'YTickLabel',{'0.7',num2str(tol_A),'1'});
end  
 
% Fig4F, bottom panel

axes9 = axes('Position',[0.4750656167979 0.112916328188465 0.293963254593176 0.114541023558083]);
hold(axes9,'on');
pp=zeros(length(h_examples_multh),1);
hold on
for k=1:length(h_examples_multh)
    x=[step:step:maxCplot];
    AA=1./(1+10.^(loga_all_models(hinde_multh(k,1),Rinde_multh(k,1))-b_all_models(hinde_multh(k,1),Rinde_multh(k,1))*log10(x)));
    x=[0,x];
    AA=[0,AA];
    pp(k)=plot(x,AA,'-','Color',Colors_fits_multh{k,1},'LineWidth',2);
    a=10.^(loga_all_models(hinde_multh(k,1),Rinde_multh(k,1)));
    b=b_all_models(hinde_multh(k,1),Rinde_multh(k,1));
    C_intercept=((a*tol_A)/(1-tol_A))^(1/b);
    plot(C_intercept,1/(1+a*C_intercept^(-b)),'*','MarkerEdgeColor',Colors_fits_multh{k,1},'MarkerFaceColor',Colors_fits_multh{k,1},'MarkerSize',10,'LineWidth',3)
    plot([0,max(x)],[tol_A,tol_A],'-k','LineWidth',0.5)
end
hold off
xlabel('Complexity','FontSize',22)
ylabel('Accuracy','FontSize',22)
xlim(axes9,[0 maxCplot]);
ylim(axes9,[0.7 1]);
box(axes9,'on');
if(NSP)
    set(axes9,'FontSize',18,'XTick',[0 100 200 300],'XTickLabel',{'0' '100' '200' '300'},'YTick',[0.7 tol_A 1],'YTickLabel',{'0.7',num2str(tol_A),'1'});
else
    set(axes9,'FontSize',18,'XTick',[0 50 100 150 200],'XTickLabel',{'0' '50' '100' '150' '200'},'YTick',[0.7 tol_A 1],'YTickLabel',{'0.7',num2str(tol_A),'1'});
end

% Supplementary Fig3

figure2=figure;
axes6 = axes('Position',[0.100738740631241 0.5466632191173 0.200290938351194 0.196634912483025]);
hold(axes6,'on');
xlim([-2,2])
ylim([0 1600])
histogram(Resid_all);
xlabel('Residuals','FontSize',22)
box(axes6,'on');
set(axes6,'FontSize',18);
title('A','FontSize',22)
axes7 = axes('Position',[0.369269883663103 0.546709991876523 0.200687717790559 0.196438587216734]);
hold(axes7,'on');
semilogx(Loss_fit_all,Resid_all,'MarkerFaceColor','b','Marker','.','LineStyle','none','Color','b');
ylabel('Residuals','FontSize',22)
xlabel('Fitted values','FontSize',22)
xlim(axes7,[0 100]);
ylim(axes7,[-4 4]);
box(axes7,'on');
set(axes7,'FontSize',18);
set(axes7,'XMinorTick','on','XScale','log','XTick',[0.0001 0.01 1 100],'XTickLabel',{'10^{-4}','10^{-2}','10^{0}','10^{2}'});
title('B','FontSize',22)
axes8 = axes('Position',[0.648797747399609 0.546709991876523 0.200990259868711 0.197019441250235]);
hold(axes8,'on');
histogram(pValueModelvsConst(~isnan(pValueModelvsConst)),'Normalization','cdf');
xlabel('p-values','FontSize',22)
ylabel('Cumulative pdf','FontSize',22)
xlim([-0.01,0.25]);
ylim([0,1.05]);
percentage_10percentSignifLevel=100*(length(find(pValueModelvsConst(:)<=0.1))/length(pValueModelvsConst(~isnan(pValueModelvsConst))));
disp(percentage_10percentSignifLevel);
box(axes8,'on')
set(axes8,'FontSize',18);
title('C','FontSize',22)

end


