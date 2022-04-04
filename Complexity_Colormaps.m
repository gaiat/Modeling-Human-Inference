
% With paper_fig=1, it plots Fig6.
% With paper_fig=2, it plots Supplementary Fig5.

paper_fig=2; 

sigma0=1;
if(paper_fig==1)
    reward_funct={'Gaussian'}; 
    fitness_funct={'Ratio'}; 
elseif(paper_fig==2)
    reward_funct={'Gaussian','Exponential','Linear'}; 
    fitness_funct={'Ratio','Difference'}; 
end
sigmaR=[0.1];
w1=1;
w2=0.05;
time_action_pred_fit=0.15;
time_action_pred=0.15; 
R_all_sim=[0,0.1:0.1:6];
h_all_sim=[0.02:0.02:1]; 
h_maxplot=0.86; 
titles={'A','B','C','D','E','F'};
if(paper_fig==1)
    pos{1,1}= [0.22489579257636 0.694486889658082 0.205707019328586 0.234026514078716];
    pos{2,1}=[0.524164935610611 0.691307879772521 0.206119023397757 0.235347342232337];
else
    pos{1,1}=[0.17647544114277 0.726196376852794 0.272382552898938 0.192443482535637];
    pos{2,1}=[0.536246276067527 0.720977939544457 0.270109235352532 0.192022490272563];
    pos{3,1}=[0.177996376083078 0.446764091858026 0.273840763936783 0.192202514000352];
    pos{4,1}=[0.535253227408143 0.442095588235294 0.270109235352532 0.194294615691958];
    pos{5,1}=[0.177996376083078 0.165514091858026 0.271854666618014 0.192202514000352];
    pos{6,1}=[0.533267130089374 0.16447025052192 0.269116186693148 0.191589071052391];
end

figure1=figure;
if(paper_fig==1)
    fig_panels={{'prediction','rel',sigmaR(1,1),fitness_funct{1,1},reward_funct{1,1}};{'estimation','rel',sigmaR(1,1),fitness_funct{1,1},reward_funct{1,1}}};
else
   fig_panels={{'prediction','rel',sigmaR(1,1),fitness_funct{1,1},reward_funct{1,2}};{'estimation','rel',sigmaR(1,1),fitness_funct{1,1},reward_funct{1,2}};...
       {'prediction','rel',sigmaR(1,1),fitness_funct{1,1},reward_funct{1,3}};{'estimation','rel',sigmaR(1,1),fitness_funct{1,1},reward_funct{1,3}};...
         {'prediction','rel',sigmaR(1,1),fitness_funct{1,2},reward_funct{1,3}};{'estimation','rel',sigmaR(1,1),fitness_funct{1,2},reward_funct{1,3}}};
end

for p=1:size(fig_panels,1)
    load(['DiminReturns_',fig_panels{p,1}{1,1},'/b_all_models'],'b_all_models') 
    load(['DiminReturns_',fig_panels{p,1}{1,1},'/loga_all_models'],'loga_all_models')   
    loga_all_models=loga_all_models([1:length(h_all_sim)],[1:length(R_all_sim)]);
    b_all_models=b_all_models([1:length(h_all_sim)],[1:length(R_all_sim)]);
    a_all_models=10.^loga_all_models;
    if(strcmp(fig_panels{p,1}{1,4},'Ratio'))
        if(strcmp(fig_panels{p,1}{1,5},'Gaussian'))
            Copt=(time_action_pred-time_action_pred_fit) + ((a_all_models.*sqrt(b_all_models))./(fig_panels{p,1}{1,3})).^(1./b_all_models); % Copt(hind,Rind)
        elseif(strcmp(fig_panels{p,1}{1,5},'Exponential'))
            Copt=(time_action_pred-time_action_pred_fit) + ((a_all_models.*b_all_models)./(fig_panels{p,1}{1,3})).^(1./b_all_models); % Copt(hind,Rind)
        elseif(strcmp(fig_panels{p,1}{1,5},'Linear'))
            Copt=(time_action_pred-time_action_pred_fit) + ((a_all_models.*(b_all_models+1))./(sqrt(2*fig_panels{p,1}{1,3}))).^(1./b_all_models); % Copt(hind,Rind)
        end
    elseif(strcmp(fig_panels{p,1}{1,4},'Difference'))
        Cost_all=[time_action_pred:0.01:350];
        if(strcmp(fig_panels{p,1}{1,5},'Gaussian'))
            Copt=zeros(length(h_all_sim),length(R_all_sim));
            for h=1:length(h_all_sim)
                for R=1:length(R_all_sim)
                    Reward=(2./(fig_panels{p,1}{1,3}*sqrt(2*pi))).*exp(-(((a_all_models(h,R).*Cost_all.^(-b_all_models(h,R))).^2)./(2*(fig_panels{p,1}{1,3}^2))));
                    Fitness=w1*Reward-w2*Cost_all;
                    [~,Copt_ind]=max(Fitness);
                    Copt(h,R)=(time_action_pred-time_action_pred_fit) + Cost_all(1,Copt_ind);
                end
            end
        elseif(strcmp(fig_panels{p,1}{1,5},'Exponential'))
            Copt=zeros(length(h_all_sim),length(R_all_sim));
            for h=1:length(h_all_sim)
                for R=1:length(R_all_sim)
                    Reward=(1./fig_panels{p,1}{1,3}).*exp(-(a_all_models(h,R).*Cost_all.^(-b_all_models(h,R))./fig_panels{p,1}{1,3}));
                    Fitness=w1*Reward-w2*Cost_all;
                    [~,Copt_ind]=max(Fitness);
                    Copt(h,R)=(time_action_pred-time_action_pred_fit) + Cost_all(1,Copt_ind);
                end
            end              
        elseif(strcmp(fig_panels{p,1}{1,5},'Linear'))
            Copt=(time_action_pred-time_action_pred_fit) + (a_all_models.^(-1).*b_all_models.^(-1).*fig_panels{p,1}{1,3}.*w1.^(-1).*w2).^(((-1)+(-1).*b_all_models).^(-1)); 
        end       
    end
    if(strcmp(fig_panels{p,1}{1,1},'estimation'))
        Copt(:,1)=time_action_pred*ones(length(h_all_sim),1); % when R=0 in the estimation task, the Evidence is the model with min inaccuracy (=0, i.e. max reward) and min complexity (=time_action_pred)
    end
    Copt(Copt<time_action_pred)=time_action_pred;
    
    if(strcmp(fig_panels{p,1}{1,1},'prediction'))
        subploth{p}=subplot('Position',pos{p,1},'FontSize',18,'XTick',[0.02 0.2 0.4 0.6 0.8],'XTickLabel',{'0.02','0.2','0.4','0.6','0.8'},'YTick',[0.1 2 4 6],'YTickLabel',{'0.1','2','4','6'});
    else
        subploth{p}=subplot('Position',pos{p,1},'FontSize',18,'XTick',[0.02 0.2 0.4 0.6 0.8],'XTickLabel',{'0.02','0.2','0.4','0.6','0.8'},'YTick',[0 2 4 6],'YTickLabel',{'0','2','4','6'});
    end
    hold(subploth{p},'on'); 
    surf(h_all_sim,R_all_sim,log10(Copt'), 'LineStyle', 'none');
    view(0,90);
    grid off
    colormap('parula')
    xlabel('Volatility','FontSize',22);
    ylabel('Noise','FontSize',22);
    caxis([log10(time_action_pred),2.5]);
    set(gca,'YDir','normal')
    if(strcmp(fig_panels{p,1}{1,1},'prediction'))
        ylim([R_all_sim(1,2),R_all_sim(1,end)])
    else
        ylim([R_all_sim(1,1),R_all_sim(1,end)])
    end
    xlim([h_all_sim(1,1),h_maxplot])
    if(paper_fig==1)
        if(p==1)
            colorbar('southoutside','Position',[0.222446008140986 0.575954508529651 0.506891531291191 0.0250744661039213],'FontSize',18,...
                'Ticks',[-0.5 0 0.5 1 1.5 2 2.5],'TickLabels',{'-0.5','0','0.5','1','1.5','2','\geq 2.5'});
        end
    else
        if(p==1)
            colorbar('southoutside','Position',[0.174779672490539 0.0450367647058824 0.623631449654447 0.0273451741903312],'FontSize',18,...
                'Ticks',[-0.5 0 0.5 1 1.5 2 2.5],'TickLabels',{'-0.5','0','0.5','1','1.5','2','\geq 2.5'});
        end
    end 
end

if(paper_fig==1)
    annotation(figure1,'textbox',[0.270031545741326 0.934824533747419 0.133123024611819 0.0365556449690648],...
        'String',{'Prediction'},'LineStyle','none','FontWeight','bold','FontSize',22);
    annotation(figure1,'textbox',[0.570347003154578 0.931575143008169 0.137539428265689 0.0365556449690648],...
        'String',{'Estimation'},'LineStyle','none','FontWeight','bold','FontSize',22);
    annotation(figure1,'textbox',[0.220189274447951 0.935636881529059 0.0447949515344217 0.0398050356114618],...
        'String',{'A'},'LineStyle','none','FontWeight','bold','FontSize',25);
    annotation(figure1,'textbox',[0.519242902208216 0.935636881529059 0.0447949515344217 0.0398050356114618],...
        'String',{'B'},'LineStyle','none','FontWeight','bold','FontSize',25);
else
    annotation(figure1,'textbox',[0.173515987456888 0.919092763881999 0.0447949515344217 0.0398050356114625],...
        'String',{'A'},'LineStyle','none','FontWeight','bold','FontSize',25,'FitBoxToText','off');
    annotation(figure1,'textbox',[0.532152534780212 0.913578057999646 0.0447949515344217 0.0398050356114625],...
        'String',{'B'},'LineStyle','none','FontWeight','bold','FontSize',25,'FitBoxToText','off');
    annotation(figure1,'textbox',[0.175648066061244 0.639680999176113 0.0447949515344216 0.0398050356114624],...
        'String','C','LineStyle','none','FontWeight','bold','FontSize',25,'FitBoxToText','off');
    annotation(figure1,'textbox',[0.533145583439596 0.634166293293748 0.0447949515344217 0.0398050356114624],...
        'String','D','LineStyle','none','FontWeight','bold','FontSize',25,'FitBoxToText','off');
    annotation(figure1,'textbox',[0.533145583439596 0.354754528587852 0.0447949515344217 0.039805035611462],...
        'String','F','LineStyle','none','FontWeight','bold','FontSize',25,'FitBoxToText','off');
    annotation(figure1,'textbox',[0.176641114720625 0.35567364623491 0.0447949515344216 0.039805035611462],...
        'String','E','LineStyle','none','FontWeight','bold','FontSize',25,'FitBoxToText','off');
    annotation(figure1,'textbox',[0.577298343770267 0.95455308418464 0.137539428265689 0.0365556449690648],...
        'String',{'Estimation'},'LineStyle','none','FontWeight','bold','FontSize',22,'FitBoxToText','off');
    annotation(figure1,'textbox',[0.224351307409647 0.955964239629772 0.133123024611819 0.036555644969065],...
        'String',{'Prediction'},'LineStyle','none','FontWeight','bold','FontSize',22,'FitBoxToText','off'); 
end


        