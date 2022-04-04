
% With KalmanHGF=1, it plots Supplementary Fig9A-C.
% With KalmanHGF=0, paper_fig=1, it plots Fig5.
% With KalmanHGF=0, paper_fig=2, it plots Supplementary Fig4A-F.
% With KalmanHGF=0, paper_fig=3, it plots Supplementary Fig4G-N.

KalmanHGF=0;
paper_fig=1; 
mark_tested=0; % 1--> the 6 conditions tested in the psychophysics experiment are marked on the color maps; 0 --> experimental conditions are not marked

if(KalmanHGF)
    Nlevels_HGF=2; % number of levels in the Hierarchical Gaussian filter
    nb_models=10;
else
    nb_models=8;
end
R_all_sim=[0,[0.1:0.1:6]];
h_all_sim=[0.02:0.02:1];  
sigma0=1;
color_Evidence=[0.68,0.92,1];
color_Prior=[0,0,0];
color_ML=[51/255,51/255,1];
color_DR1=[0,204/255,0];
color_SA1=[1,0,0];
color_Kalman=[0.91,0.53,0.05]; 
color_HGF=[0.494,0.184,0.556]; 
color_DR2=[ 1,160/255,1];
color_SA2=[1,0.97,0];
color_Bayes=[0.8,0.8,0.8];
if(KalmanHGF)
    myMap=[color_Evidence; color_Prior; color_ML; color_DR1; color_SA1; color_Kalman; color_HGF; color_DR2; color_SA2; color_Bayes]; 
else
    myMap=[color_Evidence; color_Prior; color_ML; color_DR1; color_SA1; color_DR2; color_SA2; color_Bayes];
end
titles={'A','B','C','D','E','F'};
if(KalmanHGF||(paper_fig==2)||(paper_fig==3))   
    tolerance=[0.2,0.1,0.02];
    pos{1,1}=[0.0851548699870835 0.740426750433459 0.170965114999447 0.151376146788991];
    pos{2,1}=[0.300689164458411 0.738948410209271 0.170908678231641 0.152905198776758];
    pos{3,1}=[0.512830610127788 0.738964853413032 0.172153259234803 0.152992904507355];
    pos{4,1}=[0.0834113484203331 0.510684838406009 0.171451598724491 0.155073584219512];
    pos{5,1}=[0.299528897425891 0.508492060297534 0.170722644070282 0.155963302752294];
    pos{6,1}=[0.514356170378579 0.507204195179124 0.171011205651764 0.155963302752295];    
else
    tolerance=0.1;   
    pos{1,1}= [0.22489579257636 0.694486889658082 0.205707019328586 0.234026514078716];
    pos{2,1}= [0.524164935610611 0.691307879772521 0.206119023397757 0.235347342232337]; 
end
if(KalmanHGF)
    fig_panels={{'estimation','rel',tolerance(1,1)};{'estimation','rel',tolerance(1,2)};{'estimation','rel',tolerance(1,3)}};
else
    if(paper_fig==1)
        fig_panels={{'prediction','rel',tolerance(1,1)};{'estimation','rel',tolerance(1,1)}};
    elseif(paper_fig==2)
        fig_panels={{'prediction','rel',tolerance(1,1)};{'prediction','rel',tolerance(1,2)};{'prediction','rel',tolerance(1,3)};...
            {'prediction','abs',tolerance(1,1)};{'prediction','abs',tolerance(1,2)};{'prediction','abs',tolerance(1,3)}};
    elseif(paper_fig==3)
        fig_panels={{'estimation','rel',tolerance(1,1)};{'estimation','rel',tolerance(1,2)};{'estimation','rel',tolerance(1,3)};...
            {'estimation','abs',tolerance(1,1)};{'estimation','abs',tolerance(1,2)};{'estimation','abs',tolerance(1,3)}};
    end
end

figure1=figure;

for p=1:size(fig_panels,1)
    if(KalmanHGF)
        load(['DiminReturns_',fig_panels{p,1}{1,1},'/Complexity_all_models_KalmanHGF'],'Complexity_all_models')
        load(['DiminReturns_',fig_panels{p,1}{1,1},'/AccuracyLoss_all_models_KalmanHGF'],'AccuracyLoss_all_models') 
    else
        if(strcmp(fig_panels{p,1}{1,2},'rel'))
            load(['DiminReturns_',fig_panels{p,1}{1,1},'/Complexity_all_models'],'Complexity_all_models')
            load(['DiminReturns_',fig_panels{p,1}{1,1},'/AccuracyLoss_all_models'],'AccuracyLoss_all_models') 
        else
            load(['DiminReturns_',fig_panels{p,1}{1,1},'/Complexity_all_models_abs'],'Complexity_all_models')
            load(['DiminReturns_',fig_panels{p,1}{1,1},'/AccuracyLoss_all_models_abs'],'AccuracyLoss_all_models') 
            load(['out_optim_',fig_panels{p,1}{1,1},'/MSErr_Bayes'],'MSErr')
            MSErrBayes=MSErr;
            MSEB=[];
            for hind=1:length(h_all_sim)
                for Rind=1:length(R_all_sim)
                    MSEB=[MSEB;mean(MSErrBayes{hind,Rind})];
                end
            end
            fig_panels{p,1}{1,3}=nanmean(MSEB)*fig_panels{p,1}{1,3}/(sigma0^2);
        end
    end
    BestModel=zeros(length(h_all_sim),length(R_all_sim));
    for hind=1:length(h_all_sim)
        for Rind=1:length(R_all_sim)
            Models_below_thr=find(AccuracyLoss_all_models{hind,Rind}<fig_panels{p,1}{1,3});
            if(isempty(Models_below_thr))
                BestModel(hind,Rind)=nb_models;
            else
                BestModel(hind,Rind)=Models_below_thr(1,1);
            end
        end
    end
    subploth{p}=subplot('Position',pos{p,1},'FontSize',12,'XTick',[0.02 0.2 0.4 0.6 0.8 1],'XTickLabel',{'0.02','0.2','0.4','0.6','0.8','1'},'YTick',[0 2 4 6],'YTickLabel',{'0','2','4','6'});
    hold(subploth{p},'on');
    hold on
    xlim([h_all_sim(1,1),h_all_sim(1,end)])
    ylim([R_all_sim(1,1),R_all_sim(1,end)])
    colormap(gca,myMap);
    caxis([0.5,nb_models+0.5])
    if(p==1)
        if(KalmanHGF||(paper_fig==2)||(paper_fig==3))
            set(subploth{1},'CLim',[0.5 nb_models+0.5]);
            c=colorbar('peer',subploth{p});
            c.Ticks=[];
            c.Location='SouthOutside';
            c.Position= [0.048810250152532 0.384240454914703 0.772422208663817 0.0333062550771731];
        else
            set(subploth{1},'CLim',[0.5 nb_models+0.5]);
            c=colorbar('peer',subploth{p});
            c.Ticks=[];
            c.Location='SouthOutside';
            c.Position= [0.127914303717707 0.575954508529651 0.67989918084436 0.0276198212835089];
        end
    end
    surf(h_all_sim,R_all_sim,BestModel','LineStyle','none');
    view(0,90);
    grid off
    set(gca,'YDir','normal')
    % marking conditions tested in the psychophysics experiment
    hR_sorted{1,1}=[0.1,0.01;0.1,0.8;0.1,4];
    Markers_hRcond{1,1}={'o';'o';'o'};
    hR_sorted{1,2}=[0.08,1;0.38,1;0.8,1];
    Markers_hRcond{1,2}={'d';'d';'d'};
    Marker_Sizes_hRcond=sqrt([15;60;105]);
    MarkerLineWidth=1.5;          
    if((mark_tested)&&strcmp(fig_panels{p,1}{1,1},'estimation')&&(p==2))
        for dataset_option=1:2
            % pl_legend{1,dataset_option}=zeros(1,size(hR_sorted{1,1},1));
            % labels_legend{1,dataset_option}=cell(1,size(hR_sorted{1,1},1));
            for hRind=1:size(hR_sorted{1,dataset_option},1)
                plot3(hR_sorted{1,dataset_option}(hRind,1),hR_sorted{1,dataset_option}(hRind,2),nb_models+1,...
                    'Marker',Markers_hRcond{1,dataset_option}{hRind,1},'MarkerSize',Marker_Sizes_hRcond(hRind,1),'Color','k','LineStyle','none','LineWidth',MarkerLineWidth);
            end
        end
    end
    hold off
    xlabel('Volatility','FontSize',12)
    ylabel('Noise','FontSize',12)
end

if(~KalmanHGF)
    if(paper_fig==1)
        annotation(figure1,'textbox',[0.454746601135585 0.513403736799348 0.104584527220629 0.0564551522832771],...
            'String',{'Sliding','Window'},'LineStyle','none','HorizontalAlignment','center','FontSize',12,'FitBoxToText','off');
        annotation(figure1,'textbox',[0.385769438089917 0.515028432168967 0.0831925924290513 0.0553273626999846],...
            'String',{'Delta','Rule'},'LineStyle','none','HorizontalAlignment','center','FontSize',12,'FitBoxToText','off');
        annotation(figure1,'textbox',[0.271524728209604 0.522581230528667 0.138975502782403 0.0473454955213762],...
            'String',{'Memoryless','Model'},'LineStyle','none','HorizontalAlignment','center','FontSize',12,'FitBoxToText','off');
        annotation(figure1,'textbox',[0.216489931315465 0.515028432168967 0.0792639024317876 0.0553994394778637],...
            'String',{'Prior','Model'},'LineStyle','none','HorizontalAlignment','center','FontSize',12,'FitBoxToText','off');
        annotation(figure1,'textbox',[0.130298988597782 0.526949614528908 0.0789014973751493 0.0428134548372453],...
            'String',{'Evidence'},'LineStyle','none','HorizontalAlignment','center','FontSize',12,'FitBoxToText','off');
        annotation(figure1,'textbox',[0.622267187174609 0.51177904142973 0.119736107685302 0.0577187921360295],...
            'String',{'Mixture','of Sliding','Windows'},'LineStyle','none','HorizontalAlignment','center','FontSize',12,'FitBoxToText','off');
        annotation(figure1,'textbox',[0.53405042598998 0.510966693744921 0.119736107685302 0.0577187921360295],...
            'String',{'Mixture','of Delta','Rules'},'LineStyle','none','HorizontalAlignment','center','FontSize',12,'FitBoxToText','off');
        annotation(figure1,'textbox',[0.687490998863395 0.51096669374492 0.159719213082884 0.0580959589326526],...
            'String',{'Exact','Bayesian','Model'},'LineStyle','none','HorizontalAlignment','center','FontSize',12,'FitBoxToText','off');
        annotation(figure1,'textbox',[0.270031545741326 0.934824533747419 0.133123024611819 0.0365556449690648],...
            'String',{'Prediction'},'LineStyle','none','FontWeight','bold','FontSize',18);
        annotation(figure1,'textbox',[0.570347003154578 0.931575143008169 0.137539428265689 0.0365556449690648],...
            'String',{'Estimation'},'LineStyle','none','FontWeight','bold','FontSize',18);
        annotation(figure1,'textbox',[0.220189274447951 0.935636881529059 0.0447949515344217 0.0398050356114618],...
            'String',{'A'},'LineStyle','none','FontWeight','bold','FontSize',20);
        annotation(figure1,'textbox',[0.519242902208216 0.935636881529059 0.0447949515344217 0.0398050356114618],...
            'String',{'B'},'LineStyle','none','FontWeight','bold','FontSize',20);
    else
        annotation(figure1,'textbox',[0.278212399912622 0.950259139758787 0.246497789227098 0.0365556449690647],...
        'String',{'Decreasing tolerance'},'LineStyle','none','FontSize',18,'FitBoxToText','off');
        annotation(figure1,'arrow',[0.0811470408785829 0.691656686175699],[0.945572705117779 0.94476035743297],'LineWidth',2,'HeadWidth',18,'HeadLength',18);
        annotation(figure1,'textbox',[0.687929255267061 0.771730300568638 0.124151281645688 0.0818594638505279],...
        'String',{'Normalized inaccuracy'},'LineStyle','none','HorizontalAlignment','center','FontSize',18,'FitBoxToText','off');
        annotation(figure1,'textbox',[0.687614399023794 0.536961819658808 0.16595485051861 0.0818594638505278],...
        'String',{'Non-normalized','inaccuracy'},'LineStyle','none','HorizontalAlignment','center','FontSize',18,'FitBoxToText','off');
    end
end

   
    