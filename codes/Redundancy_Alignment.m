
% This code calculates Redundancy and Alignment and plots Fig3

h_slice=0.32; 
task='prediction';
N=2; 
nb_restarts=10; 
nb_iterat=10;
h_all_sim=[0.02:0.02:0.98]; 
R_all_sim=[0,0.1:0.1:6]; 
sigma0=1;
model_families={'SA','DR'};
% Parameters for Fig3A
xmixt=0.85;
ymixt=0.15;
xirr=0.77;
yirr=0.55;
theta1=[acos(((xmixt+xirr)/2-xmixt)/norm([xmixt-xirr,ymixt-yirr]./2)):0.0000001:acos(-(xmixt-xmixt*cos(pi/4))/(norm([xmixt-xmixt*cos(pi/4),ymixt-xmixt*sin(pi/4)])))];
xrel=xmixt+norm([xmixt-xirr,ymixt-yirr])*cos(theta1(1,1)+pi/2);
yrel=ymixt+norm([xmixt-xirr,ymixt-yirr])*sin(theta1(1,1)+pi/2);

% Calculation of Redundancy and Alignment 
meanLogEigRatio=cell(length(model_families),1);
meanAlign_diffLRs_eigvmin=cell(length(model_families),1);
meanLogEigmax=cell(length(model_families),1);
meanLogEigmin=cell(length(model_families),1);
for mm=1:length(model_families)   
    model=model_families{1,mm};   
    load(['out_optim_',task,'/MSErr_',model,'_N1'],'MSErr')
    MSErr1=MSErr;
    load(['out_optim_',task,'/LearningRates_',model,'_N1'],'LearningRates')  
    LearningRate1=LearningRates;
    load(['out_optim_',task,'/Hes_',model,'_N1'],'Hes')
    Hes1=Hes;
    load(['out_optim_',task,'/MSErr_',model,'_N2'],'MSErr')
    MSErr2=MSErr;
    load(['out_optim_',task,'/LearningRates_',model,'_N2'],'LearningRates')  
    LearningRate2=LearningRates;
    load(['out_optim_',task,'/Hes_',model,'_N2'],'Hes')
    Hes2=Hes;
    meanLogEigRatio{mm,1}=zeros(length(R_all_sim),length(h_all_sim));
    meanAlign_diffLRs_eigvmin{mm,1}=zeros(length(R_all_sim),length(h_all_sim));
    meanLogEigmax{mm,1}=zeros(length(R_all_sim),length(h_all_sim));
    meanLogEigmin{mm,1}=zeros(length(R_all_sim),length(h_all_sim));    
    for Rind=1:length(R_all_sim)
       for hind=1:length(h_all_sim)
           [~,indminMSE1]=min(MSErr1{hind,Rind},[],2);
           [~,indminMSE2]=min(MSErr2{hind,Rind},[],2);
           q=0;
           logEigRatio=[];
           align_diffLRs_eigvmin=[];
           logEigmax_v=[];
           logEigmin_v=[];
           for iterat=1:nb_iterat
               hs1=Hes1{hind,Rind}{iterat,indminMSE1(iterat)};
               hs2=Hes2{hind,Rind}{iterat,indminMSE2(iterat)};
               [V,D] = eig(hs2);             
               min_reached=(~any(diag(D)<=0))&&(hs1>0); 
               if(min_reached)                   
                   q=q+1;                
                   LR_1=LearningRate1{1,1}{hind,Rind}(iterat,indminMSE1(iterat));
                   LR_2=zeros(N,1);
                   for nb_node=1:N
                       LR_2(nb_node,1)=LearningRate2{nb_node,1}{hind,Rind}(iterat,indminMSE2(iterat));
                   end
                   [~,Is]=sort(diag(D),'ascend');
                   eigmin=D(Is(1),Is(1)); 
                   eigvmin=V(:,Is(1));
                   eigmax=D(Is(2),Is(2)); 
                   eigvmax=V(:,Is(2));
                   logEigRatio=[logEigRatio;log10(eigmax/eigmin)];
                   diffLRs_1=repmat(LR_1,2,1)-LR_2;  
                   align_1=abs((pi/2)-abs(acos(dot(diffLRs_1,eigvmin)/(norm(diffLRs_1)*norm(eigvmin)))))/(pi/2); 
                   align_diffLRs_eigvmin = [align_diffLRs_eigvmin; align_1]; 
                   logEigmax_v=[logEigmax_v;log10(eigmax)];
                   logEigmin_v=[logEigmin_v;log10(eigmin)];
               end
           end
           meanLogEigRatio{mm,1}(Rind,hind)=mean(logEigRatio);          
           meanAlign_diffLRs_eigvmin{mm,1}(Rind,hind)=mean(align_diffLRs_eigvmin);
           meanLogEigmax{mm,1}(Rind,hind)=mean(logEigmax_v);
           meanLogEigmin{mm,1}(Rind,hind)=mean(logEigmin_v);
       end
    end
end


% Figure 3
   
h_slice_ind=find(round(100*h_all_sim)==round(100*h_slice));
figure1=figure;
axis_font_size=16;
label_font_size=16;
title_font_size=19;
title_font_size_large=20;
title_font_size_huge=23;

% Fig3B

% colormaps

axes1 = axes('Position',[0.270487706375434 0.500406173842404 0.13330175519937 0.146222583265629],'FontSize',axis_font_size);  
hold(axes1,'on');
surf(h_all_sim,R_all_sim,meanLogEigRatio{1,1},'FaceColor','interp', 'LineStyle', 'none');
plot3(h_all_sim(1,h_slice_ind)*ones(1,length(R_all_sim)),R_all_sim,meanLogEigRatio{1,1}(:,h_slice_ind), '-r','LineWidth',2);
plot3(h_all_sim(1,h_slice_ind+1)*ones(1,length(R_all_sim)),R_all_sim,meanLogEigRatio{1,1}(:,h_slice_ind+1), '-r','LineWidth',2);
plot3([h_all_sim(1,h_slice_ind),h_all_sim(1,h_slice_ind+1)],[R_all_sim(1,1),R_all_sim(1,1)],[meanLogEigRatio{1,1}(1,h_slice_ind),meanLogEigRatio{1,1}(1,h_slice_ind+1)], '-r','LineWidth',2);
plot3([h_all_sim(1,h_slice_ind),h_all_sim(1,h_slice_ind+1)],[R_all_sim(1,end),R_all_sim(1,end)],[meanLogEigRatio{1,1}(end,h_slice_ind),meanLogEigRatio{1,1}(end,h_slice_ind+1)], '-r','LineWidth',2);
xlim([h_all_sim(1,1),h_all_sim(1,end)]);
ylim([R_all_sim(1,1),R_all_sim(1,end)]);
view(0,90);
grid off
set(axes1,'CLim',[0 8],'XTick',[0.02 0.3 0.6 0.9],'YTick',[0 2 4 6]);
colormap('parula');
colorbar('peer',axes1,'Position',[0.434346639429612 0.499593826157595 0.0102869520229256 0.146222583265638],'Ticks',[2 4 6 8],'TickLabels',{'2','4','6','\geq 8'});
xlabel('Volatility','FontSize',label_font_size,'FontName','Helvetica');
ylabel('Noise','FontSize',label_font_size,'FontName','Helvetica');

axes2 = axes('Parent',figure1,'Position',[0.486129000051307 0.500406173842404 0.13330175519937 0.146222583265629],'FontSize',axis_font_size);
hold(axes2,'on');
surf(h_all_sim,R_all_sim,meanLogEigRatio{2,1},'FaceColor','interp', 'LineStyle', 'none');
plot3(h_all_sim(1,h_slice_ind)*ones(1,length(R_all_sim)),R_all_sim,meanLogEigRatio{2,1}(:,h_slice_ind), '-r','LineWidth',2);
plot3(h_all_sim(1,h_slice_ind+1)*ones(1,length(R_all_sim)),R_all_sim,meanLogEigRatio{2,1}(:,h_slice_ind+1), '-r','LineWidth',2);
plot3([h_all_sim(1,h_slice_ind),h_all_sim(1,h_slice_ind+1)],[R_all_sim(1,1),R_all_sim(1,1)],[meanLogEigRatio{2,1}(1,h_slice_ind),meanLogEigRatio{2,1}(1,h_slice_ind+1)], '-r','LineWidth',2);
plot3([h_all_sim(1,h_slice_ind),h_all_sim(1,h_slice_ind+1)],[R_all_sim(1,end),R_all_sim(1,end)],[meanLogEigRatio{2,1}(end,h_slice_ind),meanLogEigRatio{2,1}(end,h_slice_ind+1)], '-r','LineWidth',2);
xlim([h_all_sim(1,1),h_all_sim(1,end)]);
ylim([R_all_sim(1,1),R_all_sim(1,end)]);
view(0,90);
grid off
set(axes2,'CLim',[0 8],'XTick',[0.02 0.3 0.6 0.9],'YTick',[0 2 4 6],'YAxisLocation','right');
colormap('parula');
xlabel('Volatility','FontSize',label_font_size,'FontName','Helvetica');
ylabel('Noise','FontSize',label_font_size,'FontName','Helvetica');

% Slices:

axes5 = axes('Position',[0.180789862256061 0.500406173842404 0.0481445174330471 0.140243704305425],'FontSize',axis_font_size);  
hold(axes5,'on');
axis(axes5,[min(R_all_sim),max(R_all_sim),0,0.15+max([max(meanLogEigRatio{1,1}(:,h_slice_ind)),max(meanLogEigRatio{2,1}(:,h_slice_ind))])]);
plot(R_all_sim,meanLogEigRatio{1,1}(:,h_slice_ind),'ok','MarkerFaceColor','k');
axesLimits1 = xlim(axes5);
xplot1 = linspace(axesLimits1(1),axesLimits1(2));
fitResults1 = polyfit(R_all_sim',meanLogEigRatio{1,1}(:,h_slice_ind),4);
yplot1 = polyval(fitResults1,xplot1);
view([90 -90]); 
plot(xplot1,yplot1,'DisplayName','4th degree','Tag','4th degree','Parent',axes5,'Color','r','LineWidth',3);
set(axes5,'XAxisLocation','top','XTick',zeros(1,0));
ylabel('Redundancy','FontSize',label_font_size)

axes6 = axes('Position',[0.661423945301417 0.500406173842404 0.0481445174330471 0.140243704305425],'FontSize',axis_font_size);   
hold(axes6,'on');
axis(axes6,[min(R_all_sim),max(R_all_sim),0,0.15+max([max(meanLogEigRatio{1,1}(:,h_slice_ind)),max(meanLogEigRatio{2,1}(:,h_slice_ind))])]); 
plot(R_all_sim,meanLogEigRatio{2,1}(:,h_slice_ind),'ok','MarkerFaceColor','k');
axesLimits1 = xlim(axes6);
xplot1 = linspace(axesLimits1(1),axesLimits1(2));
fitResults1 = polyfit(R_all_sim',meanLogEigRatio{2,1}(:,h_slice_ind),4);
yplot1 = polyval(fitResults1,xplot1);
view([90 -90]); 
plot(xplot1,yplot1,'DisplayName','4th degree','Tag','4th degree','Parent',axes6,'Color','r','LineWidth',3);
set(axes6,'XAxisLocation','bottom','XTick',zeros(1,0));
ylabel('Redundancy','FontSize',label_font_size)

% Fig3C

% colormaps

axes3 = axes('Parent',figure1,'Position',[0.270487706375434 0.272948822095857 0.132895401718021 0.141348497156772],'FontSize',axis_font_size);
hold(axes3,'on');
surf(h_all_sim,R_all_sim,meanAlign_diffLRs_eigvmin{1,1},'FaceColor','interp', 'LineStyle', 'none');
plot3(h_all_sim(1,h_slice_ind)*ones(1,length(R_all_sim)),R_all_sim,meanAlign_diffLRs_eigvmin{1,1}(:,h_slice_ind), '-r','LineWidth',2);
plot3(h_all_sim(1,h_slice_ind+1)*ones(1,length(R_all_sim)),R_all_sim,meanAlign_diffLRs_eigvmin{1,1}(:,h_slice_ind+1), '-r','LineWidth',2);
plot3([h_all_sim(1,h_slice_ind),h_all_sim(1,h_slice_ind+1)],[R_all_sim(1,1),R_all_sim(1,1)],[meanAlign_diffLRs_eigvmin{1,1}(1,h_slice_ind),meanAlign_diffLRs_eigvmin{1,1}(1,h_slice_ind+1)], '-r','LineWidth',2);
plot3([h_all_sim(1,h_slice_ind),h_all_sim(1,h_slice_ind+1)],[R_all_sim(1,end),R_all_sim(1,end)],[meanAlign_diffLRs_eigvmin{1,1}(end,h_slice_ind),meanAlign_diffLRs_eigvmin{1,1}(end,h_slice_ind+1)], '-r','LineWidth',2);
xlim([h_all_sim(1,1),h_all_sim(1,end)]);
ylim([R_all_sim(1,1),R_all_sim(1,end)]);
view(0,90);
grid off
set(axes3,'CLim',[0 1],'XTick',[0.02 0.3 0.6 0.9],'YTick',[0 2 4 6]);
colormap('parula')
colorbar('peer',axes3,'Position',[0.434922004939554 0.272136474411048 0.0102328101701734 0.141348497156783],'Ticks',[0 0.2 0.4 0.6 0.8 1],'TickLabels',{'0','0.2','0.4','0.6','0.8','1'});
xlabel('Volatility','FontSize',label_font_size,'FontName','Helvetica')
ylabel('Noise','FontSize',label_font_size,'FontName','Helvetica')
    
axes4 = axes('Parent',figure1,'Position',[0.486129000051307 0.272948822095857 0.132895401718021 0.141348497156772],'FontSize',axis_font_size);
hold(axes4,'on');
surf(h_all_sim,R_all_sim,meanAlign_diffLRs_eigvmin{2,1},'FaceColor','interp', 'LineStyle', 'none');
plot3(h_all_sim(1,h_slice_ind)*ones(1,length(R_all_sim)),R_all_sim,meanAlign_diffLRs_eigvmin{2,1}(:,h_slice_ind), '-r','LineWidth',2);
plot3(h_all_sim(1,h_slice_ind+1)*ones(1,length(R_all_sim)),R_all_sim,meanAlign_diffLRs_eigvmin{2,1}(:,h_slice_ind+1), '-r','LineWidth',2);
plot3([h_all_sim(1,h_slice_ind),h_all_sim(1,h_slice_ind+1)],[R_all_sim(1,1),R_all_sim(1,1)],[meanAlign_diffLRs_eigvmin{2,1}(1,h_slice_ind),meanAlign_diffLRs_eigvmin{2,1}(1,h_slice_ind+1)], '-r','LineWidth',2);
plot3([h_all_sim(1,h_slice_ind),h_all_sim(1,h_slice_ind+1)],[R_all_sim(1,end),R_all_sim(1,end)],[meanAlign_diffLRs_eigvmin{2,1}(end,h_slice_ind),meanAlign_diffLRs_eigvmin{2,1}(end,h_slice_ind+1)], '-r','LineWidth',2);
xlim([h_all_sim(1,1),h_all_sim(1,end)]);
ylim([R_all_sim(1,1),R_all_sim(1,end)]);
view(0,90);
grid off
set(axes4,'CLim',[0 1],'XTick',[0.02 0.3 0.6 0.9],'YTick',[0 2 4 6],'YAxisLocation','right');
colormap('parula')
xlabel('Volatility','FontSize',label_font_size,'FontName','Helvetica')
ylabel('Noise','FontSize',label_font_size,'FontName','Helvetica')

% Slices

axes7 = axes('Parent',figure1,'Position',[0.180789862256061 0.272948822095857 0.0481445174330471 0.140243704305425],'FontSize',axis_font_size);
hold(axes7,'on');
axis(axes7,[min(R_all_sim),max(R_all_sim),0,0.15+max([max(meanAlign_diffLRs_eigvmin{1,1}(:,h_slice_ind)),max(meanAlign_diffLRs_eigvmin{2,1}(:,h_slice_ind))])]); 
plot(R_all_sim,meanAlign_diffLRs_eigvmin{1,1}(:,h_slice_ind),'ok','MarkerFaceColor','k');
axesLimits1 = xlim(axes7);
xplot1 = linspace(axesLimits1(1),axesLimits1(2));
fitResults1 = polyfit(R_all_sim',meanAlign_diffLRs_eigvmin{1,1}(:,h_slice_ind),4);
yplot1 = polyval(fitResults1,xplot1);
view([90 -90]); 
plot(xplot1,yplot1,'DisplayName','4th degree','Tag','4th degree','Parent',axes7,'Color','r','LineWidth',3);
set(axes7,'XAxisLocation','top','XTick',zeros(1,0),'YTick',[0 1],'YTickLabel',[0,1]);
ylabel('Alignment','FontSize',label_font_size)

axes8 = axes('Parent',figure1,'Position',[0.661423945301417 0.272948822095857 0.0481445174330471 0.140243704305425],'FontSize',axis_font_size);
hold(axes8,'on');
axis(axes8,[min(R_all_sim),max(R_all_sim),0,0.15+max([max(meanAlign_diffLRs_eigvmin{1,1}(:,h_slice_ind)),max(meanAlign_diffLRs_eigvmin{2,1}(:,h_slice_ind))])]); 
plot1=plot(R_all_sim,meanAlign_diffLRs_eigvmin{2,1}(:,h_slice_ind),'ok','MarkerFaceColor','k');
axesLimits1 = xlim(axes8);
xplot1 = linspace(axesLimits1(1),axesLimits1(2));
fitResults1 = polyfit(R_all_sim',meanAlign_diffLRs_eigvmin{2,1}(:,h_slice_ind),4);
yplot1 = polyval(fitResults1,xplot1);
view([90 -90]); 
plot(xplot1,yplot1,'DisplayName','4th degree','Tag','4th degree','Parent',axes8,'Color','r','LineWidth',3);
set(axes8,'XAxisLocation','bottom','XTick',zeros(1,0),'YTick',[0 1],'YTickLabel',[0,1]);
ylabel('Alignment','FontSize',label_font_size)
  

% Fig3A
 
Reduction_rotation=-acosd((xmixt-xmixt*cos(pi/4))/norm([xmixt-xmixt*cos(pi/4),ymixt-xmixt*sin(pi/4)]));
Irr_rotation=-acosd((xmixt-xirr)/norm([xmixt-xirr,ymixt-yirr]));
Rel_rotation=Irr_rotation+90;
Color_model_red=[102,102,102]./255;
Color_irr=[0 0.6 1];
Color_rel=[0,204,51]./255;

% Fig3A, left panel

axes7 = axes('Parent',figure1,...
    'Position',[0.270487706375434 0.761462225832639 0.131510650132896 0.201169780666143],'FontSize',axis_font_size);
hold(axes7,'on');
set(axes7,'XTick',[0 1],'YTick',[0 1]);
plot([1,1],[0,1],'-','Color',[0.8,0.8,0.8],'LineWidth',1)  
plot([0,1],[1,1],'-','Color',[0.8,0.8,0.8],'LineWidth',1)  
plot([xmixt,xmixt*cos(pi/4)],[ymixt,xmixt*sin(pi/4)],'-','Color',Color_model_red,'LineWidth',8) 
plot([0,1],[0,1],'-k','LineWidth',3)
plot([0,1],[0,0],'-k','LineWidth',3)
plot([0,0],[0,1],'-k','LineWidth',3)
plot(xmixt,ymixt,'Marker','o','MarkerSize',25,'MarkerEdgeColor','k','MarkerFaceColor','k')
plot(xmixt*cos(pi/4),xmixt*sin(pi/4),'Marker','o','MarkerSize',25,'MarkerEdgeColor','k','MarkerFaceColor','k')
xlabel('\alpha_1','Color','r','FontSize',18,'FontWeight','bold');
ylabel('\alpha_2','Color','r','FontSize',18,'FontWeight','bold');
text('FontWeight','bold','FontSize',18,'Rotation',45,'String','\alpha_1=\alpha_2','Position',[0.0843940904327378 0.197100428019626 0],'Color','r');
text('FontWeight','bold','FontSize',16,'String',{'Mixture','Model'},'Position',[0.936896935908019 0.138357901789988 0]);
text('HorizontalAlignment','center','FontWeight','bold','FontSize',16,'String',{'Single-Unit','Model'},'Position',[0.369908035400289 0.661403035138375 0]);
text('HorizontalAlignment','center','FontWeight','bold','FontSize',13,'Rotation',Reduction_rotation,'String',{'Model','reduction'},'Position',[0.634795597086021 0.316524966675681 0],'Color',Color_model_red);

% Fig3A, right panel

axes8 = axes('Parent',figure1,'Position',[0.486129000051307 0.761462225832639 0.131510650132896 0.201169780666143],'FontSize',axis_font_size);
hold(axes8,'on');
set(axes8,'XTick',[0 1],'YTick',[0 1]);
plot(xmixt+norm([xmixt-xirr,ymixt-yirr].*(3/4))*cos(theta1),ymixt+norm([xmixt-xirr,ymixt-yirr].*(3/4))*sin(theta1),'-','Color','r','LineWidth',1)
patch([xmixt+norm([xmixt-xirr,ymixt-yirr].*(3/4))*cos(theta1),xmixt],[ymixt+norm([xmixt-xirr,ymixt-yirr].*(3/4))*sin(theta1),ymixt],'r')  
plot(xrel,yrel,'.k')
plot([xmixt,xrel],[ymixt,yrel],'-','Color',Color_rel,'LineWidth',8) 
plot(xirr,yirr,'.k')
plot([xmixt,xirr],[ymixt,yirr],'-','Color',Color_irr,'LineWidth',8)
plot([1,1],[0,1],'-','Color',[0.8,0.8,0.8],'LineWidth',1)  
plot([0,1],[1,1],'-','Color',[0.8,0.8,0.8],'LineWidth',1)  
plot([xmixt,xmixt*cos(pi/4)],[ymixt,xmixt*sin(pi/4)],'-','Color',Color_model_red,'LineWidth',8) 
plot([0,1],[0,1],'-k','LineWidth',3)
plot([0,1],[0,0],'-k','LineWidth',3)
plot([0,0],[0,1],'-k','LineWidth',3)
plot(xmixt,ymixt,'Marker','o','MarkerSize',25,'MarkerEdgeColor','k','MarkerFaceColor','k')
plot(xmixt*cos(pi/4),xmixt*sin(pi/4),'Marker','o','MarkerSize',25,'MarkerEdgeColor','k','MarkerFaceColor','k')
xlabel('\alpha_1','Color',[1 0 0],'FontSize',18,'FontWeight','bold');
ylabel('\alpha_2','Color',[1 0 0],'FontSize',18,'FontWeight','bold');
text('FontWeight','bold','FontSize',18,'Rotation',45,'String','\alpha_1=\alpha_2','Position',[0.0885093167701864 0.201132686084142 0],'Color','r');
text('FontWeight','bold','FontSize',16,'String',{'Mixture','Model'},'Position',[0.937311777272682 0.130081286279954 0]);
text('HorizontalAlignment','center','FontWeight','bold','FontSize',16,'String',{'Single-Unit','Model'},'Position',[0.369410225762695 0.661370647523399 0]);
text('HorizontalAlignment','center','FontWeight','bold','FontSize',13,'Rotation',Reduction_rotation,'String',{'Model','reduction'},'Position',[0.63253758362729 0.328556175697449 0],'Color',Color_model_red);
text('FontWeight','bold','FontSize',22,'String','\theta','Position',[0.679885057471264 0.523086053412463 0],'Color','r');
text('FontWeight','bold','FontSize',13,'Rotation',Irr_rotation,'String','Irrelevant','Position',[0.822126273822412 0.612173854589804 0],'Color',Color_irr);
text('FontWeight','bold','FontSize',13,'Rotation',Rel_rotation,'String','Relevant','Position',[0.46325977217488 0.0218522992190878 0],'Color',Color_rel);

% Fig3D

nb_bins_2Dhist=10;
min_nb_datapoints_2Dhist=150; 

% Fig3D, left

axes11 = axes('Parent',figure1,'Position',[0.270487706375434 0.0495532087733533 0.132895401718021 0.141348497156772],'FontSize',axis_font_size);
hold(axes11,'on');
[N_2Dhist,c_2Dhist]=hist3([meanLogEigRatio{1,1}(:),meanAlign_diffLRs_eigvmin{1,1}(:)],[nb_bins_2Dhist,nb_bins_2Dhist]); 
tot_count_per_Redund=sum(N_2Dhist,2);
Excluded_Redund_Ind=find(tot_count_per_Redund<min_nb_datapoints_2Dhist);
N_2Dhist(Excluded_Redund_Ind,:)=NaN;
min_Excluded_Redund_Ind=min(Excluded_Redund_Ind);
if(isempty(min_Excluded_Redund_Ind))
    min_Excluded_Redund_Ind=length(tot_count_per_Redund)+1;
end
for ind=1:length(c_2Dhist{1,1})
    N_2Dhist(ind,:)=N_2Dhist(ind,:)./tot_count_per_Redund(ind,1);
end
surf(c_2Dhist{1,1}(1,[1:min_Excluded_Redund_Ind-1]),c_2Dhist{1,2},N_2Dhist([1:min_Excluded_Redund_Ind-1],:)','LineStyle', 'none');
xlim([min(c_2Dhist{1,1}(1,[1:min_Excluded_Redund_Ind-1])),max(c_2Dhist{1,1}(1,[1:min_Excluded_Redund_Ind-1]))]);
ylim([min(c_2Dhist{1,2}),max(c_2Dhist{1,2})]);
clim_min=0;
clim_max=0.1*ceil(10*max(max(N_2Dhist([1:min_Excluded_Redund_Ind-1],:),[],1),[],2));
view(0,90);
grid off
caxis([clim_min,clim_max]);
colormap('parula')
colorbar('peer',axes11,'Position',[0.434158415841584 0.0454914703493095 0.0105683948318591 0.146222583265638],'Ticks',[clim_min 0.1 0.2 clim_max],'TickLabels',{num2str(clim_min),'0.1','0.2',['\geq',num2str(clim_max)]});
set(axes11,'XTick',[2 4 6 8],'XTickLabels',{'2','4','6','8'},'YTick',[0.3 0.6 0.9]);
xlabel('Redundancy','FontSize',label_font_size,'FontName','Helvetica')
ylabel('Alignment','FontSize',label_font_size,'FontName','Helvetica')

%Fig3D, right

axes12 = axes('Parent',figure1,'Position',[0.486129000051307 0.0503655564581624 0.132895401718021 0.141348497156772],'FontSize',axis_font_size);
hold(axes12,'on');
[N_2Dhist,c_2Dhist]=hist3([meanLogEigRatio{2,1}(:),meanAlign_diffLRs_eigvmin{2,1}(:)],[nb_bins_2Dhist,nb_bins_2Dhist]); 
tot_count_per_Redund=sum(N_2Dhist,2);
Excluded_Redund_Ind=find(tot_count_per_Redund<min_nb_datapoints_2Dhist);
N_2Dhist(Excluded_Redund_Ind,:)=NaN;
min_Excluded_Redund_Ind=min(Excluded_Redund_Ind);
if(isempty(min_Excluded_Redund_Ind)==1)
    min_Excluded_Redund_Ind=length(tot_count_per_Redund)+1;
end
for ind=1:length(c_2Dhist{1,1})
    N_2Dhist(ind,:)=N_2Dhist(ind,:)./tot_count_per_Redund(ind,1);
end
surf(c_2Dhist{1,1}(1,[1:min_Excluded_Redund_Ind-1]),c_2Dhist{1,2},N_2Dhist([1:min_Excluded_Redund_Ind-1],:)','LineStyle', 'none');
xlim([min(c_2Dhist{1,1}(1,[1:min_Excluded_Redund_Ind-1])),max(c_2Dhist{1,1}(1,[1:min_Excluded_Redund_Ind-1]))]);
ylim([min(c_2Dhist{1,2}),max(c_2Dhist{1,2})]);
view(0,90);
grid off
caxis([clim_min,clim_max]);
colormap('parula')
set(axes12,'YAxisLocation','right');
set(axes12,'XTick',[1 2 3 4 5 6 7 8],'XTickLabels',{'1','2','3','4','5','6','7','8'},'YTick',[0.3 0.6 0.9]);
xlabel('Redundancy','FontSize',label_font_size,'FontName','Helvetica')
ylabel('Alignment','FontSize',label_font_size,'FontName','Helvetica')

% Titles

annotation(figure1,'textbox',[0.361542649090529 0.959194964461158 0.178987031729579 0.0422420785932595],...
    'String',{'Alignment definition'},'LineStyle','none','HorizontalAlignment','center','FontSize',title_font_size,'FitBoxToText','off');
annotation(figure1,'textbox',[0.39177272845191 0.64400406275522 0.164461242923514 0.0422420785932596],...
    'String','Redundancy','LineStyle','none','FontSize',title_font_size,'FitBoxToText','off');
annotation(figure1,'textbox',[0.402877871216184 0.413297320269446 0.139886574468751 0.0422420785932597],...
    'String',{'Alignment'},'LineStyle','none','FontSize',title_font_size,'FitBoxToText','off');
annotation(figure1,'textbox',[0.356701108554812 0.188277011577325 0.139886574468752 0.0422420785932597],...
    'String','p(Alignment|Redundancy)','LineStyle','none','FontSize',19,'FitBoxToText','off');
annotation(figure1,'textbox',[0.245794899620185 0.648065802196084 0.181921081513191 0.0763606803384281],'Color',[0.34901961684227 0.200000002980232 0.329411774873734],...
    'String',['Sliding','newline','Windows'],'LineStyle','none','HorizontalAlignment','center','FontSize',title_font_size_large,'FitBoxToText','off');
annotation(figure1,'textbox',[0.471437054793628 0.647253454511272 0.160684486172188 0.076360680338428],'Color',[0.34901961684227 0.200000002980232 0.329411774873734],...
    'String',['Delta','newline','Rules'],'LineStyle','none','HorizontalAlignment','center','FontSize',title_font_size_large,'FitBoxToText','off');
annotation(figure1,'textbox',[0.270959623642457 0.948009748172217 0.037807183364839 0.0534272948821931],...
    'String','A','LineStyle','none','FontWeight','bold','FontSize',title_font_size_huge,'FitBoxToText','off');
annotation(figure1,'textbox',[0.177090931293027 0.631194151096666 0.037807183364839 0.053427294882193],...
    'String','B','LineStyle','none','FontWeight','bold','FontSize',title_font_size_huge,'FitBoxToText','off');
annotation(figure1,'textbox',[0.177985942118038 0.398050365556453 0.037807183364839 0.0534272948821929],...
    'String','C','LineStyle','none','FontWeight','bold','FontSize',23,'FitBoxToText','off');
annotation(figure1,'textbox',[0.269941589120612 0.177091795288378 0.037807183364839 0.0534272948821928],...
    'String','D','LineStyle','none','FontWeight','bold','FontSize',23,'FitBoxToText','off');


    