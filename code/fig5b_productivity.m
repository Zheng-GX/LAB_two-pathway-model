close all;clear;clc
set(0,'defaultAxesFontName', 'Arial','defaultAxesFontsize',17);%
set(0,'defaultTextFontName', 'Arial','defaultTextFontsize',17);%
set(0,'defaultfigurecolor','w')   
set(0,'defaultAxeslinewidth',1.5)
%% exp
ref_data=[1.16 	1.58 	2.39 	2.05 	1.67 	1.29 ];
%%
t       = 30;
dt      = 0.01;

pfit =[0.676531208214299
23.4834371824889
0.277324834016154
2.37005936793198
33543.1381801989
0.697347295185046
4.89460275931591
0.604350944322394];
var_names = {'mu_max','K_I','alp','be', 'gama_XS','gama_AS','K_S','Q0'};
for i = 1:length(pfit)
    eval([var_names{i} ' = pfit(i);']);
end
K_Q     = K_S;
v=mu_max;
m       = 5   ;
nS      = 1   ;
nQ      = 2   ;
%% figure1
r=[0.91,0.3,0.06];
b=[0,0.39,0.65];
g=[0,0,0];
%
X0(1)   = 1.17  ;S0(1)   = 9;      A0(1)   = 4.27  ;           
X0(2)   = 1.33  ;S0(2)   = 21.4;   A0(2)   = 5.44  ;
X0(3)   = 1.12  ;S0(3)   = 35.5;   A0(3)   = 3.0   ;  
X0(4)   = 0.96  ;S0(4)   = 48.1;   A0(4)   = 6.0   ;  
X0(5)   = 0.93  ;S0(5)   = 61.2;   A0(5)   = 5.9   ;  
X0(6)   = 1.10  ;S0(6)   = 77.1;   A0(6)   = 7.89  ;
%
imax=t/dt;
tspan=0:dt:t;
X=zeros(1,imax+1);
Q=zeros(1,imax+1);
A=zeros(1,imax+1);
S=zeros(1,imax+1);
tP=[];
for j=1:6
    X(1)=X0(j) ; A(1)=A0(j) ; S(1) =S0(j) ;  %
    Q(1)=Q0;
    for i=1:imax
        %
        funcQ=Q(i)./(1+Q(i));
        funcS=((S(i)/K_S).^nS)./(1+((S(i)/K_S).^nS));
        f_A=((A(i)/K_I).^m)./(1+((A(i)/K_I).^m));
        funcV=1-((S(i)/K_Q).^nQ)./(1+((S(i)/K_Q).^nQ));
        %
        X(i+1)=X(i)+dt.*mu_max.*X(i).*funcQ.*funcS.*(1-f_A);
        Q(i+1)=Q(i)+dt.* funcV.*v.*Q(i);
        if Q(i+1)==inf
            Q(i+1)=Q(i);
        end
        A(i+1)=A(i)+dt.*alp.*X(i)+be.*(X(i+1)-X(i));
        S(i+1)=S(i)-(X(i+1)-X(i))./gama_XS- ...
            (A(i+1)-A(i))./gama_AS;
        if S(i+1)<0
            S(i+1)=0;
            % if length(tP)==j-1
            %     tP(j)=i*dt;
            % end
        end
    end
x{j}=X;
a{j}=A;
s{j}=S;
end
%%
% figure
% set(gcf,'Position',[308.2,348.2,1439.2,588.8])
% for j=1:6
%     %作图
%     subplot(2,3,j)
%     yyaxis left
%     plot(tspan,s{j},'-','color',g,'LineWidth',1.5)
%     hold on
%     plot(tspan,a{j},'-','color',b,'LineWidth',1.5)
%     xlabel('Time (h)');
%     ylabel({'Lactose and lactic acid','(g/dm^3)'} );
% 
%     yyaxis right
%     plot(tspan,x{j},'-','color',r,'LineWidth',1.5)
%     ylabel('Biomass (g/dm^3)');
%     ylim([0 8])
%     %%%%%%%%%%%%%
%     hold on
%     y_limits = ylim; 
%     plot([tP(j) tP(j)],y_limits,'r--','LineWidth',1)
% end
%%
% figure
% set(gcf,'Position',[308.2,348.2,1439.2,588.8])
% for j=1:6
%     %作图
%     subplot(2,3,j)
%     plot(tspan,log(x{j}),'-','color',r,'LineWidth',1.5)
%     xlabel('Time (h)');
%     ylabel('Ln(Biomass)');
%     ylim([0 8])
%     %%%%%%%%%%%%%
%     hold on
%     y_limits = ylim; 
%     plot([tP(j) tP(j)],y_limits,'r--','LineWidth',1)
% end
%%
% for j=1:6
% productivity(j)=(a{j}(round(tP(j)/dt))-a{j}(1))/tP(j);
% end
% S_0= categorical(S0);  
% figure
% 
% bar(S_0, productivity,'FaceColor',[0.90,0.90,0.90]);
% hold on
% bar(S_0, ref_data,'FaceColor',[0.7,0.70,0.70]);
% ylim([1 3.1])
% xlabel('S_0 (g/dm^3)');
% ylabel('Productivity (g/dm^3/h)');
% legend('this work','ref.2006')
%% cut off productivity
tP=[4.5 9 11 14 12 15];
for j=1:6
productivity(j)=(a{j}(round(tP(j)/dt))-a{j}(1))/tP(j);
end
figure
set(gcf,'Position',[979,431,300,227])
plot(S0, productivity,'ko-','LineWidth',1);
ylim([0.5 2.8])
xlabel('S_0 (g/dm^3)','FontSize',14);
ylabel({'Productivity','(g/dm^3/h)'},'FontSize',14);
set(gca,'fontsize',14,'LineWidth',1)
%% productivity with t
figure
set(gca,'Position',[0.21,0.24,0.71,0.73])
hold on;
% linsty={':','--','-','-.','-','--'};
linsty={'-o','-s','-d','-v','-*','-^'};
colors = [49,98,167
    51,137,188
    111,170,146
    133,193,125
    199,129,75
    164,49,70]/255;
for j=1:6
    acid=a{j}(1:(round(tP(j)/dt))+1);
    acid0=a{j}(1);
    timeP{j}=0:dt:tP(j);
T_productivity{j}=(acid-acid0)./timeP{j};
plot(timeP{j},T_productivity{j},linsty{j},'MarkerSize',4,'MarkerIndices',1:30:length(timeP{j}),'color',colors(j,:),'LineWidth',1.5);
end
set(gca,'YTick',0:3)
box on
le=legend({'S_0=9 g/dm^3','S_0=21.4 g/dm^3','S_0=35.5 g/dm^3','S_0=48.1 g/dm^3','S_0=61.2 g/dm^3','S_0=77.1 g/dm^3'},'FontSize',13,'LineWidth',0.5);
set(le,'box','off')
xlabel('Time (h)');
ylabel('Productivity (g/dm^3/h)');
ylim([0.5 2.5])