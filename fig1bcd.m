close all ; clear;clc;
%%
K_S = 10; 
K_Q = 10;
ns  = 1;
nQ  = 1;
%
dS=0.0001;
S0=100;
S=0:dS:S0;
Q=0:0.0001:10;
A=0:0.0001:100;
m=5;
K_I=19.9;
%
muS=((S/K_S).^ns)./(1+((S/K_S).^ns));
fQS=1./(1+((S/K_Q).^nQ));
muQ=Q./(1+Q);
muA=1./(1+((A/K_I).^m));
%
fontsize=17;
set(0,'defaultAxesFontName', 'Arial','defaultAxesFontsize',fontsize);
set(0,'defaultTextFontName', 'Arial','defaultTextFontsize',fontsize);
set(0,'defaultfigurecolor','w')
%%fig1B
Orange=[0.88,0.51,0.13];
Blue=[0.49,0.69,0.81];
Green=[0.56,0.68,0.45];
textsize=20;
%fig1B
figure
plot(S,muS,'color',Blue, 'LineWidth',3);
ylim([0 1.1])
xlabel('S (mg/dm^3)','fontsize',textsize);ylabel('\mu_S','fontsize',textsize);
set(gca,'linewidth',1)
set(gcf,'position',[538.6,453.8,381.6,276.8])
%%fig1C
% figure
% plot(A,muA,'color',Blue, 'LineWidth',3);
% ylim([0 1.1])
% xlim([0 50])
% xlabel('A (mg/dm^3)','fontsize',textsize);ylabel('\mu_A','fontsize',textsize);
%%fig1D
figure
plot(S,fQS,'color',Orange, 'LineWidth',3)
ylim([0 1.1])
xlabel('S (mg/dm^3)','fontsize',textsize);ylabel('f_Q(S)','fontsize',textsize);
set(gca,'linewidth',1)
set(gcf,'position',[538.6,453.8,381.6,276.8])
%%fig1E
figure
plot(Q,muQ,'color',Orange, 'LineWidth',3)
ylim([0 1.1])
xlabel('Q','fontsize',textsize);ylabel('\mu_Q','fontsize',textsize);
set(gca,'linewidth',1)
set(gcf,'position',[538.6,453.8,381.6,276.8])