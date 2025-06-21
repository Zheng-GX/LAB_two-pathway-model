close all;clear;clc
set(0,'defaultAxesFontName', 'Arial','defaultAxesFontsize',14);%坐标轴
set(0,'defaultTextFontName', 'Arial','defaultTextFontsize',14);%文字
set(0,'defaultfigurecolor','w')   
set(0,'defaultLegendInterpreter','tex');
set(0,'defaultTextInterpreter','tex');
set(0,'defaultAxeslinewidth',1)
set(0,'defaultfigureposition',[500,300,450,320]);
%% exp
%fig1
time= [0.00 	2.00 	4.00 	6.00 	8.00 	15.00 	17.00 	19.00 	21.00 	23.00 	25.00 	27.00 	29.00 	31.00 	39.00 	41.00 	43.00 	45.00 	47.00 	49.00 ]';
XE1 = [0.15 	0.23 	0.38 	0.49 	0.75 	1.68 	1.70 	1.79 	1.94 	1.99 	1.93 	1.98 	2.01 	2.02 	1.97 	2.05 	2.07 	2.07 	2.08 	2.10 ];
SE1 = [10.47 	9.81 	9.58 	8.18 	7.40 	2.71 	1.86 	1.48 	1.55 	1.44 	1.17 	1.40 	1.47 	1.32 	1.08 	1.06 	1.53 	1.24 	1.55 	1.08 ];
AE1 = [1.09 	1.77 	2.08 	3.49 	4.33 	9.26 	10.14 	10.56 	10.48 	10.59 	10.95 	10.64 	10.57 	10.72 	10.92 	11.03 	10.49 	10.81 	10.46 	10.96 ];
XE2 = [0.15 	0.22 	0.32 	0.54 	0.68 	1.15 	1.48 	1.53 	1.47 	1.78 	1.81 	2.07 	2.06 	2.04 	2.15 	2.07 	2.26 	2.22 	2.22 	2.29 ];
SE2 = [21.00 	20.10 	19.15 	18.05 	17.28 	12.49 	12.09 	11.43 	10.67 	10.31 	9.58 	9.95 	9.73 	9.32 	8.19 	8.55 	7.46 	7.62 	7.37 	7.25 ];
AE2 = [1.50 	2.36 	3.44 	4.56 	5.37 	10.42 	10.82 	11.54 	12.39 	12.70 	13.49 	13.04 	13.33 	13.76 	14.98 	14.57 	15.74 	15.61 	15.80 	15.97 ];
XE = [XE1';XE2'];
SE = [SE1';SE2'];
AE = [AE1';AE2'];
% XE = [XE1'];
% SE = [SE1'];
% AE = [AE1'];
%% fit
global nA nQ dt t Q0 var_names nS abc
t       = 50  ;
dt      = 0.01;
KS      = 8   ;
Q0      = 1  ;
nA      = 5  ;
nQ      = 2 ;
nS      = 1;
abc=[1 1 1];
y_data = [XE*abc(1); AE*abc(2); SE*abc(3)];
y_mean = mean(y_data);
SST = sum((y_data - y_mean).^2);% 
%
p0 =[0.583253564063580
6.20368715117026
0
7.93466081597068
1.05745537098682
7.11125084570527
9.53556909659832
0.935147932615253];
% p0 =[0.767143092904400
% 6.11693342834676
% 0
% 8.06781541959975
% 1.06228092342874
% 5.14964422722209
% 6.84035358009619
% 8.35484934214176
% 0.596221409911512];
%%%%%%%%%
% var_names = {'mu_max','alpha0', 'beta0',  'gama_XS',  'gama_AS','KQ','KI','KS','Q0'};
var_names = {'mu_max','alpha0', 'beta0',  'gama_XS',  'gama_AS','KI','KS','Q0'};
%%%%%%%%% 
options = optimoptions('lsqcurvefit', 'Algorithm', 'levenberg-marquardt');
% options = optimset('tolx',1e-16);
pfit=p0;
p00=zeros(length(p0),1);
while ~all(round(p00, 4) == round(p0, 4), 'all')
    p00=p0;
    [pfit,resnorm] = lsqcurvefit(@Lplantarum_FitFunc,p0,time,y_data,zeros(length(p0),1), [],options);
    p0=pfit;
    resnorm
    % 计算R²
    R_squared = 1 - resnorm / SST;
    R_squared
end
%%
for i = 1:length(pfit)
    eval([var_names{i} ' = pfit(i);']);
end
%
r=[0.91,0.3,0.06];
b=[0,0.39,0.65];
g=[0,0,0];
%
v=mu_max;
KQ=KS;
imax=t/dt;
tspan=0:dt:t;
X=zeros(1,imax+1);
Q=zeros(1,imax+1);
A=zeros(1,imax+1);
S=zeros(1,imax+1);
%%
X(1)=0.15 ;
A(1)=1.09 ; 
S(1)=10.47; 
Q(1)=Q0;
for i=1:imax
    
    funcQ=Q(i)./(1+Q(i));
    funcS=S(i)^nS./(KS^nS+S(i)^nS);
    f_A=1./(1+((A(i)/KI).^nA));
    funcV=1./(1+((S(i)/KQ).^nQ));
    
    X(i+1)=X(i)+dt.*mu_max.*X(i).*funcQ.*funcS.*f_A;
    Q(i+1)=Q(i)+dt.* funcV.*v.*Q(i);
    if Q(i+1)==inf
        Q(i+1)=Q(i);
    end
    A(i+1)=A(i)+alpha0.*(X(i+1)-X(i))+dt.*beta0.*X(i);
    S(i+1)=S(i)-(X(i+1)-X(i))./gama_XS- ...
        (A(i+1)-A(i))./gama_AS;
    if S(i+1)<0
        S(i+1)=0;
    end
end

%% 11g/l
figure
set(gca,'Position',[0.132,0.165,0.715,0.760])
yyaxis left
plot(tspan,A,'-','color',g,'LineWidth',1.5);
hold on
plot(tspan,S,'-','color',b,'LineWidth',1.5);
xlabel('Time (h)');
ylabel('Lactose and LA (g/dm^3)');
ylim([0 12])
ax=gca;
ax.YColor = 'k';

yyaxis right
plot(tspan,X,'-','color',r,'LineWidth',1.5);
ylabel('Biomass (g/dm^3)');
ylim([0 2.5])
ax.YColor = 'k';

%exp
hold on
yyaxis left
p(1)=plot(time,AE1,'s','color',g);
hold on
p(2)=plot(time,SE1,'+','color',b);
hold on
yyaxis right
p(3)=plot(time,XE1,'o','color',r,'MarkerFaceColor',r);
xlim([0 50])
le=legend(p,'LA','lactose','biomass');
set(le,'box','off')
%% 22g/l
X(1)=0.15 ;
A(1)=1.50 ; 
S(1)=21.00; 
Q(1)=Q0;
for i=1:imax
    funcQ=Q(i)./(1+Q(i));
    funcS=S(i)^nS./(KS^nS+S(i)^nS);
    f_A=1./(1+((A(i)/KI).^nA));
    funcV=1./(1+((S(i)/KQ).^nQ));
    X(i+1)=X(i)+dt.*mu_max.*X(i).*funcQ.*funcS.*f_A;
    Q(i+1)=Q(i)+dt.* funcV.*v.*Q(i);
    if Q(i+1)==inf
        Q(i+1)=Q(i);
    end
    A(i+1)=A(i)+alpha0.*(X(i+1)-X(i))+dt.*beta0.*X(i);
    S(i+1)=S(i)-(X(i+1)-X(i))./gama_XS- ...
        (A(i+1)-A(i))./gama_AS;
    if S(i+1)<0
        S(i+1)=0;
    end
end
%
figure
set(gca,'Position',[0.132,0.165,0.715,0.760])
yyaxis left
plot(tspan,A,'-','color',g,'LineWidth',1.5);
hold on
plot(tspan,S,'-','color',b,'LineWidth',1.5);
xlabel('Time (h)');
ylabel('Lactose and LA (g/dm^3)');
ylim([0 25])
ax=gca;
ax.YColor = 'k';

yyaxis right
plot(tspan,X,'-','color',r,'LineWidth',1.5);
ylabel('Biomass (g/dm^3)');
ylim([0 2.5])
ax.YColor = 'k';

%exp
hold on
yyaxis left
p(1)=plot(time,AE2,'s','color',g);
hold on
p(2)=plot(time,SE2,'+','color',b);
hold on
yyaxis right
p(3)=plot(time,XE2,'o','color',r,'MarkerFaceColor',r);
xlim([0 50])
le=legend(p,'LA','lactose','biomass');
set(le,'box','off')
%%
%
data=[mu_max
    alpha0
    beta0
    gama_XS
    gama_AS
    KI
    KS
    KQ
    Q0];
column_name=('parameter');
row_name={'μ_0';'alpha';'beta';'gama_XS';'gama_AS';'KI';'KS';'KQ';'Q0'};
%
set(figure,'position',[200 200 450 330]);
uitable(gcf,'Data',data,'Position',[20 20 400 315],'FontSize',15,'ColumnWidth', ...
    {120},'Columnname',column_name,'Rowname',row_name);
%%
% run('prediction_lactic_acid.m')