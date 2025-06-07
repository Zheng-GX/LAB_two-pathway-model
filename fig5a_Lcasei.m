close all;clear;clc
set(0,'defaultAxesFontName', 'Arial','defaultAxesFontsize',17);%
set(0,'defaultTextFontName', 'Arial','defaultTextFontsize',17);%
set(0,'defaultfigurecolor','w')   
set(0,'defaultAxeslinewidth',1.5)
%% 
%
time1= [0.0     0.5 	1.0 	1.5     2.0     2.5     3.0     3.5     4.0     4.5 ]';
X_E1 = [1.17 	1.34 	1.56 	1.56 	1.88 	1.97 	2.21 	2.67 	3.26 	3.49 ]';
S_E1 = [9 	8.11 	7.51 	6.91 	5.84 	4.79 	3.54 	2.40 	0.17 	0.03 ]';
A_E1 = [4.27 	4.40 	4.83 	5.49 	6.03 	6.73 	7.50 	8.60 	9.69 	9.49 ]';
%
time2= [0.0 	1.0 	2.0 	3.0 	4.0 	5.0 	6.0 	7.0     8.0 	9.0 ]' ;
X_E2 = [1.33 	1.64 	2.07 	2.37 	3.08 	3.76 	4.61 	5.07    4.80 	4.63 ]';
S_E2 = [21.4 	19.09 	16.58 	13.98 	10.51 	6.85 	3.31 	0.05    0.02 	0.05 ]';
A_E2 = [5.44 	5.83 	6.41 	8.22 	10.13 	12.06 	16.2 	18.29   19.18 	19.47]';
%
time3= [0.0 	1.0 	2.0 	3.0 	4.0 	5.0 	6.0 	7.0 	8.0 	9.0 	10.0 	11.0 ]';
X_E3 = [1.12 	1.45 	1.94 	2.25 	2.94 	3.63 	4.56 	5.39 	6.23 	6.83 	6.95 	7.06 ]';
S_E3 = [35.5 	34.44 	34.49 	29.78 	25.34 	21.53 	20.25 	16.84 	12.09 	6.35 	0.53 	0.16 ]';
A_E3 = [3.0 	5.12 	5.62 	6.33 	7.57 	9.51 	13.65 	17.39 	21.60 	25.03 	29.41 	31.65]';
%
time4= [0.0 	1.0 	2.0 	3.0 	4.0 	5.0 	6.0 	7.0 	8.0 	9.0 	10.0 	11.0 	12.0 	13.0 	14.0 ]';
X_E4 = [0.96 	1.22 	1.55 	1.95 	2.53 	3.17 	3.63 	4.30 	4.96 	5.89 	6.08 	6.45 	6.48 	6.66 	6.84  ]';
S_E4 = [48.1 	47.05 	46.88 	43.62 	40.43 	36.26 	33.93 	30.10 	25.95 	19.88 	15.71 	10.36 	8.11 	6.09 	3.64 ]';
A_E4 = [6.0 	6.54 	7.88 	9.20 	10.90 	12.89 	15.86 	20.80 	22.97 	26.95 	29.16 	26.95 	29.70 	35.00 	36.65 ]';
%
time5= [0.0 	1.0 	2.0 	3.0 	4.0 	5.2 	6.0 	7.0 	8.0 	9.0 	10.4 	11.0 ]';
X_E5 = [0.93 	1.18 	1.49 	2.00 	2.54 	3.29 	3.83 	4.51 	5.00 	5.37 	6.14 	6.15  ]';
S_E5 = [61.2 	59.54 	58.59 	54.52 	47.33 	46.96 	43.85 	33.53 	33.01 	31.13 	24.32 	23.57 ]';
A_E5 = [5.9 	6.92 	7.92 	9.85 	10.68 	14.18 	16.66 	16.31 	21.44 	27.97 	28.13 	32.98 ]';
%
time6= [0.0 	1.0 	2.0 	3.0 	4.0 	5.0 	6.0 	7.0 	8.0 	9.0 	10.0 	11.0 	12.0 	13.0 	14.0 ]';
X_E6 = [1.10 	1.29 	1.57 	2.30 	2.83 	3.59 	3.81 	4.43 	4.85 	5.43 	5.91 	6.12 	6.66 	6.92 	6.94 ]';
S_E6 = [77.1 	70.05 	69.79 	70.08 	65.80 	63.35 	60.24 	54.21 	45.62 	38.97 	34.65 	30.89 	26.94 	23.06 	20.65 ]';
A_E6 = [7.89 	7.85 	9.79 	11.91 	13.51 	16.76 	20.05 	23.41 	26.36 	28.92 	33.08 	35.89 	37.64 	40.46 	44.91 ]';
%
X_E  = [X_E3 ;X_E4 ];
S_E  = [S_E3 ;S_E4 ];
A_E  = [A_E3 ;A_E4 ];
time = [time3;time4];
%
global L1 L2 m nS nQ K_S K_Q Q0 dt t
L1=length(time3);L2=length(time4)+L1;
%%
%
p0 =[1.75543842251338
24.3641340477288
0.245168198927857
2.71055352508494
31420.8928852871
0.696543588186903];
m       = 5   ;
nS      = 1   ;
nQ      = 2   ;
K_S     = 3; 
K_Q     = K_S;
Q0      = 0.15;
t       = 25;
dt      = 0.01;
%%%%%%%%
y_data = [X_E; A_E; S_E];
y_mean = mean(y_data);
SST = sum((y_data - y_mean).^2);% 
%%%%%%%% 
options = optimset('tolx',1e-16);
pfit=p0;
p00=zeros(length(p0),1);
while ~all(round(p00, 4) == round(p0, 4), 'all')
    p00=p0;
    [pfit,resnorm] = lsqcurvefit(@Lcasei_FitFunc,p0,time,y_data,zeros(length(p0),1), [],options);
    p0=pfit;
    resnorm
    % 
    R_squared = 1 - resnorm / SST;
    R_squared
end
%%
r=[0.91,0.3,0.06];
b=[0,0.39,0.65];
g=[0,0,0];
%
%
mu_max  = pfit(1); 
K_I     = pfit(2); 
alpha   = pfit(3);
beta    = pfit(4);
gama_XS = pfit(5);
gama_AS = pfit(6);

X0(1)   = 1.17  ;S0(1)   = S_E1(1);      A0(1)   = 4.27  ;           
X0(2)   = 1.33  ;S0(2)   = S_E2(1);      A0(2)   = 5.44  ;
X0(3)   = 1.12  ;S0(3)   = S_E3(1);     A0(3)   = 3.0   ;  
X0(4)   = 0.96  ;S0(4)   = S_E4(1);      A0(4)   = 6.0   ;  
X0(5)   = 0.93  ;S0(5)   = S_E5(1);      A0(5)   = 5.9   ;  
X0(6)   = 1.10  ;S0(6)   = S_E6(1);      A0(6)   = 7.89  ;
%
v=mu_max;
%
imax=t/dt;
tspan=0:dt:t;
X=zeros(1,imax+1);
Q=zeros(1,imax+1);
A=zeros(1,imax+1);
S=zeros(1,imax+1);
for j=1:6
    X(1)=X0(j) ; A(1)=A0(j) ; S(1) =S0(j) ; 
    Q(1)=Q0;
    %
    T{1}=time1;X_plot{1}=X_E1;A_plot{1}=A_E1;S_plot{1}=S_E1;
    T{2}=time2;X_plot{2}=X_E2;A_plot{2}=A_E2;S_plot{2}=S_E2;
    T{3}=time3;X_plot{3}=X_E3;A_plot{3}=A_E3;S_plot{3}=S_E3;
    T{4}=time4;X_plot{4}=X_E4;A_plot{4}=A_E4;S_plot{4}=S_E4;
    T{5}=time5;X_plot{5}=X_E5;A_plot{5}=A_E5;S_plot{5}=S_E5;
    T{6}=time6;X_plot{6}=X_E6;A_plot{6}=A_E6;S_plot{6}=S_E6;

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
        A(i+1)=A(i)+dt.*alpha.*X(i)+beta.*(X(i+1)-X(i));
        S(i+1)=S(i)-(X(i+1)-X(i))./gama_XS- ...
            (A(i+1)-A(i))./gama_AS;
        if S(i+1)<0
            S(i+1)=0;
        end
    end
    %
    tspan=0:dt:t;
    %
    a=[X',A',S'];
    %
    subplot(2,3,j)
    % subplot(1,3,j)
    %
    yyaxis left
    plot(tspan,a(:,2),'-','color',g,'LineWidth',1.5)
    hold on
    plot(tspan,a(:,3),'-','color',b,'LineWidth',1.5)
    xlabel('Time (h)');
    ylabel({'Lactose and lactic acid','(g/dm^3)'} );
    %
    xlim([0 4.5])
    ylim([0 11])
    if j==2
        xlim([0 9])
        ylim([0 25])
    end
    if j==3
        xlim([0 11])
        ylim([0 45])
    end
    if j==4
        xlim([0 14])
        ylim([0 60])
    end
    if j==5
        xlim([0 12])
        ylim([0 70])
    end
    if j==6
        xlim([0 15])
        ylim([0 90])
    end
    %
    yyaxis right
    plot(tspan,a(:,1),'-','color',r,'LineWidth',1.5)
    ylabel('Biomass (g/dm^3)');
    %）
    ylim([0 8])
    %
    hold on
    yyaxis left
    plot(T{j},A_plot{j},'s','color',g)
    hold on
    plot(T{j},S_plot{j},'+','color',b)
    yyaxis right
    plot(T{j},X_plot{j},'o','color',r,'MarkerFaceColor',r)
    set(0,'defaultfigurecolor','w')
end
%%
data_row(1)=mu_max;
data_row(2)=K_I    ;
data_row(3)=alpha  ;
data_row(4)=beta   ;
data_row(5)=gama_XS;
data_row(6)=gama_AS;
data_row(7)=K_S    ;
data_row(8)=Q0    ;
data=data_row';

%
m=length(data_row);n=1;
column_name=('parameter');
row_name=([' μ_max ';'  K_I  ';' alpha ';' beta  ';'gama_XS';'gama_AS';'  K_S  ';'  Q_0  ']);
% 
set(figure(2),'position',[200 200 450 330]);
uitable(gcf,'Data',data,'Position',[20 20 400 315],'FontSize',15,'ColumnWidth', ...
        {120,120,120,120,120,120,120,120}, ...
       'Columnname',column_name,'Rowname',row_name);
%%
X=zeros(1,imax+1);
A=zeros(1,imax+1);
S=zeros(1,imax+1);
Q=zeros(1,imax+1);
X(1)=1;A(1)=5;Q(1)=Q0;

a=1;
S_0=0.5:0.1:100;
for i_0=S_0
    S(1)=i_0;
    for i=1:imax
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
        A(i+1)=A(i)+dt.*alpha.*X(i)+beta.*(X(i+1)-X(i));
        S(i+1)=S(i)-(X(i+1)-X(i))./gama_XS- ...
            (A(i+1)-A(i))./gama_AS;
        if S(i+1)<0
            S(i+1)=0;
        end
        %
        mu(i)=(log(X(i+1))-log(X(i)))./dt;
    end
    [MU ,TIME]=max(mu);
    Mu_max(a)=max(mu);
    a=a+1;
end
%%
[aa bb]=max(Mu_max);
S_0(bb);
figure
set(gca,'Position',[0.21,0.24,0.71,0.73])
plot(S_0,Mu_max,'k-', 'LineWidth',1.5);
xlabel('S_0 (g/dm^3)','FontSize',17);ylabel('Specific growth rate (h^-^1)','FontSize',17);
xlim([0 80])
ylim([0.15 0.31])
xticks(0:40:80)
yticks(0.1:0.1:0.3)
hold on
S_E  =[9	21.4	35.5	48.1	61.2	77.1];
Mu_E=[0.27875	0.22059	0.23062	0.23815	0.23276	0.22142 ];
plot(S_E,Mu_E,'o','Color','r','MarkerFaceColor','r')
h1=legend('Model predition','Experiental data');
set(h1,'Box','off')