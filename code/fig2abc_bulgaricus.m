clear ;close all;clc
%% fit
%exp
time1= [0       0.5 	1.0 	1.5 	2.0 	2.5 	3.0 	3.5     4.0 	4.5 	5.0 	5.5  ...
        6.0 	6.5 	7.0     ]';
time2= [0       0.5 	1.0 	1.5 	2.0 	2.5 	3.0 	3.5     4.0 	4.5 	5.0 	5.5  ...
        6.0 	6.5 	7.0 	7.5 	8.0 	8.5 	9.0 	9.5 	10.0 	10.5 	11.0   ]';
X_E1 = [0.01 	0.02 	0.02 	0.04 	0.05 	0.11 	0.17 	0.26 	0.44 	0.72 	1.12 	1.69 ...
        2.12 	2.50 	2.50   	]'; 
S_E1 = [26.19 	26.08 	25.94 	25.77 	25.58 	25.03 	24.51 	23.43 	21.30 	18.50 	14.90 	8.93 ...
        4.82 	1.00 	0.00    ]';
A_E1 = [0.19 	0.38 	0.48 	0.64 	0.83 	1.33 	1.84 	2.78 	4.63 	7.39 	10.48 	15.93...
        19.49 	22.99 	24.06   ]';
X_E2 = [0.01 	0.01 	0.02 	0.04 	0.06 	0.11 	0.18 	0.28 	0.45 	0.72 	1.18 	1.80 ...
        2.16 	2.53 	2.97 	3.07 	3.16 	3.17 	3.16 	3.16 	3.16 	3.16 	3.16   ]';
S_E2 = [43.01 	42.88 	42.71 	42.49 	42.31 	41.71 	41.18 	39.97 	38.26 	35.36 	31.01 	26.01...
        20.68 	15.82 	12.42 	10.61 	8.78 	6.71 	4.53 	2.91 	1.34 	1.14 	0.87   ]';
A_E2 = [0.00 	0.30 	0.31 	0.51 	0.79 	1.26 	1.90 	2.92 	4.46 	7.09 	11.06 	15.44...
        20.49 	24.70 	27.86 	29.57 	31.41 	33.50 	35.55 	37.14 	38.76 	38.98 	39.28  ]';
%
time = [time1;time2];X_E=[X_E1;X_E2];A_E=[A_E1;A_E2];S_E=[S_E1;S_E2];
%%
p0 =[2.172206652339938;21.074450083594467;0.460922666523526;10.564536590959085;0.335801086892000;1.190828228293918];
global m nS nQ K_S K_Q Q0 dt t
m       = 5   ;
nS      = 1   ;
nQ      = 2   ;
K_S     = 8; 
K_Q     = K_S;
Q0      = 0.6;
t       = 12;
dt      = 0.01;
%
options = optimset('tolx',1e-16);
% [pfit,resnorm] = lsqcurvefit(@t2_lactis_FitFunc,p0,time,[X_E,A_E,S_E],zeros(length(p0),1),[],options);
pfit=p0;
p00=zeros(length(p0),1);
while ~all(round(p00, 4) == round(p0, 4), 'all')
    p00=p0;
    [pfit,resnorm] = lsqcurvefit(@Lbulgaricus_FitFunc,p0,time,[X_E,A_E,S_E],zeros(length(p0),1),[],options);
    p0=pfit;
end
%%%%%%%%%%%%%%%%%
y_data = [X_E; A_E; S_E];
y_mean = mean(y_data);
SST = sum((y_data - y_mean).^2);% 
R_squared = 1 - resnorm / SST;
R_squared
%%
set(0,'defaultAxesFontName', 'Arial','defaultAxesFontsize',17);
set(0,'defaultTextFontName', 'Arial','defaultTextFontsize',17);
set(0,'defaultTextInterpreter','tex'); 
set(0,'defaultLegendInterpreter','tex');
set(0,'defaultfigurecolor','w')
set(0,'defaultAxeslinewidth',1.5)
set(0,'defaultfigureposition',[500,300,450,320]);
%配色
r=[0.91,0.3,0.06];
b=[0,0.39,0.65];
g=[0,0,0];
%作图
%参数设置
miu_max = pfit(1); 
K_I     = pfit(2); 
alpha   = pfit(3);
beta    = pfit(4);
gama_XS = pfit(5);
gama_AS = pfit(6);

X0(1)   = 0.02;                           
A0(1)   = 0.19;                         
S0(1)   = 26.19;
X0(2)   = 0.02;                            
A0(2)   = 0;
S0(2)   = 43.01;
v=miu_max;
imax=t/dt;
tspan=0:dt:t;
X=zeros(1,imax+1);
Q=zeros(1,imax+1);
A=zeros(1,imax+1);
S=zeros(1,imax+1);
%%%%%%%%%%%%%%%%
for j=1:2
    X(1)=X0(j) ; A(1)=A0(j) ; S(1) =S0(j) ; 
    Q(1)=Q0;
    T=time1;X_plot=X_E1;A_plot=A_E1;S_plot=S_E1;
    if j==2
        T=time2;X_plot=X_E2;A_plot=A_E2;S_plot=S_E2;
    end
    for i=1:imax
        funcQ=Q(i)./(1+Q(i));
        funcS=((S(i)/K_S).^nS)./(1+((S(i)/K_S).^nS));
        f_A=((A(i)/K_I).^m)./(1+((A(i)/K_I).^m));
        funcV=1-((S(i)/K_Q).^nQ)./(1+((S(i)/K_Q).^nQ));
        %
        X(i+1)=X(i)+dt.*miu_max.*X(i).*funcQ.*funcS.*(1-f_A);
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
    a=[X',A',S'];
    %
    figure
    hold on
    %
    yyaxis left
    plot(tspan,a(:,2),'-','color',g,'LineWidth',1.5)
    plot(tspan,a(:,3),'-','color',b,'LineWidth',1.5)
    xlabel('Time (h)');
    ylabel('Lactose and LA (g/dm^3)');
    yyaxis right
    plot(tspan,a(:,1),'-','color',r,'LineWidth',1.5)
    ylabel('Biomass (g/dm^3)');
    ylim([0 5])
    xlim([0 7])
    if j==2
        ylim([0 8])
        xlim([0 11])
    end
    %
    yyaxis left
    p(1)=plot(T,A_plot,'s','color',g);
    if j==2
        yticks(0:20:40)
    end
    p(2)=plot(T,S_plot,'+','color',b);
    yyaxis right
    p(3)=plot(T,X_plot,'o','color',r,'MarkerFaceColor',r);
    le=legend(p,'Lactic acid','Lactose','Biomass');
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'k';
    box on
    set(le,'box','off')
end
hold off
%%

data_row(1)=miu_max;
data_row(2)=K_I    ;
data_row(3)=alpha  ;
data_row(4)=beta   ;
data_row(5)=gama_XS;
data_row(6)=gama_AS;
data_row(7)=K_S    ;
data_row(8)=K_Q    ;
data_row(9)=Q0    ;
data=data_row';

%
column_name=('parameter');
row_name={'μ_max';'KI';'alpha';'beta';'gama_XS';'gama_AS';'KS';'KQ';'Q0'};
% 表格作图
set(figure,'position',[200 200 450 330]);
uitable(gcf,'Data',data,'Position',[20 20 400 315],'FontSize',15,'ColumnWidth', ...
        {120,120,120,120,120,120,120,120}, ...
       'Columnname',column_name,'Rowname',row_name);
%%
X=zeros(1,imax+1);
A=zeros(1,imax+1);
S=zeros(1,imax+1);
Q=zeros(1,imax+1);
X(1)=0.02;A(1)=0;Q(1)=Q0;

a=1;
S_0=1:0.1:100;
for i_0=S_0
    S(1)=i_0;
    for i=1:imax
        funcQ=Q(i)./(1+Q(i));
        funcS=((S(i)/K_S).^nS)./(1+((S(i)/K_S).^nS));
        f_A=((A(i)/K_I).^m)./(1+((A(i)/K_I).^m));
        funcV=1-((S(i)/K_Q).^nQ)./(1+((S(i)/K_Q).^nQ));
        %
        X(i+1)=X(i)+dt.*miu_max.*X(i).*funcQ.*funcS.*(1-f_A);
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
        miu(i)=(log(X(i+1))-log(X(i)))./dt;
    end
    [MIU ,TIME]=max(miu);
    Miu_max(a)=max(miu);
    a=a+1;
end
%%
[mu ,Time]=max(Miu_max);
S_0(Time);
figure
set(gca,'Position',[0.18,0.24,0.71,0.67])
plot(S_0,Miu_max,'k-', 'LineWidth',1.5);
xlabel('S_0 (g/dm^3)');
ylabel('Specific growth rate (h^{-1})');
axis([ 0 70 0.4 1.2])
set(gca,'Ytick',0.4:0.2:1.2)
hold on
S_E  =[5.19 	10.46 	20.09 	40.73 	61.16 ];
Miu_E=[0.69 	0.86 	0.99 	0.88 	0.76  ];
plot(S_E,Miu_E,'o','Color','k','MarkerFaceColor','k')
h1=legend('Model prediction','Experimental data');
set(h1,'box','off')
set(gca,'Xminortick','on')
set(gca,'Yminortick','on')