close all;clear;clc
set(0,'defaultAxesFontName', 'Arial','defaultAxesFontsize',14);
set(0,'defaultTextFontName', 'Arial','defaultTextFontsize',14);
set(0,'defaultfigurecolor','w')   
set(0,'defaultLegendInterpreter','tex');
set(0,'defaultTextInterpreter','tex');
set(0,'defaultAxeslinewidth',1)
set(0,'defaultfigureposition',[500,300,450,320]);
%%
time=[0.00 	2.00 	4.00 	6.00 	8.00 	15.00 	17.00 	19.00 	21.00 	23.00 	25.00 	27.00 	29.00 	31.00 	39.00 	41.00 	43.00 	45.00 	47.00 	49.00 ];
Exp_X=[0.15 	0.22 	0.32 	0.54 	0.68 	1.15 	1.48 	1.53 	1.47 	1.78 	1.81 	2.07 	2.06 	2.04 	2.15 	2.07 	2.26 	2.22 	2.22 	2.29 ];
%%
t       = 50  ;
dt      = 0.01;
KS      = 8   ;
Q0      = 1  ;
nA      = 5  ;
nQ      = 2 ;
nS      = 1;
pfit =[0.767143092904400
6.11693342834676
0
8.06781541959975
1.06228092342874
5.14964422722209
6.84035358009619
8.35484934214176
0.596221409911512];
var_names = {'mu_max','alpha0', 'beta0',  'gama_XS',  'gama_AS','KQ','KI','KS','Q0'};
for i = 1:length(pfit)
    eval([var_names{i} ' = pfit(i);']);
end
%%
%
r=[0.91,0.3,0.06];
b=[0,0.39,0.65];
g=[0,0,0];
%
v=mu_max;
%
imax=t/dt;
tspan=0:dt:t;
X=zeros(1,imax+1);
Q=zeros(1,imax+1);
A=zeros(1,imax+1);
S=zeros(1,imax+1);
X0=[0.15 0.015 0.0015] ;
A(1)=1.5 ; 
Q(1)=Q0;
S0=1:0.1:50;
%% S0-mu with dif X0
for k=1:length(X0)
    X(1)=X0(k) ;
    for j=1:length(S0)
        S(1)=S0(j);
        for i=1:imax
            %
            funcQ=Q(i)./(1+Q(i));
            funcS=S(i)^nS./(KS^nS+S(i)^nS);
            f_A=1./(1+((A(i)/KI).^nA));
            funcV=1./(1+((S(i)/KQ).^nQ));
            %
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
            mu(i)=(log(X(i+1))-log(X(i)))./dt;
        end
        Mumax{k}(j)=max(mu);
    end
end
%% time-X with dif X0
S(1)=22;
for k=1:length(X0)
    X(1)=X0(k) ;
    for i=1:imax
        %
        funcQ=Q(i)./(1+Q(i));
        funcS=S(i)^nS./(KS^nS+S(i)^nS);
        f_A=1./(1+((A(i)/KI).^nA));
        funcV=1./(1+((S(i)/KQ).^nQ));
        %
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
x{k}=X;
end
%%
figure
set(gca,'Position',[0.16,0.21,0.78,0.73])
linestyle={'-',':','-.'};
colors=[0.00,0.45,0.74
    0.93,0.69,0.13
    0.85,0.33,0.10
    ];
for a=1:3
    plot(S0,Mumax{a},linestyle{a},'color',colors(a,:),'Linewidth',2)
    hold on
end
xlabel('S_0 (g/dm^3)')
ylabel('Specific growth rate (1/h)')
% ylim([0.05 0.31])
ylim([0.05 0.38])
num = X0(2:3);
coeff = num ./ 10.^floor(log10(num)); 
exponent = floor(log10(num));       
legend(['X_0=',num2str(X0(1)),' g/dm^3'],['X_0=',num2str(X0(2)),' g/dm^3'],['X_0=',num2str(X0(3)),' g/dm^3'])
% str{1} = sprintf('%.1f × 10^{%d}', coeff(1), exponent(1));
% str{2} = sprintf('%.1f × 10^{%d}', coeff(2), exponent(2));
% legend(['X_0=',num2str(X0(1)),' (g/dm^3)'],['X_0=',str{1},' (g/dm^3)'],['X_0=',str{2},' (g/dm^3)'])
%%
figure
set(gca,'Position',[0.16,0.21,0.78,0.73])
linestyle={'-',':','-.'};
for a=1:3
    plot(tspan,x{a},linestyle{a},'color',colors(a,:),'Linewidth',2)
    hold on
end
p=scatter(time,Exp_X,'filled','k');
xlabel('Time (h)')
ylabel('Biomass conentration (g/dm^3)')
ylim([0 2.7])
h=legend(p,'0.15 g/dm^3 (exp.data)');
% set(h,'box','off')
% legend(['X_0=',num2str(X0(1)),'g/dm^3'],['X_0=',str{1},'g/dm^3'],['X_0=',str{2},'g/dm^3'])