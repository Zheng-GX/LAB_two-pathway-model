close all;clear;clc
set(0,'defaultAxesFontName', 'Arial','defaultAxesFontsize',14);
set(0,'defaultTextFontName', 'Arial','defaultTextFontsize',14);
set(0,'defaultfigurecolor','w')   
set(0,'defaultLegendInterpreter','tex');
set(0,'defaultTextInterpreter','tex');
set(0,'defaultAxeslinewidth',1)
set(0,'defaultfigureposition',[500,300,450,320]);
%%
time=[0 	2 	4 	6 	8 	15 	17 	19 	21 	23 	24 ];
Exp_X=[0.15 	0.23 	0.38 	0.50 	0.67 	1.59 	1.73 	1.84 	1.94 	1.99 	2.00 ];
%%
t       = 50  ;
dt      = 0.01;
nA      = 5  ;
nQ      = 2   ;
KS      = 8   ;
Q0      = 2   ;
nS      = 1;
pfit =[0.489598365575430
6.14257649036447
3.78945025389706e-14
0.619351654370761
1.31281159251083
5.44557622047313
7.49465385946851];
var_names = {'mu_max','alpha0', 'beta0',  'gama_XS',  'gama_AS','KQ','KI'};
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
X0=[0.15 0.001 0.00001] ;
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
S(1)=11;
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
linestyle={'-',':','-.'};
for a=1:3
    plot(S0,Mumax{a},linestyle{a},'Linewidth',2)
    hold on
end
xlabel('S_0 (g/dm^3)')
ylabel('Specific growth rate (1/h)')
ylim([0.05 0.31])
num = X0(2:3);
coeff = num ./ 10.^floor(log10(num)); 
exponent = floor(log10(num));       
str{1} = sprintf('%.1f × 10^{%d}', coeff(1), exponent(1));
str{2} = sprintf('%.1f × 10^{%d}', coeff(2), exponent(2));
legend(['X_0=',num2str(X0(1)),' (g/dm^3)'],['X_0=',str{1},' (g/dm^3)'],['X_0=',str{2},' (g/dm^3)'])
%%
figure
linestyle={'-',':','-.'};
for a=1:3
    plot(tspan,x{a},linestyle{a},'Linewidth',2)
    hold on
end
p=scatter(time,Exp_X,'filled','k');
xlabel('Time (h)')
ylabel('Biomass conentration (g/dm^3)')
ylim([0 2.2])
legend(p,'Exp.data')
% legend(['X_0=',num2str(X0(1)),'g/dm^3'],['X_0=',str{1},'g/dm^3'],['X_0=',str{2},'g/dm^3'])