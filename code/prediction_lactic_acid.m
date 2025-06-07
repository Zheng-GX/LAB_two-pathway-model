% close all;clear;clc
% set(0,'defaultAxesFontName', 'Arial','defaultAxesFontsize',14);
% set(0,'defaultTextFontName', 'Arial','defaultTextFontsize',14);
% set(0,'defaultfigurecolor','w')   
% set(0,'defaultLegendInterpreter','tex');
% set(0,'defaultTextInterpreter','tex');
% set(0,'defaultAxeslinewidth',1)
% set(0,'defaultfigureposition',[500,300,450,320]);
%%
time=[0.00 	2.00 	4.00 	6.00 	8.00 	15.00 	17.00 	19.00 	21.00 	23.00 	25.00 	27.00 	29.00 	31.00 	39.00 	41.00 	43.00 	45.00 	47.00 	49.00 ];
L11=[1.11 	1.23 	1.80 	2.45 	3.46 	9.31 	10.17 	10.55 	10.47 	10.63 	10.97 	10.65 	10.59 	10.75 	11.00 	11.02 	10.55 	10.84 	10.51 	11.05 ];
L22=[2.33 	3.07 	3.14 	3.55 	5.34 	9.96 	10.67 	10.83 	11.84 	12.60 	12.54 	13.02 	13.70 	14.46 	15.10 	14.95 	15.84 	15.18 	15.53 	15.63 ];
L33=[2.19 	2.53 	3.40 	4.55 	5.42 	10.40 	10.84 	11.58 	12.41 	12.72 	13.50 	13.08 	13.89 	13.70 	14.95 	14.99 	15.66 	15.60 	15.83 	15.96 ];
L44=[2.30 	2.52 	3.40 	4.81 	4.84 	9.37 	10.96 	10.92 	11.85 	12.36 	13.08 	13.95 	14.47 	13.85 	16.21 	16.35 	16.18 	16.54 	16.58 	16.52 ];
L55=[2.85 	2.96 	3.72 	4.80 	5.27 	9.31 	11.27 	11.38 	12.26 	12.85 	12.26 	13.47 	15.33 	15.29 	15.73 	16.89 	15.95 	17.64 	17.55 	17.39 ];
%%
% nA      = 5   ;
% nQ      = 2   ;
% Q0      = 2 ;
% KS      = 8;
% t       = 50;
% dt      = 0.01;
% %%%%%%%%%%%%
% var_names = {'mu_max','alpha0', 'beta0',  'gama_XS',  'gama_AS','KQ','KI'};
% pfit=[0.482388638515921
% 6.11055284325908
% 4.37201185678646e-14
% 1.08054051172755
% 1.34818570136841
% 6.04023916000857
% 7.16685754853557];
% for i = 1:length(pfit)
%     eval([var_names{i} ' = pfit(i);']);
% end
% %%%%%%%%%%%%
% v=mu_max;
% imax=t/dt;
% tspan=0:dt:t;
% X=zeros(1,imax+1);
% Q=zeros(1,imax+1);
% A=zeros(1,imax+1);
% S=zeros(1,imax+1);
%%
X(1)=0.15 ;
Q(1)=Q0;
A0=[1.12 2.19 2.3 2.3 2.85];
S0=[11 22 33 44 55];
alpha_0=[alpha0 alpha0 alpha0 alpha0 alpha0];
% alpha_0=[alpha0 alpha0 7 10 20];
for j=1:length(S0)
    S(1)=S0(j);
    A(1)=A0(j) ;
    alpha0=alpha_0(j);
    for i=1:imax
        %
        funcQ=Q(i)./(1+Q(i));
        funcS=S(i)./(KS+S(i));
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
lactic_acid(j,:)=A;
end
%%
%
colors = [
    0.16  0.44  0.56 
    0.30  0.64  0.61  
    0.70  0.80  0.38  
    0.94  0.72  0.26  
    0.86  0.37  0.34  
];
figure;
hold on
for i = 1:5
    plot(tspan, lactic_acid(i,:), 'Color', colors(i,:), 'LineWidth', 1.5);
end
box on
ylim([0 20])
xlabel('Time (h)');
ylabel('Lactic acid (g/dm^3)');
set(gca, 'TickDir', 'in'); % 
%%%%%
hold on
L=[L11;L22;L33;L44;L55];
M={'o','^','d','s','*'};
for i=1:5
scatter(time,L(i,:),'MarkerEdgeColor',colors(i,:),'Marker',M{i})
end
legend({'11 g/l', '22 g/l', '33 g/l', '44 g/l', '55 g/l', ...
    '11 g/l', '22 g/l','33 g/l', '44 g/l','55 g/l'},'NumColumns', 2);
set(gca,'Position',[0.132,0.165,0.715,0.760])