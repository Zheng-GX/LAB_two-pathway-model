close all;clear ;clc
%%
set(0,'defaultAxesFontName', 'Arial','defaultAxesFontsize',17);
set(0,'defaultTextFontName', 'Arial','defaultTextFontsize',17);
set(0,'defaultfigurecolor','w') 
set(0,'defaultAxeslinewidth',1.5)
set(0,'defaultTextInterpreter','tex');
set(0,'defaultLegendInterpreter','tex');
%%
%
X=zeros(1,1001);A=zeros(1,1001);S=zeros(1,1001);Q=zeros(1,1001);logX=zeros(1,1001);
Sum_S_X=zeros(1,1001);Sum_S_A=zeros(1,1001);
S_0=zeros(1,951);mu=zeros(1,1001);Miu_max=zeros(1,951);A_max=zeros(1,951);X_max=zeros(1,951);
% 
pfit=[2.172206652339938;21.074450083594467;0.460922666523526;10.564536590959085;0.335801086892000;1.190828228293918];
miu_max = pfit(1); 
K_I     = pfit(2); 
alpha   = pfit(3);
beta    = pfit(4);
gama_XS = pfit(5);
gama_AS = pfit(6);
%%%%%%%
m       = 5   ;
nS      = 1   ;
nQ      = 2   ;
K_S     = 8; 
K_Q     = K_S;
t       = 100;
dt      = 0.01;
v=miu_max;
%%%%%%%
imax=t/dt;
trspan=0:dt:t;
X(1)=0.02;A(1)=0  ; 
%%%%%%%
Q0=[2 1 0.6 0.2];
j_S0=2:0.1:70;
a=1;
%%%%%%
black=[0,0,0];
Red=[0.77,0.14,0.14];
Blue=[0.03,0.44,0.71];
Green=[0.2,0.4,0];
Orange=[0.86,0.52,0.01];
h=[Red;Orange;Green;Blue];  

for j_0=Q0
    Q(1)=j_0;
    e=1;
    for S0=j_S0
        S(1)=S0;
        for i=1:imax
            %
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
            mu(i)=(log(X(i+1))-log(X(i)))./dt;
        end
        Mumax(e)=max(mu);
        e=e+1;
    end
    hold on
    linestyle={'-',':','-.','--'};
    p(a)=plot(j_S0,Mumax,linestyle{a},'Color',h(a,:),'Linewidth',3);
    a=a+1;
end
%%
hold on
range=0:0.1:50;
plot(range,0.535*ones(length(range)),'k--')
hold on
size=80;
color1=[51,255,255]/255;
color2=[31,20,121]/255;
scatter(3.1, 0.535,size,"filled","MarkerFaceColor",color1)
hold on
scatter(40, 0.535,size,"filled","MarkerFaceColor",color2)

xlabel('S_0 (g/dm^3)');
ylabel('Specific growth rate (1/h)');
xlim([0 70])
ylim([0.2 1.6])
yticks(0.2:0.4:1.4)
box on                                                     
text(5,0.46,'L','Color','k')
text(39,0.46,'H','Color','k')
h1=legend({['Q_0=',num2str(Q0(1))],['Q_0=',num2str(Q0(2))], ...
    ['Q_0=',num2str(Q0(3))],['Q_0=',num2str(Q0(4))]}, 'NumColumns', 2,'fontsize',15);  
set(h1,'box','off')
set(gcf,"Position",[717.8,417,503.2,391.2])
set(gca,'Position',[0.16,0.21,0.78,0.73])
print('figure.pdf', '-dpdf', '-vector', '-r300');
