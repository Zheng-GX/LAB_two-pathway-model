close all;clear ;clc
%%
set(0,'defaultAxesFontName', 'Arial','defaultAxesFontsize',17);%
set(0,'defaultTextFontName', 'Arial','defaultTextFontsize',17);%
set(0,'defaultfigurecolor','w')   
set(0,'defaultAxeslinewidth',1.5)
%%
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
Q0=0.6;
v=miu_max;
imax=t/dt;
%
r=[0.91,0.3,0.06];
y=[0.95,0.82,0.00];
g=[0.2,0.8,0.14];
b=[0,0.39,0.65];
%
X=zeros(1,imax+1);A=zeros(1,imax+1);S=zeros(1,imax+1);Q=zeros(1,imax+1);
mu=zeros(1,imax+1);
%
A(1)=0  ;
X(1)=0.02;
Q(1)=Q0;
S_0=2:0.1:70;  
a=1;
for i_0=S_0
    S(1)=i_0;
    %
    e=1;
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

    [Mu_max(a),time(a)]=max(mu);
    growth_yield(a)=(X(imax)-X(1))/S(1);
    A_max(a)=A(imax);
    a=a+1;
end
[aa, bb]=max(Mu_max);
S_0(bb)
TIME=(time-1)*dt;
%%
% 1. λ-μ
figure
size=80;
color1=[26,121,20]/255;
scatter(Mu_max(bb),1./TIME(bb),size,"filled","MarkerFaceColor",color1)
hold on 
plot(Mu_max(1:bb),1./TIME(1:bb),'.','Color',[0.00,0.45,0.74]);
hold on 
plot(Mu_max((bb+1):(a-1)),1./TIME((bb+1):(a-1)),'.r')
xlabel('\mu (1/h)');
ylabel('1/\lambda (1/h)');
set(gca,'Position',[0.18,0.21,0.78,0.73])
xlim([0.2 1.2])
xticks(0.4:0.4:1.2)
ylim([0.15 0.6])
yticks(0.2:0.2:0.6)
box on
annotation('arrow', [0.431,0.616], [0.753,0.655]);
text(0.582,0.4,'$S_0$','Interpreter','latex')
text(0.5,0.23,'$S_{c}=14.3\,\rm g/dm^3$','Interpreter','latex')
set(gcf,"Position",[717.8,417,503.2,391.2])
set(gca,'Position',[0.16,0.21,0.78,0.73])
%%%%%%%%%
figure
scatter(1./Mu_max(bb),1./growth_yield(bb),size,"filled","MarkerFaceColor",color1)
hold on
plot(1./Mu_max(1:bb),1./growth_yield(1:bb),'.','Color',[0.00,0.45,0.74]);
hold on 
plot(1./Mu_max((bb+1):(a-1)),1./growth_yield((bb+1):(a-1)),'.r')
ylabel('1/Y (S_0/X)');
xlabel('1/\mu (h)');
set(gca,'Position',[0.18,0.21,0.78,0.73])
xlim([0.8 2.7])
ylim([11 19.5])
xticks(1:1:2)
yticks(12:3:18)
box on
annotation('arrow', [0.77,0.57], [0.45,0.4]);
text(1.83,15,'$S_0$','Interpreter','latex')
set(gcf,"Position",[717.8,417,503.2,391.2])
set(gca,'Position',[0.16,0.21,0.78,0.73])