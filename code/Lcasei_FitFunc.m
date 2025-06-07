function ypred = Lcasei_FitFunc(pfit,time)
%
global L1 L2  m nS nQ K_S K_Q Q0 dt t
mu_max  = pfit(1);
K_I     = pfit(2);
alpha   = pfit(3);
beta    = pfit(4);
gama_XS = pfit(5);
gama_AS = pfit(6);

Q(1)=Q0;v=mu_max;
%
imax=t/dt;

time3=time(1:L1);
time4=time(L1+1:L2);
T={time3,time4};
%
X0(1)   = 1.12  ;S0(1)   = 35.63;     A0(1)   = 3.0   ;
X0(2)   = 0.96  ;S0(2)   = 48.1;      A0(2)   = 6.0   ;
%
for n=1:2
    X(1)=X0(n);
    A(1)=A0(n);
    S(1)=S0(n);
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
    km=length(T{n});
    X_=[];A_=[];S_=[];
    for k=1:km
        X_=[X_,X(round(T{n}(k))/dt+1)];
        A_=[A_,A(round(T{n}(k))/dt+1)];
        S_=[S_,S(round(T{n}(k))/dt+1)];
    end
    Xc{n}=X_;
    Ac{n}=A_;
    Sc{n}=S_;
end
Xe=[[Xc{1}]';[Xc{2}]'];
Ae=[[Ac{1}]';[Ac{2}]'];
Se=[[Sc{1}]';[Sc{2}]'];
ypred=[ Xe; Ae; Se ];
end