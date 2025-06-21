function ypred = Lbulgaricus_FitFunc(pfit,time)

%
global m nS nQ K_S K_Q Q0 dt t
miu_max = pfit(1);
K_I     = pfit(2);
alpha   = pfit(3);
beta    = pfit(4);
gama_XS = pfit(5);
gama_AS = pfit(6);

Q(1)=Q0;v=miu_max;
X(1)=0.02;A(1)=0.19 ;S(1)=26.19 ;
%
imax=t/dt;
%
b=1;X_=[];A_=[];S_=[];
%
time1=time(1:15);
time2=time(16:38);
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
    a=time1/dt;

    if ismember(i,a+1)
        X_(b)=X(i);
        A_(b)=A(i);
        S_(b)=S(i);
        b=b+1;
    end

end

%
X(1)=0.02;A(1)=0 ;S(1)=43.01 ;

b=1;
for i=1:imax
    funcQ=Q(i)./(1+Q(i));
    funcS=((S(i)/K_S).^nS)./(1+((S(i)/K_S).^nS));
    f_A=((A(i)/K_I).^m)./(1+((A(i)/K_I).^m));
    funcV=1-((S(i)/K_Q).^nQ)./(1+((S(i)/K_Q).^nQ));
    %
    X(i+1)=X(i)+dt.*miu_max.*X(i).*funcQ.*funcS.*(1-f_A);
    Q(i+1)=Q(i)+dt.*funcV.*v.*Q(i);
    if Q(i+1)==inf
        Q(i+1)=Q(i);
    end
    A(i+1)=A(i)+dt.*alpha.*X(i)+beta.*(X(i+1)-X(i));
    S(i+1)=S(i)-(X(i+1)-X(i))./gama_XS- ...
        (A(i+1)-A(i))./gama_AS;
    if S(i+1)<0
        S(i+1)=0;
    end
    a=time2/dt;


    if ismember(i,a+1)
        X_1(b)=X(i);
        A_1(b)=A(i);
        S_1(b)=S(i);
        b=b+1;
    end

end
Xe=[X_' ;X_1'];Ae=[A_' ; A_1'];Se=[S_' ; S_1'];
ypred=[ Xe, Ae, Se ];
end