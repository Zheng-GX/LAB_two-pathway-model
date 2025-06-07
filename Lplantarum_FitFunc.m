function ypred = Lplantarum_FitFunc(pfit,time)
%
global nA nQ Q0 dt t KS var_names nS abc

for i = 1:length(pfit)
    eval([var_names{i} ' = pfit(i);']);
end

X(1)=0.15 ;
A(1)=1.09 ; 
S(1) =10.47; 
Q(1)=Q0;
v=mu_max;
%
imax=t/dt;
%
T=time;
%
for n=1:2
    if n==2
        X(1)=0.15 ;
        A(1)=1.50 ;
        S(1)=21.00;
    end
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
    X_=[];A_=[];S_=[];
    for k=1:length(T)
        X_=[X_,X(round(T(k)/dt+1))];
        A_=[A_,A(round(T(k)/dt+1))];
        S_=[S_,S(round(T(k)/dt+1))];
    end
    Xc{n}=X_;
    Ac{n}=A_;
    Sc{n}=S_;
end
Xe=[Xc{1}';Xc{2}'] ;
Ae=[Ac{1}';Ac{2}'];
Se=[Sc{1}';Sc{2}'];
% Xe=[Xc{1}'] ;
% Ae=[Ac{1}'];
% Se=[Sc{1}'];
ypred=[ Xe*abc(1); Ae*abc(2); Se*abc(3) ];
end