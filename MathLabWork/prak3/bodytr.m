%{
\dot x1 =x2
\dot x2=-(c*x1 + k*x2)/m

Fтр = 
{- mu*N*sign(\dot x1-v),\dot x1 nq v
R*sign(\dot x1 - v),\dot x1 = v
}
%}

function varargout=bodytr(t,x,flag)
    switch flag
    case ''
        varargout{1}=ft(t,x);
    case 'events'
        [varargout{1:3}]=Tr_event(t,x);
    end
end

function f=ft(t,x)
global  m c k mu g sgn v
f(1,1)=x(2);
f(2,1)=-(c*x(1) + k*x(2))/m-mu*g*sgn;
end

function [val,ist,dir]=Tr_event(t,x)
global v
val = x(2)-v;
ist=1;
dir=0;
end
