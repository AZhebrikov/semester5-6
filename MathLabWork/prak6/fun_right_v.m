
function varargout = fun_right_v(t,x,flag)
    switch flag
    case ''
        varargout{1}=f_right(t,x);
    case 'events'
        [varargout{1:3}] = Tr_events(t,x);
    end
end

function f_out = f_right(t,x)
global u
f_out(1,1)=x(2);
f_out(2,1)=u + 0.01*x(1)*cos(t);
end

function f_out = f_switching_line(t,x)
    if x(2)>=0
        f_out=x(1)+x(2)^2/2;
    else
        f_out=x(1)-x(2)^2/4;
    end
end

function [val,ist,dir]=Tr_events(t,x)
    val = f_switching_line(t,x);
    ist=1;
    dir=0;
end