function [Cost] = F1(R)
%%  MOP1 function (F1)
f1 = R(:,1);
m = size(R,2);
g = 1+9*sum(R(:,2:end),2)/(m-1);
h = 1-sqrt(f1./g);
f2 = g.*h;
Tx = [f1 f2];
Cost =Tx';
end