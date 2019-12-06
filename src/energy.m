function [e,ee,ek,ep] = energy(p,v,m,k,g)
ee=0.5*k*p*p;
ek=0.5*m*v*v;
ep=-m*g*p;
e=ee+ek+ep;
endfunction