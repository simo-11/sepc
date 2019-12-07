%
%
function i=ifs(k,m,p,v,dt,w0,pg0,pa,type)
  switch(type)
  case 1
% spring force at start of step for whole step
% - causes too large impulse if start position
%   is near maximum and step is large    
    i=-k*p*dt;
  case 2
% integrate -k*position
% position is taken as pg0+sin(x+pf)*pa
% int(sin(x))=-cos(x)
% pf (phase fix) so that initial position matches
% 
    pf=asin((p-pg0)/pa);
    i0=pg0*dt;
    ie=-cos(w0*dt+pf);
    is=-cos(pf);
    id=ie-is;
    i=-k*(i0+pa*id);
  otherwise
    error(["type " num2str(type) " not supported in ifs"]);
  endswitch
endfunction
