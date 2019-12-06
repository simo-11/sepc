%
%
function v=ifs(k,mass,position,velocity,dt,w0,vmax,pmax,pmin,type)
  switch(type)
  case 1
% spring force at start of step for whole step
% - causes too large impulse if start position
%   is near maximum and step is large    
    v=-k*position*dt;
case 2
% integrate -k*position
% position is taken as position+sin(x)*vmax*dt
% int(sin(x))=-cos(x)
% so that initial position and velocity match
% 
    v=-k*position*dt;
    if(vmax>sqrt(eps))
      xs=asin(velocity/vmax);
      xe=xs+dt*w0;
      is=-cos(xs);
      ie=-cos(xe);
      d=vmax*(ie-is);
      r=abs(d/v);
      v+=d;
    endif
  otherwise
  error(["type " num2str(type) " not supported in ifs"]);
  endswitch
endfunction
