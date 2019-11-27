%
% type 0:
% spring force at start of step for whole step
% - causes too large impulse if start position
%   is near maximum and step is large
% type 1:
%
%
function v=ifs(k,mass,position,velocity,dt,w0,type)
  switch(type)
  case 0
    v=-k*position*dt;
  otherwise
  error(["type " num2str(type) " not supported in ifs"]);
  endswitch
endfunction
