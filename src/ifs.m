%
% 
function v=ifs(k,mass,position,velocity,dt,type)
  switch(type)
  case 0
    v=-k*position*dt;
  otherwise
  error(["type " num2str(type) " not supported in ifs"]);
  endswitch
endfunction
