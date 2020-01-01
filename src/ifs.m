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
    dp=p-pg0;
    pf=real(asin(dp/pa));
    i0=pg0*dt;
    ie=-cos(w0*dt+pf);
    is=-cos(pf);
    id=ie-is;
    i=-k*(i0+pa*id);
    if(w0*dt>1.57)
      ii=i;
      imin=(pg0-pa-p-v*dt)*m/dt;
      imax=(pg0+pa-p-v*dt)*m/dt;
      if(i<imin) 
        i=imin;
      elseif(i>imax)
        i=imax;
      endif
      if(i!=ii)
        ii,imax,imin,i
      endif
    endif
   case 3
% spring force at start of step for whole step
% but check limit
    i=-k*p*dt;
    if(w0*dt>1.98)
      ii=i;
      imin=(pg0-pa-p-v*dt)*m/dt;
      imax=(pg0+pa-p-v*dt)*m/dt;
      if(i<imin) 
        i=imin;
      elseif(i>imax)
        i=imax;
      endif
      if(i!=ii)
        ii,imax,imin,i
      endif
    endif   
  otherwise
    error(["type " num2str(type) " not supported in ifs"]);
  endswitch
endfunction
