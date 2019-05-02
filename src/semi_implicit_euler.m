t = 0.0;
dt = 0.06;
velocity = 0.0;
position = 1;
k=1000;
force = 0;
mass = 1.0;
fprintf('%-12s%-12s%-12s%-12s\n','time','velocity','position','energy');
i=1;
e0=0.5*k*position*position;
clear v[tpve];
while t <= 10.0
  e=0.5*k*position*position+0.5*mass*velocity*velocity;
  force=-k*position;  
  velocity +=( force / mass ) * dt;
  position += velocity * dt;
  fprintf('%-12.2f%-12.3f%-12.3f%-12.3f\n',t,velocity,position,e);
  vt(i)=t;
  vp(i)=position;
  vv(i)=velocity;
  ve(i)=e;
  i++;
  t += dt;
 end 
 clf
 plot(vt,vp,'k',vt,vv,'r',vt,ve,'b');