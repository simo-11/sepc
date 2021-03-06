%
% vertical spring position at start of step
% gravity and potential energy are taken into account
t = 0.0;
dt = 0.1;
velocity = 0.0;
position = 0;
k=1;
g=-10;
force = 0;
mass = 1.0;
fprintf('%-12s%-12s%-12s%-12s\n','time','position','velocity','energy');
i=1;
e0=0.5*k*position*position-mass*g*position;
clear v[tpve];
while t <= 1.1*2*pi*sqrt(mass/k)
  vt(i)=t; %#ok<*SAGROW>
  vp(i)=position;
  vv(i)=velocity;
  e=0.5*k*position*position+0.5*mass*velocity*velocity-mass*g*position;
  ve(i)=e;
  force=-k*position+mass*g;
  velocity = velocity +( force / mass ) * dt;
  position = position+ velocity * dt;
  fprintf('%-12.3f%-12.3f%-12.3f%-12.3f\n',t,position,velocity,e);
  i=i+1;
  t = t+ dt;
 end 
 clf
 subplot(1,3,1);
 plot(vt,vp,'k');
 title(["wpe - position [m]  " datestr(now())]);
 subplot(1,3,2);
 plot(vt,vv,'r');
 title('velocity [m/s]');
 subplot(1,3,3);
 plot(vt,ve,'b');
 title('energy [J]');
