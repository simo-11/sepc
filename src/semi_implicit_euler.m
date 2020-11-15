%
% horizontal spring position at start of step
t = 0.0;
dt = 0.1;
velocity = 0.0;
position = 1;
k=1;
force = 0;
mass = 1.0;
fprintf('%-12s%-12s%-12s%-12s\n','time','position','velocity','energy');
i=1;
e0=0.5*k*position*position;
clear v[tpve];
tmax=1.1*2*pi*sqrt(mass/k);
size=ceil(tmax/dt);
vt=zeros(size,1);
vp=zeros(size,1);
vv=zeros(size,1);
ve=zeros(size,1);
while t <= tmax
  vt(i)=t;
  vp(i)=position;
  vv(i)=velocity;
  e=0.5*k*position*position+0.5*mass*velocity*velocity;
  ve(i)=e;
  force=-k*position;
  velocity = velocity +( force / mass ) * dt;
  position = position + velocity * dt;
  fprintf('%-12.3f%-12.3f%-12.3f%-12.3f\n',t,position,velocity,e);
  i=i+1;
  t = t + dt;
 end 
 clf
 subplot(1,3,1);
 plot(vt,vp,'k');
 title(["position [m]  " datestr(now())]);
 subplot(1,3,2);
 plot(vt,vv,'r');
 title('velocity [m/s]');
 subplot(1,3,3);
 plot(vt,ve,'b');
 title('energy [J]');
