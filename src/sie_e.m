%
% vertical energy based v. 1
% gravity and potential energy are taken into account
% 
t = 0.0;
dt = 0.1;
velocity = 0.0;
position = 0;
k=1;
g=-10;
force = 0;
mass = 1.0;
fprintf('%-12s%-12s%-12s%-12s%-12s%-12s%-12s%-12s\n',
  'time','position','ap','velocity','energy',
  'elastic-e','kinetic-e','potential-e');
i=1;
e0=0.5*k*position*position-mass*g*position;
clear v[tpve];
while t <= 1.1*2*pi*sqrt(mass/k)
  % store and print results 
  vt(i)=t;
  vp(i)=position;
  vv(i)=velocity;
  ee=0.5*k*position*position;
  ek=0.5*mass*velocity*velocity;
  ep=-mass*g*position;
  e=ee+ek+ep;
  ve(i)=e;
  ap=(position+velocity*dt/2);
  fprintf('%-12.3f%-12.3f%-12.3f%-12.3f%-12.3f%-12.3f%-12.3f%-12.3f\n',
    t,position,ap,velocity,e,ee,ek,ep);
  % Is impulse from spring - better than -k*position*dt
  % - position at middle of timestep creates damping and is not stable
  %
  Is=ifs(k,mass,position,velocity,dt,0); 
  Ig=mass*g*dt;
  I=Ig+Is;
  velocity += I/mass;
  position += velocity * dt;
  i++;
  t += dt;
 end 
 clf
 subplot(1,3,1);
 plot(vt,vp,'k');
 title(["e - position [m]  " datestr(now())]);
 subplot(1,3,2);
 plot(vt,vv,'r');
 title(['velocity [m/s] - k=' num2str(k)]);
 subplot(1,3,3);
 plot(vt,ve,'b');
 title('energy [J]');
