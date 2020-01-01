%
% vertical energy based v. 1
% gravity and potential energy are taken into account if g is set
% 
t = 0.0;
dt = 2.1;
k=1;
g=0;
mass = 1.0;
v0=sqrt(0.01/2);
p0=sqrt(0.01/2);
pg0=-mass*g/k;
fprintf('%-11s%-11s%-11s%-11s%-11s%-11s%-11s%-11s\n',
  'time','position','velocity','energy',
  'elastic-e','kinetic-e','potential-e','ifs');
i=1;
[e0,ee,ek,ep] = energy(p0,v0,mass,k,g);
pa=sqrt(2*(ee+ek)/k)
va=sqrt(2*(ee+ek)/mass);
w0=sqrt(k/mass);
pf=asin((p0-pg0)/pa);
clear v[tpvem];
clear i[s];
clear t[vpes];
velocity=v0;
position=p0;
ifsm=3;
ifsmt={"-k*p*dt","limited -k*int(p*dt)","limited -k*p*dt"};
if(w0*dt>2 && ifsm==1)
  fprintf('Integration will be unstable w0*dt= %f, i.e.>2\n',w0*dt);
endif
while t <= 3.1*2*pi/w0
  % store for printing results 
  vt(i)=t;
  vp(i)=position;
  tp(i)=pg0+sin(w0*t+pf)*pa;
  ts(i)=-k*pg0*dt-k*pa*(cos(w0*(t-dt/2)+pf)-cos(w0*(t+dt/2)+pf));
  tv(i)=cos(w0*t+pf)*va;
  vv(i)=velocity;
  [ve(i),ee,ek,ep]=energy(position,velocity,mass,k,g);
  te(i)=energy(tp(i),tv(i),mass,k,g);
  Is=ifs(k,mass,position,velocity,dt,w0,pg0,pa,ifsm); 
  fprintf('%-11.3g%-11.3g%-11.3g%-11.3g%-11.3g%-11.3g%-11.3g%-11.3g\n',
    t,position,velocity,ve(i),ee,ek,ep,Is);
  % Is impulse from spring - better than -k*position*dt
  % - position at middle of timestep creates damping and is not stable
  %
  is(i)=Is;
  velocity += Is/mass;
  position += velocity * dt;
  i++;
  t += dt;
 end 
 clf
 subplot(1,4,1);
 plot(vt,vp,'k',vt,tp,'r');
 title(["sie: p [m]  " datestr(now())]);
 subplot(1,4,2);
 plot(vt,vv,'k',vt,tv,'r');
 title(['v [m/s], k=' num2str(k) ', m=' num2str(mass) ', dt=' num2str(dt)]);
 subplot(1,4,3);
 plot(vt,ve,'k',vt,te,'r');
 title(['e [J] - \omega\Delta_t=' num2str(w0*dt)]);
 subplot(1,4,4);
 plot(vt,is,'k',vt,ts,'r');
 title(['ifs [Ns] = ' ifsmt{ifsm}]);
