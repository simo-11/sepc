%
% vertical energy based v. 1
% gravity and potential energy are taken into account if g is set
% 
t = 0.0;
dt = 0.1;
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
pa=sqrt(2*(ee+ek)/k);
vmax=sqrt(2*(ee+ek)/mass);
w0=sqrt(k/mass);
pf=asin(p0/w0);
clear v[tpvem];
clear i[s];
clear t[vp];
velocity=v0;
position=p0;
tvs=1;
vmax=velocity;
pmax=position;
pmin=position;
if(w0*dt>2)
  fprintf('Integration will be unstable w0*dt= %f, i.e.>2\n',w0*dt);
endif
ifsm=1;
ifsmt={"-k*p*dt","-k*int(p*dt)"};
while t <= 1.1*2*pi/w0
  vmax=max(vmax,velocity);
  pmax=max(pmax,position);
  pmin=min(pmin,position);
  % store for printing results 
  vt(i)=t;
  vp(i)=position;
  tp(i)=pg0+sin(w0*t+pf)*pa;
  tv(i)=cos(w0*t+pf)*vmax;
  vv(i)=velocity;
  ve(i)=energy(position,velocity,mass,k,g);
  te(i)=energy(tp(i),tv(i),mass,k,g);
  Is=ifs(k,mass,position,velocity,dt,w0,vmax,pmax,pmin,ifsm); 
  fprintf('%-11.3g%-11.3g%-11.3g%-11.3g%-11.3g%-11.3g%-11.3g%-11.3g\n',
    t,position,velocity,e,ee,ek,ep,Is);
  % Is impulse from spring - better than -k*position*dt
  % - position at middle of timestep creates damping and is not stable
  %
  is(i)=Is;
  Ig=mass*g*dt;
  I=Ig+Is;
  velocity += I/mass;
  position += velocity * dt;
  i++;
  t += dt;
 end 
 clf
 subplot(1,5,1);
 plot(vt,vp,'k',vt,tp,'r');
 title(["sie: p [m]  " datestr(now())]);
 subplot(1,5,2);
 plot(vt,vv,'k',vt,tv,'r');
 title(['v [m/s] - k=' num2str(k)]);
 subplot(1,5,3);
 plot(vt,ve,'k',vt,te,'r');
 title(['e [J] - \omega\Delta_t=' num2str(w0*dt)]);
 subplot(1,5,4);
 plot(vt,is,'b');
 title(['ifs [Ns] = ' ifsmt{ifsm}]);
% subplot(1,5,5);
% plot(vt,vm,'b');
% title(['v_{max} [m/s]']);
