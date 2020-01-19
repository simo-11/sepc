%
% vertical energy based v. 1
% gravity and potential energy are taken into account if g is set
% 
t = 0.0;
dt = 0.1;
k=1;
g=10;
mass = 1.0;
pg0=-mass*g/k
v0=sqrt(0.01/2);
p0=pg0+sqrt(0.01/2);
v0=-10;
p0=0;
i=1;
[e0,ee,ek,ep] = energy(p0,v0,mass,k,g)
[e0_pg0,ee_pg0,ek_pg0,ep_pg0]=energy(pg0,0,mass,k,g)
ea=e0-e0_pg0
pa=sqrt(2*ea/k)
va=sqrt(2*ea/mass)
w0=sqrt(k/mass)
if(pa>eps)
  pf=asin((p0-pg0)/pa)
  if(v0<0)
    pf=pi-pf
  endif
else
  pf=0
endif  
% vectors for time, position, velocity, energy
clear v[tpve];
clear i[s];
% vectors for theoretical values for velocity, position, energy, spring impulse
clear t[vpes];
% energy from computation for elastic, kinetic and potential
clear ec[ekp];
% theoretical energies
clear et[ekp];
velocity=v0;
position=p0;
ifsm=1;
ifsmt={"-k*p*dt","limited -k*int(p*dt)","limited -k*p*dt"};
if(w0*dt>2 && ifsm==1)
  fprintf('Integration will be unstable w0*dt= %f, i.e.>2\n',w0*dt);
endif
ig=-mass*g*dt % impulse from gravity
fprintf('%-11s%-11s%-11s%-11s%-11s%-11s%-11s%-11s\n',
  'time','position','velocity','energy',
  'elastic-e','kinetic-e','gravity-e','ifs');
while t <= 1.1*2*pi/w0
  % store for printing results 
  vt(i)=t;
  vp(i)=position;
  tp(i)=pg0+sin(w0*t+pf)*pa;
  ts(i)=-k*pg0*dt-k*pa*(cos(w0*(t-dt/2)+pf)-cos(w0*(t+dt/2)+pf));
  tv(i)=cos(w0*t+pf)*va;
  vv(i)=velocity;
  [ve(i),ee,ek,ep]=energy(position,velocity,mass,k,g);
  ece(i)=ee;eck(i)=ek;ecp(i)=ep;
  [te(i),ete(i),etk(i),etp(i)]=energy(tp(i),tv(i),mass,k,g);
  Is=ifs(k,mass,position,velocity,dt,w0,pg0,pa,ifsm); 
  fprintf('%-11.4g%-11.4g%-11.4g%-11.4g%-11.4g%-11.4g%-11.4g%-11.4g\n',
    t,position,velocity,ve(i),ee,ek,ep,Is);
  % Is impulse from spring - hopefully better than -k*position*dt
  is(i)=Is;
  velocity += (Is+ig)/mass;
  position += velocity * dt;
  i++;
  t += dt;
 end 
 clf
 subplot(2,4,1);
 plot(vt,vp,'k',vt,tp,'r');
 title(["sie: p [m]  "
   datestr(now())]);
 subplot(2,4,2);
 plot(vt,vv,'k',vt,tv,'r');
 title(['v [m/s], k=' num2str(k) ', m=' num2str(mass) ',' 
 ' dt=' num2str(dt)  ', g=' num2str(g)]);
 subplot(2,4,3);
 plot(vt,ve,'k',vt,te,'r');
 title(['e [J] - \omega\Delta_t=' num2str(w0*dt)
 '\omega=' num2str(w0) ]);
 subplot(2,4,4);
 plot(vt,is,'k',vt,ts,'r');
 title(['ifs [Ns] = ' ifsmt{ifsm}]);
 subplot(2,4,5);
 plot(vt,ece,'k', vt,ete,'r', vt,ece-ete,'b');
 title(['elastic energy, diff in blue [J]']);
 subplot(2,4,6);
 plot(vt,eck,'k', vt,etk,'r', vt,eck-etk,'b');
 title(['kinetic energy, diff in blue [J]']);
 subplot(2,4,7);
 plot(vt,ecp,'k', vt,etp,'r', vt,ecp-etp,'b');
 title(['potential energy, diff in blue[J]']);
