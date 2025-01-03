function dy = shapef(u,y,para,fp,ka,mu,Ga1,Ga2,uc,r0c,z0c,psi0c,h0,dt)
psi = y(1,:);
dpsi = y(2,:);
r = y(3,:);
z = y(4,:);
la = y(6,:);
al = y(7,:);
be = y(8,:);
vu = y(9,:);
r0 = spline(uc,r0c,u);
z0 = spline(uc,z0c,u);
psi0 = spline(uc,psi0c,u);

c0 = 0;
h = para;

fu = 0;
fn = fp;

dy(1,:)=dpsi;

dy(2,:)=(-1).*dpsi.*h.*r.^(-1).*cos(psi)+ be.*h.^2.*ka.^(-1).*r.^(-1).* ...
  cos(psi)+ al.*h.^2.*ka.^(-1).*r.^(-1).*sin(psi)+ h.^2.*r.^(-2).*cos( ...
  psi).*sin(psi);

dy(3,:)=h.*cos(psi);

dy(4,:)=(-1).*h.*sin(psi);

dy(5,:)=h.*r;

dy(6,:)=(-1).*fu.*h+ Ga2.*h.*vu+ (-2).*dpsi.*mu.*r.^(-1).*vu.*sin(psi)+ (-2) ...
  .*dpsi.*dt.^(-1).*mu.*r.^(-1).*z.*sin(psi).^2+ 2.*dpsi.*dt.^(-1).* ...
  mu.*r.^(-1).*z0.*sin(psi).^2+ dpsi.*dt.^(-1).*mu.*sin(2.*psi)+ (-1) ...
  .*dpsi.*dt.^(-1).*mu.*r.^(-1).*r0.*sin(2.*psi)+ (-2).*dt.^(-1).* ...
  h0.*mu.*r.^(-1).*sin(psi).*sin(psi+ (-1).*psi0);

dy(7,:)=(-1).*c0.*dpsi.*ka+ (1/2).*dpsi.^2.*h.^(-1).*ka+ (1/2).*c0.^2.*h.* ...
  ka+ (-1).*fn.*h.*r.*sin(psi)+ dpsi.*la.*r.*sin(psi)+ (-2).*dpsi.*mu.* ...
  vu.*cos(psi).*sin(psi)+ dt.^(-1).*Ga1.*h.*r.*z.*cos(psi).*sin(psi)+  ...
  (-1).*dt.^(-1).*Ga1.*h.*r.*z0.*cos(psi).*sin(psi)+ h.*la.*sin(psi) ...
  .^2+ (-1/2).*h.*ka.*r.^(-2).*sin(psi).^2+ dt.^(-1).*Ga1.*h.*r.^2.* ...
  sin(psi).^2+ (-1).*dt.^(-1).*Ga1.*h.*r.*r0.*sin(psi).^2+ 2.*h.*mu.* ...
  r.^(-1).*vu.*cos(psi).*sin(psi).^2+ (-2).*dpsi.*dt.^(-1).*mu.*z.* ...
  cos(psi).*sin(psi).^2+ 2.*dpsi.*dt.^(-1).*mu.*z0.*cos(psi).*sin( ...
  psi).^2+ (-2).*dpsi.*dt.^(-1).*mu.*r.*sin(psi).^3+ 2.*dpsi.*dt.^(-1) ...
  .*mu.*r0.*sin(psi).^3+ 2.*dt.^(-1).*h.*mu.*r.^(-1).*z.*cos(psi).* ...
  sin(psi).^3+ (-2).*dt.^(-1).*h.*mu.*r.^(-1).*z0.*cos(psi).*sin(psi) ...
  .^3+ 2.*dt.^(-1).*h.*mu.*sin(psi).^4+ (-2).*dt.^(-1).*h.*mu.*r.^(-1) ...
  .*r0.*sin(psi).^4;

dy(8,:)=(-1).*fn.*h.*r.*cos(psi)+ dpsi.*la.*r.*cos(psi)+ (-2).*dpsi.*mu.* ...
  vu.*cos(psi).^2+ dt.^(-1).*Ga1.*h.*r.*z.*cos(psi).^2+ (-1).*dt.^(-1) ...
  .*Ga1.*h.*r.*z0.*cos(psi).^2+ h.*la.*cos(psi).*sin(psi)+ dt.^(-1).* ...
  Ga1.*h.*r.^2.*cos(psi).*sin(psi)+ (-1).*dt.^(-1).*Ga1.*h.*r.*r0.* ...
  cos(psi).*sin(psi)+ 2.*h.*mu.*r.^(-1).*vu.*cos(psi).^2.*sin(psi)+ ( ...
  -2).*dpsi.*dt.^(-1).*mu.*z.*cos(psi).^2.*sin(psi)+ 2.*dpsi.*dt.^( ...
  -1).*mu.*z0.*cos(psi).^2.*sin(psi)+ (-2).*dpsi.*dt.^(-1).*mu.*r.* ...
  cos(psi).*sin(psi).^2+ 2.*dpsi.*dt.^(-1).*mu.*r0.*cos(psi).*sin( ...
  psi).^2+ 2.*dt.^(-1).*h.*mu.*r.^(-1).*z.*cos(psi).^2.*sin(psi).^2+ ( ...
  -2).*dt.^(-1).*h.*mu.*r.^(-1).*z0.*cos(psi).^2.*sin(psi).^2+ 2.* ...
  dt.^(-1).*h.*mu.*cos(psi).*sin(psi).^3+ (-2).*dt.^(-1).*h.*mu.*r.^( ...
  -1).*r0.*cos(psi).*sin(psi).^3;

dy(9,:)=(-1).*h.*r.^(-1).*vu.*cos(psi)+ (-1).*dpsi.*dt.^(-1).*z.*cos(psi)+  ...
  dpsi.*dt.^(-1).*z0.*cos(psi)+ (-1).*dpsi.*dt.^(-1).*r.*sin(psi)+  ...
  dpsi.*dt.^(-1).*r0.*sin(psi)+ (-1).*dt.^(-1).*h.*r.^(-1).*z.*cos( ...
  psi).*sin(psi)+ dt.^(-1).*h.*r.^(-1).*z0.*cos(psi).*sin(psi)+ (-1).* ...
  dt.^(-1).*h.*sin(psi).^2+ dt.^(-1).*h.*r.^(-1).*r0.*sin(psi).^2;
