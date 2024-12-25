function res = twobcf(ya,yb,para,vu,ka,rb,ri)
  psi0 = 0.000;
  h0 = para;
  c0 = 0;
  alpha_i = (ka*ri*(c0^2*h0^2 - ya(2)^2))/(2*h0^2);

res = [ ya(1)-psi0;... % psi(0) = 0;
    ya(3)-ri;...   % r(0) = r0;
    ya(5);...      % area(0) = 0;
    ya(7)-alpha_i;...
    ya(8)-0;       % beta(0) = 0;
    ya(9)-0;      % vu(0) = 0;
    yb(3)-rb;...   % r(1) = rb;
    yb(4)-0;...    % z(1) = 0;
    yb(1)-0; ...      % fixed-hinge;
    yb(9)-vu;];