function dx = my_pendcart(x,m,M,L,g,d,u)

Sx = sin(x(2));
Cx = cos(x(2));
D = m*L*L*(M+m*(1-Cx^2));

dx(1,1) = x(3);
dx(2,1) = x(4);
dx(3,1) = (1/D)*(-m^2*L^2*g*Cx*Sx + m*L^2*(m*L*x(4)^2*Sx - d*x(3))) + m*L*L*(1/D)*u;
dx(4,1) = (1/D)*((m+M)*m*g*L*Sx - m*L*Cx*(m*L*x(4)^2*Sx - d*x(3))) - m*L*Cx*(1/D)*u;

end