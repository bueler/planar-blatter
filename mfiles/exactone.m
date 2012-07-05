function uexact = exactone(x,z,prm)
% EXACTONE Return exact solution described in exercise B.1 (a):
%    u(x,z) = A cosh(lambda (x-L)) sin(2 lambda z)

lambda = 3 * pi / (4 * prm.H0); % k=1 case
uexact = cosh(lambda * (x - prm.L))  / cosh(lambda * prm.L);
uexact = prm.u0 .* uexact .* sin(2 * lambda * z);

