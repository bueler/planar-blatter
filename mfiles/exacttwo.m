function uexact = exacttwo(x,z,prm)
% EXACTTWO Return exact solution described in exercise B.1 (b):
%    u(x,z) = B sinh(lambda x) sin(2 lambda z)

lambda = 3 * pi / (4 * prm.H0); % k=1 case
uexact = sinh(lambda * x)  / (prm.nu0 * lambda * cosh(lambda * prm.L));
uexact = prm.sigma0 .* uexact .* sin(2 * lambda * z);

