function prm = getparams;
% GETPARAMS Puts the major parameters for Blatter-related codes in one place.

spera = 31556926.0;
n = 3.0;
Dtyp = 1e-3/spera;
B = 1.9e7;
nu0 = (B / 2) * Dtyp^((1/n)-1);
Lmargin = 11000.0;
L = 10000.0;
u0 = 20.0 / spera;

prm = struct('secpera',spera,...  % seconds in one year
             'n',n,...            % Glen flow law exponent
             'g',9.81,...         % m s-2
             'rho',910.0,...      % kg m-3
             'H0',1000.0,...      % m
             'Lmargin',Lmargin,...% m; location of margin in SIA exact soln used to generate geometry
             'L',L,...            % m; horizontal length of computational domain
             'B',B,...            % Pa s^{1/3}; 1/10th of EISMINT-Ross value
             'Dtyp',Dtyp,...      % s-1; typical strain-rate of 1 m a-1 in 1 km
             'nu0',nu0,...        % Pa s; viscosity corresponding to B,Dtyp,n; in TEST1 and TEST2
             'u0',u0,...          % s-1; 20 m/a; velocity used in TEST1 b.c.
             'sigma0',1.0e5,...   % Pa; = 1 bar; deviatoric stress for TEST2 b.c.
             'bedamp',100.0,...   % m; amplitude of bed bumpiness when used to generate geometry
             'bedwlen',4400.0,... % m; wavelength of bed bumpiness when ...
             'testcase',0,...     % default only
             'deltax',-1,...      % default only
             'J',-1);             % default only

