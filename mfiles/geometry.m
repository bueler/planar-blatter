function [h,b,dhdx] = geometry(x,prm)
% GEOMETRY Return requested information on geometry, as a function of x.  The
%   input x may be a vector, in which case each output is a vector of the same
%   size.

if (prm.testcase == 1) || (prm.testcase == 2)
  h = prm.H0 * ones(size(x));
  b = zeros(size(x));
  dhdx = zeros(size(x));
  return
end

% see formula (21) in Bueler et al (2005)
n = prm.n;  L = prm.Lmargin;
pp = 1 + 1/n;  pow = n/(2*n+2);
C = prm.H0 / (1-1/n)^pow;

ind = (x < L);
s = x(ind) / L;
chi = zeros(size(x));
chi(ind) = pp * s - 1/n + (1 - s).^pp - s.^pp;

h = C * chi.^pow;

if nargout > 1
  b = prm.bedamp * cos(2 * pi * x / prm.bedwlen);
end

if nargout > 2 
  dC = (pp / L) * pow * prm.H0 / (1-1/n)^pow;
  dpow = - (n+2) / (2*n+2);
  dchi = 1 - (1 - s).^(pp-1) - x.^(pp-1);
  dhdx = dC * chi.^dpow .* dchi;
end

