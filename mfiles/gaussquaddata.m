function [xipq, zetapq, xpq, zpq, Jpq, Bpq] = gaussquaddata(xvert,zvert,prm)
% GAUSSQUADDATA Provide Gaussian quadrature-related numbers for quadrilateral
%   with vertices xvert,zvert.  Inputs xvert, zvert and outputs xipq, zetapq,
%   xpq, zpq, Jpq are all 4x1 vectors.  Output Bpq is a list of four 2x2
%   matrices.  Specifically, Bpq(:,:,p) is 2x2 for each p=1,2,3,4.

om = 1/sqrt(3);
xipq   = [-om  om  om -om]';
zetapq = [-om -om  om  om]';

dz12 = zvert(1)-zvert(2);
dz14 = zvert(1)-zvert(4);
dz32 = zvert(3)-zvert(2);
dz34 = zvert(3)-zvert(4);

xpq = zeros(4,1);  zpq = xpq;  Jpq = xpq;
Bpq = zeros(2,2,4);

for k=1:4               % index of Gauss point at which we want outputs
  xi   = xipq(k);
  zeta = zetapq(k);

  % get all \chi_r involved in computing x,z for s Gauss point
  chi = zeros(4,1);
  for r=1:4
    chi(r) = chifcn(r,xi,zeta);
  end

  % x,z coordinates of Gauss points
  xpq(k) = xvert' * chi;
  zpq(k) = zvert' * chi;

  % Jacobian determinant of reference map at Gauss points
  tmp = dz14 * (xi-1) + dz32 * (xi+1);
  Jpq(k) = (prm.deltax / 8) * tmp;
  
  % build quadratic form matrix B at Gauss points
  E = zeros(2,2);
  E(1,1) = tmp;
  E(1,2) = - dz12 * (zeta-1) + dz34 * (zeta+1);
  E(2,2) = 2 * prm.deltax;
  E = E / 4;
  Bpq(:,:,k) = E' * diag([4 1]) * E;
end

