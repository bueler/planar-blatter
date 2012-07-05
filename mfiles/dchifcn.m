function dchi = dchifcn(r,xi,zeta)
% DCHIFCN  Compute two partial derivatives
%   \partial \chi_r(\xi,\zeta) / \partial \xi
%   \partial \chi_r(\xi,\zeta) / \partial \zeta
% on the reference quad [-1,1]x[-1,1], returned as a 2x1 vector.

switch r
  case 1
    dchi = [  (zeta-1) / 4;   (xi-1) / 4];
  case 2
    dchi = [- (zeta-1) / 4; - (xi+1) / 4];
  case 3
    dchi = [  (zeta+1) / 4;   (xi+1) / 4];
  case 4
    dchi = [- (zeta+1) / 4; - (xi-1) / 4];
  otherwise
    error('only r=1,2,3,4 are allowed');
end

