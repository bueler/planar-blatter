function chi = chifcn(r,xi,zeta)
% CHIFCN  Compute \chi_r(\xi,\zeta) on the reference quad [-1,1]x[-1,1].

if r==1
  chi =   (xi-1) * (zeta-1) / 4;
elseif r==2
  chi = - (xi+1) * (zeta-1) / 4;
elseif r==3
  chi =   (xi+1) * (zeta+1) / 4;
elseif r==4
  chi = - (xi-1) * (zeta+1) / 4;
else
  error('only r=1,2,3,4 are allowed');
end

