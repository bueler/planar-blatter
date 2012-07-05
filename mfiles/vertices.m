function [xx, zz] = vertices(i,j,prm)
% VERTICES  Return the (x,z) coordinates of the vertices of quadrilateral with
%   lower-left indices (i,j).  Returned vectors xx, zz are each 4x1 vectors
%   giving these coordinates in lower-left-start and counter-clockwise order.

xx = zeros(4,1); zz = xx;
xx(1) = i * prm.deltax;
xx(2) = xx(1) + prm.deltax;
xx(3) = xx(2);
xx(4) = xx(1);

[h,b] = geometry(xx(1:2),prm);
H = h - b;
lower = (j / prm.J)^1.5;
upper = ((j+1) / prm.J)^1.5;
zz(1) = b(1) + lower * H(1);
zz(2) = b(2) + lower * H(2);
zz(3) = b(2) + upper * H(2);
zz(4) = b(1) + upper * H(1);

