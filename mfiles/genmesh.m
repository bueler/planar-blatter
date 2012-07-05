function [xx,zz] = genmesh(I,J,x,h,b,doshow)
% GENMESH  Generate and show the mesh from input vectors x,h,b defining the
% geometry.  Vectors x,h,b must be row vectors of length I+1.  Output mesh
% arrays for plot are (J+1) x (I+1) arrays.

q = 1.5;
xx = repmat(x(:)',J+1,1);  % x(:) is a column vector
zz = zeros(size(xx));
if doshow > 0
  for i=1:I+1
    plot([x(i) x(i)],[b(i) h(i)],'k','linewidth',3.0)
  end
end
for j=0:J
  zz(j+1,:) = b + (j^q/J^q) * (h-b);
  if doshow > 0
    plot(x,zz(j+1,:),'k','linewidth',3.0)
    plot(x,zz(j+1,:),'ko','markersize',6.0,'linewidth',3.0,'markerfacecolor','k')
  end
end

