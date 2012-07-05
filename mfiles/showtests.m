function showtests(I,J)
%SHOWTESTS  Color plot of numerical solns to TEST1 and TEST2 on I x J mesh.
%  (figures "tuned" to Octave's stupidities ...)

if nargin < 1, I = 20; end
if nargin < 2, J = I; end
fprintf('showing TEST1 and TEST2 on %d x %d mesh\n',I,J)

prm = getparams;
prm.J = J;
prm.deltax = prm.L / I;  dx = prm.deltax;

figure(1),clf
for testcase = [1 2]
  fprintf('TEST%d:\n',testcase)
  prm.testcase = testcase;
  [Ucol, A, b] = linearfem(I,J,prm.testcase);
  fprintf('    done\n')

  x = linspace(0,prm.L,I+1);
  [h,bed] = geometry(x,prm);
  [xx,zz] = genmesh(I,J,x,h,bed,0);
  U = reshape(Ucol',J+1,I+1);

  subplot(2,1,testcase)
  imagesc(xx,zz,flipud(U) * prm.secpera,[-20,20])
  colorbar
  hold off, view(2), xlabel x, ylabel z
  set(gca,'yticklabel',{'1000','800','600','400','200','0'})
end

