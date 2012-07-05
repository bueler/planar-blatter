function [dx,err] = testsolns(I,J,testcase,doplot);
% TESTSOLNS  Generate color plot of numerical solutions 
% to TEST1 or TEST2 on mesh with I quads in horizontal and J in vertical.
% Also generate verification data if testcase >= 1.  Runs TEST1 by default.

if nargin < 3, testcase = 1; end

prm = getparams;
prm.testcase = testcase;
prm.J = J;
prm.deltax = prm.L / I;  dx = prm.deltax;
fprintf('  delta x               = %.2f  m\n', dx )

% NUMERICAL SOLN
[Ucol, A, b] = linearfem(I,J,prm.testcase);
N = length(Ucol);
fprintf('  N, sparsity (nnz/N^2) = %d, %.3f\n', N, nnz(A)/(N*N) )
fprintf('  est. of cond(A)       = %.2f\n', condest(A) )
fprintf('  |A*U-b|_2 / |b|_2     = %.2e\n', norm(A*Ucol-b) / norm(b) )

% generically, need this stuff below
x = linspace(0,prm.L,I+1);
[h,bed] = geometry(x,prm);
[xx,zz] = genmesh(I,J,x,h,bed,0);
U = reshape(Ucol',J+1,I+1);

if doplot > 0
  figure(1), clf, subplot(1,2,2)
  surf(xx,zz,U * prm.secpera,'edgecolor','none')
  hold off, view(2), xlabel x, ylabel z, colorbar('North')
  title('numerical')
end

if testcase > 0
  % EXACT SOLN
  switch testcase
    case 1
      uu = exactone(xx,zz,prm);
    case 2
      uu = exacttwo(xx,zz,prm);
    otherwise
      error('testcase not implemented')
  end
  fprintf('  ||uexact||_infty      = %.3f  m/a\n', max(max(abs(uu))) * prm.secpera )
  err = max(max(abs(U-uu))) * prm.secpera;
  fprintf('  ||U-uexact||_infty    = %.3f  m/a\n', err )

  if doplot > 0
    subplot(1,2,1)
    surf(xx,zz,uu * prm.secpera,'edgecolor','none')
    hold off, view(2), xlabel x, ylabel z
    title('exact')
  end
else
  err = NaN;
end

