function [U,A,b] = linearfem(I,J,testcase)
% LINEARFEM  Solve the linearized Blatter equations using Q1 finite elements
%   and m=2 point Gaussian quadrature.  Mesh is an IxJ array of quadrilaterals.

global prm
prm = getparams;
if nargin < 3, testcase = 0; end
prm.testcase = testcase;
prm.J = J;
prm.deltax = prm.L / I;

N = (I+1)*(J+1);
A = sparse(N,N,9*N);                         % empty sparse format but with preallocation
b = zeros(N,1);

tic                                          % will time assembly and solve

% easy equations for Dirichlet nodes
dirscale = prm.nu0;                          % scale for better conditioning
for i=0:I                                    % bdry I
   m = nfcn(i,0,J);
   A(m,m) = dirscale;  b(m) = 0.0;
end
for j=1:J-1                                  % bdry II
   m = nfcn(0,j,J);
   [discard, zz] = vertices(0,j,prm);
   A(m,m) = dirscale;  b(m) = dirscale*alpha(zz(1)); % lower-left z
end
m = nfcn(0,J,J);
A(m,m) = dirscale;  b(m) = dirscale*alpha(zz(4)); % use upper-left of prev

% assembly of contributions from each quad_ij
for i = 0:I-1
   for j = 0:J-1
      ii = [i   i+1 i+1 i  ];                % traverse vertices of quad_ij
      jj = [j   j   j+1 j+1];                %    in counter-clockwise order
      % generate geometry and gauss point data for current quad:
      [xvert, zvert] = vertices(i,j,prm);
      [xipq, zetapq, xpq, zpq, Jpq, Bpq] = gaussquaddata(xvert,zvert,prm);
      viscpq = viscosity(xpq,zpq,prm.testcase,prm.nu0);
      for r = 1:4
         if ~isdirichlet(ii(r),jj(r))
            row = nfcn(ii(r),jj(r),J);       % row = index of equation
            b(row) = b(row) + YY(r,xipq,zetapq,xpq,Jpq);
            if ii(r)==I                      % test for Neumann bdry IV
              b(row) = b(row) + ZZ(r,zvert,zetapq,zpq);
            end
            for s = 1:4
               col = nfcn(ii(s),jj(s),J);    % col = index of unknown
               X = XX(r,s,xipq,zetapq,viscpq,Jpq,Bpq);
               if isdirichlet(ii(s),jj(s))   % symmetric Dirichlet
                  b(row) = b(row) - X * alpha(zvert(s));
               else
                  A(row,col) = A(row,col) + X;
               end 
            end % for s
         end
      end % for r
   end % for j
end % for i

fprintf('    [system assembly time = %.2f seconds]\n',toc)

% solve the linear system
tic
U = A \ b;
fprintf('    [linear solve time    = %.2f seconds]\n',toc)

end % function linearfem

% GRID HELPER FUNCTIONS

  function n = nfcn(i,j,J)
    % NFCN compute global node number from local indices
    n = i * (J+1) + j + 1;
  end

  function z = isdirichlet(i,j)
    % ISDIRICHLET  is node (i,j) on the Dirichlet part of the boundary?
    z = (i == 0) || (j == 0);
  end

% VISCOSITY and BOUNDARY CONDITION HELPER FUNCTIONS

  function nu = viscosity(xx,zz,testcase,nu0)
    switch testcase
      case 0
        error('viscosity not implemented')
      case 1
        nu = repmat(nu0,size(xx));
      case 2
        nu = repmat(nu0,size(xx));
      otherwise
        error('unknown testcase')
   end
  end

  function u = alpha(z)
    global prm
    if prm.testcase == 1
      u = exactone(0.0,z,prm);
    elseif prm.testcase == 2
      u = 0;
    else
      u = 0;  % left boundary is ice divide case
    end
  end

  function sigma = beta(z)
    global prm
    if prm.testcase == 1
      sigma = 0;
    elseif prm.testcase == 2
      lambda = 3 * pi / (4 * prm.H0); % k=1 case
      sigma = 2 * prm.sigma0 * sin(2 * lambda * z);
    else
      sigma = 0;  % left boundary is ice divide case
    end
    
  end


% MATRIX ASSEMBLY (QUAD INTEGRAL) HELPER FUNCTIONS

  function X = XX(r,s,xipq,zetapq,nupq,Jpq,Bpq)
    X = 0;
    for k=1:4
      dchis = dchifcn(s,xipq(k),zetapq(k));
      dchir = dchifcn(r,xipq(k),zetapq(k));
      % weights w_pq all = 1
      X = X + nupq(k) * (dchir' * Bpq(:,:,k) * dchis) / abs(Jpq(k));
    end
  end

  function Y = YY(r,xipq,zetapq,xpq,Jpq)
    global prm
    % get h'(x) for current Gauss points
    [discard1,discard2,dhdx] = geometry(xpq,prm);
    Y = 0.0;
    for k=1:4
      % weights w_pq all = 1
      Y = Y + dhdx(k) * chifcn(r,xipq(k),zetapq(k)) * abs(Jpq(k));
    end
    Y = - prm.rho * prm.g * Y;
  end

  function Z = ZZ(r,zvert,zetapq,zpq)
    % only evaluates correctly if along right boundary component IV
    Z = 0;
    for p=1:2
      Z = Z + beta(zpq(p+1)) * chifcn(r,1,zetapq(p+1));
    end
    Z = (zvert(3)-zvert(2)) * Z;
  end

