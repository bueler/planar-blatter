function [A,b] = sparsepattern(I,J)
% SPARSEPATTERN  Generate the sparsity pattern, though not the correct entries,
% for the Q_1 finite element matrix for the Blatter equations.

N = (I+1)*(J+1);
A = sparse(N,N);                             % empty sparse format
b = zeros(N,1);

% easy equations for Dirichlet nodes
for i=0:I                                    % bdry I
   m = nfcn(i,0,J);  A(m,m) = 1.0;  b(m) = 0.0;
end
for j=1:J                                    % bdry II
   m = nfcn(0,j,J);  A(m,m) = 1.0;  b(m) = alpha(j);
end

% assembly of contributions from each quad_ij
for i = 0:I-1
   for j = 0:J-1
      ii = [i   i+1 i+1 i  ];                % traverse vertices of quad_ij
      jj = [j   j   j+1 j+1];                %    in counter-clockwise order
      for r = 1:4
         if ~isdirichlet(ii(r),jj(r))
            row = nfcn(ii(r),jj(r),J);            % row = index of equation
            b(row) = b(row) + YY(i,j,r);
            if isneumannIV(ii(r),jj(r),I)
              b(row) = b(row) + ZZ(i,j,r);
            end
            for s = 1:4
               col = nfcn(ii(s),jj(s),J);         % col = index of unknown
               X = XX(i,j,r,s);
               if isdirichlet(ii(s),jj(s))   % symmetric Dirichlet
                  b(row) = b(row) - X * alpha(jj(s));
               else
                  A(row,col) += X;
               end 
            end % for s
         end
      end % for r
   end % for j
end % for i

  function n = nfcn(i,j,JJ)
    % NN compute global node number from local indices
    n = i * (JJ+1) + j + 1;
  end

  function z = isdirichlet(i,j)
    z = (i == 0) || (j == 0);
  end

  function z = isneumannIV(i,j,II)
    z = (i == II);
  end

  function u = alpha(j)
    u = 0;  % left boundary is ice divide case
  end

  function x = XX(i,j,r,s)
    % FAKE
    if (r>=s)
      x = i+j + r-s;
    else
      x = XX(i,j,s,r);
    end
  end

  function y = YY(i,j,r)
    % FAKE
    y = i+j+r;
  end

  function z = ZZ(i,j,r)
    % FAKE
    z = -i-j-r;
  end

