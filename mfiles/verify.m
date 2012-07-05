function verify(N)
% VERIFY Trivial script to use TESTSOLNS to verify TEST1 and TEST2 and generate
% figure for the paper.  Note 'verify(N)' does N levels of refinement.

if nargin < 1, N = 3; end

JJ = [10 20 40 80 160];
JJ = JJ(1:N);

dx = zeros(size(JJ));
err = zeros(2,length(JJ));

for k = 1:length(JJ)
  j = JJ(k);
  for testcase = [1 2]
    fprintf('TEST%d, case I = J = %d:\n',testcase,j)
    [dx(k),err(testcase,k)] = testsolns(j,j,testcase,0);
  end
end

for testcase = [1 2]
  p = polyfit(log(dx),log(err(testcase,:)),1);
  fprintf('\nTEST%d VERIFY conclusion:  convergence rate O(dx^%.3f)\n',...
          testcase,p(1));
end

figure
loglog(dx,err','o','markersize',12,'linewidth',3.0)
axis([30 2000 0.005 5])
legend('TEST1','TEST2','location','southeast')

set(gca,'xtick',fliplr(dx))
xl = {'62.5','125','250','500','1000'};
set(gca,'xticklabel',xl(end-length(dx)+1:end))

set(gca,'ytick',[0.01 0.03 0.1 0.3 1 3])
set(gca,'yticklabel',{'0.01', '0.03', '0.1', '0.3', '1', '3'})

set(gca,'fontsize',14.0)
xlabel('\Delta x   (m)','fontsize',12.0)
ylabel('numerical error   (m/a)','fontsize',12.0)
grid on

print('linearverify.pdf','-dpdf')

