function figgen(mlist,makepdfs);
% FIGGEN  Plot the profile of an ice sheet in several ways, to help illustrate
% the notes about the flow line Blatter implementation.  In the case of a flat
% bed the profile would be an exact solution to the steady SIA, but we are
% also using that profile with a nonflat bed.

if nargin < 1, mlist = 1:4; end
if nargin < 2, makepdfs = 0; end

prm = getparams;
fignames = {'profile.pdf','bdryblatter.pdf','conform.pdf','q1mesh.pdf'};

xmin = 0.0;  xmax = prm.L;
xmid = (xmax + xmin)/2;
pts = 500;
prm.deltax = prm.L / (pts-1);
prm.J = 10; % dummy only
x = linspace(xmin,xmax,pts);
[h,b] = geometry(x,prm);

for m = mlist
  figure(m), clf

  plot(x,h,'k','linewidth',5.0), hold on
  plot(x,b,'k','linewidth',5.0)
  plot([xmin xmin],[b(1) h(1)],'k','linewidth',5.0)
  plot([xmax xmax],[b(end) h(end)],'k','linewidth',5.0)

  if m==1
    % profile overview with velocity arrows
    text(xmid,1.1*h(floor(pts/2)),'z = h(x)','fontsize',18.0)
    text(0.65*xmid,1.1*min(b),'z = b(x)','fontsize',18.0)
    plot_arrow(0.2*xmax,0.5*(b(1)+h(1)), 0.35*xmax,0.48*(b(1)+h(1)),...
               'linewidth', 3, 'headwidth', 0.02 );
    plot_arrow(0.2*xmax,0.2*(b(1)+h(1)), 0.31*xmax,0.19*(b(1)+h(1)),...
              'linewidth', 3, 'headwidth', 0.02 );
    plot_arrow(0.4*xmax,0.5*(b(1)+h(1)), 0.56*xmax,0.51*(b(1)+h(1)),...
               'linewidth', 3, 'headwidth', 0.02 );
    plot_arrow(0.4*xmax,0.2*(b(1)+h(1)), 0.49*xmax,0.19*(b(1)+h(1)),...
               'linewidth', 3, 'headwidth', 0.02 );
    plot_arrow(0.7*xmax,0.5*(b(1)+h(1)), 0.85*xmax,0.45*(b(1)+h(1)),...
               'linewidth', 3, 'headwidth', 0.02 );
    plot_arrow(0.7*xmax,0.2*(b(1)+h(1)), 0.83*xmax,0.19*(b(1)+h(1)),...
               'linewidth', 3, 'headwidth', 0.02 );
    text(0.5*xmax,0.35*(b(1)+h(1)),'velocity  (u,w)','fontsize',18.0)

  elseif m==2
    % boundary condition figure
    text(0.6*xmid,1.3*min(b),'I:    u = 0','fontsize',18.0)
    text(-0.1*xmid,0.3*(b(1)+h(1)),'II:    u = \alpha(z)','rotation',90.0,'fontsize',18.0)
    text(0.15*xmid,1.35*h(floor(pts/2)),'III:   \sigma_{ij} n_j = 0  \Leftrightarrow  (4 \nu u_x, \nu u_z) \cdot n = 0','fontsize',18.0,'rotation',-20.0)
    text(1.05*xmax,1.25*h(end),'IV:  2 \nu u_x = \beta(z)','rotation',-90.0,'fontsize',18.0)

  elseif m==3
    % conforming curves
    q = 1.5;
    J = 5;
    for j=1:J-1
      z = b + (j^q/J^q) * (h-b);
      plot(x,z,'k','linewidth',4.0)
    end

  elseif m==4
    % Q1 grid
    clf
    I = 9;
    J = 5;
    xc = linspace(xmin,xmax,I+1);
    [hc,bc] = geometry(xc,prm);
    plot(xc,bc), hold on  % get plot started
    genmesh(I,J,xc,hc,bc,1);

  else,  error('how did I get here?'),  end

  hold off
  if m<3, axis tight, end
  axis off
  if makepdfs>0, print(fignames{m},'-dpdf'), end
end

