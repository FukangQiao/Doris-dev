function fractopo(D2)
%
% fractopo(D2)
%
% script to create fractal topography
% based on fractal dimension D2
% 
% example: fractopo(2.2)
% 
% Ramon Hanssen, April 2000

if nargin==0,help fractopo;break;end

% stel sampling interval = 20 m
  dx   = 0.02;
  N    = 128;
  % Fractal dimension
  %D2   =2.2;
  beta = 7-2*D2; %8/3;
  fprintf(1,'Fractal dimension = %3.1f\n',D2);
  fprintf(1,'Used 1D beta      = %3.3f\n',beta);

% simulate fractal surface
  doplot = 'n';
  [s] = fracsurf(N,beta,doplot);
  % Scale to another range
  s = 10000 *(s-0.5);  
  %s = s - 0.5;  
%    figure(1);imagesc(s);colorbar
%    xlabel(['dx = ',num2str(dx),' km']);ylabel(['dy = ',num2str(dx),' km']);

% Make a new colormap;
  j = jet(64);
  topomap = [(j(1,:));j(64:-1:15,:)];

% Show 3D simulated topography  
  % Sealevel = 1 standard deviation below mean
  fprintf(1,'Sealevel = 1 standard deviation below mean\n');
  standdev = std(s(:));
  sea = mean(s(:)) - standdev; 
  figure;colormap(topomap)
    s(find(s<sea))=sea*ones(size(find(s<sea)));
    s = s -sea;                % sea level = 0 m
    mesh(s); set(gca,'Zlim',[0 8*standdev]) 
    view(-15,60)



