% anafracdemo
% Demonstration script to play around with fractal stuff
%
% Ramon Hanssen, May 2000


fprintf(1,' First simulate a 2D fractal surface with beta = 8/3, 50 by 50 km, 200 m sampling \n');

  N    = 256; 
  dx   = 0.2;
  beta = 8/3; 

[fsurf] = fracsurf(N,beta);
 
fprintf(1,' Done.\n');
fprintf(1,' Now show and analyze this surface\n');

 fignr=1
 rotational = 'n'
 anafrac2D(fsurf,dx,fignr,rotational)

fprintf(1,' Simulate Same area with 3 regimes:\n');
  w1 = 2;  % km
  w2 = 1;  % km
  beta   = ['5/3'; '8/3';'2/3'];

fprintf(1,' a)       w > %3.1f km \t(beta = %s) \n',w1,beta(1,:));
fprintf(1,' b) %3.1f < w < %3.1f km \t(beta = %s) \n',w2,w1,beta(2,:));
fprintf(1,' c)       w < %3.1f km \t(beta = %s) \n',w2,beta(3,:));
  
  df = 1/(N*dx);
  P1 = 100 * ((2/w1) / df)/N;
  P2 = 100 * ((2/w2) / df)/N;

  betanum   = [eval(beta(1,:)),eval(beta(2,:)),eval(beta(3,:))];
  regime = [P1,P1+P2,100];

for ii = 1:3,
  PowerCoefficient = [800,1000,2000,3000];
  [fsurf1] = fracsurfatmo(N,betanum,regime,'n',PowerCoefficient(ii));
  [freq,powerv,powerh] =  anafrac2D(fsurf1,dx,3,rotational);
  lustrumboekversie = 1
  if lustrumboekversie == 1, % zwart-wit versie NL talig voor lustrumboek
    figure(5);loglog(freq,powerh,'k'); hold on
    xlabel(['Golfgetal [cycli/km]']);ylabel(['Energie [mm^2]']);
  else
    figure(5);loglog(freq,powerh,'b'); hold on
    figure(5);loglog(freq,powerv,'r');hold on
    xlabel(['Cycles per km']);ylabel(['Power [mm^2]']);
  end
end
   % fix the axis limits
     XLIM = get(gca,'xlim');
     YLIM = get(gca,'ylim');
     set(gca,'xlim',XLIM)
     set(gca,'ylim',YLIM)
   % plot power lines
     hold on
     for i = -20:25
       % turb = 10^i*XLIM.^(-5/3); % voor regime III
         turb = 10^i*XLIM.^(-8/3); % voor regime II
       h=plot(XLIM,turb);
       set(h,'linestyle',':');
     end
     hold off




fprintf(1,' \n');

%figure(1);print dummy1.eps -depsc
%figure(2);print dummy2.eps -depsc
%figure(3);print dummy3.eps -depsc
%figure(4);print dummy4.eps -depsc

% figure(3);print simu3regimes.eps -depsc
% figure(5);print simu3regimesb.eps -depsc
% !mv simu3regimes.eps ~/sar/pres/lustrumboek/.
% !mv simu3regimesb.eps ~/sar/pres/lustrumboek/.

if lustrumboekversie == 1
  figure(5);print simu3regimesb.jpg -djpeg100
  !mv simu3regimesb.jpg ~/sar/pres/lustrumboek/.
end
