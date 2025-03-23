% checkfrac
%
%
% script to check what the 1D and 2D power spectra of a 
% simulated fractal surface yield

clear 


% stel sampling interval = 20 m
  dx   = 0.02;
  N    = 128;
  beta = 8/3;
  fprintf(1,'Input beta = %6.3f\n',beta);

% simulate fractal surface
  doplot = 'n';
  [s] = fracsurf(N,beta,doplot);
  % Scale to another range
  %s = 100 *2 .* s;  
  %s = s - 0.5;  
    figure(1);imagesc(s);colorbar
    xlabel(['dx = ',num2str(dx),' km']);ylabel(['dy = ',num2str(dx),' km']);

%%%% 1D spectra %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1D spectra ( horizontal )
  fsh = fft(s')';
  % Periodogram spectral estimator
  powerh = (1/length(fsh)) * mean(( fsh .* conj(fsh) ));

% 1 enkele rij voor de vergelijking
  fs1 = fft(s(34,:));
  % Periodogram spectral estimator
  p1  = (1/length(fs1)) * ( fs1 .* conj(fs1) );

% 1D spectra (vertical)
  fsv = (fft(s))';
  % Periodogram spectral estimator
  powerv = (1/length(fsv)) * mean( fsv .* conj(fsv) );

% The frequency coordinate
  Nyquist = 1/(2*dx);
  df      = 1/(N*dx);
  freq  = [-Nyquist : df : ( Nyquist - df )];
  freq  = fftshift(freq);

% Use only positive part of spectrum
  powerh = powerh(1: N/2);
  powerv = powerv(1: N/2);
  power1 = p1    (1: N/2);
  freq   = freq  (1: N/2);

  figure(2); loglog(freq,powerh);
  hold on;
    loglog(freq,powerv,'r');
%    loglog(freq,power1,'g');
  hold off
  xlabel(['Cycles per km']);ylabel(['Power']);

% Calculate slopes from spectra
  [C0h,betah] = pslope(freq,powerh);
  D2h         = (7-betah)/2;
    fprintf(1,'Horizontal: C0 = %6.6f, beta = %6.3f, Fractal dim = %3.1f\n',C0h,betah,D2h);
  [C0v,betav] = pslope(freq,powerv);
  D2v         = (7-betav)/2;
    fprintf(1,'Vertical  : C0 = %6.6f, beta = %6.3f, Fractal dim = %3.1f\n',C0v,betav,D2v);
  [C01,beta1] = pslope(freq,power1);
  D21         = (7-beta1)/2;
    fprintf(1,'Sample    : C0 = %6.6f, beta = %6.3f, Fractal dim = %3.1f\n',C01,beta1,D21);


%%%% 2D spectra (rotationally averaged) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
  fprintf(1,'Calculating rotionally averaged spectrum\n');

  fs2D    = fftshift(fft2(s));
  freq    = [-Nyquist : df : ( Nyquist - df )];
  [FX,FY] = meshgrid(freq,freq);
  k       = sqrt( FX.^2 + FY.^2 );
  pow2D   = ( 1/N^2 ) .* (fs2D .* conj(fs2D) );
  pow2D1D = zeros(N/2-1,1);
  for i = 1 : N/2-1
     sel        = find( k >= i * df & k < (i+1) * df ); 
     pow2D1D(i) = mean(pow2D(sel));
  end
  freqk = [1:N/2-1] *df;
  figure(2);hold on;loglog(freqk, pow2D1D,'m');hold off

  [C0r,betar] = pslope(freqk,pow2D1D);
  D2r         = (7-betar+1)/2;
    fprintf(1,'Rotational: C0 = %6.6f, beta = %6.3f, Fractal dim = %3.1f\n',C0r,betar,D2r);
  legend(['Hor, (',num2str(betah),')'],...
         ['Ver, (',num2str(betav),')'],...
         ['Rot, (',num2str(betar),')']);

% Calculate variace fractal surface
  fprintf(1,'Variance fractal surface               : %f\n',cov(s(:)) );
  fprintf(1,'Variance cross-section fractal surface : %f\n',cov(s(34,:)) );
  %figure(3);hist(s(:),20)

% Show 3D simulated topography  
  % Sealevel = 1 standard deviation below mean
  standdev = std(s(:));
  sea = mean(s(:)) - standdev; 
  figure(3);colormap jet
    s(find(s<sea))=sea*ones(size(find(s<sea)));
    s = s -sea;                % sea level = 0 m
    mesh(s); set(gca,'Zlim',[0 8*standdev]) 
    view(-15,60)

