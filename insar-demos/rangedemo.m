% RANGEDEMO demonstrates adaptive filter in range.
%   show range spectrum of 2 images
%   and of interferogram;
%   test peak estimation of FFT of (non oversampled!) interferogram.
%   test for range filtering in Doris software.
%   Change parameters in this script for more control.
%
%   See also FILTRANGE, FILTERBLOCK, INSAR toolbox.
%

% $Revision: 1.2 $  $Date: 2000/03/29 18:45:43 $
% Bert Kampes, 03-Mar-2000

%%% Initialize
if (~exist('m','var') & ~exist('s','var')) clear; end; close all force; more off;
showspectrummaster=1;
showspectrumslave=1;
dohamming=1;%			apply hamming 2 times


%%% Initial main parameters of images.
%imagedir ='/home/fmr/d1/scripts/matlab/insar/';
imagedir = [];
master   = 'master2_128_512';%		SHORT
slave    = 'slave2_128_512';%		FLOAT
L        = 128;%			number of lines
P        = 512;

%%% Read images from disk if not present.
if (~exist('m','var')) m=freadbk([imagedir, master],L,'cmplxint16'); end
if (~exist('s','var')) s=freadbk([imagedir, slave], L,'cmplxfloat32'); end


%%% Plot variables.
sub = 10;%			do a plot per sub
subplots = floor(P/sub);%	number of subplots
xaxis = 0:P-1;
xaxis = xaxis - P/2;


%%% Show spectrum in range for master if requested.
if (showspectrummaster==1)
  total = zeros(1,P);
  plotnumber = 0;
  for i=1:L
    Mrange = fft(m(i,:));
    Mrange = abs(fftshift(Mrange));
    if (rem(i,sub)==0)
      plotnumber = plotnumber + 1;
      figure(plotnumber);
      %subplot(subplots,1,plotnumber),
      plot(xaxis,Mrange);
      eval(['title (''range spectrum for master line ',num2str(i),''')']);
    end
    total = total + Mrange;
  end
  %%% Averaged spectrum
  figure(plotnumber+1);
  total = total / L;
  plot(xaxis,total);
  eval(['title (''average range spectrum for all ',num2str(L),' lines (master)'')']);

  %%% Plot for all 
  figure(plotnumber+2);
  imagesc(abs(fftshift(fft(m,[],2)))); colorbar
  eval(['title (''range spectrum for all ',num2str(L),' lines (master)'')']);
end


%%% Show spectrum in range for slave if requested.
if (showspectrumslave==1)
  disp ('press key to continue with slave');
  pause;
  close all force;
  total = zeros(1,P);
  plotnumber = 0;
  for i=1:L
    Srange = fft(s(i,:));
    Srange = abs(fftshift(Srange));
    if (rem(i,sub)==0)
      plotnumber = plotnumber + 1;
      figure(plotnumber);
      plot(xaxis,Srange);
      eval(['title (''range spectrum for slave line ',num2str(i),''')']);
    end
    total = total + Srange;
  end
  %%% Averaged spectrum
  figure(plotnumber+1);
  total = total / L;
  plot(xaxis,total);
  eval(['title (''average range spectrum for all ',num2str(L),' lines (slave)'')']);
  
  %%% Plot for all 
  figure(plotnumber+2);
  imagesc(abs(fftshift(fft(m,[],2)))); colorbar
  eval(['title (''range spectrum for all ',num2str(L),' lines (slave)'')']);
end 


%%% Test interferogram gen. by pointwize mult.
disp ('press key to continue with complex interferogram');
pause;
close all force;
cint = m .* conj(s);%		complex interferogram

% Show phase map
figure (1);
imagesc(angle(cint)); colorbar
title ('phase of complex interferogram.');

% Show range spectrum for some (total) lines
% use mean over small block in azimuth (8 lines?)
plotnumber = 1;
fftlength    = P;
meanoverrows =   8;
for i=1:L
  if (rem(i,sub)==0)
    Irange = fft(cint(i-7:i,:),[],2);%			fft over rows
    Irange = mean(Irange,1);%				mean over rows
    Irange = Irange .* conj(Irange);
    [maxvalue,ii]=max(Irange);
    SNR = maxvalue ./ ((sum(Irange)-maxvalue)/P);

    shifts(plotnumber) = ii;%				check this by swapping m/s
    %shifts(plotnumber) = fftlength-ii%			check this by swapping m/s
    snrs(plotnumber) = SNR;
    plotnumber = plotnumber + 1;
    figure(plotnumber);
    subplot(2,1,1), plot(xaxis,Irange,'r');
    eval(['title (''range spectrum for interferogram line ',num2str(i),''')']);
    Irange = fft(cint(i,:));
    Irange = abs(Irange);
    hold on
    subplot(2,1,2), plot(xaxis,Irange,'b');
    title('fft only over 1 line');
    hold off
  end
end

disp('maxima complex interferogram spectrum');
shifts=shifts
snrs=snrs



%%% EOF
more on
