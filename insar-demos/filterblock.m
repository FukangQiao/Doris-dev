function [masterfiltered,slavefiltered]=filterblock(master,slave,nlmean,SNRthreshold,RSR,RBW,dooversample,dohamming,debug)
% filterblock filters master,slave based on mean fft inrange (kinda periodogram). 
% [mf,sf] = filterblock(master,slave,nlmean,snrthreshold,rsr,rbw,dohamming)
% shifted spectrum computed by complex interferogram
% mean over nlmean lines, filtering based on rsr,rbw. optionally do haming filtering twice
% result is a matrix size(l-nlmean+1,:)
% intented to be used in sequential calls
% nlmean is odd!
% fftlength = power of 2
%
% first oversample data! (otherwise no diff in estimate +11 or -4 for N=15.)
% and do it with fr, not with indices.


%%% Handle input
if (nargin < 9 ) debug        = 3; end;
if (nargin < 8 ) dohamming    = 0; end;
if (nargin < 7 ) dooversample = 0; end;
%debug = 0;%			0: none; 1: little; 2: more; 3: most;


%%% Handle input (II)
if (rem(nlmean,2)==0) error('only odd nlmean');   end;
if (isreal(master))   error('only complex data'); end;
if (isreal(slave))    error('only complex data'); end;

%
[lines, fftlength] = size(master);
outputlines = lines - nlmean + 1;
startline   = ((nlmean-1)/2) + 1;%		first output line in master/slave
endline     = lines - ((nlmean-1)/2);%		last output line in master/slave
if (~ispow2(fftlength)) error('only pow2 fftlengths.'); end;

%
if (outputlines < 1)
  warning('for this block not possible, continuing');
  masterfiltered = master(startline:endline,:);
  slavefiltered  = slave (startline:endline,:);
  return;
end;

% Compute interferogram
if (dooversample==1)
  cint    = master .* conj(slave);%		complex interferogram (not ued.)
  master2 = oversample(master,1,2);%		oversample in range
  slave2  = oversample(slave,1,2);%		oversample in range
  cint2   = master2 .* conj(slave2);%		oversampled complex interferogram
  CINT2   = fft(cint2,[],2);%                   	FFT over whole range (rows)
  CINT2   = CINT2 .* conj(CINT2);%		power
else % do not oversamle
  cint   = master .* conj(slave);%		complex interferogram
  CINT   = fft(cint,[],2);%                   	FFT over whole range (rows)
  CINT2  = CINT .* conj(CINT);%			power
end


% get shift
MASTER   = fft(master,[],2);%                   	FFT over range (rows)
SLAVE    = fft(slave,[],2);%                    	FFT over range (rows)

% new way of dealing with shft etc, not as index but as frequency
% this is checked...
deltaf   = RSR/fftlength;%				interval
freqaxis = -RSR/2:deltaf:RSR/2-deltaf;%			fftlength

% compute (inverse) weighting hamming here?
if (dohamming)
  hammingwindow = myhamming(freqaxis,RBW,RSR,0.75); 
  hammingwindow(find(hammingwindow~=0))=1./hammingwindow(find(hammingwindow~=0));
end;

total  = sum(CINT2(1:nlmean,:),1);%		sum over rows
if (debug>=1)disp('matrices SHIFTS, SNRS');SHIFTS=[];SNRS=[];end;
for ii = 1:outputlines
  totalsum = sum(total);
  % get shift, snr
  [maxvalue,shift]=max(total);
  if (dooversample==1)
    SNR   = 2*fftlength * maxvalue/(totalsum-maxvalue);%
  else
    SNR   = fftlength * maxvalue/(totalsum-maxvalue);%
  end
  if (debug>=3) disp([maxvalue, (totalsum-maxvalue)/fftlength]); end;

  shift = shift - 1;%					array index starts at 1 in matlab
  if (debug>=1) SHIFTS=[SHIFTS,shift]; SNRS=[SNRS,SNR];end;
  if (debug>=2) disp([shift, SNR]); end;
% not required?.
  negshift = 0;
  % oversample
  if (dooversample==1)
    if (shift > fftlength)
      negshift = 1;%					indicate neg. shift
      shift = abs(shift-2*fftlength);
    end;
  else
    if (shift > fftlength/2) negshift = 1; shift = abs(shift-fftlength); end;
  end

  % do actual filtering
  if (SNR>SNRthreshold)
    linenumber = ii+floor(nlmean/2);%			output line in master/slave

    if (debug>=3)
      % plot all spectra, and filtered
      figure(1);
      subplot(3,2,1), plot(abs(MASTER(linenumber,:))); title('original spectrum master');
      subplot(3,2,2), plot(abs(SLAVE(linenumber,:)));  title('original spectrum slave');
      subplot(3,2,3),
	%%plot(total); title('peak estimation, red is weighted cov.');
	% compute also weigthed corr. estimation, corrected for #points.
	% this might be important for large shifts, but not normally.
        %%hold on;
	% weighting should be different, this is for sqrt of powerspectrum??
	% BK
        ll       = .5*length(total);
        %xxx      = abs(-ll/2:ll/2-1)+1;
        %xxx      = (-ll/2:ll/2-1).^2 +1;
        xxx      = (-ll:ll-1).^2 +1;
	weigcov  = total./(xxx);
	%weigcov = ll.*total./(xxx);
        plot(weigcov,'r+');
	hold off;
        [mmax,msh]=max(weigcov);
	mSNR=2*fftlength*mmax./(sum(weigcov)-mmax);
        if (dooversample==0) mSNR=mSNR./2.; end;
        msh=msh-1;
	disp([msh, mSNR]);
      %subplot(3,2,5), plot(abs(MASTER(linenumber,:)));title('filtered spectrum master');
      %subplot(3,2,6), plot(abs(SLAVE(linenumber,:))); title('filtered spectrum slave');
    end;


    if (dohamming)
    % weighting with hamming over new area
      %-shift?
      hammingwindow2 = myhamming(freqaxis-.5*shift*deltaf,RBW-abs(shift)*deltaf,RSR,0.75); 

      % make filter in 1 go
      filter = hammingwindow .* hammingwindow2;%	composite filter
    else
      % shift new window so that it is centered around new mean and widt ok.
      filter = myrect((freqaxis-.5*shift*deltaf) ./ (RBW-abs(shift)*deltaf));
    end;

    if (debug>=3)
    % plot filters
      if (dohamming)
        subplot(3,2,4), plot(fliplr(fftshift(hammingwindow)),'g'); 
	hold on;
        subplot(3,2,4), plot(fliplr(fftshift(hammingwindow2)),'b'); 
      end;
      subplot(3,2,4), plot(fliplr(fftshift(filter)),'r--');%
      hold off;
      title('filters: hamming or rect, striped is composite.');
      %disp('press key'); pause;
    end;

    filter = ifftshift(filter);%

    % which side to filter?
    % this seems to be ok from exp. with data
    %if (negshift==1)
    if (negshift==0)
      MASTER(linenumber,:) = MASTER(linenumber,:) .* filter;
      SLAVE(linenumber, :) = SLAVE(linenumber, :) .* fliplr(filter);
    else
      MASTER(linenumber,:) = MASTER(linenumber,:) .* fliplr(filter);
      SLAVE(linenumber, :) = SLAVE(linenumber, :) .* filter;
    end;


    if (debug>=3)
      % plot all spectra, and filtered
      figure(1);
      %subplot(3,2,1), plot(abs(MASTER(linenumber,:)));title('original spectrum master');
      %subplot(3,2,2), plot(abs(SLAVE(linenumber,:))); title('original spectrum slave');
      %subplot(3,2,3), plot(total); title('peak estimation');
      subplot(3,2,5), plot(abs(MASTER(linenumber,:))); title('filtered spectrum master');
      subplot(3,2,6), plot(abs(SLAVE(linenumber,:)));  title('filtered spectrum slave');
      disp('press key'); pause;
    end;

  else% threshold
    disp('SNR < threshold; no filtering');
    disp([SNR SNRthreshold]);
  end;

  % update sums over blocks, do this with swap function...(?)
  % update total
  if (ii~=outputlines)%					if this then break.
    line1 = CINT2(ii,:);
    lineN = CINT2(nlmean+ii,:);%			next one
    %diff = lineN-line1;
    %total = total + diff;
    total = total - line1 + lineN;%			why divide...
  end;
end;

% Inverse FFT over filtered part
masterfiltered = ifft(MASTER(startline:endline,:),[],2);
slavefiltered  = ifft(SLAVE (startline:endline,:),[],2);



if (debug >= 1)
  figure(2);
    cint2 = ifft(MASTER,[],2).*conj(ifft(SLAVE,[],2));
    subplot(2,2,1), imagesc(angle(cint)); colorbar;
    title('phase map of complex interferogram.');
    subplot(2,2,2), imagesc(angle(cint2)); colorbar;
    title('phase map of filtered complex interferogram.');
    subplot(2,2,3), imagesc(abs(fft(master,[],2))); colorbar;
    title('range spectra of master.');
    subplot(2,2,4), imagesc(abs(MASTER)); colorbar;
    title('range spectra of filtered master.');
%  figure(3);
%    imagesc(CINT2); colorbar;
%    title('powerspectrum of complex interferogram.');
  disp('Press key to continue');
  pause;
end;


%%% EOF
