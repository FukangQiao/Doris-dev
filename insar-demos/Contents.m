% InSAR matlab demos
%
% MATLAB demos for Radar Interferometry Delft University of Technology.
% 
% 
% Images:
%	master2_128_512	- raster format binary SHORT complex data.
%			  ERS2, orbit 01393, ESA SLC cut out, Italy(?), Bp=35m.
%	slave2_128_512	- raster format binary float complex data.
%			  ERS1, orbit 21066, ESA SLC cut out, resampled, Italy(?)
%	master3_128_512	- raster format binary float complex data.
%			  ERS1, orbit ?????, ESA SLC cut out, Plymouth(?), Bp=200m.
%	slave3_128_512	- raster format binary float complex data.
%			  ERS1, orbit ?????, ESA SLC cut out, resampled, Plymouth(?)
%
%
% Demos:
%       rangedemo.m	- show range spectrum master, slave, interferogram.
%	filterblock.m	- function called by filtrange for blockwize range filtering.
%	filtperiodo.m	- compare periodogram with fft over length.
%	filtrange.m	- adaptive range filter (simulation for Doris software development).
%
%
% See also the InSAR toolbox.
%

%%% EOF
%%% $Revision: 1.2 $  $Date: 2000/03/29 18:45:17 $
%%% Bert Kampes, 03-Mar-2000
