function [rttime] = Pix2Time(rtpix,Mod,Freq,isot,varargin)
% From rtpix a retention vector (in pixels), a point per line 
% (output in (s;min), input in (pixel;pixel) (with first value at (1;1)
% BEWARE, this is not the same as GC Image-way I believe.)
% Mod is modulation rate, in [s], Freq is sampling frequency, in [Hz], isot
% = if you want to suppress some time at the beginning of the chromatogram,
% in [min].
% (by default, isot = 0 if not given)
% One option:
% -dim: if set to 1 or 2, the unit conversion will be done assuming the
% data is only first (1) or second (2) dimension. Can also be a vector
% (e.g. [2,1]), in case the data order is not the one expected by default.
% 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% © All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
% Switzerland, Environmental Chemistry Modeling Laboratory, 2014. 
% Authors : Jonas Gros and J. Samuel Arey.
% See academic license terms stated in file:
% LICENSE.txt
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% First, if ever the user used some options, give their values to the
% corresponding variables:
for i=1:2:length(varargin)
    if strcmpi(varargin{i},'dim')
        dim = varargin{i+1};
    else
        error(['The option ''' varargin{i} ''' is unknown'])
    end
end

if ~exist('isot','var')
    isot = 0;
end
if ~exist('dim','var')
    rttime(:,2)=rtpix(:,2)/Freq;
    rttime(:,1)=(rtpix(:,1))*Mod/60+isot;
elseif length(dim)>1
    rttime = zeros(size(rtpix));
    for k = 1:length(rtpix(1,:))
        if dim(k)==2
            rttime(:,k) = rtpix(:,k)/Freq;
        else
            rttime(:,k) = (rtpix(:,k))*Mod/60+isot;
        end
    end
elseif length(dim)==1
    if dim==2
            rttime = rtpix/Freq;
    else
            rttime = (rtpix)*Mod/60+isot;
    end
end

