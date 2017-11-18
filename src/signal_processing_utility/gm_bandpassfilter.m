function [  dataFiltered, filter,nzero, amp,f  ] = ...
    gm_bandpassfilter( fq,val,dt,n, interp_method, data )
% This file is part of GeoMat.
% 
% GeoMat is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% Foobar is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
% Contact: Chen Qi Email: Geophysics.Chen.Qi@gmail.com
% =========================================================================
% Frequency domain genwav equation similar to Dr. Hilterman's
%
% Input
% =====
% fq (vector) = query frequency in Hz
% val (vector) = amplitude for the query frequency fq (val length should equal to fq length)
% dt (integer) = sampling rate in second
% n (integer) = sample number in the total filter (need to be an odd number)
% interp_method (string) = valid interpolation method allowed, refer to interp built-in function
%               = 'linear', [default]
% data (matrix)[optional] = column vectors for filtering
%
% Output
% ======
% filter (vector) = the designed filter with time zero at middle
% nzero (integer) = index in the filter for time zero
% amp (vector) = designed amplitude spectrum
% f (vector) = frequency axis for amp

df = 1/dt/n;
f = 0:df:1/2/dt;
f = f(:);
fq = fq(:);
val = val(:);

amp = interp1(fq,val,f,interp_method,0);
amp_ifft = [amp;flipud(amp(2:end-1))];

filter = [ifftshift(real(ifft(amp_ifft)));0];

nzero = ceil((n+1)/2);
dataFiltered = [];
if(exist('data','var'))
    % convolve filter with each trace in trMat
    dataFiltered = data;
    for k = 1:size(dataFiltered,2)
        dataFiltered(:,k) = gm_conv(data(:,k),filter,nzero);
    end
end

end

