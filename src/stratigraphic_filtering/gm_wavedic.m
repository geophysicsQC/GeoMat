function [ wavedic] = gm_wavedic(rc,N,rc0 )
% This file is part of GeoMat.
% 
% GeoMat is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% GeoMat is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with GeoMat.  If not, see <http://www.gnu.org/licenses/>.
% Contact: Chen Qi Email: Geophysics.Chen.Qi@gmail.com
% =========================================================================
% the function calculate the wavefield dictionary
%
% Input
% =====
% rc - reflection coefficients of the earth
% N - length of sample points of each wavelets
% rc0 - reflection coefficient of the surface reflector
%
% Output
% ======
% wavedic - matrix of wavefield dictionary
%
% Dependencies
% =============
% gm_reflectivity_multiples

wavedic = zeros(N,length(rc));
rc = reshape(rc,length(rc),1);
inverse_rc = -flipud(rc);
h = waitbar(0,'Initializing waitbar...');
for n = 1:length(rc)
    % Calculate the downward transmitted wave
    if n~=1
        down_transmitted = gm_reflectivity_multiples(rc(1:n-1),N,rc0,'trans');
    else
        down_transmitted = zeros(N,1);
        down_transmitted(1) = 1;
    end
    % Calculate the upward transmitted wave
    if n~=1
        up_transmitted = gm_reflectivity_multiples([inverse_rc(end-n+2:end);-rc0],N,inverse_rc(end-n+1),'trans');
    else
        up_transmitted = gm_reflectivity_multiples(-rc0,N,inverse_rc(end),'trans');
    end
        % Calculate the wave at the layer uppon its transmission
        up_transmitted = up_transmitted/(1+rc0);
    % Convolve the downgoing and upgoing
    tempt = conv(up_transmitted,down_transmitted);
    wavedic(:,n) = tempt(1:N);
    waitbar(n/length(rc),h,['Wavefield dictionary ',num2str(floor(n/length(rc)*100)),'% finished']);
end
close(h);

end

