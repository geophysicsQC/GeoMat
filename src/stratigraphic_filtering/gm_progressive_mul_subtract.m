function [ subtraction_panel ] = gm_progressive_mul_subtract( ref_coef, flood_level )
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
%
% Input
% =====
% ref_coef - 1D array represents the primary-only reflection coefficients.
% flood_level - stabalization factor that will increase the stability of
%               the frequency domain inversion. This value should be
%               between [0, 1).
%
% Output
% ======
% subtraction_panel - a matrix with dimension
%                     [length(ref_coef),length(ref_coef)]
%
% Note
% ====
% 1. ref_coef should be calculated using gm_normal_incident_rc() function.
% 2. suface reflection coefficient assumes to be 0.
%
% Reference
% =========
% 1. Chen Qi and Fred Hilterman (2017). ”Well ties for seismic with severe 
% stratigraphic filtering.” GEOPHYSICS, 82(5), IM31-IM39.
% https://doi.org/10.1190/geo2016-0695.1
%
% Dependency
% ==========
% gm_reflectivity_multiples
% gm_wavefield_dictionary

% convert ref_coef to column vector
ref_coef = ref_coef(:);

% calculate the primary plus multiples
pm = gm_reflectivity_multiples(ref_coef, length(ref_coef), 0);

% initialize output matrix
subtraction_panel = repmat(ref_coef, 1, length(ref_coef));

% calculate wavefield dictionary
wd = gm_wavefield_dictionary(ref_coef, length(ref_coef), 0);

% calculate frequency domain wavefield dictionary
nfft = 2^(nextpow2(length(ref_coef)));
WD = zeros(nfft, length(ref_coef));
for n = 1:size(WD,2)
    WD(:,n) = fft(wd(:,n), nfft);
end

% samples less than flood_amplitude will be set to it
abs_WD = abs(WD);
flood_amplitude = max(max(abs_WD)) * flood_level;
WD(abs_WD < flood_amplitude) = flood_amplitude; 

% initialize the remained trace, which will be updated and assigned to
% subtraction panel in the loop. we only need to loop 1 less time than the
% sample number since the last column of the subtraction_panel is the
% primary reflection coefficient series.
remain = pm;
for n = 1:size(subtraction_panel,2) - 1
    % subtract multiples after sample n
    remain(n:end) = ...
        remain(n:end) - ref_coef(n) * wd(1:end-n+1, n);
    % get trace to be deconvolved
    decon_trace = remain(n+1:end);
    subtraction_panel(n+1:end,n) = real(ifft(fft(decon_trace, nfft)./WD(:,n)));
end


end

