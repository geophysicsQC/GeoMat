function [ timecur, twt ] = gm_depth2time( ...
    depcur, velocity, depstep, dt)
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
% along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
% Contact: Chen Qi Email: Geophysics.Chen.Qi@gmail.com
% =========================================================================
% 
% Input
% =====
% depcur - matrix of depth domain each column represents one log curve
% velocity - vector of velocity used for converting
% depstep - depth sampling interval (1-way)
% dt - two-way sampling interval for the time domain curve
%
% Output
% ======
% timecur - matrix of time domain log curves correponding to depcur matrix
% twt - vector of two-way sampling time (convenient for plotting)

velocity = velocity(:);
time_consumed = cumsum(2*depstep./velocity);
time_consumed = [0;time_consumed(1:end-1)];
twt = 0:dt:max(time_consumed);
timecur = zeros(length(twt), size(depcur,2));
for n = 1:size(depcur,2)
    timecur(:,n) = interp1(time_consumed, depcur(:,n), twt);
end

end

