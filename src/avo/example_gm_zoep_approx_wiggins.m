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
% along with GeoMat If not, see <http://www.gnu.org/licenses/>.
% Contact: Chen Qi Email: Geophysics.Chen.Qi@gmail.com
% =========================================================================
% Show the example of avo plotting and visually correlate to the reference
% book's plotting (Hilterman, 2001, page 3-31, Figure 3.C.14)

% The result shows different than the plotting of the book because they are
% using exact zoeppriz equation but we are using approximation.
incident_angle = 0:40;
rc = gm_zoep_approx_wiggins([10000,13000], [2.4,2.64], [5800,7500], deg2rad(incident_angle));
plot(incident_angle, rc);
xlabel('Incident Angle Degree');
ylabel('Amplitude');
ylim([-0.2,0.3])
xlim([0,40])