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
% This is the script file for testing the progressive polynomial division
%% Create the cyclically stratified sequence
N = 15; % number of cycles
a = 12000 * 2.3; % impedance for shale
b = 7700  * 1.3; % impedance for coal
coal_thick = [a;b;b;a]; % thick coal is with two unit thickness
imp = ... % impedance log is not connected
    [a * ones(30,1); repmat(coal_thick,N,1);...
    a * ones(200,1)];
dt = 0.001; % ms
% calculate reflection coefficients out of impedance
c = gm_imp2rc(imp);

% apply progressive deconvolution to generate multiple subtraction panel
mul_subtraction_panel = gm_progressive_mul_subtract(c, 0.001);

% filter the result to the seismic bandwidth
freq_corner = [1,10,100,120]; % Hz
mul_subtraction_panel = gm_bandpassfilter(freq_corner, ...
    [0,1,1,0], dt, 501, 'linear', mul_subtraction_panel);

%% generate report figure
figure;

subplot(4, 4, 1:3);
time_axis = 0:dt:(length(c)-1)*dt;
plot(time_axis, imp(2:end), 'linewidth', 2); 
xlabel('Time (s)'); ylabel('Impedance'); grid on;
xlim([0, 0.15]);
title('Progressive Multiple Subtraction Demo')

subplot(4, 4, [5:7, 9:11, 13:15]);
imagesc(time_axis, time_axis, mul_subtraction_panel); colormap jet;
xlabel('Time (s)'); ylabel('Time (s)'); axis tight; grid on;
gm_wva(mul_subtraction_panel, time_axis, time_axis, 10, [], [], [], [], 'hold');
xlim([0, 0.15]); ylim([0, 0.15]);

subplot(4, 4, [8, 12, 16])
plot(imp(2:end), time_axis, 'linewidth', 2);
xlabel('Impedance'); ylabel('Time (s)'); grid on; axis ij;
ylim([0, 0.15]);
set(gca, 'yaxislocation', 'right', 'xaxislocation', 'top');
set(gcf, 'position', [680   309   612   669]);
