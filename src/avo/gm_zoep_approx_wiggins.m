function [ rc, A, B, C ] = gm_zoep_approx_wiggins( ...
    alpha, rho, beta, theta)
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
% calculate reflectivities based on the linear approximation of zoeppritz
% equation. This function is based on Wiggins et al.,1983
% 
% Input
% =====
% alpha - p-wave velocity vector
% rho - density
% beta - s-wave velocity vector
% theta - incident angles
%
% Output
% ======
% rc - matrix of reflectivities, the column number equals to length(theta)
% A, B, C - three terms in the AVO equation, they are defined as:
%   A = (E + F)/2
%   B = E/2 - 2*H^2*D
%   C = E/2
%   where
%   E = delta_alpha/mean_alpha
%   F = delta_rho/mean_rho
%   G = delta_beta/mean_beta
%   H = mean_beta/mean_alpha
%   D = 2*G + F
%   delta_alpha = transmitted alpha - incident alpha
%   mean_alpha = (transmitted alpha + incident alpha)/2
%   (same logic applis to rho and beta)
% 
%   The output reflectivity can be expressed by A, B, C and theta
%   rc = A + B*sin(theta)^2 + C*sin(theta)^2*tan(theta)^2
% Reference
% =========
% Fred Hilterman, Seismic Amplitude Interpretation, SEG, 2001, 6-19 and
% 3-15

% convert all input to column vectors
alpha = alpha(:);
rho = rho(:);
beta = beta(:);
theta = theta(:);

% initialize output
rc = zeros(length(alpha)-1, length(theta));

% temporary parameters
delta_alpha = alpha(2:end) - alpha(1:end-1);
delta_rho = rho(2:end) - rho(1:end-1);
delta_beta = beta(2:end) - beta(1:end-1);
mean_alpha = (alpha(2:end) + alpha(1:end-1))/2;
mean_rho = (rho(2:end) + rho(1:end-1))/2;
mean_beta = (beta(2:end) + beta(1:end-1))/2;

E = delta_alpha./mean_alpha;
F = delta_rho./mean_rho;
G = delta_beta./mean_beta;
H = mean_beta./mean_alpha;
D = 2*G + F;

% Generate output
A = (E + F)/2;
B = E/2 - 2*H.^2.*D;
C = E/2;
for n = 1:length(theta)
    rc(:,n) = A + B.*sin(theta).^2 + C.*sin(theta).^2.*tan(theta).^2;
end

end

