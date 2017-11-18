function [ data_output ] = gm_reflectivity_multiples( rc,N,r0,wavetype )
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
%==========================================================================
% This is a basic function calculating 1D synthetic seismogram. According to Rosa's book, JUL Verion2
% The program will calculate upgoing wavefield at the first layer as Si (with transmission loss, internal
% multiples, primary reflections) and the measured seismogram as S (with transmission loss, internal
% multiples, primary reflections and free surface multiples). The Si and S have a simple relationship in
% Z domain: S = Si/(1+Si).S and Si formula:
% S = Z^k G*/(H+r0 Z^k G*); Si = Z^k G*/H. Where G* = G(Z^(-1)) and H* = H(Z^(-1));
% G and H could be derived by this recursive formula system:
% [H;G]_k = [1,r_k Z;r_k,Z][H;G]_(k-1) and H_1 = 1, G_1 = r_1
%
% Input
% =====
% rc - reflection coefficients from interface 1 to interface n
% N - the number of samples in the synthetic seismogram. By default, the length of the seismogram is
% length(rc)+1 where the extra 1 is the very first output sample which equals to r0.
% r0 - reflection coefficient for the free surface. By default, r0 = 0
%
% Output
% ======
% Depend on wavetype
% wavetype == 'ref' - Measured wavefield in the first layer (includes free surface multiples. If r0 = 0 then S = Si
% which means no free surface multiples).
%
% wavetype == 'trans' - The escaping wave at the bottom of all the layers (downgoing wavefield at bottom)
% The program will evaluate number of the output variables to decide whether or not to proceed to
% calculate S and E. This is useful for frequent calling of this function to get Si only.
%
% wavetype == 'mulfil' - multiple filter without free surface multiples.
%
% Dependencies
% ============
% gm_invzdeconv

if nargin < 2
    N = length(rc);
end
if nargin < 3
    r0 = 0;
end
if nargin < 4
    wavetype = 'ref';
end
H = zeros(length(rc),1);
G = H;
H(1) = 1;
G(1) = rc(1);
tc = 1 - rc; % Transmission coefficients
TC = tc(1); % Result of continues multiplication
% Use recursive formula to calculate H and G
for n = 2:length(H)
    H_temp = H(1:n);
    H(1:n) = H(1:n) + [0;rc(n)*G(1:n-1)];
    G(1:n) = rc(n)*H_temp + [0;G(1:n-1)];
end

HZG = [0;r0*flipud(G)]+[H;0]; % Save time

if strcmp(wavetype,'ref')
    S = gm_invzdeconv(flipud(G),HZG,N);
    S = reshape(S,length(S),1);
    data_output = S;
end
if strcmp(wavetype,'trans')
    for n = 2:length(rc)
        TC = TC*tc(n);
    end
    E = TC*gm_invzdeconv(1,HZG,N);
    E = reshape(E,length(E),1);
    data_output = E;
end
if strcmp(wavetype,'mulfil')
    TC2 = 1-rc(1)^2;
    for n = 2:length(rc)
        TC2 = TC2*(1-rc(n)^2);
    end
    HH = conv(H,flipud(H));
    HH = HH(1:length(H));
    E = TC2*gm_invzdeconv(1,HH,N);
    E = reshape(E,length(E),1);
    data_output = E;
end
end

