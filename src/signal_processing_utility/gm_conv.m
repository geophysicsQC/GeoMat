function [ b ] = gm_conv( a,sig,nzero )
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
% The function's premitive is the matlab built-in function conv. The difference
% is that this function only returns the b with size(b) = size(a) and a is
% filtered by sig whose time zero is at nzero sample point. 
% By default, nzero is at mid of the sig, let nzero = 1 for a causal signal.

if nargin < 3
    nzero = ceil((length(a) + 1)/2);
end
if length(a)==1
    a = a(:);
end
sig = sig(:);
b = zeros(size(a));
for k = 1:size(a,2)
    c = conv(a(:,k),sig);
    b(:,k) = c(nzero:nzero + length(a) - 1);
end

end

