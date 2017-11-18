function [ c,r ] = gm_invzdeconv( a,b,n )
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
% inverse z transfrom based on matlab deconv function
%
% Input
% =====
% a - numerator (will be converted to column)
% b - denomitor (will be converted to column)
% n - terms of output (default = length(b))
%
% (NOTE: the order of a and b for the z-transform depends on the definition of
% the z-transfrom. In geophysical aspects we uses z^k and in other cases they
% use z^-k so in geophysics a and b should be in the ascending order of z and
% in other cases should be in descending order)
%
% Output
% ======
% c - output resulted quotient (in column)
% r - remaining (in column)
% 
% by Chen Qi
% Univ of Houston

if nargin < 3
     n = length(b);
end

a = reshape(a,length(a),1);
b = reshape(b,length(b),1);

zero_needed = length(b) + n - 1 - length(a);
a = [a;zeros(zero_needed,1)];
[c,r] = deconv(a,b);

end

