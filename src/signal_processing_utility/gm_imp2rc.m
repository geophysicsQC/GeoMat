function [ rc ] = gm_imp2rc( imp,convention )
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
% calculate reflection coefficient from impedance
% convention: if convention is <'US'>,then the rc will be positive for harder
% rock reflection. else if convention is 'UK', then the rc will be negative
% for harder rock reflection.

if nargin < 2
    convention = 'US';
end
if strcmp(convention,'US')
    rc = (imp(2:end,:) - imp(1:end-1,:))./(imp(2:end,:) + imp(1:end-1,:));
elseif strcmp(convention,'UK')
    rc = -(imp(2:end,:) - imp(1:end-1,:))./(imp(2:end,:) + imp(1:end-1,:));
end


end

