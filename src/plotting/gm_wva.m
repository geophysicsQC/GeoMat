function gm_wva( ...
    d,y,x,index_incre,lvl,clip,line_color,face_color,mode , trace_balance )
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
% wiggle variable area plot of seismic data
%
% Input
% =====
% d = matrix of data
% y = [] = vector y means axis tick in y, scalar y means sampling rate in y
% x = [] = same for y but for horizontal direction (x direction)
% index_incre = 1 = index increment. 1 stands for no
%                  skipping. This is useful when there are too many 
%                  traces on screen.
% lvl = 1 = positive = level of gain & fill peaks
%         = negative = level of gain & fill troughs
% clip = 3 = clip level
% line_color = 'k' = line color
% face_color = 'k' = face (patch) color
% mode = 'new'  = open a new figure and new axis
%        'hold' = plot on the current axis (default)
%        'wipe' = clear current axis and draw a new one
% trace_balance = [] = don't balance trace
%               = 'max' = balance according to maximum values
%               = 'rms' = using rms to balance trace

% avoid nan values
d(isnan(d)) = 0;

% avoid bugs in plotting
d(end,:) = nan;

% input check
if ~exist('index_incre','var')||isempty(index_incre)
    index_incre = 1;
end
dori = d;
d = d(:,1:index_incre:size(d,2));
if ~exist('trace_balance','var')||isempty(trace_balance)
    trace_balance = [];
end
if ~exist('x','var')||isempty(x)
    x = 1:index_incre:size(dori,2);
else
    x = x(1:index_incre:length(x));
end

if ~exist('clip','var')||isempty(clip)
    clip = 3;
end
if ~exist('line_color','var')||isempty(line_color)
    line_color = 'k';
end
if ~exist('face_color','var')||isempty(face_color)
    face_color = 'k';
end
if ~exist('y','var')||isempty(y)
    y = 1:size(d,1);
end
if ~exist('lvl','var')||isempty(lvl)
    lvl = 1;
end
if ~exist('mode','var')||isempty(mode)
    mode = 'hold';
end
if strcmp(mode,'new')
    figure;
    ax =  gca;
elseif strcmp(mode,'hold')
    ax = gca;
elseif strcmp(mode,'wipe')
    cla;
    ax = gca;
else
    error('Invalid mode!');
end

if isscalar(y)
    y = 0:y:(size(d,1)-1)*y;
end

ns = size(d,1);

% draw line
% normalize d
if length(x)>1
    dx = mean(diff(x));
else
    dx = 1;
end

% ---- trace balance ----
if strcmp(trace_balance,'max')
    max_val = max(d);
    for k = 1:size(d,2)
        d(:,k) = d(:,k)/max_val(k);
    end
elseif strcmp(trace_balance,'rms')
    rms_val = rms(d);
    for k = 1:size(d,2)
        d(:,k) = d(:,k)/rms_val(k);
    end
end

fct = (1/max(max(abs(d)))) * (dx*abs(lvl));
clip = dx*clip;
dn = d*fct;
dn(dn>clip) = clip;
dn(dn<-clip) = -clip;

x = x(:)';
dn = dn + repmat(x,ns,1);
y = y(:); y = repmat([y;nan],length(x),1);
dn = [dn;ones(1,size(dn,2))*nan];
ndata = size(dn,1)*size(dn,2);
dns = reshape(dn,ndata,1);
line(dns,y,'color',line_color);
axis ij;
xlim([min(x)-dx,max(x)+dx]);
dy = y(2) - y(1);
ylim([min(y),max(y)-dy]);

if ~strcmp(face_color,'none')
% fill areas
for iter = 1:size(d,2)
    nv = 1; % number of vertex
    kf = 1; % face number
    % for each trace,
    % find samples need to be filled
    tr = dn(:,iter);
    x0 = x(iter);
    if lvl >= 0
        nfil = find(tr>(x0));
    else
        nfil = find(tr<(x0));
    end
    if ~isempty(nfil)
        nbreak = find(diff(nfil)~=1);
        
        grp = zeros(length(nbreak)+1,2);
        grp(:,1) = [nfil(1);nfil(nbreak+1)]; % group startings
        grp(:,2) = [nfil(nbreak);nfil(end)]; % group endings
        
        face = nan*ones(size(grp,1),ns);
        vertex = nan*ones(ns,2);
        
        for k = 1:size(grp,1)
            n01 = grp(k,1);
            n02 = grp(k,2);
            n1 = n01-1; % sample above
            n2 = n02+1; % sample below
            nf = 1; % number of vertex in the current face
            x1 = tr(n01);
            if n1~=0
                x2 = tr(n1);
                y2 = y(n1);
            else
                x2 = x0;
                y2 = y(1);
            end
            y1 = y(n01);
            vertex(nv,:) =  [x0,intercept(x0,x1,y1,x2,y2)];
            face(kf,nf) = nv;
            % go to next vertex and face
            nv = nv + 1;
            nf = nf + 1;
            y12 = y(n01:n02);
            x12 = tr(n01:n02);
            vertex(nv:nv+n02-n01,:) = [x12,y12];
            face(kf,nf:nf+n02-n01) = nv:nv+n02-n01;
            nv = nv + n02 - n01 + 1;
            nf = nf + n02 - n01 + 1;
            x1 = tr(n02);
            if n2~=(nfil(end)+1)
                x2 = tr(n2);
                y2 = y(n2);
            else
                x2 = x0;
                % y2 = y(end-1); Chen 3/4
                if(n2 < (length(y)-1))
                    y2 = y(n2);
                else
                    y2 = y(end-1);
                end
            end
            y1 = y(n02);
            vertex(nv,:) =  [x0,intercept(x0,x1,y1,x2,y2)];
            face(kf,nf) = nv;
            % go to next vertex and face
            nv = nv + 1;
            kf = kf + 1;
        end

        patch('faces',face(1:kf-1,:),'vertices',vertex(1:nv-1,:),...
            'facecolor',face_color,'edgecolor','none','parent',ax);
    end
end
end




end

function y = intercept(x,x1,y1,x2,y2)
y = (y2-y1)/(x2-x1)*(x-x1) + y1;
end

