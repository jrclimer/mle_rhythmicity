function [f,inds,uqinds,inds2] = rhythmicity_pdf(varargin)
%RHYTHMICITY_PDF Parametric distribution of lags for rhythmic neurons
%   RHYTHMICITY_PDF generates the PMF of lags based on the parameters
%   in phat.
%
%   PARAMETERS
%       T [0.6]: The width of the anaysis window
%       t [0:0.01:0.6]: Times to analyze the liklihood
%       noskip (false): If true, skipping removed from model.
%       phat [log10(0.3) 0.2 log10(0.5) 10 0.2(Not present if noskip=true)
%           0.8]: Set of parameters for the parametric model of lags. [tau
%           b c f s(Not present if noskip=true) r].  Each an also be
%           entered individually (e.g., ...'tau',-0.5,...)
%
%           tau (log10(sec)): Exponential falloff rate of the whole
%               distribution
%           b: Baseline probability
%           c (log10(sec)): Falloff rate of the rhythmicity magnitude
%           f (Hz): Frequency of the rhymicity
%           s: Skipping(Not present if noskip=true)
%           r: Rhythmicity
%
%       phat can also be a cell array.  The values can either be scalars,
%       or have the same size as t/inds.
%
%       normalize (true): Normalizes the PDF to be a probability
%           distribution (The integral over the window is 1). If false, the
%           value at 0 is 1.
%       dt: 0.001 - Number of steps.
%       inds: Indecies of spike times. Precalculated, adds a lot of speed.
%
%       If the following are given, greatly speeds up.
%       uqinds: The indecies of unique values of the parameters (nonscalar)
%       inds2: The uniques indecies of values of the parameters (nonscalar)
%
%   RETURNS
%       f: The liklihood of each lag
%       inds: The index of each lag (from being binned by dt)
%       uqinds: The indecies of unique values of the parameters (nonscalar)
%       inds2: The uniques indecies of values of the parameters (nonscalar)
% RELEASE NOTES
%   v0.1 2014-04-30 Updated from rhythmicity_pdf
%   v0.3 2014-07-27 Release for review
%
% Copyright (c) 2014, Trustees of Boston University
% All rights reserved.
%
% This file is part of mle_rhythmicity
%
% This code has been freely distributed by the authors under the BSD 
% licence (http://opensource.org/licenses/BSD-2-Clause). If used or
% modified, we would appreciate it if you cited our paper:
%
% Climer, J. R., DiTullio, R., Newman, E. L., Hasselmo, M. E., Eden, U. T. 
% (2014), Examination of rhythmicity of extracellularly recorded neurons in
% the entorhinal cortex. Hippocampus, Epub ahead of print. doi:
% 10.1002/hipo.22383.
persistent FIELDS

if ~exist('FIELDS','var')||isempty(FIELDS)
    FIELDS = {'tau','b','c','f','s','r'};
end

%% Parse input
ip = inputParser;
ip.CaseSensitive = true;
ip.addParamValue('T',0.6,@(x)isnumeric(x)&&isscalar(x));
ip.addParamValue('t',linspace(0,0.6,1000),@(x)isnumeric(x));
ip.addParamValue('phat',{log10(0.3) 0.2 log10(0.5) 10 0.2 0.8});
ip.addParamValue('normalize',true);
ip.addParamValue('noskip',false);
ip.addParamValue('dt',0.001,@(x)isnumeric(x)&&isscalar(x));
ip.addParamValue('inds',[]);
ip.addParamValue('inds2',[]);
ip.addParamValue('uqinds',[]);

for k = FIELDS
    ip.addParamValue(k{1},NaN);
end
ip.parse(varargin{:});
T = ip.Results.T;
t = ip.Results.t;
noskip = ip.Results.noskip;
phat = ip.Results.phat;
normalize = ip.Results.normalize;
dt = ip.Results.dt;
inds = ip.Results.inds;
inds2 = ip.Results.inds2;
uqinds = ip.Results.uqinds;

t0 = 0:dt:T;

if isempty(inds), [~,inds] = histc(t,t0); end;
inds(inds==0) = 1;

for i=1:numel(FIELDS)
    if ~isnan(ip.Results.(FIELDS{i}))
        phat{i} = ip.Results.(FIELDS{i});
        eval([FIELDS{i} ' = ip.Results.' FIELDS{i} ';']);
    else
        if iscell(phat)
            if ~noskip||i<5
                eval([FIELDS{i} ' = phat{i};']);
            elseif i==5
                s = 0;
            elseif i==6
                r = phat{end};
            end
        else
            if ~noskip||i<5
                eval([FIELDS{i} ' = phat(i);']);
            elseif i==5
                s = 0;
            elseif i==6
                r = phat(end);
            end
        end
    end
end

eval(['phat={' strjoin(FIELDS,' ') '};']);

b(b<0) = 0;


%%
if iscell(phat)&&all(cellfun(@isscalar,phat)) % All scalars, cheap
    % Parametric distribution
    if s==0
    fx = @(t)(1-b).*exp(-t.*10.^(-tau)).*(r.*exp(-t.*10.^(-c)).*...
        cos(2*pi*f*t)...
        +1)+b;
    else
        fx = @(t)(1-b).*exp(-t.*10.^(-tau)).*(r.*exp(-t.*10.^(-c)).*...
        ((2+2*sqrt(1-s)-s)*cos(2*pi*f*t)+4*s*cos(pi*f*t)+2-2*sqrt(1-s)-3*s)/4 ...
        +1)+b;
    end
    
    phat = cell2mat(phat);
    f = fx(t0);
    if normalize
        f = f/sum(f);
    end
    f = f(inds);
    
    f(f<=0|isnan(f)|isinf(f)) = realmin;
    
    if (ismember('t',ip.UsingDefaults)&&~isequal(size(inds),size(f)))||...
            (~isequal(size(t),size(f)))
        f = f';
    end
else % Not all scalars  - conditional parametarization
    if isempty(inds2) % Calculate unique sets if not done already
        [vals,uqinds,inds2] = unique(cat(2,phat{~cellfun(@isscalar,phat)}),'rows');        
    else
        vals = cat(2,phat{~cellfun(@isscalar,phat)});
        vals = vals(uqinds,:);
    end
    
    % Reshape for full distributions
    temp = FIELDS(~cellfun(@isscalar,phat));
    for j=1:numel(temp)
       eval([temp{j} '=repmat(vals(:,j),[1 numel(t0)]);']); 
    end
    t = repmat(t0,[size(vals,1) 1]);
    
    if ~isscalar(s)
	fx = zeros(size(s));
    % Find full distributions
    fx(s==0) = (1-b(s==0)).*exp(-t(s==0).*10.^(-tau(s==0))).*(r(s==0).*exp(-t(s==0).*10.^(-c(s==0))).*...
        cos(2*pi*f(s==0).*t(s==0))...
        +1)+b(s==0);
    fx(s>0) = (1-b(s>0)).*exp(-t(s>0).*10.^(-tau(s>0))).*(r(s>0).*exp(-t(s>0).*10.^(-c(s>0))).*...
        ((2+2*sqrt(1-s(s>0))-s(s>0))*cos(2*pi*f(s>0)*t)+4*s(s>0)*cos(pi*f(s>0)*t)+2-2*sqrt(1-s(s>0))-3*s(s>0))/4 ...
        +1)+b(s>0);
    else
       if s==0
           fx = (1-b).*exp(-t.*10.^(-tau)).*(r.*exp(-t.*10.^(-c)).*...
        cos(2*pi*f.*t)...
        +1)+b;
       else
           fx = (1-b).*exp(-t.*10.^(-tau)).*(r.*exp(-t.*10.^(-c)).*...
        ((2+2*sqrt(1-s)-s)*cos(2*pi*f.*t)+4*s*cos(pi*f.*t)+2-2*sqrt(1-s)-3*s)/4 ...
        +1)+b;
       end
    end
    
    norms = sum(fx,2);
    
    % find values
    f = fx(sub2ind(size(fx),inds2,inds))./norms(inds2);
    f = f(:);
    
    f(f<=0) = realmin;
    
end
end

