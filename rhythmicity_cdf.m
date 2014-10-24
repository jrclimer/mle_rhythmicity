function [ F ] = rhythmicity_cdf( varargin )
%RHYTHMICITY_CDF Parametric cumulative distribution of lags for rhythmic neurons
%   RHYTHMICITY_PDF generates the cumulative liklihood of lags based on the 
%   parameters in phat.
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
%       normalize (true): Normalizes the PDF to be a probability
%           distribution (The integral over the window is 1). If false, the
%           value at 0 is 1.
%
% RELEASE NOTES
%   v0.1 2014-04-30 Updated from rhythmicity_cdf
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
%% Parse input
ip = inputParser;
ip.CaseSensitive = true;
ip.KeepUnmatched = true;
ip.addParamValue('t',linspace(0,0.6,1000),@(x)isnumeric(x)&&all(x(:)>=0));
ip.addParamValue('T',0.6,@(x)isscalar(x)&&isnumeric(x));
ip.addParamValue('dt',0.001,@(x)isscalar(x)&&isnumeric(x));
ip.addParamValue('inds',[]);
ip.parse(varargin{:});

t = ip.Results.t;
T = ip.Results.T;
dt = ip.Results.dt;

if ~ismember('t',ip.UsingDefaults)
    temp = true(size(varargin));
    temp(find(cellfun(@(x)isequal(x,'t'),varargin))) = false;
    temp(find(~temp)+1) = false;
    varargin = varargin(temp);
    clear temp;
end

t0 = 0:dt:T;

inds = ip.Results.inds;
    
if ~ismember('inds',ip.UsingDefaults)
    temp = true(size(varargin));
    temp(find(cellfun(@(x)isequal(x,'inds'),varargin))) = false;
    temp(find(~temp)+1) = false;
    varargin = varargin(temp);
    clear temp;    
end

if isempty(inds)
    [~,inds] = histc(t,t0);
end

F = rhythmicity_pdf('inds',1:floor(T/dt)+1,varargin{:});

% Integrate CDF
F=arrayfun(@(i)sum(F(1:i)),inds);

end

