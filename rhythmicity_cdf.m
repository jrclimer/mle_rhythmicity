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
% This code has been freely distributed by the authors. If used or
% modified, we would appreciate it if you cited our paper, <paper
% information here>
%
% Copyright (c) 2014, Trustees of Boston University
% All rights reserved.
%
% This file is part of mle_rhythmicity
%
% This version of mle_rhythmicity is solely for the purposes of review and
% demonstration of its functionality by editors and reviewers.
% Redistribution to others and other uses in source and binary forms, with
% or without modification, is prohibited.  Upon acceptance for publication,
% this code and future versions modified by the authors will be made
% permanently available on GitHub under the BSD
% license (Available at http://opensource.org/licenses/bsd-license.php)
% allowing future users to freely distribute and modify the code.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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

