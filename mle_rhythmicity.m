function [ stats ] = mle_rhythmicity( x,session_dur,varargin )
%MLE_RHYTHMICITY Use maximum liklihood estimation to find rhythmicity
%parameters
%
% INPUT
%   x: Data. Can either be spike timestamps, or a cell array of lags
%   session_dur: Duration of the session
%
% PARAMETERS
%   T (0.6): Examinination window
%   plotit (true): If true, plots histogram and distribution estimation
%   noskip (false): If true, skipping removed from model.
%   dt (0.001): dt for model
%
% OUTPUT
%   stats: Output structure of the fit
%       Fs: Firing frequency (Hz)
%       FsX: Multiplier of the firing rate in the examination window
%       Fsci: 95% confidence intervals for the firing frequency
%       FsX: 95% confidence intervals for the multiplier
%       LL: Log liklihood of fit
%       LLs: Log l
%       a: Maximum liklihood estimator for the magnitude of the rhythmicity
%       aci: 95% confidence intervals for the rhythmicity
%       ci: 95% confidence intervals for the parameters
%       flat_LL: Log liklihood of the arrhythmic fit
%       noskip_LL: Log liklihood of the non-skipping fit
%       noskip_phat: Maximum liklihood estimators for the non skipping fit
%           [tau b c f r] (Not present if noskip=true)
%       D_rhyth: Deviance of rhythmic versus arhythmic fit
%       D_sk: Deviance of the non-skipping versus full fit. (Not present if
%           noskip=true)
%       p_rhyth: p that the cell is rhythmic (chi2-test).
%       p_sk: p that the cell is skipping (chi2-test) (Not present if
%           noskip=true)
%       phat: Maximum liklihood estimators [tau b c f s(Not present if
%           noskip=true) r]
%       x: Cell array of the lags following each spike
%       x_: All lags
%
% RELEASE NOTES
%   v0.1 2014-04-30 Updated from mle_rhythmicity
%   v0.2 2014-05-02 Beta release
%   v0.21 2014-05-06 Bug fixes with noskip=false
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

%% Warnings
warning('off','stats:mle:IterLimit');


%% Parse inputs
ip = inputParser;
ip.addParamValue('T',0.6);
ip.addParamValue('plotit',true);
ip.addParamValue('noskip',false);
ip.addParamValue('dt',0.001);
ip.addParamValue('stats0',[]);
ip.parse(varargin{:});
for j = fields(ip.Results)'
    eval([j{1} ' = ip.Results.' j{1} ';']);
end

% If timestamps, find lags
if ~iscell(x)
    x_ = x(:)';
    x = cell(size(x_));
    for i=1:numel(x)
        x{i} = x_(x_>x_(i)&x_<=x_(i)+T)-x_(i);
    end
end

[Fs,Fsci]=poissfit(numel(x));
Fs = Fs/session_dur;
Fsci = Fsci/session_dur;
[FsX,FsX_ci] = poissfit(cellfun(@numel,x));
FsX = (FsX/T)/Fs;
FsX_ci = (FsX_ci/T)/Fs;
try
x_ = [x{:}];
catch err
    x_ = cat(1,x{:});
end
x_ = x_(:);

t0 = 0:dt:T;
[~,inds] = histc(x_,t0);


%% Initial guess using particle swarm
PopulationSize = 75;
lbnd = [-1 0 -1 1 0];
ubnd = [1 1 1 13 1];
InitialPopulation = rand(PopulationSize,numel(lbnd)).*repmat(ubnd-lbnd,[PopulationSize 1])+repmat(lbnd,[PopulationSize 1]);
%
% close all;
phat0 = pso(@(phat)-sum(log(rhythmicity_pdf('inds',inds,'phat',phat,'noskip',true))),size(lbnd,2) ...
    ,[],[],[],[]...
    ,lbnd...
    ,ubnd...
    ,[] ...
    ,psooptimset('Generations',100,'InitialPopulation',InitialPopulation,'PopulationSize',PopulationSize,'ConstrBoundary','reflect' ...
    ...,'PlotFcns',{@psoplotswarm}...
    ,'Display','off'...
    ) ...
    );

% Solve
[noskip_phat,ci]  = mle(x_,'pdf',@(x,tau,b,c,f,r,varargin)rhythmicity_pdf('inds',inds,'phat',[tau b c f r],'noskip',true),....
        'start',phat0,....
        'lowerbound',[-inf 0 -inf 1 0],....
        'upperbound',[inf 1 inf 13 1],....
        'alpha',0.05);
noskip_LL = sum(log(rhythmicity_pdf('inds',inds,'phat',noskip_phat,'noskip',true)));

if noskip
    phat = noskip_phat;
    LL = noskip_LL;
    clear noskip_phat;
    clear noskip_LL;
else % Need to do full fit
    % Initial guess using particle swarm
    lbnd = [lbnd 0];
    ubnd = [ubnd 1];
    InitialPopulation = [InitialPopulation rand(PopulationSize,1)];
    
%     close all;
    phat0 = pso(@(phat)-sum(log(rhythmicity_pdf('inds',inds,'phat',phat))),size(lbnd,2) ...
        ,[],[],[],[]...
        ,lbnd...
        ,ubnd...
        ,[] ...
        ,psooptimset('Generations',100,'InitialPopulation',InitialPopulation,'PopulationSize',PopulationSize,'ConstrBoundary','reflect' ...
        ....,'PlotFcns',{@psoplotswarm}...
        ,'Display','off'...
        ) ...
        );
    
    % fit
    [phat,ci] = mle(x_,'pdf',@(x,tau,b,c,f,s,r,varargin)rhythmicity_pdf('inds',inds,'phat',[tau b c f s r]),....
        'start',phat0,....
        'lowerbound',[-inf 0 -inf 1 0 0],....
        'upperbound',[inf 1 inf 13 1 1],....
        'alpha',0.05);
    
    % Testing
    LL = sum(log(rhythmicity_pdf('inds',inds,'phat',phat)));
    
    D_sk = 2*(LL-noskip_LL);
    p_sk = 1-chi2cdf(D_sk,1);
end

% Arrhythmic fit
flat_phat = mle(x_,'pdf',@(x,tau,b,varargin)rhythmicity_pdf('inds',inds,'phat',[tau,b,1,1,0,0]),....
        'start',phat(1:2),....
        'lowerbound',[-inf 0],....
        'upperbound',[inf 1],....
        'optimfun','fmincon',...
        'alpha',0.05);
flat_LL = sum(log(rhythmicity_pdf('inds',inds,'phat',[flat_phat,1,1,0,0])));
%% More testing
D_rhyth = 2*(LL-flat_LL);
p_rhyth = 1-chi2cdf(D_rhyth,size(lbnd,2)-2);

%% Amplitude stuff
[a,aci] = mle(x_,'pdf',@(x,a,varargin)rhythmicity_pdf('inds',inds,'phat',[phat(1:end-1) a/(1-phat(2))],'noskip',noskip),...
    'start',(1-phat(2))*phat(end),...
    'lowerbound',0,....
    'upperbound',1-phat(2));
    
%% Plotting
if plotit
    %
    if noskip
       phat = [phat(1:end-1) 0 phat(end)]; 
    end
    
    ts = linspace(0,T,61);
    hist(x_,ts);
    xlim([min(ts) max(ts)]);
    
    hold on;
    ps = rhythmicity_pdf('t',ts,'phat',phat);
    ps = ps/sum(ps)*numel(x_);
    plot(ts,ps,'Color',[1 0 0],'LineWidth',2);
    
    ps = rhythmicity_pdf('t',ts,'phat',[phat(1:2) phat(3:end)*0]);
    ps = ps/sum(ps)*numel(x_);
    plot(ts,ps,'c--','LineWidth',2);
    hold off
    legend('data','MLE','flat');
    
    
    if noskip
        phat = [phat(1:end-2) phat(end)];
    end
        
    if noskip
        title(['a=' sprintf('%2.2g',a) ', p_{rhyth}=' sprintf('%2.2g',p_rhyth)]);
    else
        title(['a=' sprintf('%2.2g',a) ', p_{rhyth}=' sprintf('%2.2g',p_rhyth) ', p_{skip}=' sprintf('%2.2g',p_sk)]);
    end
end

se = diff(ci)/2/norminv(1-0.05/2);

%% Format output
clear c;clear i;clear j;clear LLs;clear T;clear f;clear faxis;clear ip;clear phats;
clear session_dur;clear varargin;clear plotit;clear ps;clear ans;clear InitialPopulation;clear PopulationSize;clear dt;clear lbnd;clear noskip;
clear phat0; clear t0;clear ts;clear ubnd;clear stats0;
c = who;
stats = struct;
for i = c'
    eval(['stats.' i{1} '=' i{1} ';']);
end

warning('on','stats:mle:IterLimit');
end
