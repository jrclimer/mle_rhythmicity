function [plr,pw,pa,plr_ip,plr_a,D,W,Wa,D_ip,D_a,stats1,stats2] = compare_mle_rhythmicity(stats1,dur1,stats2,dur2,varargin)
%COMPARE_MLE_RHYTHMICITY Compare MLE rhythmicity results
%
% INPUT
%   stats1: First dataset. Can either be input for mle_rhythmicity or the
%       output struct of mle_rhythmicity
%   dur1: Duration of the first session
%   stats2: Second dataset. Can either be input for mle_rhythmicity of the
%       output struct of mle_rhythmicity
%   dur2: Duration of the second session
%
% OUTPUT
%   plr = Significance of overall liklihood-ratio test
%   pw = Significance of Wald tests
%   pa = Significance of Wald test on amplitude
%   plr_ip = Significance of likelihood-ratio tests on individual parameters allowing all
%       others to vary
%   plr_a = Significance of likelihood-ratio test on amplitude allowing all other paremeters
%       to vary
%   D = Deviance (D) for of overall D-test
%   W = Wald statsitic for of Wald tests
%   Wa = Wald statsitic for Wald test on amplitude
%   D_ip = D statsitics for likelihood-ratio tests on individual parameters allowing all
%       others to vary
%   D_a = Deviance (D) for likelihood-ratio test on amplitude allowing all other paremeters
%       to vary
%   stats1 = Results of mle_rhythmicity(stats1,dur1,varargin{:}) (or stats1 if stats1
%       is a struct)
%   stats2 = Results of mle_rhythmicity(stats2,dur2,varargin{:}) (or stats2 if stats2
%       is a struct)
%
% PARAMETERS - options for mle_rhythmicity fits for stats1 and stats2 if stats1 and stats2
%   aren't structs. See doc mle_rhythmicity for more details.
%
% RELEASE NOTES
%   v0.1 2014-04-30 Updated from mle_rhythmicity
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
% This code has been freely distributed by the authors under the BSD 
% licence (http://opensource.org/licenses/BSD-2-Clause). If used or
% modified, we would appreciate it if you cited our paper:
%
% Climer, J. R., DiTullio, R., Newman, E. L., Hasselmo, M. E., Eden, U. T. 
% (2014), Examination of rhythmicity of extracellularly recorded neurons in
% the entorhinal cortex. Hippocampus, Epub ahead of print. doi:
% 10.1002/hipo.22383.

if ~isstruct(stats1)
    x1 = stats1;
    x2 = stats2;
    
    stats1 = mle_rhythmicity(x1,varargin{:},dur1);
    stats2 = mle_rhythmicity(x2,varargin{:},dur2);
end

noskip = numel(stats1.phat)==5;

%%
if noskip
    NVAR = 5;
else
    NVAR = 6;
end
%%
% keyboard
%%
PopulationSize = 51;
PopulationSize = ceil(PopulationSize/3)*3;
lbnd = [-1 0 -1 1 0];
ubnd = [1 1 1 13 1];
if ~noskip
   lbnd = [lbnd 0];
   ubnd = [ubnd 1];
end

% InitialPopulation = rand(PopulationSize,numel(lbnd)).*repmat(ubnd-lbnd,[PopulationSize 1])+repmat(lbnd,[PopulationSize 1]);

InitialPopulation = normrnd(...
    [repmat(stats1.phat,[PopulationSize/3 1]);repmat(stats2.phat,[PopulationSize/3 1]);repmat(mean([stats1.phat;stats2.phat]),[PopulationSize/3 1])],...
    repmat(ubnd-lbnd,[PopulationSize 1])/4);

temp = 2*repmat(lbnd,[PopulationSize 1])-InitialPopulation;
InitialPopulation(InitialPopulation<repmat(lbnd,[PopulationSize 1])) = ...
    temp(InitialPopulation<repmat(lbnd,[PopulationSize 1]));

temp = 2*repmat(ubnd,[PopulationSize 1])-InitialPopulation;
InitialPopulation(InitialPopulation>repmat(ubnd,[PopulationSize 1])) = ...
    temp(InitialPopulation>repmat(ubnd,[PopulationSize 1]));

%%
close all;
phat0 = pso(@(phat)-sum(log(rhythmicity_pdf('inds',[stats1.inds;stats2.inds],'phat',phat,'noskip',noskip))),size(lbnd,2) ...
    ,[],[],[],[]...
    ,lbnd...
    ,ubnd...
    ,[] ...
    ,psooptimset('Generations',100,'InitialPopulation',InitialPopulation,'PopulationSize',PopulationSize,'ConstrBoundary','reflect' ...
    ...,'PlotFcns',{@psoplotswarm}...
    ,'Display','off'...
    ) ...
    );
%%
phat = mle([stats1.x_;stats2.x_],'pdf',@(x,varargin)rhythmicity_pdf('inds',[stats1.inds(:);stats2.inds(:)],'phat',[varargin{1:numel(ubnd)}],'noskip',noskip),...
    'start',phat0,...
    'lowerbound',[-inf 0 -inf 1 0 0],....
    'upperbound',[inf 1 inf 13 1 1],....
    'alpha',0.05);

LL = sum(log(rhythmicity_pdf('inds',[stats1.inds(:);stats2.inds(:)],'phat',phat,'noskip',noskip)));
D = 2*(stats1.LL+stats2.LL-LL);
plr = 1-chi2cdf(D,NVAR);

%%
lbnd = [-inf 0 -inf 1 0 0];
ubnd = [inf 1 inf 13 1 1];

vars = {'tau','b','c','f','s','r'};

if noskip
    lbnd = [lbnd(1:end-2) lbnd(end)];
    ubnd = [ubnd(1:end-2) ubnd(end)];
    vars = [vars(1:end-2) vard(end)];
end

vars1 = cellfun(@(x)[x '1'],vars,'UniformOutput',false);
vars2 = cellfun(@(x)[x '2'],vars,'UniformOutput',false);

inds = [stats1.inds(:);stats2.inds(:)];

fun = @(x,varargin)(x<=numel(stats1.x_)).*rhythmicity_pdf('inds',inds(x),'phat',[varargin{1:NVAR-1} varargin{end}],'noskip',noskip)+...
    (x>=numel(stats1.x_)).*rhythmicity_pdf('inds',inds(x),'phat',[varargin{NVAR:end}],'noskip',noskip);

phat = mle(1:numel(inds),'pdf',fun,...
    'start',[stats1.phat(1:end-1) stats2.phat(1:end-1) mean([stats1.a stats2.a])],...
    'lowerbound',[repmat(lbnd(1:end-1),[1 2]) 0],...
    'upperbound',[repmat(ubnd(1:end-1),[1 2]) 1]);
%%
temp = num2cell(phat);
LL = sum(log(fun((1:numel(inds))',temp{:})));
D_a = 2*(stats1.LL+stats2.LL-LL);
plr_a = 1-chi2cdf(D_a,1);
%%
D_ip = zeros(1,NVAR);
plr_ip = zeros(1,NVAR);
%%
for j=1:NVAR
   fun = @(x,varargin)(x<=numel(stats1.x_)).*rhythmicity_pdf('inds',inds(x),'phat',[varargin{1:j-1} varargin{end} varargin{j:NVAR-1}],'noskip',noskip)+...
       (x>numel(stats1.x_)).*rhythmicity_pdf('inds',inds(x),'phat',[varargin{NVAR:NVAR+j-2} varargin{end} varargin{NVAR+j-1:end-1}],'noskip',noskip);

   phat = mle(1:numel(inds),'pdf',fun,...
    'start',[stats1.phat([1:j-1 j+1:end]) stats2.phat([1:j-1 j+1:end]) mean([stats1.phat(j) stats2.phat(j)])],...
    'lowerbound',[repmat(lbnd([1:j-1 j+1:end]),[1 2]) lbnd(j)],...
    'upperbound',[repmat(ubnd([1:j-1 j+1:end]),[1 2]) ubnd(j)]);

    temp = num2cell(phat);
    LL = sum(log(fun((1:numel(inds))',temp{:})));
    D_ip(j) = 2*(stats1.LL+stats2.LL-LL);
    plr_ip(j) = 1-chi2cdf(D_ip(j),1);
    
end

W = (stats2.phat-stats1.phat)./...
    sqrt(mean(([diff(stats1.ci);diff(stats2.ci)]/2/normcdf(1-0.05/2)).^2));
pw = normcdf(-abs(W))*2;

Wa = (stats2.a-stats1.a)./...
    sqrt(mean(([diff(stats1.aci);diff(stats2.aci)]/2/normcdf(1-0.05/2)).^2));
pa = normcdf(-abs(Wa))*2;
end