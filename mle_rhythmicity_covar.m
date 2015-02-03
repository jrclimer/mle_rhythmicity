function [ stats ] = mle_rhythmicity_covar( x,session_dur,y_ts,y0,varargin )
%MLE_RHYTHMICITY Use maximum liklihood estimation to find rhythmicity
%parameters and how much the frequency coveries with a parameter y
%
% INPUT
%   x: Timestamps of spikes
%   y_ts: Tiemstamps of covarying parameter
%   y: The covarying parameter
%
% PARAMETERS
%   T (0.6): Examinination window
%   plotit (true): If true, plots histogram and distribution estimation
%   noskip (false): If true, skipping removed from model.
%   dt (0.001): dt for model
%   stats0 ([]): Output of mle_rhythmicity_fast
%   precision (1e-2): Precision of coviariate (at higher values (less
%       precise, but runs much faster)
%
% OUTPUT
%   stats: Output structure of the fit
%       D: The deviance of the fit with frequency covarying with y
%       LL: Log liklihood of fit with with frequency covarying with y
%       ci: 95% confidence intervals for the parameters
%       f0: The mean frequency (from stats0.phat(4))
%       f_b: The frequency y-intercept
%       f_m: The frequency slope
%       p: The significance of the deviance
%       phat:  Maximum liklihood estimators [tau b c fb fm s(Not present if
%           noskip=true) r].  Note that f is now split into fb and fm,
%           where fb is the Y-intercept for the correlation between speed
%           and rhythmicity and fm is the slope.
%       x: Cell array of the lags following each spike
%       x_: All lags
%       y: Cell array of the corvariate for each lag
%       y_: Covariate for all lags
%
% RELEASE NOTES
%   v0.1 2014-04-30 Updated from mle_rhythmicity
%   v0.2 2014-05-02 Beta release
%   v0.21 2014-05-06 Updated bug fixes with noskip=false and plotting
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
% permanently available at
% https://mind4.bu.edu/distributed_code/mle_rhythmicity_fast/ under the BSD
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

%% Parse inputs
ip = inputParser;
ip.addParamValue('T',0.6);
ip.addParamValue('plotit',true);
ip.addParamValue('noskip',false);
ip.addParamValue('precision',1e-2);
ip.addParamValue('dt',0.001);
ip.addParamValue('stats0',[]);
ip.addParamValue('ts2',[]);
ip.parse(varargin{:});
for j = fields(ip.Results)'
    eval([j{1} ' = ip.Results.' j{1} ';']);
end

% Do initial fit if necessary
if isempty(stats0), stats0 = mle_rhythmicity_fast(x,session_dur,varargin{:}); end

% Format data
y = interp1(y_ts,y0,x,'nearest','extrap');
if ~isempty(ts2)
    yts2 = interp1(y_ts,y0,ts2,'nearest','extrap');
    yts2 = yts2(:)';
end
    
spk_y = y;
x_ = x(:)';
y_ = y(:)';

x = cell(size(x_));
y = x;

if isempty(ts2)
    for i=1:numel(x)
        y{i} = y_(x_>x_(i)&x_<=x_(i)+T);
        x{i} = x_(x_>x_(i)&x_<=x_(i)+T)-x_(i);
    end
else
    for i=1:numel(x)
        y{i} = yts2(ts2>x_(i)&ts2<=x_(i)+T);
        %x{i} = x_(x_>x_(i)&x_<=x_(i)+T)-x_(i);
        x{i} = ts2(ts2>x_(i)&ts2<x_(i)+T)-x_(i);
    end
end
y(cellfun(@(x) length(x)==0, x))= [];
x(cellfun(@(x) length(x)==0, x)) = [];

x_ = cell2mat(cellfun(@(x) x(:)', x,'UniformOutput',0));
x_ = x_(:);

y_ = [y{:}];
y_ = y_(:);
y1 = y_;
y_ = round(y_/precision)*precision;

t0 = 0:dt:T;
[~,inds] = histc(x_,t0);

%% Initial guess from quintiles of lags
n = 5;

f = zeros(n,1);
ys = n*f;

for i=1:n
   f(i) = mle(inds(y_>=quantile(y_,(i-1)/n)&y_<=quantile(y_,i/n)),...
       'pdf',@(inds,f)rhythmicity_pdf_fast('inds',inds,'phat',[stats0.phat(1:3) f stats0.phat(5:end)],'noskip',noskip,'T',T,'dt',dt),...
       'start',stats0.phat(4));
   ys(i) = mean(y_(y_>=quantile(y_,(i-1)/n)&y_<=quantile(y_,i/n)));
end
ys = [ys./ys ys];
f = inv(ys'*ys)*ys'*f;

%% Set up uqings and inds2 (makes faster)
[~,~,uqinds,inds2]=rhythmicity_pdf_fast('inds',inds,'phat',[num2cell(stats0.phat(1:3)) {f(1)+f(2)*y_} num2cell(stats0.phat(5:end))],'noskip',noskip,'T',T,'dt',dt);

%% Fit
phat0 = mle(inds,'pdf',@(inds,fb,fm)rhythmicity_pdf_fast('inds',inds,'inds2',inds2,'uqinds',uqinds,'phat',[num2cell(stats0.phat(1:3)) {fb+fm*y_} num2cell(stats0.phat(5:end))],'noskip',noskip,'T',T,'dt',dt),...
    'start',f',...
    'lowerbound',[-2 -0.5],...
    'upperbound',[15 0.5],...
    'alpha',0.05);%,'options',statset('display','iter'));

%
if noskip
   [phat,ci] = mle(inds,'pdf',@(inds,tau,b,c,fb,fm,r)rhythmicity_pdf_fast('inds',inds,'inds2',inds2,'uqinds',uqinds,'phat',{tau,b,c,fb+fm*y_,r},'noskip',true,'T',T,'dt',dt),...
       'start',[stats0.phat(1:3) phat0 stats0.phat(5:end)],...
       'lowerbound',[-inf 0 -inf -inf -inf 0],....
        'upperbound',[inf 1 inf inf inf 1],....
        'alpha',0.05);%,'options',statset('display','iter'));    
else
    [phat,ci] = mle(inds,'pdf',@(inds,tau,b,c,fb,fm,s,r)rhythmicity_pdf_fast('inds',inds,'inds2',inds2,'uqinds',uqinds,'phat',{tau,b,c,fb+fm*y_,s,r},'noskip',true,'T',T,'dt',dt),...
       'start',[stats0.phat(1:3) phat0 stats0.phat(5:end)],...
       'lowerbound',[-inf 0 -inf -inf -inf 0 0],....
        'upperbound',[inf 1 inf inf inf 1 1],....
        'alpha',0.05);%,'options',statset('display','iter'));   
end

%% Extract values and do tests
f_b = phat(4);
f_m = phat(5);
f0 = stats0.phat(4);

LL = sum(log(rhythmicity_pdf_fast('inds',inds,'inds2',inds2,'uqinds',uqinds,'phat',[num2cell(phat(1:3)) {phat(4)+phat(5)*y_} num2cell(phat(6:end))],'noskip',noskip,'T',T,'dt',dt,'T',T,'dt',dt)));
D = 2*(LL-stats0.LL);
p = 1-chi2cdf(D,1);

%% Plotting
if plotit
    %%
    clf;
        
    tbin = 0:dt:T;
    
    vbin = [floor(quantile(y_,0.03)/10)*10 ceil(quantile(y_,0.97)/10)*10];
    
    %
    for i=1:ceil(range(vbin/5))
      if median(hist(y_,vbin(1):i:vbin(2)))>250
          break
      end
    end
    
    if i>3
        i = round(i/5)*5;
    end
    
    vbin = vbin(1):i:vbin(2);
    
    % vbin = 0:10:150;
    cnts = histcn([y_ x_],vbin,tbin);
   
    
   nwaves = ceil((phat(4)+150*phat(5))*0.6);
   
   %subplot(2,1,1);
   figure
   plot(x_,y_,'ko','MarkerSize',1);
   hold on;
   for i=1:nwaves
       plot((i-1)./(phat(4)+vbin*phat(5)),vbin,'--','Color','r','LineWidth',2);
   end
   hold off;
   set(gca,'YDir','Normal');
   xlim(minmax(tbin));ylim(minmax(vbin));
   title(['f_b=' sprintf('%2.2g',phat(4)) ',f_m=' sprintf('%2.2g',phat(5)) ',p=' sprintf('%2.2g',p)])
   xlabel('Lag (s)');ylabel('Speed (cm/s)');
   
   for i=1:nwaves
      break 
   end
   
   %subplot(2,1,2);   
   figure
   myfilter = fspecial('gaussian',[3 3], 0.5);
   imagesc(tbin,vbin,imfilter(diag(sum(cnts,2).^-1)*cnts,myfilter,'replicate'));
   set(gca,'YDir','Normal');
  
   hold on;
   for i=1:nwaves
       plot((i-1)./(phat(4)+vbin*phat(5)),vbin,'--','Color','k','LineWidth',2);
   end
   xlabel('Lag (s)');ylabel('Speed (cm/s)');
   
   hold off;
end
%% Format output
clear T; clear cnts; clear dr; clear i; clear inds; clear inds2; clear ip;
clear j; clear myfilter; clear n; clear noskip; clear nwaves; clear phat0; clear plotit;
clear precision; clear spk_y; clear t0; clear tbin; clear uqinds; clear varargin;
clear vbin; clear y0; clear y1; clear y_ts; clear ys;clear dt;clear c; clear stats;clear f;

c = who;
stats = struct;
for i = c'
    eval(['stats.' i{1} '=' i{1} ';']);
end

end

