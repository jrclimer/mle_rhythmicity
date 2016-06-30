function [ cif, cif_int ] = cif_generator( funname )
%CIF_GENERATOR Generate function handles for rhythm conditional intensity 
% This helper function generates function handles for fitting rhythmic
% distributions of lags. Using functional programming in this way greatly
% speeds up evaluation during runtime.
%
% INPUT
% funname: String name of the CIF handles to be generated. Can be:
%   flat - Exponential decay with no oscillation
%   pure - Pure sinusoidal decay
%   noskip - Decaying rhythm with no skipping
%   full - Decaying rhythm with skipping
%
% OUTPUT
%   cif - Function handle for the Conditional Intensity Function (CIF).
%       with the form @(lags_list,...). lags_list is the lags at which
%       to evaluate the CIF, and it is followed by the appropriate
%       parameters
%   cif_int - The function handle for the definite integral of the CIF,
%       with the form @(max_lag,...). 
%
% Copyright 2015-2016 Trustees of Boston University
% All rights reserved.
%
% This file is part of mle_rhythmicity revision 2.0. The last committed
% version of the previous revision is the SHA starting with 93862ac...
%
% This code has been freely distributed by the authors under the BSD
% license (http://opensource.org/licenses/BSD2-Clause). If used or
% modified, we would appreciate if you cited our papers:
%
% Climer JR, DiTullio R, Newman EL, Hasselmo ME, Eden UT. (2014),
% Examination of rhythmicity of extracellularly recorded neurons in the
% entorhinal cortex. Hippocampus, 25:460-473. doi: 10.1002/hipo.22383.
%
% Hinman et al., Multiple Running Speed Signals in Medial Entorhinal
% Cortex, Neuron (2016). http://dx.doi.org/10.1016/j.neuron.2016.06.027

switch funname
    case 'flat'% Non-rhythmic
        cif = @(lags_list,tau,b)(1-b)*exp(-lags_list*10^-tau)+b;
        cif_int = @(max_lag,tau,b)10^tau*(1-b)*(1-exp(-max_lag*10^-tau))+b*max_lag;
    case 'pure'% Pure sinusoid
        cif = @(lags_list,f)cos(2*pi*f*lags_list)+1;
        cif_int = @(max_lag,f)sin(2*pi*f*max_lag)/(2*pi*f)+max_lag;
    case 'noskip'% Non-skipping
        cif = @(lags_list,tau,b,c,f,r)(1-b).*exp(-lags_list.*10.^-tau).*(r.*exp(-lags_list.*10.^-c).*...
           cos(2*pi*f.*lags_list)...
           +1)+b;
       cif_int = @(max_lag,tau,b,c,f,r)b*max_lag+(1-b).*10.^tau.*(1-exp(-max_lag*10.^-tau))+r.*(1-b).*...
           exp(-max_lag*(10.^-c+10.^-tau)).*((10.^-c+10.^-tau).*(exp(max_lag.*(10.^-c+10.^-tau))-cos(2*pi*f*max_lag))+2*pi*f.*sin(2*pi*f*max_lag))./...
           (4*pi^2*f.^2+(10.^-c+10.^-tau).^2);
    case 'full'
        fixb = @(b)max(b,realmin);% b>0
        cif = @(lags_list,tau,b,c,f,s,r)(1-b).*exp(-lags_list.*10.^-tau).*(r.*exp(-lags_list.*10.^-c).*...
           ((2+2*sqrt(1-s)-s).*cos(2*pi*f.*lags_list)+4*s.*cos(pi*f.*lags_list)+2-2*sqrt(1-s)-3*s)/4 ...
           +1)+b;
       cif = @(lags_list,tau,b,c,f,s,r)cif(lags_list,tau,fixb(b),c,f,s,r);
       cif_int = @skipping_integral;
    otherwise
end

end

% Integral of the full distribution with skipping
function D = skipping_integral(max_lag,tau,b,c,f,s,r)
    A = 10.^-c+10.^-tau;
    B = exp(max_lag*A);
    b = max(b,realmin);
    
    D=b*max_lag+...
        ...
        (1-b).*...
        (10.^tau).*...
        (1-exp(-max_lag*10.^-tau))+...
        (r.*(1-b)/4).*...
    ((2+2*sqrt(1-s)-s).*((B.^-1).*(A.*(B-cos(2*pi*max_lag*f))+2*pi*f.*sin(2*pi*max_lag*f))./...
    (4*pi^2*f.^2+A.^2))+...
    4*s.*((B.^-1).*(A.*(B-cos(pi*max_lag*f))+pi*f.*sin(pi*max_lag*f))./...
    (pi^2*f.^2+A.^2))+...
    (2-2*sqrt(1-s)-3*s).*(10.^(c+tau)).*(1-(B.^-1))./(10.^c+10.^tau));
end

