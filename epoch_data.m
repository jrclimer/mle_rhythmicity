function [ x, T ] = epoch_data( spk_ts, epochs, varargin )
%EPOCH_DATA Gets the epoched cell array of lags
%
% INPUT
%   spk_ts: The spike times
%   epochs: A mX2 matrix of the epochs, where the first column
%       is the start times and the second is the end.
% 
% PARAMETERS
%   T (0.6): Examinination window
%   drop_tails (false): Whether to drop the data at the ends of the
%       epochs. If true, will increase the falloff.
%
% Copyright (c) 2014, 2015 Trustees of Boston University
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
ip = inputParser;
ip.addParamValue('T',0.6);
ip.addParamValue('drop_tails',false);
ip.parse(varargin{:});
for j = fields(ip.Results)'
    eval([j{1} ' = ip.Results.' j{1} ';']);
end

if drop_tails
   dt = T; 
else
   dt = 0; 
end

x = {};
for i=1:size(epochs,1)
    x1 = spk_ts(spk_ts>=epochs(i,1)&spk_ts<epochs(i,2)-dt);
    x2 = spk_ts(spk_ts>=epochs(i,1)&spk_ts<epochs(i,2));
    x = [x;arrayfun(@(y)x2(x2>y&x2<=y+T)-y,x1,'UniformOutput',false)];
end

T = sum(diff(epochs,[],2));

end

