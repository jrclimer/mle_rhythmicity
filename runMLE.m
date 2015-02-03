function stats_covar = runMLE(root,varargin)

% 2014-05-06 Jason Climer - Handled NaN velocity values by interpolating velocity

%     addpath mle_rhythmicity_fast/
%     startup

if ~exist('runSpeed','var')
    runSpeed = 150;
end

startupJason

epoch0 = root.epoch;
root.epoch = [0 inf];
vid_ts = CMBHOME.Utils.ContinuizeEpochs(root.ts);
subv = min(vid_ts);
vid_ts = vid_ts - subv;

vel = CMBHOME.Utils.ContinuizeEpochs(root.vel) * root.spatial_scale;
temp = 1:numel(vel);
vel = interp1(temp(~isnan(vel)),vel(~isnan(vel)),temp);

root.epoch = epoch0;
root = CMBHOME.Utils.RunningEpochs(root,0,runSpeed,.5);

spk_ts = CMBHOME.Utils.ContinuizeEpochs(root.cel_ts);
spk_ts = spk_ts - subv;

dur = sum(root.epoch(:,2) - root.epoch(:,1));


figure
stats = ...
    mle_rhythmicity_fast(spk_ts,dur,'plotit',1,varargin{:});

figure
stats_covar = ...
    mle_rhythmicity_covar(spk_ts, dur,vid_ts,vel,'plotit',1,'stats0',stats,varargin{:});

end
