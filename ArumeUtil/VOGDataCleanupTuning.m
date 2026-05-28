
function VOGDataCleanupTuning(data)
app = InteractiveUI('VOG data cleanup',@(app) (VOGDataCleanupTuningUpdate(app)), .2);


app.Data.data = data;

app.AddSlider('HVelMax',   1000,  [0 100000])
app.AddSlider('VVelMax',   1000,  [0 100000])
app.AddSlider('TVelMax',   200,  [0 100000])

app.AddSlider('AccelMax',   50000,  [0 1000000])
app.AddSlider('TAccelMax',   50000,  [0 1000000])


app.AddSlider('HPosMaxRange',   90,  [0 10000])
app.AddSlider('VPosMaxRange',   90,  [0 10000])
app.AddSlider('TPosMaxRange',   90,  [0 10000])

app.AddSlider('BadDataPadding',   200,  [0 10000])
app.Open();
end

function [f, h] = VOGDataCleanupTuningUpdateInitPlots(app)

f = figure;
h.H = subplot(2,1,1,'nextplot','add');
h.LeftX = plot(1,1);
h.RightX = plot(1,1);
set(gca,'ylim',[-50 50])
h.V = subplot(2,1,2,'nextplot','add');
h.LeftY = plot(1,1);
h.RightY = plot(1,1);
set(gca,'ylim',[-50 50])

end

function VOGDataCleanupTuningUpdate(app)

if ( ~isfield(app.Data,'f') || isempty(app.Data.f) )
    [f,h] = VOGDataCleanupTuningUpdateInitPlots(app);
    app.Data.f = f;
    app.Data.h = h;
end

data = app.Data.data;
h = app.Data.h;


totalNumberOfFrames = data.FrameNumber(end) - data.FrameNumber(1)+1;
totalTime           = data.Time(end) - data.Time(1);
rawSampleRate       = (totalNumberOfFrames-1)/totalTime;

% find what signals are present in the data
[eyes, eyeSignals, headSignals] = VOGAnalysis.GetEyesAndSignals(data);

for i=1:length(eyes)

    % collect signals
    dt = diff(data.Time);
    x = data.([eyes{i} 'X']);
    y = data.([eyes{i} 'Y']);
    t = data.([eyes{i} 'T']);
    vx = [0;diff(x)./dt];
    vy = [0;diff(y)./dt];
    vt = [0;diff(t)./dt];
    accx = [0;diff(vx)./dt];
    accy = [0;diff(vy)./dt];
    acct = [0;diff(vt)./dt];
    acc = sqrt(accx.^2+accy.^2);

    % sometimes we may also want to save the gaze position
    % in raw pixel coordinates
    if any(strcmp([eyes{i} 'PixX'],data.Properties.VariableNames))
        xpx = data.([eyes{i} 'PixX']);
        ypx = data.([eyes{i} 'PixY']);
        % tpx = cleanedData.([eyes{i} 'PixT']);
        vxpx = [0;diff(xpx)./dt];
        vypx = [0;diff(ypx)./dt];
        % vtpx = [0;diff(tpx)./dt];
        accxpx = [0;diff(vxpx)./dt];
        accypx = [0;diff(vypx)./dt];
        % acctpx = [0;diff(vtpx)./dt];
        accpx = sqrt(accxpx.^2+accypx.^2);
    end




    % find blinks and other abnormal pupil sizes or eye movements

    positionOutOfRange = abs(x) > app.Values.HPosMaxRange ...	% Horizontal eye position out of range
        | abs(y) > app.Values.VPosMaxRange;         	% Vertical eye position out of range

    velocityOutOfRange = abs(vx) > app.Values.HVelMax ...	% Horizontal eye velocity out of range
        | abs(vy) > app.Values.VVelMax;             % Vertical eye velocity out of range

    accelerationOutOfRange = acc>app.Values.AccelMax;

    badData = positionOutOfRange | velocityOutOfRange | accelerationOutOfRange;

    badData = imclose(badData,ones(10));
    badData  = boxcar( badData , round(app.Values.BadDataPadding/1000*rawSampleRate))>0;

    x(badData) = nan;
    y(badData) = nan;
    t(badData) = nan;
    set(h.([eyes{i} 'X']), 'xdata', data.Time,  'ydata', x )
    set(h.([eyes{i} 'Y']), 'xdata', data.Time,  'ydata', y )
end


end
% 
% opt =  VOGAnalysis.GetParameters
% 
%        smoothRloessSpan: 5
%                           BadDataPadding: 200
%                              pupilSizeTh: 10
%                        pupilSizeChangeTh: 10000
%                             HPosMaxRange: 1000
%                             VPosMaxRange: 1000
%                             TPosMaxRange: 20
%                                  HVelMax: 1000
%                                  VVelMax: 1000
%                                  TVelMax: 200
%                                 AccelMax: 50000
%                                TAccelMax: 50000
%                      DETECT_FLAT_PERIODS: 0
%                          Remove_Bad_Data: 1
%           Interpolate_Spikes_of_Bad_Data: 1
%     Interpolate_Pupil_Spikes_of_Bad_Data: 1
%                                    windw: 0.2000
% 
% 
% 
% 
% 
%         function [cleanedData] = CleanData(calibratedData, params)
%             % CLEAN DATA Cleans all the data that may be blinks or bad tracking
%             %
%             %   [cleanedData] = CleanData(calibratedData, params)
%             %
%             %   Inputs:
%             %       - calibratedData: calibrated data
%             %       - params: parameters for the processing.
%             %
%             %   Outputs:
%             %       - cleanedData: cleaned data
% 
%             try
%                 tic
% 
% 
%                 % find what signals are present in the data
%                 [eyes, eyeSignals, headSignals] = VOGAnalysis.GetEyesAndSignals(calibratedData);
% 
% 
%                 % ---------------------------------------------------------
%                 % Interpolate missing frames
%                 %----------------------------------------------------------
%                 % Find missing frames and intorpolate them
%                 % It is possible that some frames were dropped during the
%                 % recording. We will interpolate them. But only if they are
%                 % just a few in a row. If there are many we will fill with
%                 % NaNs. The fram numbers and the timestamps will be
%                 % interpolated regardless. From now on frame numbers and
%                 % timestamps cannot be NaN and they must follow a continued
%                 % growing interval
% 
%                 cleanedData = table;    % cleaned data
% 
%                 % calcualte the samplerate
%                 totalNumberOfFrames = calibratedData.FrameNumber(end) - calibratedData.FrameNumber(1)+1;
%                 totalTime           = calibratedData.Time(end) - calibratedData.Time(1);
%                 rawSampleRate       = (totalNumberOfFrames-1)/totalTime;
% 
%                 % find dropped and not dropped frames
%                 notDroppedFrames = calibratedData.FrameNumber - calibratedData.FrameNumber(1) + 1;
%                 droppedFrames = ones(max(notDroppedFrames),1);
%                 droppedFrames(notDroppedFrames) = 0;
%                 interpolableFrames = droppedFrames-imopen(droppedFrames,ones(3)); % 1 or 2 frames in a row, not more
% 
%                 % TODO: deal with concatenated files that may have large
%                 % gaps on timestamps. Maybe interpolate only up to some
%                 % certain duration. Otherwise leave the gap. But then maybe
%                 % all cleanup needs to be done in chunks to not contaminate
%                 % discontinuos recordings. Or maybe they just need to be
%                 % concatenated after cleanup... not sure!
% 
%                 % create the new continuos FrameNumber and Time variables
%                 % but also save the original raw frame numbers and time
%                 % stamps with NaNs in the dropped frames.
%                 cleanedData.FrameNumber                 = (1:max(notDroppedFrames))';
%                 cleanedData.Time                        = (cleanedData.FrameNumber-1)/rawSampleRate;
%                 cleanedData.RawFrameNumber              = nan(height(cleanedData), 1);
%                 cleanedData.LeftCameraRawFrameNumber  	= nan(height(cleanedData), 1);
%                 cleanedData.RightCameraRawFrameNumber 	= nan(height(cleanedData), 1);
%                 cleanedData.RawTime                     = nan(height(cleanedData), 1);
%                 cleanedData.RawFrameNumber(notDroppedFrames)            = calibratedData.FrameNumber;
%                 cleanedData.LeftCameraRawFrameNumber(notDroppedFrames)  = calibratedData.LeftFrameNumberRaw;
%                 cleanedData.RightCameraRawFrameNumber(notDroppedFrames) = calibratedData.RightFrameNumberRaw;
%                 cleanedData.RawTime(notDroppedFrames)                   = calibratedData.Time;
%                 cleanedData.DroppedFrame                                = droppedFrames;
% 
%                 % interpolate signals
%                 signalsToInterpolate = {};
%                 rawSignalsToInterpolate = {};
%                 rawIntSignalsToInterpolate = {};
%                 for i=1:length(eyes)
%                     for j=1:length(eyeSignals)
%                         signalsToInterpolate{end+1}         = [eyes{i} eyeSignals{j}];
%                         rawSignalsToInterpolate{end+1}      = [eyes{i} 'Raw' eyeSignals{j}];
%                         rawIntSignalsToInterpolate{end+1}   = [eyes{i} 'RawInt' eyeSignals{j}];
%                     end
%                 end
%                 for j=1:length(headSignals)
%                     signalsToInterpolate{end+1} = ['Head' headSignals{j}];
%                     rawSignalsToInterpolate{end+1}      = [eyes{i} 'Raw' 'Head' headSignals{j}];
%                     rawIntSignalsToInterpolate{end+1}   = [eyes{i} 'RawInt' 'Head' headSignals{j}];
%                 end
% 
%                 for i=1:length(signalsToInterpolate)
%                     signalName = signalsToInterpolate{i};
%                 for i=1:length(eyes)
% 
%                     rawSignalName = rawSignalsToInterpolate{i};
%                     rawIntSignalName = rawIntSignalsToInterpolate{i}; 
% 
%                     cleanedData.(signalName)       = nan(height(cleanedData), 1); % signal that will be cleaned
%                     cleanedData.(rawSignalName)    = nan(height(cleanedData), 1); % raw signal with nans in dropped frames
%                     cleanedData.(rawIntSignalName) = nan(height(cleanedData), 1); % almost raw signal with some interpolated dropped frames
% 
%                     cleanedData.(rawSignalName)(notDroppedFrames)  = calibratedData.(signalName);
% 
%                     % interpolate missing frames but only if they are
%                     % 2 or less in a row. Otherwise put nans in there.
%                     datInterp = interp1(notDroppedFrames, cleanedData.(rawSignalName)(notDroppedFrames),  cleanedData.FrameNumber );
%                     datInterp(droppedFrames & ~interpolableFrames) = nan;
%                     cleanedData.(rawIntSignalName) = datInterp;
% 
%                     cleanedData.(signalName) = datInterp;
%                 end
% 
%                 % ---------------------------------------------------------
%                 % End interpolate missing samples
%                 %----------------------------------------------------------
% 
% 
%                 % ---------------------------------------------------------
%                 % Find bad samples
%                 %----------------------------------------------------------
%                 % We will use multiple heuristics to determine portions of
%                 % data that may not be good. Then we will interpolate short
%                 % spikes of bad data while removing everything else bad
%                 % plus some padding around
%                 % Find bad samples
%                 for i=1:length(eyes)
% 
%                     badData = isnan(cleanedData.([eyes{i} 'X'])) | isnan(cleanedData.([eyes{i} 'Y']));
%                     pupilSizeChangeOutOfrange = nan(size(badData));
% 
%                     % Calculate a smooth version of the pupil size to detect changes in
%                     % pupil size that are not normal. Thus, must be blinks or errors in
%                     % tracking. Downsample the signal to speed up the smoothing.
%                     if ( ismember('Pupil', eyeSignals) && length( cleanedData.([eyes{i} 'Pupil'])) > 200)
%                         pupil = cleanedData.([eyes{i} 'Pupil']);
%                         pupilDecimated = pupil(1:25:end); %decimate the pupil signal
%                         % if ( exist('smooth','file') )
%                             % pupilSmooth = smooth(pupilDecimated,params.CleanUp.smoothRloessSpan*rawSampleRate/25/length(pupilDecimated),'rloess');
%                             % if all(pupilSmooth == 0)
%                                 pupilSmooth = nanmedfilt(pupilDecimated,round(params.CleanUp.smoothRloessSpan*rawSampleRate/25));
%                             % end
%                         % else
%                             pupilSmooth = nanmedfilt(pupilDecimated,round(params.CleanUp.smoothRloessSpan*rawSampleRate/25));
%                         % end
%                         pupilSmooth = interp1((1:25:length(pupil))',pupilSmooth,(1:length(pupil))');
% 
% %                         cleanedData.([eyes{i} 'Pupil']) = pupilSmooth;
% 
%                         % find blinks and other abnormal pupil sizes or eye movements
%                         pth = std(pupilSmooth,'omitnan')*params.CleanUp.pupilSizeTh; %pth = mean(pupilSmooth,'omitnan')*params.CleanUp.pupilSizeTh/100;
%                         pupilSizeChangeOutOfrange = abs(pupilSmooth-pupil) > pth ...                 % pupil size far from smooth pupil size
%                             | abs([0;diff(pupil)*rawSampleRate]) > params.CleanUp.pupilSizeChangeTh;        % pupil size changes too suddenly from sample to sample
% 
%                         % find spikes (Single bad samples surrounded by at least 1 good sample to each side) to interpolate
%                         if ( params.CleanUp.Interpolate_Pupil_Spikes_of_Bad_Data )
%                             pupilSpikes  = pupilSizeChangeOutOfrange & ( boxcar(pupilSizeChangeOutOfrange, 3)*3 == 1 );
%                             cleanedData.([eyes{i} 'PupilSpikes']) = pupilSpikes;
%                             for j=1:length(eyeSignals)
%                                 cleanedData.([eyes{i} eyeSignals{j}])(pupilSpikes) = interp1(find(~pupilSpikes),cleanedData.([eyes{i}  eyeSignals{j}])(~pupilSpikes),  find(pupilSpikes));
%                             end
% 
%                             badData = badData | (pupilSizeChangeOutOfrange & ~pupilSpikes);
%                         else
%                             badData = badData | pupilSizeChangeOutOfrange;
%                         end
%                     end
% 
%                     % collect signals
%                     dt = diff(cleanedData.Time);
%                     x = cleanedData.([eyes{i} 'X']);
%                     y = cleanedData.([eyes{i} 'Y']);
%                     t = cleanedData.([eyes{i} 'T']);
%                     vx = [0;diff(x)./dt];
%                     vy = [0;diff(y)./dt];
%                     vt = [0;diff(t)./dt];
%                     accx = [0;diff(vx)./dt];
%                     accy = [0;diff(vy)./dt];
%                     acct = [0;diff(vt)./dt];
%                     acc = sqrt(accx.^2+accy.^2);
% 
%                     % sometimes we may also want to save the gaze position
%                     % in raw pixel coordinates
%                     if any(strcmp([eyes{i} 'PixX'],cleanedData.Properties.VariableNames))
%                         xpx = cleanedData.([eyes{i} 'PixX']);
%                         ypx = cleanedData.([eyes{i} 'PixY']);
%                         % tpx = cleanedData.([eyes{i} 'PixT']);
%                         vxpx = [0;diff(xpx)./dt];
%                         vypx = [0;diff(ypx)./dt];
%                         % vtpx = [0;diff(tpx)./dt];
%                         accxpx = [0;diff(vxpx)./dt];
%                         accypx = [0;diff(vypx)./dt];
%                         % acctpx = [0;diff(vtpx)./dt];
%                         accpx = sqrt(accxpx.^2+accypx.^2);
%                     end
% 
%                     % find blinks and other abnormal pupil sizes or eye movements
% 
%                     positionOutOfRange = abs(x) > params.CleanUp.HPosMaxRange ...	% Horizontal eye position out of range
%                         | abs(y) > params.CleanUp.VPosMaxRange;         	% Vertical eye position out of range
% 
%                     velocityOutOfRange = abs(vx) > params.CleanUp.HVelMax ...	% Horizontal eye velocity out of range
%                         | abs(vy) > params.CleanUp.VVelMax;             % Vertical eye velocity out of range
% 
%                     accelerationOutOfRange = acc>params.CleanUp.AccelMax;
% 
%                     badData = badData | positionOutOfRange | velocityOutOfRange | accelerationOutOfRange;
% 
%                     badFlatPeriods = nan(size(badData));
%                     if ( params.CleanUp.DETECT_FLAT_PERIODS )
%                         % if three consecutive samples are the same value this main they
%                         % are interpolated
%                         badFlatPeriods =  boxcar([nan;abs(diff(x))],2) == 0 ...
%                             | boxcar([nan;abs(diff(y))],2) == 0;
%                         badData = badData | badFlatPeriods;
%                     end
% 
% %                      badDataSpikes = VOGAnalysis.FindSpikyNoisePeriods(x,params.CleanUp.windw,rawSampleRate) ...
% %                          | VOGAnalysis.FindSpikyNoisePeriods(y,params.CleanUp.windw,rawSampleRate);
% 
% %                      badData = badData | badDataSpikes;
% 
%                     % spikes of good data in between bad data are probably bad
%                     badData = imclose(badData,ones(10));
%                     badDataT = badData | abs(t) > params.CleanUp.TPosMaxRange | abs(vt) > params.CleanUp.TVelMax | abs(acct) > params.CleanUp.TAccelMax;
%                     badDataT = imclose(badDataT,ones(10));
% 
%                     % but spikes of bad data in between good data can be
%                     % interpolated
%                     % find spikes of bad data. Single bad samples surrounded by at least 2
%                     % good samples to each side
%                     spikes  = badData & ( boxcar(~badData, 3)*3 >= 2 );
%                     spikest = badDataT & ( boxcar(~badDataT, 3)*3 >= 2 );
% 
%                     % TODO: maybe better than blink span find the first N samples
%                     % around the blink that are within a more stringent criteria
%                     if ( params.CleanUp.BadDataPadding > 0 )
%                         if ( params.CleanUp.Interpolate_Spikes_of_Bad_Data)
%                             badData  = boxcar( badData  & ~spikes, round(params.CleanUp.BadDataPadding/1000*rawSampleRate))>0;
%                             badDataT = boxcar( badDataT & ~spikest, round(params.CleanUp.BadDataPadding/1000*rawSampleRate))>0;
%                         else
%                             badData  = boxcar( badData, round(params.CleanUp.BadDataPadding/1000*rawSampleRate))>0;
%                             badDataT = boxcar( badDataT, round(params.CleanUp.BadDataPadding/1000*rawSampleRate))>0;
%                         end
%                     end
% 
%                     cleanedData.([eyes{i} 'Spikes']) = spikes;
%                     cleanedData.([eyes{i} 'BadData']) = badData;
% %                     cleanedData.([eyes{i} 'BadDataSpikes']) = badDataSpikes;
%                     cleanedData.([eyes{i} 'SpikesT']) = spikest;
%                     cleanedData.([eyes{i} 'BadDataT']) = badDataT;
% 
%                     cleanedData.([eyes{i} 'BadPupil'])          = pupilSizeChangeOutOfrange;
%                     cleanedData.([eyes{i} 'BadPosition'])       = positionOutOfRange;
%                     cleanedData.([eyes{i} 'BadVelocity'])       = velocityOutOfRange;
%                     cleanedData.([eyes{i} 'BadAcceleration'])   = accelerationOutOfRange;
%                     cleanedData.([eyes{i} 'BadFlatPeriods'])    = badFlatPeriods;
% 
%                     % Clean up data
%                     for j=1:length(eyeSignals)
%                         if ( ~strcmp(eyeSignals{j},'T') )
%                             if ( params.CleanUp.Remove_Bad_Data )
%                                 badData = cleanedData.([eyes{i} 'BadData']);
%                                 % put nan on bad samples of data (blinks)
%                                 cleanedData.([eyes{i} eyeSignals{j}])(badData) = nan;
%                             end
%                             if ( params.CleanUp.Interpolate_Spikes_of_Bad_Data )
%                                 spikes = cleanedData.([eyes{i} 'Spikes']);
%                                 % interpolate single spikes of bad data
%                                 cleanedData.([eyes{i} eyeSignals{j}])(spikes)  = interp1(find(~spikes),cleanedData.([eyes{i} eyeSignals{j}])(~spikes),  find(spikes));
%                             end
%                         else
%                             if ( params.CleanUp.Remove_Bad_Data )
%                                 badDataT = cleanedData.([eyes{i} 'BadDataT']);
%                                 % put nan on bad samples of data (blinks)
%                                 cleanedData.([eyes{i} eyeSignals{j}])(badDataT) = nan;
%                             end
%                             if ( params.CleanUp.Interpolate_Spikes_of_Bad_Data )
%                                 spikest = cleanedData.([eyes{i} 'SpikesT']);
%                                 % interpolate single spikes of bad data
%                                 cleanedData.([eyes{i} eyeSignals{j}])(spikest)  = interp1(find(~spikest),cleanedData.([eyes{i} eyeSignals{j}])(~spikest),  find(spikest));
%                             end
%                         end
%                     end
% 
%                 end
% 
%                 timeCleaning = toc;
% 
%                 cprintf('blue', sprintf('++ VOGAnalysis :: Data has %d dropped frames, %d were interpolated.\n', ...
%                     sum(cleanedData.DroppedFrame), sum(interpolableFrames)) );
% 
% 
%                 Lbad = nan;
%                 Rbad = nan;
%                 LbadT = nan;
%                 RbadT = nan;
% 
%                 if ( any(contains(eyes,'Left') ))
%                     Lbad = round(mean(~cleanedData.LeftBadData)*100);
%                     LbadT = round(mean(~cleanedData.LeftBadDataT)*100);
%                 end
%                 if ( any(contains(eyes,'Right') ))
%                     Rbad = round(mean(~cleanedData.RightBadData)*100);
%                     RbadT = round(mean(~cleanedData.RightBadDataT)*100);
%                 end
%                 cprintf('blue', sprintf('++ VOGAnalysis :: Data cleaned in %0.1f s: LXY %d%%%% RXY %d%%%% LT %d%%%% RT %d%%%% is good data.\n', ...
%                     timeCleaning, Lbad, Rbad, LbadT, RbadT ));
% 
% 
%                 cleanedData.Properties.UserData.Eyes = eyes;
%                 cleanedData.Properties.UserData.Signals = eyeSignals;
%                 cleanedData.Properties.UserData.EyeSignals = intersect(eyeSignals, {'X', 'Y','T'});
%                 cleanedData.Properties.UserData.LEFT = any(strcmp(eyes,'Left'));
%                 cleanedData.Properties.UserData.RIGHT = any(strcmp(eyes,'Right'));
% 
%                 cleanedData.Properties.UserData.params = params;
% 
%             catch ex
%                 getReport(ex)
%             end
%         end