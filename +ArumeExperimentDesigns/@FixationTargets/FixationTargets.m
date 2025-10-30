classdef FixationTargets < ArumeExperimentDesigns.EyeTracking
    %OPTOKINETICTORSION Summary of this class goes here
    %   Detailed explanation goes here

    properties
        fixRad = 20;
        fixColor = [255 0 0];
        targetPositions =[];
    end

    % ---------------------------------------------------------------------
    % Experiment design methods
    % ---------------------------------------------------------------------
    methods ( Access = protected )
        function dlg = GetOptionsDialog( this, importing )
            if( ~exist( 'importing', 'var' ) )
                importing = 0;
            end
            dlg = GetOptionsDialog@ArumeExperimentDesigns.EyeTracking(this, importing);
            dlg.Debug.DisplayVariableSelection = 'TrialNumber TrialResult TargetPosition'; % which variables to display every trial in the command line separated by spaces

            dlg.DisplayOptions.ScreenWidth = { 55 '* (cm)' [1 3000] };
            dlg.DisplayOptions.ScreenHeight = { 31 '* (cm)' [1 3000] };
            dlg.DisplayOptions.ScreenDistance = { 67 '* (cm)' [1 3000] };

            dlg.TrialDuration =  { 10 '* (s)' [1 100] };
            dlg.NumberRepetitions = 10;            
            dlg.DisplayOptions.SelectedScreen = { 1 '* (screen)' [0 5] };


            dlg.TargetSize = 1;
            dlg.Calibration_Type = { {'Center dot' '5 dots' '{9 dots}' '13 dots' '17 dots'} };
            dlg.Calibration_Distance_H = { 10 '* (deg)' [1 3000] };
            dlg.Calibration_Distance_V = { 10 '* (deg)' [1 3000] };

            dlg.BackgroundBrightness = 255/2;
        end


        function trialTable = SetUpTrialTable(this)
            Screen('Preference', 'SkipSyncTests', 0);                  % DO NOT skip sync tests (set to 1 or 2 for debugging only)
            Screen('Preference', 'VisualDebugLevel', 1);               % Minimize startup splash
            Screen('Preference', 'SuppressAllWarnings', 0);            % Show all warnings
            Screen('Preference', 'ConserveVRAM', 4096);

            h = this.ExperimentOptions.Calibration_Distance_H;
            v = this.ExperimentOptions.Calibration_Distance_V;
            temp = 1.5;

            switch(this.ExperimentOptions.Calibration_Type)
                case 'Center dot'
                    this.targetPositions = {[0,0]};
                case '5 dots'
                    this.targetPositions = {[0,0],[h,v],[h,-v],[-h,v],[-h,-v]};
                case '9 dots'
                    this.targetPositions = {[0,0],[h,0],[-h,0],[0,v],[0,-v],[h,v],[h,-v],[-h,v],[-h,-v]};
                case '13 dots'
                    this.targetPositions = {[0,0],[h,0],[-h,0],[0,v],[0,-v],[h,v],[h,-v],[-h,v],[-h,-v],...
                        [h/temp,v/temp],[h/temp,-v/temp],[-h/temp,v/temp],[-h/temp,-v/temp]};
                case '17 dots'
                    this.targetPositions = {[0,0],[h,0],[-h,0],[0,v],[0,-v],[h,v],[h,-v],[-h,v],[-h,-v],...
                        [h/temp,0],[-h/temp,0],[0,v/temp],[0,-v/temp],[h/temp,v/temp],[h/temp,-v/temp],[-h/temp,v/temp],[-h/temp,-v/temp]};
            end
            t = ArumeCore.TrialTableBuilder();

            t.AddConditionVariable("TargetPosition", { ...
                [0 0], [0 10], [10 0], [10 10], [-10 10], [-10 -10], [10 -10] [-10 0], [0 -10] ...
                [0 2], [0 4], [0 6], [0 8], [0 -2], [0 -4], [0 -6], [0 -8], ...
                [2 0], [4 0], [6 0], [8 0], [-2 0], [-4 0], [-6 0], [-8 0], ...
                [2 2], [4 4], [6 6], [8 8], [-2 2], [-4 4], [-6 6], [-8 8], ...
                [2 -2], [4 -4], [6 -6], [8 -8], [-2 -2], [-4 -4], [-6 -6], [-8 -8] ...
                });

            % Add all conditions to a single block with N repetitions
            nReps = this.ExperimentOptions.NumberRepetitions;
            t.AddBlock(1:height(t.ConditionTable), 1);

            trialSequence = 'Random';
            blockSequence = 'Sequential';
            blockSequenceRepeatitions = nReps;
            abortAction = 'Repeat';
            trialsPerSession = 1000;  % You can use a large number if not splitting

            trialTable = t.GenerateTrialTable(trialSequence, blockSequence, blockSequenceRepeatitions, abortAction, trialsPerSession);
        end


        function [trialResult, thisTrialData] = runTrial(this, thisTrialData)
            try
                Enum = ArumeCore.ExperimentDesign.getEnum();
                graph = this.Graph;
                trialResult = Enum.trialResult.CORRECT;

                lastFlipTime        = GetSecs;
                secondsRemaining    = this.ExperimentOptions.TrialDuration;
                thisTrialData.TimeStartLoop = lastFlipTime;

                if (~isempty(this.eyeTracker))
                    thisTrialData.EyeTrackerFrameStartLoop = this.eyeTracker.RecordEvent( ...
                        sprintf('TRIAL_START_LOOP %d %d [%d %d]', ...
                        thisTrialData.TrialNumber, ...
                        thisTrialData.Condition) );
                end

                while secondsRemaining > 0
                    secondsElapsed      = GetSecs - thisTrialData.TimeStartLoop;
                    secondsRemaining    = this.ExperimentOptions.TrialDuration - secondsElapsed;

                    if secondsRemaining > 0
                        %-- Draw fixation spot as a white cross (+)
                        Screen('FillRect', graph.window, this.ExperimentOptions.BackgroundBrightness);

                        [mx, my] = RectCenter(this.Graph.wRect);

                        % Get the stimulus position in degrees from the trial table
                        targetDeg = thisTrialData.TargetPosition{1};
                        dx = targetDeg(1);
                        dy = targetDeg(2);

                        % === LOG POSITION INTO TRIAL DATA ===
                        thisTrialData.StimulusPosition_X = targetDeg(1);
                        thisTrialData.StimulusPosition_Y = targetDeg(2);

                        % Convert visual degrees to pixels
                        pixelsPerDegree = this.Graph.pxWidth / this.ExperimentOptions.DisplayOptions.ScreenWidth * ...
                            this.ExperimentOptions.DisplayOptions.ScreenDistance;

                        targetHPix = pixelsPerDegree * tand(dx);
                        targetYPix = pixelsPerDegree * tand(dy);

                        % Fixation location
                        fixX = mx + targetHPix / 2;
                        fixY = my + targetYPix / 2;

                        % Fixation cross parameters
                        targetSizeDeg = this.ExperimentOptions.TargetSize;
                        crossLength = pixelsPerDegree * tand(targetSizeDeg); % in pixels
                        crossThickness = 2;
                        crossColor = [0, 0, 0];

                        % Define cross lines centered on fixX, fixY
                        crossCoords = [ ...
                            -crossLength/2, 0; ...
                            crossLength/2, 0; ...
                            0, -crossLength/2; ...
                            0,  crossLength/2 ...
                            ]';

                        Screen('DrawLines', graph.window, crossCoords, crossThickness, crossColor, [fixX, fixY], 2);

                        Screen('DrawingFinished', graph.window);
                    else
                        Screen('FillRect', graph.window, this.ExperimentOptions.BackgroundBrightness);
                        Screen('DrawingFinished', graph.window);
                    end

                    % Flip screen buffer
                    this.Graph.Flip();
                end

            catch ex
                rethrow(ex)
            end
        end

    end

    methods ( Access = public )
    end

    % ---------------------------------------------------------------------
    % Plot methods
    % ---------------------------------------------------------------------
    methods ( Access = public )
        function [out, options] = PlotAggregate_MicroCorrelations(this, sessions, options)


            out = [];
            if ( nargin == 1 )
                options = this.PlotAggregate_MicroCorrelations('get_defaults');
            end

            if ( ischar(sessions) )
                command = sessions;
                switch( command)
                    case 'get_options'
                        options =[];
                        return
                    case 'get_defaults'
                        optionsDlg = VOGAnalysis.PlotAggregate_MicroCorrelations('get_options');
                        options = StructDlg(optionsDlg,'',[],[],'off');
                        return
                end
            end


            data = [];
            % for each subj, ignoring 000
            for subj = 1:length(sessions)
                %load data
                AnalysisResults_QuickPhases = sessions(subj).analysisResults.QuickPhases;

                % pick the center positions
                center_indx = AnalysisResults_QuickPhases.Left_X_MeanPosition <= 2 &...
                    AnalysisResults_QuickPhases.Left_X_MeanPosition >= -2 &...
                    AnalysisResults_QuickPhases.Left_Y_MeanPosition <= 2 &...
                    AnalysisResults_QuickPhases.Left_Y_MeanPosition >= -2;

                d = AnalysisResults_QuickPhases(center_indx,:);

                % create the table
                data(subj).X_Vergence = d.Left_X_Displacement - d.Right_X_Displacement;
                data(subj).Y_Vergence = d.Left_Y_Displacement - d.Right_Y_Displacement;
                data(subj).T_Vergence = d.Left_T_Displacement - d.Right_T_Displacement;

                data(subj).X_Version = (d.Left_X_Displacement + d.Right_X_Displacement)./2;
                data(subj).Y_Version = (d.Left_Y_Displacement + d.Right_Y_Displacement)./2;
                data(subj).T_Version = (d.Left_T_Displacement + d.Right_T_Displacement)./2;

                data(subj).Left_X = d.Left_X_Displacement;
                data(subj).Left_Y = d.Left_Y_Displacement;
                data(subj).Left_T = d.Left_T_Displacement;

                data(subj).Right_X = d.Right_X_Displacement;
                data(subj).Right_Y = d.Right_Y_Displacement;
                data(subj).Right_T = d.Right_T_Displacement;


            end


            %% create the corr table

            % combine all data
            dataAll.X_Vergence = [];
            dataAll.Y_Vergence = [];
            dataAll.T_Vergence = [];

            dataAll.X_Version = [];
            dataAll.Y_Version = [];
            dataAll.T_Version = [];

            dataAll.Left_X = [];
            dataAll.Left_Y = [];
            dataAll.Left_T = [];

            dataAll.Right_X = [];
            dataAll.Right_Y = [];
            dataAll.Right_T = [];

            for subj = 1:length(data)
                dataAll.X_Vergence = [dataAll.X_Vergence; data(subj).X_Vergence];
                dataAll.Y_Vergence = [dataAll.Y_Vergence; data(subj).Y_Vergence];
                dataAll.T_Vergence = [dataAll.T_Vergence; data(subj).T_Vergence];

                dataAll.Y_Version = [dataAll.Y_Version; data(subj).Y_Version];
                dataAll.X_Version = [dataAll.X_Version; data(subj).X_Version];
                dataAll.T_Version = [dataAll.T_Version; data(subj).T_Version];

                dataAll.Left_X = [dataAll.Left_X; data(subj).Left_X];
                dataAll.Left_Y = [dataAll.Left_Y; data(subj).Left_Y];
                dataAll.Left_T = [dataAll.Left_T; data(subj).Left_T];

                dataAll.Right_X = [dataAll.Right_X; data(subj).Right_X];
                dataAll.Right_Y = [dataAll.Right_Y; data(subj).Right_Y];
                dataAll.Right_T = [dataAll.Right_T; data(subj).Right_T];
            end

            dataAll_VV = [...
                dataAll.X_Version,...
                dataAll.Y_Version,...
                dataAll.X_Vergence,...
                dataAll.Y_Vergence,...
                dataAll.T_Vergence,...
                dataAll.T_Version...
                ];

            dataAll_LR = [...
                dataAll.Left_X,...
                dataAll.Left_Y,...
                dataAll.Left_T,...
                dataAll.Right_X,...
                dataAll.Right_Y,...
                dataAll.Right_T...
                ];

            [R,P] = corr(dataAll_VV,'rows','complete');

            x = [-1:0.5:1];
            c = 0;
            figure,
            for i = 1:6
                for j = 1:6
                    c = c + 1;
                    h(i,j) = subplot(6,6,c);
                    if i <= j, continue; end

                    %indx = (dataAll_(:,j) < 1 & dataAll_(:,j) > -1) & (dataAll_(:,i) < 1 & dataAll_(:,i) > -1);

                    fitlm_{i,j} = fitlm(dataAll_VV(:,j),dataAll_VV(:,i),'RobustOpts','on');
                    %fitlm_{i,j}.Coefficients = robustfit(dataAll_VV(:,j),dataAll_VV(:,i));

                    plot(x,x*fitlm_{i,j}.Coefficients{2,1} + fitlm_{i,j}.Coefficients{1,1},'r--','LineWidth',2),
                    %plot(fitlm_{i,j}),legend off
                    hold on
                    plot(dataAll_VV(:,j),dataAll_VV(:,i),'b.');
                    xlim([-1,1])
                    ylim([-1,1])
                    title(strcat("R = ",num2str(round(R(c),2)),", P = ",num2str(P(c)) ,", Slope = ",num2str( fitlm_{i,j}.Coefficients{2,1} ) ))


                end
            end

            linkaxes(h(:))

            n = [...
                "X Version",...
                "Y Version",...
                "X Vergence",...
                "Y Vergence",...
                "T Vergence",...
                "T Version"...
                ];

            c=1;
            for i = 1:6:36
                subplot(6,6,i)
                ylabel(n(c))
                c=c+1;
            end


            for i = 31:36
                subplot(6,6,i)
                xlabel(n(i-30))
            end

            sgtitle("Center Fixational Micro Saccades Positional Displacement (degree)")


            figure,
            subplot(1,2,1)
            heatmap(R)
            title("R")

            subplot(1,2,2)
            heatmap(P)
            title("P")
            sgtitle("Center Fixational Micro Saccades Correlation parameters for Vergence and Version")

            figure,
            bar(mean(dataAll_VV,'omitnan'))
            hold on
            errorbar(mean(dataAll_VV,'omitnan'),std(dataAll_VV,'omitnan'),'o')
            xticklabels(n)
            ylabel('Amplitude (degree)')


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%% LEFT RIGHT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % [R,P] = corr(dataAll_LR,'rows','complete');
            %
            %             x = [-1:0.5:1];
            %             c = 0;
            %             figure,
            %             for i = 1:6
            %                 for j = 1:6
            %                     c = c + 1;
            %                     h2(i,j) = subplot(6,6,c);
            %                     if i >= j, continue; end
            %
            %                     %indx = (dataAll_(:,j) < 1 & dataAll_(:,j) > -1) & (dataAll_(:,i) < 1 & dataAll_(:,i) > -1);
            %
            %                     fitlm_{i,j} = fitlm(dataAll_LR(:,j),dataAll_LR(:,i));
            %
            %                     plot(x,x*fitlm_{i,j}.Coefficients{2,1} + fitlm_{i,j}.Coefficients{1,1},'r--','LineWidth',2),
            %                     %plot(fitlm_{i,j}),legend off
            %                     hold on
            %                     plot(dataAll_LR(:,j),dataAll_LR(:,i),'b.');
            %                     xlim([-1,1])
            %                     ylim([-1,1])
            %                     title(strcat("R = ",num2str(round(R(c),2)), ", P = ",num2str(P(c)) ,", Slope = ",num2str( fitlm_{i,j}.Coefficients{2,1} ) ))
            %
            %
            %                 end
            %             end
            %
            %             linkaxes(h2(:))
            %
            %             n = [...
            %                 "Left X",...
            %                 "Left Y",...
            %                 "Left T",...
            %                 "Right X",...
            %                 "Right Y",...
            %                 "Right T"...
            %                 ];
            %
            %             c=1;
            %             for i = 1:6:36
            %                 subplot(6,6,i)
            %                 ylabel(n(c))
            %                 c=c+1;
            %             end
            %
            %
            %             for i = 31:36
            %                 subplot(6,6,i)
            %                 xlabel(n(i-30))
            %             end
            %
            %             sgtitle("Center Fixational Micro Saccades Positional Displacement (degree)")
            %
            %
            %             figure,
            %             subplot(1,2,1)
            %             heatmap(R)
            %             title("R")
            %
            %             subplot(1,2,2)
            %             heatmap(P)
            %             title("P")
            %             sgtitle("Center Fixational Micro Saccades Correlation parameters for Left and Right")
            %




        end

        function Plot_Listings(this)
            figure
            s = this.Session.samplesDataTable;
            plot(s.LeftX, s.LeftT,'o')
            hold
            plot(s.RightX, s.RightT,'o')
            set(gca,'xlim',[-20 20],'ylim',[-20 20])
            xlabel('Horizontal')
            ylabel('Torsion')
        end

        function [out, options] = Plot_XY(this)

            s = this.Session.samplesDataTable;

            figure
            plot(s.LeftX, s.LeftY,'o','MarkerSize',3);
            hold
            plot(s.RightX, s.RightY,'o','MarkerSize',3);

            set(gca,'xlim',[-40 40], 'ylim', [-30 30])
        end

        function [out, options] = Plot_PositionVelocity(this)

            p = this.Session.analysisResults.SlowPhases.X_MeanPosition;
            pp = round(p/2.5)*2.5;
            [means,pred,grp,sem] = grpstats(this.Session.analysisResults.SlowPhases.X_MeanVelocity,pp, ["median","predci","gname", "sem"],"Alpha",0.1);
            figure
            subplot(2,2,1)
            errorbar(str2double(grp), means, sem,'o')
            set(gca,'xlim',[-15 15], 'ylim', [-3 3])
            xlabel('Horizontal position (deg)')
            ylabel('Horizontal drift velocity (deg/s)')
            title("fixation data")

            p = this.Session.analysisResults.SlowPhases.Y_MeanPosition;
            pp = round(p/2.5)*2.5;
            [means,pred,grp,sem] = grpstats(this.Session.analysisResults.SlowPhases.Y_MeanVelocity,pp, ["median","predci","gname", "sem"],"Alpha",0.1);

            subplot(2,2,2)
            errorbar(str2double(grp)*-1, means*-1, sem,'o')
            xlabel('Vertical position (deg)')
            ylabel({'Vertical drift velocity (deg/s)' ; '*-(y) to correct for flipped sign in data collection'})
            set(gca,'xlim',[-15 15], 'ylim', [-3 3])
            title("fixation data")


            p = this.Session.analysisResults.SlowPhases.Left_X_MeanPosition - this.Session.analysisResults.SlowPhases.Right_X_MeanPosition;
            v = this.Session.analysisResults.SlowPhases.Left_X_MeanVelocity - this.Session.analysisResults.SlowPhases.Right_X_MeanVelocity;
            pp = round(p*2.5)/2.5;
            [means,pred,grp,sem] = grpstats(v,pp, ["median","predci","gname", "sem"],"Alpha",0.1);

            subplot(2,2,3)
            errorbar(str2double(grp), means, sem,'o')
            set(gca,'xlim',[-5 5], 'ylim', [-3 3])
            xlabel('Vergence position (deg)')
            ylabel('Vergence drift velocity (deg/s)')
            title("fixation data")

            % p = this.Session.analysisResults.SlowPhases.Y_MeanPosition;
            % pp = round(p/2.5)*2.5;
            % [means,pred,grp] = grpstats(this.Session.analysisResults.SlowPhases.Y_MeanVelocity,pp, ["mean","predci","gname"],"Alpha",0.1);
            %
            % subplot(2,2,4)
            % plot(str2double(grp), means,'o')

            % set(gca,'xlim',[-15 15], 'ylim', [-3 3])
        end

        function [out, options] = Plot_DriftDiffusionAnalysis(this)
            slowPhases = this.Session.analysisResults.SlowPhases;
            data = this.Session.samplesDataTable;
            time = data.Time;

            condition = 'fix';
            eyes = ["right", "left"];

            for jk = 1:length(eyes)
                eye = eyes(jk);
                %% Define Colors
                colors = defineColors();

                %% Plot Legend (only if condition == 'fix')
                if strcmp(condition, 'fix')
                    plotLegend(colors);
                end

                %% Initialize variables for analysis
                D_est_all = [];
                timeCoord = [];
                coord_labels = {};
                axis_labels = {};
                avg_taus = [];
                all_meanDispY = [];
                all_meanDispX = [];
                all_dur = [];
                num = 0;

                samplerate = data.Properties.UserData.sampleRate;
                dt = 1 / samplerate;          % Sampling interval

                minsamplelag = round(50/1000*samplerate);
                maxsamplelag = round(250/1000*samplerate);
                lags = minsamplelag:maxsamplelag;          % Corresponding to 50 ms to 250 ms
                all_msd_raw = [];
                all_msd_det = [];
                all_n_raw = [];
                all_n_det = [];
                trend_slopes_X = [];
                trend_slopes_Y = [];
                trend_intercepts_X = [];
                trend_intercepts_Y = [];

                slowPhases = slowPhases(slowPhases.DurationMs >= 200,:); % keep only long fixations (>100ms)

                %% Main analysis loop over slow phases
                for i = 1:size(slowPhases,1)
                    % Extract data for current slow phase
                    startInd = slowPhases.StartIndex(i) + 50;
                    endInd = slowPhases.EndIndex(i);

                    if strcmp(eye, 'right')
                        meanDispX = slowPhases.Right_X_MeanVelocity(i)*(slowPhases.DurationMs(i)/1000);
                        meanDispY = slowPhases.Right_Y_MeanVelocity(i)*(slowPhases.DurationMs(i)/1000);
                        eyeX = data.RightX(startInd:endInd) * 60; % degrees to arcminutes
                        eyeY = data.RightY(startInd:endInd) * -60;
                    else
                        meanDispX = slowPhases.Left_X_MeanVelocity(i)*(slowPhases.DurationMs(i)/1000);
                        meanDispY = slowPhases.Left_Y_MeanVelocity(i)*(slowPhases.DurationMs(i)/1000);
                        eyeX = data.LeftX(startInd:endInd) * 60; % degrees to arcminutes
                        eyeY = data.LeftY(startInd:endInd) * -60;
                    end

                    % tic
                    badSamples = eyeX/60 > 15 | eyeX/60 < -15 | eyeY/60 > 15 | eyeY/60 < -15;
                    if ( sum(badSamples)>0)
                        disp('baaaaaaaaaad');
                        continue;
                    end
                    % eyeX(badSamples) = nan;
                    % eyeY(badSamples) = nan;
                    % toc

                    % Skip if too many NaNs
                    if mean(isnan(eyeX) | isnan(eyeY)) > 0.1
                        disp('baaaaaaaaaad');
                        continue
                    end


                    % Build time vector
                    N = length(eyeX);
                    t = (0:N-1)' * dt;

                    %% --- Fit and remove linear trend ---
                    p_x = polyfit(t, eyeX, 1);
                    p_y = polyfit(t, eyeY, 1);

                    trend_slopes_X(i,1) = p_x(1);
                    trend_slopes_Y(i,1) = p_y(1);
                    trend_intercepts_X(i,1) = p_x(2);
                    trend_intercepts_Y(i,1) = p_y(2);

                    eyeX_detrended = eyeX - polyval(p_x, t);
                    eyeY_detrended = eyeY - polyval(p_y, t);

                    %% --- Compute MSDs for raw and detrended traces ---
                    [taus_raw, msd_raw, n_raw] = compute_msd(eyeX, eyeY, dt, lags);
                    [taus_det, msd_det, n_det] = compute_msd(eyeX_detrended, eyeY_detrended, dt, lags);

                    if isempty(avg_taus)
                        avg_taus = taus_raw;
                    end

                    all_msd_raw(i,:) = msd_raw;
                    all_msd_det(i,:) = msd_det;
                    all_n_raw(i,:) = n_raw;
                    all_n_det(i,:) = n_det;

                    all_meanDispY = [all_meanDispY; meanDispY*-1];
                    all_meanDispX = [all_meanDispX; meanDispX];
                    all_dur = [all_dur; slowPhases.DurationMs(i)/1000];


                    % Plot and label based on condition
                    if strcmp(condition, 'fix')
                        coord = slowPhases.TargetPosition{i};
                        if ~isempty(coord)
                            plot(taus_raw, msd_raw, '-', 'Color', getLineColor(coord, colors), 'HandleVisibility', 'off');
                            num = num + 1;
                            xlabel(['Number of phases: ' num2str(num)]);
                            title(['Eye: ' num2str(eye) ' Diffusion per Slow Phase (not detrended)']);

                            % Estimate diffusion constant from linear fit slope
                            validIdx = ~isnan(msd_raw);
                            if sum(validIdx) < 2
                                disp('baaaaaaaaaad');
                                continue;  % Not enough points to fit a line
                            end

                            p = polyfit(taus_raw(validIdx), msd_raw(validIdx), 1);
                            D_est = p(1) / 4; % arcmin^2/sec
                            D_est_all = [D_est_all; D_est];
                            timeCoord = [timeCoord; time(startInd)];

                            % Label coordinate and axis/quadrant
                            coord_labels{end+1,1} = mat2str(coord);
                            axis_labels{end+1,1} = getAxisLabel(coord);
                        end
                    else
                        plot(taus_raw, msd_raw);
                        title('Diffusion per Slow Phase');
                        xlabel(['Number of phases: ' num2str(num)]);
                        num = num + 1;
                    end
                    hold on
                end

                %% Plot boxplots if condition == 'fix'
                if strcmp(condition, 'fix')
                    plotBoxplots(D_est_all, coord_labels, axis_labels, eye);
                    title(['Eye: ' num2str(eye) ' Average MSD by coord/axis across included slow phases: ' num2str(num) ' (not detrended)']);
                end

                %% Plot average MSD across all slow phases
                if ~isempty(all_msd_raw)
                    figure;
                    subplot(2,1,1);
                    msd_raw_avg = mean(all_msd_raw, 'omitnan', Weights=all_n_raw);
                    msd_det_avg = mean(all_msd_det, 'omitnan', Weights=all_n_det);

                    %% --- Plot comparison ---
                    
                    plot(avg_taus, msd_raw_avg, 'b', 'LineWidth', 1.5); hold on;
                    plot(avg_taus, msd_det_avg, 'k--', 'LineWidth', 1.5);

                    xlabel('Time lag (s)');
                    ylabel('Average MSD (arcmin^2)');
                    title(['Eye: ' num2str(eye) ' Weighted average MSD across ' num2str(num) ' slow phases']);
                    grid on;

                    % --- Fit and overlay linear trends for diffusion constants ---
                    p_raw = polyfit(avg_taus, msd_raw_avg, 1);
                    p_det = polyfit(avg_taus, msd_det_avg, 1);
                    plot(avg_taus, polyval(p_raw, avg_taus), 'b:','LineWidth', 1);
                    plot(avg_taus, polyval(p_det, avg_taus), 'k:','LineWidth', 1);
                    legend('Raw MSD (systematic + random)', 'Detrended MSD (random only)', 'raw trend linear fit', 'detrended linear fit');

                    D_raw = p_raw(1)/4;
                    D_det = p_det(1)/4;

                    text(avg_taus(end)*0.6, max(msd_raw_avg)*0.8, ...
                        {['D_{raw} = ' num2str(D_raw, '%.3f') ' arcmin^2/s'], ...
                        ['D_{detrended} = ' num2str(D_det, '%.3f') ' arcmin^2/s'], ...
                        ['Reduction = ' num2str(100*(1 - D_det/D_raw), '%.1f') '%']});
                end


                %% Plot average MSD across all slow phases
                if ~isempty(all_msd_raw)
                    subplot(2,1,2);
                    msd_raw_avg = mean(all_msd_raw, 'omitnan');
                    msd_det_avg = mean(all_msd_det, 'omitnan');

                    %% --- Plot comparison ---
                    
                    plot(avg_taus, msd_raw_avg, 'b', 'LineWidth', 1.5); hold on;
                    plot(avg_taus, msd_det_avg, 'k--', 'LineWidth', 1.5);

                    xlabel('Time lag (s)');
                    ylabel('Average MSD (arcmin^2)');
                    title(['Eye: ' num2str(eye) ' Unweighted average MSD across ' num2str(num) ' slow phases']);
                    grid on;

                    % --- Fit and overlay linear trends for diffusion constants ---
                    p_raw = polyfit(avg_taus, msd_raw_avg, 1);
                    p_det = polyfit(avg_taus, msd_det_avg, 1);
                    plot(avg_taus, polyval(p_raw, avg_taus), 'b:','LineWidth', 1);
                    plot(avg_taus, polyval(p_det, avg_taus), 'k:','LineWidth', 1);
                    legend('Raw MSD (systematic + random)', 'Detrended MSD (random only)', 'raw linear fit', 'detrended linear fit');

                    D_raw = p_raw(1)/4;
                    D_det = p_det(1)/4;

                    text(avg_taus(end)*0.6, max(msd_raw_avg)*0.8, ...
                        {['D_{raw} = ' num2str(D_raw, '%.3f') ' arcmin^2/s'], ...
                        ['D_{detrended} = ' num2str(D_det, '%.3f') ' arcmin^2/s'], ...
                        ['Reduction = ' num2str(100*(1 - D_det/D_raw), '%.1f') '%']});
                end

                if exist('trend_slopes_X', 'var')
                    figure;
                    subplot(2,2,1);
                    histogram(trend_slopes_X, 30);
                    xlabel('Systematic drift velocity X (arcmin/s)');
                    ylabel('Count');
                    title(['Eye: ' num2str(eye) ' Distribution of X drift slopes']);

                    subplot(2,2,2);
                    histogram(trend_slopes_Y, 30);
                    xlabel('Systematic drift velocity Y (arcmin/s)');
                    ylabel('Count');
                    title(['Eye: ' num2str(eye) ' Distribution of Y drift slopes']);

                    subplot(2,2,[3 4]);
                    hold on;
                    quiver(zeros(size(trend_slopes_X)), zeros(size(trend_slopes_Y)), ...
                        trend_slopes_X, trend_slopes_Y, 0, ...
                        'Color', [0.7 0.7 0.7], 'MaxHeadSize', 0.1);

                    % Mean drift vector in red
                    quiver(0, 0, mean(trend_slopes_X, 'omitnan'), mean(trend_slopes_Y, 'omitnan'), ...
                        0, 'r', 'LineWidth', 2, 'MaxHeadSize', 0.5);

                    axis equal; grid on;
                    xlabel('X drift velocity (arcmin/s)');
                    ylabel('Y drift velocity (arcmin/s)');
                    title(['Eye: ' num2str(eye) ' Systematic drift directions (gray) and mean (red)']);
                end


                %% weighted average displacement
                dispX = mean(all_meanDispX, 'omitnan', Weights=all_dur);
                dispY = mean(all_meanDispY, 'omitnan', Weights=all_dur);

                %% --- Helper Functions --- %%
            end
            function colors = defineColors()
                colors.posX = [1 0 0];      % Red
                colors.negX = [0 0 1];      % Blue
                colors.posY = [1 0.5 0];    % Orange+
                colors.negY = [0 1 0];      % Green

                colors.Q1 = (colors.posX + colors.posY) / 2;
                colors.Q2 = (colors.negX + colors.posY) / 2;
                colors.Q3 = (colors.negX + colors.negY) / 2;
                colors.Q4 = (colors.posX + colors.negY) / 2;

                colors.origin = [0 0 0];    % Black
                colors.gray = [0.5 0.5 0.5];
            end

            function plotLegend(colors)
                figure; hold on;
                % Dummy points for legend
                scatter(NaN, NaN, 100, colors.posX, 'filled');
                scatter(NaN, NaN, 100, colors.negX, 'filled');
                scatter(NaN, NaN, 100, colors.posY, 'filled');
                scatter(NaN, NaN, 100, colors.negY, 'filled');
                scatter(NaN, NaN, 100, colors.Q1, 'filled');
                scatter(NaN, NaN, 100, colors.Q2, 'filled');
                scatter(NaN, NaN, 100, colors.Q3, 'filled');
                scatter(NaN, NaN, 100, colors.Q4, 'filled');
                scatter(NaN, NaN, 100, colors.origin, 'filled');
                legend({'+X Axis', '-X Axis', '+Y Axis', '-Y Axis', ...
                    'Q1 (+X,+Y)', 'Q2 (-X,+Y)', 'Q3 (-X,-Y)', 'Q4 (+X,-Y)', 'Origin'}, ...
                    'Location', 'northeastoutside', 'Box', 'on');
                hold off;
            end

            function [taus, msd, n] = compute_msd(x, y, dt, lags)  % edited to collect systematic drift
                N = length(x);
                taus = lags * dt;
                msd = nan(size(lags));
                n = nan(size(lags));
                for idx = 1:length(lags)
                    lag = lags(idx);
                    dx = x(1+lag:N) - x(1:N-lag);
                    dy = y(1+lag:N) - y(1:N-lag);
                    msd(idx) = mean(dx.^2 + dy.^2, 'omitnan');
                    n(idx) = sum(~isnan(dx.^2 + dy.^2));
                end
            end

            function [taus, msd, n] = compute_msd_detrended(x, y, dt, lags)
                % Remove overall linear trend before computing MSD
                N = length(x);

                % Build time vector
                t = (0:N-1)' * dt;

                % Fit linear trend to x and y
                p_x = polyfit(t, x, 1);  % p_x(1) = slope (velocity)
                p_y = polyfit(t, y, 1);

                % Subtract trend from x and y (i.e., detrend)
                x_detrended = x - polyval(p_x, t);
                y_detrended = y - polyval(p_y, t);

                % Compute MSD on detrended signal
                taus = lags * dt;
                msd = nan(size(lags));
                n = nan(size(lags));
                for idx = 1:length(lags)
                    lag = lags(idx);
                    dx = x_detrended(1+lag:N) - x_detrended(1:N-lag);
                    dy = y_detrended(1+lag:N) - y_detrended(1:N-lag);
                    msd(idx) = mean(dx.^2 + dy.^2, 'omitnan');
                    n(idx) = sum(~isnan(dx.^2 + dy.^2));
                end
            end


            function color = getLineColor(coord, colors)
                if isempty(coord)
                    color = colors.gray;
                    return;
                end

                x = coord(1);
                y = coord(2);
                if x > 0 && y == 0
                    color = colors.posX;
                elseif x < 0 && y == 0
                    color = colors.negX;
                elseif y > 0 && x == 0
                    color = colors.posY;
                elseif y < 0 && x == 0
                    color = colors.negY;
                elseif x > 0 && y > 0
                    color = colors.Q1;
                elseif x < 0 && y > 0
                    color = colors.Q2;
                elseif x < 0 && y < 0
                    color = colors.Q3;
                elseif x > 0 && y < 0
                    color = colors.Q4;
                elseif x == 0 && y == 0
                    color = colors.origin;
                else
                    color = colors.gray;
                end
            end

            function label = getAxisLabel(coord)
                x = coord(1);
                y = coord(2);
                if x > 0 && y == 0
                    label = '+X Axis';
                elseif x < 0 && y == 0
                    label = '-X Axis';
                elseif y > 0 && x == 0
                    label = '+Y Axis';
                elseif y < 0 && x == 0
                    label = '-Y Axis';
                elseif x > 0 && y > 0
                    label = 'Q1 (+X,+Y)';
                elseif x < 0 && y > 0
                    label = 'Q2 (-X,+Y)';
                elseif x < 0 && y < 0
                    label = 'Q3 (-X,-Y)';
                elseif x > 0 && y < 0
                    label = 'Q4 (+X,-Y)';
                elseif x == 0 && y == 0
                    label = 'Origin';
                else
                    label = 'Other';
                end
            end

            function plotBoxplots(D_est_all, coord_labels, axis_labels, eye)
                figure;
                subplot(2,1,1);
                boxplot(D_est_all, categorical(coord_labels));
                xlabel('Target Coordinate');
                ylabel('Estimated Diffusion Constant (arcmin^2/sec)');
                title(['Eye: ' num2str(eye) ' By Exact Coordinate (not detrended)']);
                ylim([-100, 500]);
                xtickangle(45);

                subplot(2,1,2);
                boxplot(D_est_all, categorical(axis_labels));
                xlabel('Axis / Quadrant');
                ylabel('Estimated Diffusion Constant (arcmin^2/sec)');
                title(['Eye: ' num2str(eye) ' Grouped by Axis / Quadrant (not detrended)']);
                ylim([-100, 500]);
                xtickangle(45);
            end
        end
    end
end