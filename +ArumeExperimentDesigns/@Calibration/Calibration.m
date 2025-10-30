classdef Calibration < ArumeExperimentDesigns.EyeTracking
    %Calibration Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
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
              
            %% Change defauls
            dlg.DisplayOptions.ScreenWidth = { 142.8 '* (cm)' [1 3000] }; % these settings are for the main eye tracking room
            dlg.DisplayOptions.ScreenHeight = { 80 '* (cm)' [1 3000] }; 
            dlg.DisplayOptions.ScreenDistance = { 85 '* (cm)' [1 3000] }; 
            dlg.TrialDuration =  { 3 '* (s)' [1 100] };
            dlg.DisplayOptions.SelectedScreen = { 1 '* (screen)' [0 5] };
            dlg.DisplayOptions.StereoMode = { 0 '* (mode)' [0 9] };  % 0=no stereo, 4=stereo
            dlg.NumberOfRepetitions = 1;
            
            %% Add new options
            dlg.TargetSize = 3;
            dlg.Calibration_Type = { {'5 dots' '9 dots' '13 dots' '17 dots' '{Stereo}'} };
            dlg.Calibration_Distance_H = { 5 '* (deg)' [1 3000] };
            dlg.Calibration_Distance_V = { 5 '* (deg)' [1 3000] };
            dlg.BackgroundBrightness = 255/2; % scalar between 0 (black) and 255 (white)
        end


        function trialTable = SetUpTrialTable(this)
            
            % Update the screen parameters for stereo automatically if
            % you're doing a stereo calibration
            if (this.ExperimentOptions.Calibration_Type == "Stereo")
                this.ExperimentOptions.DisplayOptions.ScreenWidth = 60;
                this.ExperimentOptions.DisplayOptions.ScreenHeight = 33.5;
                this.ExperimentOptions.DisplayOptions.ScreenDistance = 57;
            end
            
            % Prep the target positions depending on the calibration type
            h = this.ExperimentOptions.Calibration_Distance_H;
            v = this.ExperimentOptions.Calibration_Distance_V;
            temp = 1.5;
            switch(this.ExperimentOptions.Calibration_Type)
                
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
                case 'Stereo'
                    this.targetPositions = {[0,0],[h,0],[-h,0],[0,v],[0,-v],[h,v],[h,-v],[-h,v],[-h,-v]};
            end
            
            % Build the trial table using updated methods 
            t = ArumeCore.TrialTableBuilder();
            t.AddConditionVariable( 'TargetPosition', 1:length(this.targetPositions));
            if this.ExperimentOptions.Calibration_Type == "Stereo" % if the calibration is a stereo calibration, add another column for eye
                t.AddConditionVariable( 'Eye', ["right" "left"] ); 
                t.AddBlock(find(t.ConditionTable.Eye=="right"), 1);
                t.AddBlock(find(t.ConditionTable.Eye=="left"), 1);
            end
            trialSequence = 'Random';
            blockSequence =  'Random';
            blockSequenceRepeatitions = this.ExperimentOptions.NumberOfRepetitions;
            abortAction = 'Repeat';
            trialsPerSession = 100000;
            trialTable = t.GenerateTrialTable(trialSequence, blockSequence,blockSequenceRepeatitions, abortAction,trialsPerSession);
            
        end

        function [trialResult, thisTrialData] = runTrial( this, thisTrialData )
            try
                
                Enum = ArumeCore.ExperimentDesign.getEnum();
                graph = this.Graph;
                trialResult = Enum.trialResult.CORRECT;
                
                lastFlipTime        = GetSecs;
                secondsRemaining    = this.ExperimentOptions.TrialDuration;
                thisTrialData.TimeStartLoop = lastFlipTime;
                
                if ( ~isempty(this.eyeTracker) )
                    thisTrialData.EyeTrackerFrameStartLoop = this.eyeTracker.RecordEvent(sprintf('TRIAL_START_LOOP %d %d %d', thisTrialData.TrialNumber, thisTrialData.Condition, thisTrialData.TargetPosition) );
                end
                
                while secondsRemaining > 0

                    secondsElapsed      = GetSecs - thisTrialData.TimeStartLoop;
                    secondsRemaining    = this.ExperimentOptions.TrialDuration - secondsElapsed;
                    
                    
                    % -----------------------------------------------------------------
                    % --- Drawing of stimulus -----------------------------------------
                    % -----------------------------------------------------------------
                    
                    % If the calibration is a stereo calibration
                    if this.ExperimentOptions.Calibration_Type == "Stereo"
                        ok=cell2mat(this.targetPositions(thisTrialData.TargetPosition));
                        targetHPix = this.Graph.pxWidth/(this.ExperimentOptions.DisplayOptions.ScreenWidth/2) * this.ExperimentOptions.DisplayOptions.ScreenDistance * tand(ok(1));
                        targetVPix = this.Graph.pxWidth/(this.ExperimentOptions.DisplayOptions.ScreenWidth/2) * this.ExperimentOptions.DisplayOptions.ScreenDistance * tand(ok(2));
                        if thisTrialData.Eye == "left"
                            % Draw left stim:
                            Screen('SelectStereoDrawBuffer', this.Graph.window, 0);
                            Screen('FillRect', this.Graph.window, this.ExperimentOptions.BackgroundBrightness);
                            Screen('DrawDots', this.Graph.window, [targetHPix; targetVPix], this.ExperimentOptions.TargetSize, this.fixColor, this.Graph.wRect(3:4)/2, 1); % fixation spot
                            %Screen('FrameRect', this.Graph.window, [1 0 0], [], 5);
                            
                        elseif thisTrialData.Eye == "right"
                            % Draw right stim:
                            Screen('SelectStereoDrawBuffer', this.Graph.window, 1);
                            Screen('FillRect', this.Graph.window, this.ExperimentOptions.BackgroundBrightness);
                            Screen('DrawDots', this.Graph.window, [targetHPix; targetVPix], this.ExperimentOptions.TargetSize, this.fixColor, this.Graph.wRect(3:4)/2, 1); % fixation spot
                            %Screen('FrameRect', this.Graph.window, [0 1 0], [], 5);
                        end
                    
                    % If the calibration is a non-stereo calibration
                    else
                        Screen('FillRect', graph.window, this.ExperimentOptions.BackgroundBrightness);
                        %-- Draw fixation spot
                        [mx, my] = RectCenter(graph.wRect);
                        targetPix = graph.pxWidth/this.ExperimentOptions.DisplayOptions.ScreenWidth * this.ExperimentOptions.DisplayOptions.ScreenDistance * tand(this.ExperimentOptions.TargetSize);
                        targetHPix = graph.pxWidth/this.ExperimentOptions.DisplayOptions.ScreenWidth * this.ExperimentOptions.DisplayOptions.ScreenDistance * tand(this.targetPositions{thisTrialData.TargetPosition}(1));
                        targetYPix = graph.pxWidth/this.ExperimentOptions.DisplayOptions.ScreenWidth * this.ExperimentOptions.DisplayOptions.ScreenDistance * tand(this.targetPositions{thisTrialData.TargetPosition}(2));
                        fixRect = [0 0 targetPix/2 targetPix/2];
                        fixRect = CenterRectOnPointd( fixRect, mx+targetHPix/2, my+targetYPix/2 );
                        Screen('FillOval', graph.window, this.fixColor, fixRect);
                        fixRectCenter = CenterRectOnPointd( fixRect./2, mx+targetHPix/2, my+targetYPix/2 );
                        Screen('FillOval', graph.window, [250,250,250], fixRectCenter);
                        Screen('DrawingFinished', graph.window); % Tell PTB that no further drawing commands will follow before Screen('Flip')
                    end
                    
                    % -----------------------------------------------------------------
                    % --- END Drawing of stimulus -------------------------------------
                    % -----------------------------------------------------------------
                    
                    % -----------------------------------------------------------------
                    % -- Flip buffers to refresh screen -------------------------------
                    % -----------------------------------------------------------------
                    this.Graph.Flip();
                    % -----------------------------------------------------------------
                    
                end
            catch ex
                rethrow(ex)
            end
            
        end        
    end
    
    methods ( Access = public )
         function [analysisResults, samplesDataTable, trialDataTable, sessionDataTable]  = RunDataAnalyses(this, analysisResults, samplesDataTable, trialDataTable, sessionDataTable, options)
             
             [analysisResults, samplesDataTable, trialDataTable, sessionDataTable]  = RunDataAnalyses@ArumeExperimentDesigns.EyeTracking(this, analysisResults, samplesDataTable, trialDataTable, sessionDataTable, options);

             targetPositions_ = cell2mat(this.targetPositions');

             calibrationPointsX = targetPositions_(trialDataTable.TargetPosition,1);
             calibrationPointsY = targetPositions_(trialDataTable.TargetPosition,2);

             %%%%%%%%%%% SR 3/13/25 THIS IS TO PROCESS OLDER EXPTS WITH A
             %%%%%%%%%%% CALIBRATION THAT HAD THE YS FLIPPED. THIS WILL
             %%%%%%%%%%% ULTIMATELY BE FIXED IN THE CALIBRATION DISPLAY
             %%%%%%%%%%% PART TO AVOID PROBLEMS IN THE FUTURE (will
             %%%%%%%%%%% implement that fix once OST Vergence Free Viewing
             %%%%%%%%%%% is done collecting data)
             if this.ExperimentOptions.Calibration_Type == "Stereo"
                calibrationPointsY=-calibrationPointsY;
             end
             
             % Get the indices from the start of the trial+0.5 seconds to
             % the end of the trial 
             fstart = round(trialDataTable.SampleStartTrial + 0.500*samplesDataTable.Properties.UserData.sampleRate);
             fstops = trialDataTable.SampleStopTrial;
             
             % Make a new "sample table"  (targetPosition) that has the 
             % calibration pos per sample row
             % If it's a stereo calibration, then 
             if ismember('Eye', trialDataTable.Properties.VariableNames)
                 
                 t = nan(size(samplesDataTable,1),4);
                 if trialDataTable.Eye(1) =="left"
                    calibrationPointsLeft = [[calibrationPointsX(1:length(calibrationPointsX)/2); nan(length(calibrationPointsX)/2,1)], [calibrationPointsY(1:length(calibrationPointsY)/2); nan(length(calibrationPointsY)/2,1)]];
                    calibrationPointsRight = [[nan(length(calibrationPointsX)/2,1); calibrationPointsX(length(calibrationPointsX)/2+1:end)], [nan(length(calibrationPointsY)/2,1); calibrationPointsY(length(calibrationPointsY)/2+1:end)]];
                 elseif trialDataTable.Eye(1) == "right"
                    calibrationPointsLeft = [[nan(length(calibrationPointsX)/2,1); calibrationPointsX(length(calibrationPointsX)/2+1:end)], [nan(length(calibrationPointsY)/2,1); calibrationPointsY(length(calibrationPointsY)/2+1:end)]];
                    calibrationPointsRight = [[calibrationPointsX(1:length(calibrationPointsX)/2); nan(length(calibrationPointsX)/2,1)], [calibrationPointsY(1:length(calibrationPointsY)/2); nan(length(calibrationPointsY)/2,1)]];
                 else
                     disp('Debug here')
                 end
                 for i=1:length(fstart)
                     t(fstart(i):fstops(i),1) = calibrationPointsLeft(i,1); %x
                     t(fstart(i):fstops(i),2) = calibrationPointsLeft(i,2); %y
                     t(fstart(i):fstops(i),3) = calibrationPointsRight(i,1); %x
                     t(fstart(i):fstops(i),4) = calibrationPointsRight(i,2); %y
                 end
                 samplesDataTable.Target_LeftX = t(:,1);
                 samplesDataTable.Target_LeftY = t(:,2);
                 samplesDataTable.Target_RightX = t(:,3);
                 samplesDataTable.Target_RightY = t(:,4);
                 
             % If it's not a stereo calibration, make this targetPosition
             % table like normal 
             else
                 t = nan(size(samplesDataTable,1),2);
                 for i=1:length(fstart)
                     t(fstart(i):fstops(i),1) = calibrationPointsX(i);
                     t(fstart(i):fstops(i),2) = calibrationPointsY(i);
                 end
                 
                 samplesDataTable.Target_x = t(:,1);
                 samplesDataTable.Target_y = t(:,2);
                 samplesDataTable.Target_LeftX = t(:,1);
                 samplesDataTable.Target_LeftY = t(:,2);
                 samplesDataTable.Target_RightX = t(:,1);
                 samplesDataTable.Target_RightY = t(:,2);
             end
             
             analysisResults.calibrationTableDPI = VOGAnalysis.CalculateCalibrationDPI(samplesDataTable);
             analysisResults.calibrationTableCR = VOGAnalysis.CalculateCalibrationCR(samplesDataTable);
             analysisResults.calibrationTable = VOGAnalysis.CalculateCalibration(samplesDataTable);
             
         end
    end
      
    % ---------------------------------------------------------------------
    % Plot methods
    % ---------------------------------------------------------------------
    methods ( Access = public )
        
        function Plot_CalibrationDebug(this)


            analysisResults = this.Session.analysisResults;
            samplesDataTable = this.Session.samplesDataTable;
            trialDataTable = this.Session.trialDataTable;
            sessionDataTable = this.Session.sessionDataTable;

            targetPositions_ = cell2mat(this.targetPositions');

            calibrationPointsX = targetPositions_(trialDataTable.TargetPosition,1);
            calibrationPointsY = targetPositions_(trialDataTable.TargetPosition,2);

            % Get the indices from the start of the trial+0.5 seconds to
            % the end of the trial
            fstart = round(trialDataTable.SampleStartTrial + 0.500*samplesDataTable.Properties.UserData.sampleRate);
            fstops = trialDataTable.SampleStopTrial;

            % Make a new "sample table"  (targetPosition) that has the
            % calibration pos per sample row
            % If it's a stereo calibration
            if ismember('Eye', trialDataTable.Properties.VariableNames)
                t = nan(size(samplesDataTable,1),4);
                if trialDataTable.Eye(1) =="left"
                    calibrationPointsLeft = [[calibrationPointsX(1:length(calibrationPointsX)/2); nan(length(calibrationPointsX)/2,1)], [calibrationPointsY(1:length(calibrationPointsY)/2); nan(length(calibrationPointsY)/2,1)]];
                    calibrationPointsRight = [[nan(length(calibrationPointsX)/2,1); calibrationPointsX(length(calibrationPointsX)/2+1:end)], [nan(length(calibrationPointsY)/2,1); calibrationPointsY(length(calibrationPointsY)/2+1:end)]];
                elseif trialDataTable.Eye(1) == "right"
                    calibrationPointsLeft = [[nan(length(calibrationPointsX)/2,1); calibrationPointsX(length(calibrationPointsX)/2+1:end)], [nan(length(calibrationPointsY)/2,1); calibrationPointsY(length(calibrationPointsY)/2+1:end)]];
                    calibrationPointsRight = [[calibrationPointsX(1:length(calibrationPointsX)/2); nan(length(calibrationPointsX)/2,1)], [calibrationPointsY(1:length(calibrationPointsY)/2); nan(length(calibrationPointsY)/2,1)]];
                else
                    disp('Debug here')
                end
                for i=1:length(fstart)
                    t(fstart(i):fstops(i),1) = calibrationPointsLeft(i,1); %x
                    t(fstart(i):fstops(i),2) = calibrationPointsLeft(i,2); %y
                    t(fstart(i):fstops(i),3) = calibrationPointsRight(i,1); %x
                    t(fstart(i):fstops(i),4) = calibrationPointsRight(i,2); %y
                end
                targetPosition = table();
                targetPosition.LeftX = t(:,1);
                targetPosition.LeftY = t(:,2);
                targetPosition.RightX = t(:,3);
                targetPosition.RightY = t(:,4);

                % If it's not a stereo calibration, make this targetPosition
                % table like normal
            else
                t = nan(size(samplesDataTable,1),2);
                for i=1:length(fstart)
                    t(fstart(i):fstops(i),1) = calibrationPointsX(i);
                    t(fstart(i):fstops(i),2) = calibrationPointsY(i);
                end

                targetPosition = table();
                targetPosition.x = t(:,1);
                targetPosition.y = t(:,2);
                targetPosition.LeftX = t(:,1);
                targetPosition.LeftY = t(:,2);
                targetPosition.RightX = t(:,1);
                targetPosition.RightY = t(:,2);
            end


            samplesDataTable.LeftX = samplesDataTable.LeftX_UNCALIBRATED;
            samplesDataTable.LeftY = samplesDataTable.LeftY_UNCALIBRATED;
            samplesDataTable.RightX = samplesDataTable.RightX_UNCALIBRATED;
            samplesDataTable.RightY = samplesDataTable.RightY_UNCALIBRATED;


            calibratedCalibrationData   = VOGAnalysis.CalibrateDataDPI(samplesDataTable, analysisResults.calibrationTableDPI);
            PlotCalibrationDPI(analysisResults.calibrationTableDPI, samplesDataTable, calibratedCalibrationData, targetPosition)
            title('DPI calibration')


            calibratedCalibrationData   = VOGAnalysis.CalibrateDataCR(samplesDataTable, analysisResults.calibrationTableCR);
            PlotCalibrationCR(analysisResults.calibrationTableCR, samplesDataTable, calibratedCalibrationData, targetPosition)
            title('CR calibration')

            calibratedCalibrationData   = VOGAnalysis.CalibrateData(samplesDataTable, analysisResults.calibrationTable);
            PlotCalibration(analysisResults.calibrationTable, samplesDataTable, calibratedCalibrationData, targetPosition)
            title('pupil calibration')
        end

        % needs to be modified for large saccades
        function [out, options] = PlotAggregate_MircoCorrelations(this, sessions, options)


            out = [];
            if ( nargin == 1 )
                options = this.PlotAggregate_MircoCorrelations('get_defaults');
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
            end


            %% create the corr table

            % combine all data
            dataAll.X_Vergence = [];
            dataAll.Y_Vergence = [];
            dataAll.T_Vergence = [];

            dataAll.X_Version = [];
            dataAll.Y_Version = [];
            dataAll.T_Version = [];

            for subj = 1:length(data)
                dataAll.X_Vergence = [dataAll.X_Vergence; data(subj).X_Vergence];
                dataAll.Y_Vergence = [dataAll.Y_Vergence; data(subj).Y_Vergence];
                dataAll.T_Vergence = [dataAll.T_Vergence; data(subj).T_Vergence];

                dataAll.Y_Version = [dataAll.Y_Version; data(subj).Y_Version];
                dataAll.X_Version = [dataAll.X_Version; data(subj).X_Version];
                dataAll.T_Version = [dataAll.T_Version; data(subj).T_Version];
            end

            dataAll_ = [...
                dataAll.X_Vergence,...
                dataAll.Y_Vergence,...
                dataAll.T_Vergence,...
                dataAll.X_Version,...
                dataAll.Y_Version,...
                dataAll.T_Version...
                ];

            [R,P] = corr(dataAll_,'rows','complete');

            c = 0;
            figure,
            for i = 1:6
                for j = 1:6
                    c = c + 1;
                    h(i,j) = subplot(6,6,c);
                    if i >= j, continue; end
                    plot(dataAll_(:,j),dataAll_(:,i),'.');
                    xlim([-1,1])
                    ylim([-1,1])
                    title(strcat("R = ",num2str(round(R(c),2))))
                end
            end

            linkaxes(h(:))

            n = [...
                "H Vergence",...
                "V Vergence",...
                "T Vergence",...
                "H Version",...
                "V Version",...
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

            sgtitle("Center Fixational Micro Saccades Amplitude (degree)")
        end


    end
end


        function PlotCalibration(calibrationCoefficients, rawCalibrationData, calibratedCalibrationData, targetPosition)
            
%             rawCalibrationData.LeftX = rawCalibrationData.LeftX - rawCalibrationData .LeftCR1X;
%             rawCalibrationData.LeftY = rawCalibrationData.LeftY - rawCalibrationData.LeftCR1Y;
%             rawCalibrationData.RightX = rawCalibrationData.RightX - rawCalibrationData.RightCR1X;
%             rawCalibrationData.RightY = rawCalibrationData.RightY - rawCalibrationData.RightCR1Y;
            
            
            figure
            t = calibratedCalibrationData.Time;
            
            subplot(4,3,1,'nextplot','add')
            title('Targets X/Y (deg)');
            plot(targetPosition.LeftX,targetPosition.LeftY,'.')
            
            subplot(4,3,2,'nextplot','add')
            title('Eye positions X/Y (pixels)');
            plot(rawCalibrationData.LeftX_UNCALIBRATED(~isnan(targetPosition.LeftX)),rawCalibrationData.LeftY_UNCALIBRATED(~isnan(targetPosition.LeftY)),'.');
            plot(rawCalibrationData.RightX_UNCALIBRATED(~isnan(targetPosition.RightX)),rawCalibrationData.RightY_UNCALIBRATED(~isnan(targetPosition.RightY)),'.');
            
            subplot(4,3,3,'nextplot','add')
            title('Eye positions in time');
            plot(t,rawCalibrationData.LeftX_UNCALIBRATED);
            plot(t,rawCalibrationData.RightX_UNCALIBRATED);
            plot(t,rawCalibrationData.LeftY_UNCALIBRATED);
            plot(t,rawCalibrationData.RightY_UNCALIBRATED);
            plot(t,targetPosition.LeftX);
            plot(t,targetPosition.LeftY);
            
            x = [-30 30];
            y = [-22 22];
            subplot(4,3,4,'nextplot','add')
            title('Horizontal target vs eye pos. (deg)');
            %h1 = plot(targetPosition.LeftX+randn(size(t))/5,rawCalibrationData.LeftX_UNCALIBRATED,'.');
            h1 = plot((targetPosition.LeftX(~isnan(targetPosition.LeftX)))+(randn(size(targetPosition.LeftX(~isnan(targetPosition.LeftX))))/5),rawCalibrationData.LeftX_UNCALIBRATED(~isnan(targetPosition.LeftX)),'.')
            %h2 = plot(targetPosition.RightX+randn(size(t))/5,rawCalibrationData.RightX_UNCALIBRATED,'.');
            h2 = plot((targetPosition.RightX(~isnan(targetPosition.RightX)))+(randn(size(targetPosition.RightX(~isnan(targetPosition.RightX))))/5),rawCalibrationData.RightX_UNCALIBRATED(~isnan(targetPosition.RightX)),'.')
            line(x,calibrationCoefficients.OffsetX('LeftEye')+calibrationCoefficients.GainX('LeftEye')*x,'color',get(h1,'color'));
            line(x,calibrationCoefficients.OffsetX('RightEye')+calibrationCoefficients.GainX('RightEye')*x,'color',get(h2,'color'));
            
            subplot(4,3,5,'nextplot','add')
            title('Vertical target vs eye pos. (deg)');
            %h1 = plot(targetPosition.LeftY+randn(size(t))/5,rawCalibrationData.LeftY_UNCALIBRATED,'.');
            h1 = plot((targetPosition.LeftY(~isnan(targetPosition.LeftY)))+(randn(size(targetPosition.LeftY(~isnan(targetPosition.LeftY))))/5),rawCalibrationData.LeftY_UNCALIBRATED(~isnan(targetPosition.LeftY)),'.')
            %h2 = plot(targetPosition.RightY+randn(size(t))/5,rawCalibrationData.RightY_UNCALIBRATED,'.');
            h2 = plot((targetPosition.RightY(~isnan(targetPosition.RightY)))+(randn(size(targetPosition.RightY(~isnan(targetPosition.RightY))))/5),rawCalibrationData.RightY_UNCALIBRATED(~isnan(targetPosition.RightY)),'.')
            line(y,calibrationCoefficients.OffsetY('LeftEye')+calibrationCoefficients.GainY('LeftEye')*y,'color',get(h1,'color'));
            line(y,calibrationCoefficients.OffsetY('RightEye')+calibrationCoefficients.GainY('RightEye')*y,'color',get(h2,'color'));
            
            subplot(4,3,6,'nextplot','add')
            title('Targets and eye positions X/Y (deg)');
            plot(calibratedCalibrationData.LeftX(~isnan(targetPosition.LeftX)),calibratedCalibrationData.LeftY(~isnan(targetPosition.LeftY)),'.','markersize',1)
            plot(calibratedCalibrationData.RightX(~isnan(targetPosition.RightX)),calibratedCalibrationData.RightY(~isnan(targetPosition.RightY)),'.','markersize',1)
            plot(targetPosition.LeftX,targetPosition.LeftY,'.','markersize',20,'color','k')
            set(gca,'xlim',x,'ylim',y);
            
            subplot(4,3,[7:9],'nextplot','add')
            plot(t,calibratedCalibrationData.LeftX);
            plot(t,calibratedCalibrationData.RightX);
            plot(t,targetPosition.LeftX,'color','b','linewidth',3)
            plot(t,targetPosition.RightX,'color','r','linewidth',3)
            set(gca,'ylim',x);
            ylabel('Horizontal (deg)');
            legend('Left Eye','Right Eye')
            
            subplot(4,3,[10:12],'nextplot','add')
            plot(t,calibratedCalibrationData.LeftY);
            plot(t,calibratedCalibrationData.RightY);
            plot(t,targetPosition.LeftY,'color','b','linewidth',3)
            plot(t,targetPosition.RightY,'color','r','linewidth',3)
            set(gca,'ylim',y);
            ylabel('Vertical (deg)');
            legend('Left Eye','Right Eye')
        end

          function PlotCalibrationCR(calibrationCoefficients, rawCalibrationData, calibratedCalibrationData, targetPosition)
            
             lx = rawCalibrationData.LeftX_UNCALIBRATED;
             ly = rawCalibrationData.LeftY_UNCALIBRATED;
             rx = rawCalibrationData.RightX_UNCALIBRATED;
             ry = rawCalibrationData.RightY_UNCALIBRATED;
            
            
            figure
            t = calibratedCalibrationData.Time;
            
            subplot(4,3,1,'nextplot','add')
            title('Targets X/Y (deg)');
            plot(targetPosition.LeftX,targetPosition.LeftY,'.')
            
            subplot(4,3,2,'nextplot','add')
            title('Eye positions X/Y (pixels)');
            plot(lx,ly,'.');
            plot(rx,ry,'.');
            
            subplot(4,3,3,'nextplot','add')
            title('Eye positions in time');
            plot(t,lx);
            plot(t,rx);
            plot(t,ly);
            plot(t,ry);
            plot(t,targetPosition.LeftX);
            plot(t,targetPosition.LeftY);
            x = [-30 30];
            y = [-22 22];
            
            subplot(4,3,4,'nextplot','add')
            title('Horizontal target vs eye pos. (deg)');
            h1 = plot(targetPosition.LeftX+randn(size(t))/5,lx,'.');
            h2 = plot(targetPosition.RightX+randn(size(t))/5,rx,'.');
            line(x,calibrationCoefficients.OffsetX('LeftEye')+calibrationCoefficients.GainX('LeftEye')*x,'color',get(h1,'color'));
            line(x,calibrationCoefficients.OffsetX('RightEye')+calibrationCoefficients.GainX('RightEye')*x,'color',get(h2,'color'));
            
            subplot(4,3,5,'nextplot','add')
            title('Vertical target vs eye pos. (deg)');
            h1 = plot(targetPosition.LeftY+randn(size(t))/5,ly,'.');
            h2 = plot(targetPosition.RightY+randn(size(t))/5,ry,'.');
            line(y,calibrationCoefficients.OffsetY('LeftEye')+calibrationCoefficients.GainY('LeftEye')*y,'color',get(h1,'color'));
            line(y,calibrationCoefficients.OffsetY('RightEye')+calibrationCoefficients.GainY('RightEye')*y,'color',get(h2,'color'));
            
            subplot(4,3,6,'nextplot','add')
            title('Targets and eye positions X/Y (deg)');
            plot(calibratedCalibrationData.LeftX(~isnan(targetPosition.LeftX)),calibratedCalibrationData.LeftY(~isnan(targetPosition.LeftY)),'.','markersize',1)
            plot(calibratedCalibrationData.RightX(~isnan(targetPosition.RightX)),calibratedCalibrationData.RightY(~isnan(targetPosition.RightY)),'.','markersize',1)
            plot(targetPosition.LeftX,targetPosition.LeftY,'.','markersize',20,'color','k')
            set(gca,'xlim',x,'ylim',y);
            
            subplot(4,3,[7:9],'nextplot','add')
            plot(t,calibratedCalibrationData.LeftX);
            plot(t,calibratedCalibrationData.RightX);
            plot(t,targetPosition.LeftX,'color','b','linewidth',3)
            plot(t,targetPosition.RightX,'color','r','linewidth',3)
            set(gca,'ylim',x);
            ylabel('Horizontal (deg)');
            legend('Left Eye','Right Eye')
            
            subplot(4,3,[10:12],'nextplot','add')
            plot(t,calibratedCalibrationData.LeftY);
            plot(t,calibratedCalibrationData.RightY);
            plot(t,targetPosition.LeftY,'color','b','linewidth',3)
            plot(t,targetPosition.RightY,'color','r','linewidth',3)
            set(gca,'ylim',y);
            ylabel('Vertical (deg)');
            legend('Left Eye','Right Eye')
          end



          function PlotCalibrationDPI(calibrationCoefficients, rawCalibrationData, calibratedCalibrationData, targetPosition)
            
             lx = rawCalibrationData.LeftX_UNCALIBRATED;
             ly = rawCalibrationData.LeftY_UNCALIBRATED;
             rx = rawCalibrationData.RightX_UNCALIBRATED;
             ry = rawCalibrationData.RightY_UNCALIBRATED;
            
            
            figure
            t = calibratedCalibrationData.Time;
            
            subplot(4,3,1,'nextplot','add')
            title('Targets X/Y (deg)');
            plot(targetPosition.LeftX,targetPosition.LeftY,'.')
            
            subplot(4,3,2,'nextplot','add')
            title('Eye positions X/Y (pixels)');
            plot(lx,ly,'.');
            plot(rx,ry,'.');
            
            subplot(4,3,3,'nextplot','add')
            title('Eye positions in time');
            plot(t,rawCalibrationData.LeftCR1X);
            plot(t,rawCalibrationData.RightCR1X);
            plot(t,rawCalibrationData.LeftCR4X);
            plot(t,rawCalibrationData.RightCR4X);

            plot(t,rawCalibrationData.LeftCR1Y);
            plot(t,rawCalibrationData.RightCR1Y);
            plot(t,rawCalibrationData.LeftCR4Y);
            plot(t,rawCalibrationData.RightCR4Y);

            plot(t,targetPosition.LeftX);
            plot(t,targetPosition.LeftY);
            x = [-30 30];
            y = [-22 22];
            
            subplot(4,3,4,'nextplot','add')
            title('Horizontal target vs eye pos. (deg)');
            h1 = plot(targetPosition.LeftX+randn(size(t))/5,lx,'.');
            h2 = plot(targetPosition.RightX+randn(size(t))/5,rx,'.');
            line(x,calibrationCoefficients.OffsetX('LeftEye')+calibrationCoefficients.GainX('LeftEye')*x,'color',get(h1,'color'));
            line(x,calibrationCoefficients.OffsetX('RightEye')+calibrationCoefficients.GainX('RightEye')*x,'color',get(h2,'color'));
            
            subplot(4,3,5,'nextplot','add')
            title('Vertical target vs eye pos. (deg)');
            h1 = plot(targetPosition.LeftY+randn(size(t))/5,ly,'.');
            h2 = plot(targetPosition.RightY+randn(size(t))/5,ry,'.');
            line(y,calibrationCoefficients.OffsetY('LeftEye')+calibrationCoefficients.GainY('LeftEye')*y,'color',get(h1,'color'));
            line(y,calibrationCoefficients.OffsetY('RightEye')+calibrationCoefficients.GainY('RightEye')*y,'color',get(h2,'color'));
            
            subplot(4,3,6,'nextplot','add')
            title('Targets and eye positions X/Y (deg)');
            plot(calibratedCalibrationData.LeftX(~isnan(targetPosition.LeftX)),calibratedCalibrationData.LeftY(~isnan(targetPosition.LeftY)),'.','markersize',1)
            plot(calibratedCalibrationData.RightX(~isnan(targetPosition.RightX)),calibratedCalibrationData.RightY(~isnan(targetPosition.RightY)),'.','markersize',1)
            plot(targetPosition.LeftX,targetPosition.LeftY,'.','markersize',20,'color','k')
            set(gca,'xlim',x,'ylim',y);
            
            subplot(4,3,[7:9],'nextplot','add')
            plot(t,calibratedCalibrationData.LeftX);
            plot(t,calibratedCalibrationData.RightX);
            plot(t,targetPosition.LeftX,'color','b','linewidth',3)
            plot(t,targetPosition.RightX,'color','r','linewidth',3)
            set(gca,'ylim',x);
            ylabel('Horizontal (deg)');
            legend('Left Eye','Right Eye')
            
            subplot(4,3,[10:12],'nextplot','add')
            plot(t,calibratedCalibrationData.LeftY);
            plot(t,calibratedCalibrationData.RightY);
            plot(t,targetPosition.LeftY,'color','b','linewidth',3)
            plot(t,targetPosition.RightY,'color','r','linewidth',3)
            set(gca,'ylim',y);
            ylabel('Vertical (deg)');
            legend('Left Eye','Right Eye')
        end