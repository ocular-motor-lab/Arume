classdef FixationTargets_BEAST < ArumeExperimentDesigns.EyeTracking
    %OPTOKINETICTORSION Summary of this class goes here
    %   Detailed explanation goes here

    properties
        fixRad = 20;
        fixColor = [255 0 0];
        targetPositions =[];
        newCenter = [0, 0];
        iftrial1 = 1;
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

            dlg.DisplayOptions.ScreenWidth = { 142.8 '* (cm)' [1 3000] };
            dlg.DisplayOptions.ScreenHeight = { 80 '* (cm)' [1 3000] };
            dlg.DisplayOptions.ScreenDistance = { 208.5 '* (cm)' [1 3000] };

            dlg.TrialDuration =  { 1 '* (s)' [1 100] };
            dlg.NumberRepetitions = 5;

            dlg.TargetSize = 1;
            dlg.Experiment_Type = { {'Center dot' '2 dots - Horizontal' '{2 dots - Vertical}' '4 dots - Square'} };
            dlg.Calibration_Distance_H = { 2 '* (deg)' [1 3000] };
            dlg.Calibration_Distance_V = { 2 '* (deg)' [1 3000] };

            dlg.RasterCenter_x = {0 '* (pixel)' [-3000 3000]};
            dlg.RasterCenter_y = {0 '* (pixel)' [-3000 3000]};
            
            dlg.CrossSize = {40 '* (pixel)' [1 3000]};

            dlg.BackgroundBrightness = 0;
        end


        function trialTable = SetUpTrialTable(this)
            h = this.ExperimentOptions.Calibration_Distance_H;
            v = this.ExperimentOptions.Calibration_Distance_V;

            switch(this.ExperimentOptions.Experiment_Type)
                case 'Center dot'
                    this.targetPositions = {[0,0]};
                case '2 dots - Horizontal'
                    this.targetPositions = {[h,0],[-h,0]};
                case '2 dots - Vertical'
                    this.targetPositions = {[0,v],[0,-v]};
                case '4 dots - Square'
                    this.targetPositions = {[h,0],[-h,0],[0,v],[0,-v]};
            end

            targets = 1:length(this.targetPositions);

            %-- condition variables ---------------------------------------
            i= 0;

            i = i+1;
            conditionVars(i).name   = 'TargetPosition';
            conditionVars(i).values = targets;


            trialTableOptions = this.GetDefaultTrialTableOptions();
            trialTableOptions.trialSequence = 'Sequential';
            trialTableOptions.trialAbortAction = 'Drop';
            trialTableOptions.trialsPerSession = (length(targets))*this.ExperimentOptions.NumberRepetitions;

            trialTableOptions.blockSequence       = 'Sequential';	% Sequential, Random, Random with repetition, ...numberOfTimesRepeatBlockSequence = 1;
            trialTableOptions.blocksToRun         = 1;
            trialTableOptions.blocks                = struct( 'fromCondition', 1, 'toCondition', length(targets), 'trialsToRun', length(targets));
            %trialTableOptions.blocks(2)             = struct( 'fromCondition', 2, 'toCondition', length(targets), 'trialsToRun', length(targets)-1 );
            %trialTableOptions.blocks(3)             = struct( 'fromCondition', 1, 'toCondition', 1, 'trialsToRun', 1 );
            trialTableOptions.numberOfTimesRepeatBlockSequence = this.ExperimentOptions.NumberRepetitions;
            trialTable = this.GetTrialTableFromConditions(conditionVars, trialTableOptions);

%             %-- Find the center of Raster ---------------------------------
%             gamepad = ArumeHardware.GamePad;
%             right = 0;
%             this.newCenter = [this.ExperimentOptions.RasterCenter_x, this.ExperimentOptions.RasterCenter_y];
%             move_newCenter = 5;
%             this.Graph = ArumeCore.PTB( );
%             this.Graph.Init( this );
%             
%             while(right == 0)
%                 [ ~, ~, right, a, b, x, y] = gamepad.Query;
%                 if a == 1, this.newCenter(2) = this.newCenter(2) + move_newCenter;
%                 elseif b == 1, this.newCenter(1) = this.newCenter(1) + move_newCenter;
%                 elseif x == 1, this.newCenter(1) = this.newCenter(1) - move_newCenter;
%                 elseif y == 1, this.newCenter(2) = this.newCenter(2) - move_newCenter;
%                 elseif right == 1, gamepad.Close;
%                 end
% 
%                 drawCross(this)
%             end
%             this.ExperimentOptions.RasterCenter_x = this.newCenter(1);
%             this.ExperimentOptions.RasterCenter_x = this.newCenter(2);
% 
%             this.Graph.Clear(this.Graph);
        end

        function drawCross(this)
            graph = this.Graph;
            Screen('FillRect', graph.window, this.ExperimentOptions.BackgroundBrightness);

            %-- Draw fixation spot
            [mx, my] = RectCenter(graph.wRect);
            mx = mx + this.newCenter(1);
            my = my + this.newCenter(2);

            Screen('DrawLine', graph.window, [250,250,250], mx-40, my, mx+40, my,5);
            Screen('DrawLine', graph.window, [250,250,250], mx, my-40, mx, my+40,5);
            Screen('DrawingFinished', graph.window); % Tell PTB that no further drawing commands will follow before Screen('Flip')
            this.Graph.Flip();
        end


        function [trialResult, thisTrialData] = runTrial( this, thisTrialData )

            %for the first time
            if this.iftrial1 == 1
                this.iftrial1 = 0;
                try
                    %-- Find the center of Raster ---------------------------------
                    gamepad = ArumeHardware.GamePad;
                    right = 0;
                    this.newCenter = [this.ExperimentOptions.RasterCenter_x, this.ExperimentOptions.RasterCenter_y];
                    move_newCenter = 5;
                    this.Graph = ArumeCore.PTB( );
                    this.Graph.Init( this );

                    while(right == 0)
                        [ ~, ~, right, a, b, x, y] = gamepad.Query;
                        if a == 1, this.newCenter(2) = this.newCenter(2) + move_newCenter;
                        elseif b == 1, this.newCenter(1) = this.newCenter(1) + move_newCenter;
                        elseif x == 1, this.newCenter(1) = this.newCenter(1) - move_newCenter;
                        elseif y == 1, this.newCenter(2) = this.newCenter(2) - move_newCenter;
                        %elseif right == 1, gamepad.Close;
                        end

                        drawCross(this)
                    end
                    this.ExperimentOptions.RasterCenter_x = this.newCenter(1);
                    this.ExperimentOptions.RasterCenter_y = this.newCenter(2);

                    this.Graph.Clear(this.Graph);
                catch
                end
            end

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

                    if secondsRemaining > 0
                        % -----------------------------------------------------------------
                        % --- Drawing of stimulus -----------------------------------------
                        % -----------------------------------------------------------------
                        Screen('FillRect', graph.window, this.ExperimentOptions.BackgroundBrightness);

                        %-- Draw fixation spot
                        [mx, my] = RectCenter(this.Graph.wRect);
                        mx = this.newCenter(1) + mx;
                        my = this.newCenter(2) + my;
                        % this.Graph.pxWidth
                        % targetHPix
                        targetPix = this.Graph.pxWidth/this.ExperimentOptions.DisplayOptions.ScreenWidth * this.ExperimentOptions.DisplayOptions.ScreenDistance * tand(this.ExperimentOptions.TargetSize);
                        targetHPix = this.Graph.pxWidth/this.ExperimentOptions.DisplayOptions.ScreenWidth * this.ExperimentOptions.DisplayOptions.ScreenDistance * tand(this.targetPositions{thisTrialData.TargetPosition}(1));
                        targetYPix = this.Graph.pxWidth/this.ExperimentOptions.DisplayOptions.ScreenWidth * this.ExperimentOptions.DisplayOptions.ScreenDistance * tand(this.targetPositions{thisTrialData.TargetPosition}(2));
                        fixRect = [0 0 targetPix/2 targetPix/2];
                        fixRect = CenterRectOnPointd( fixRect, mx+targetHPix/2, my+targetYPix/2 );
                        Screen('FillOval', graph.window, this.fixColor, fixRect);

                        fixRect2 = CenterRectOnPointd( fixRect./2, mx+targetHPix/2, my+targetYPix/2 );
                        Screen('FillOval', graph.window, [250,250,250], fixRect2);

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


                    % -----------------------------------------------------------------
                    % --- Collecting responses  ---------------------------------------
                    % -----------------------------------------------------------------

                    % -----------------------------------------------------------------
                    % --- END Collecting responses  -----------------------------------
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

            fstart = round(trialDataTable.SampleStartTrial + 0.500*samplesDataTable.Properties.UserData.sampleRate);
            fstops = trialDataTable.SampleStopTrial;

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

            analysisResults.calibrationTableCR = VOGAnalysis.CalculateCalibrationCR(samplesDataTable, targetPosition);
            analysisResults.calibrationTable = VOGAnalysis.CalculateCalibration(samplesDataTable, targetPosition);

        end
    end

    % ---------------------------------------------------------------------
    % Plot methods
    % ---------------------------------------------------------------------
    methods ( Access = public )


    end
end




