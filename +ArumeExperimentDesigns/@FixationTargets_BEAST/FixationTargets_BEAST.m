classdef FixationTargets_BEAST < ArumeExperimentDesigns.EyeTracking
    %OPTOKINETICTORSION Summary of this class goes here
    %   Detailed explanation goes here

    properties
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
            dlg.NumberRepetitions = 1;

            dlg.TargetSize = 1;
            dlg.Experiment_Type = { {'BeaSTCal' '{9DotsCal - Center - 9DotsCal}' 'Torsion'} };
            dlg.Calibration_Distance_H = { 2 '* (deg)' [1 3000] };
            dlg.Calibration_Distance_V = { 2 '* (deg)' [1 3000] };

            dlg.CenterFixation_Duration = 15;

            dlg.BeaSTCal_Speed_degPerSec = 5;
            dlg.BeaSTCal_Step_deg = 0.5;

            dlg.Torsion_Direction = { {'{CW}' 'CCW'} };

            dlg.RasterCenter_x = {0 '* (pixel)' [-3000 3000]};
            dlg.RasterCenter_y = {0 '* (pixel)' [-3000 3000]};

            dlg.CrossSize = {40 '* (pixel)' [1 3000]};

            dlg.CrossColor_R = 0;
            dlg.CrossColor_G = 255;
            dlg.CrossColor_B = 255;

            dlg.BackgroundBrightness = 0;
        end


        function trialTable = SetUpTrialTable(this)
            h = this.ExperimentOptions.Calibration_Distance_H;
            v = this.ExperimentOptions.Calibration_Distance_V;

            switch(this.ExperimentOptions.Experiment_Type)
                case 'BeaSTCal'
                    this.targetPositions = {[0,0]};

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
                    trialTableOptions.blocksToRun         = 3;
                    trialTableOptions.blocks                = struct( 'fromCondition', 1, 'toCondition', length(targets), 'trialsToRun', length(targets));

                    trialTableOptions.numberOfTimesRepeatBlockSequence = this.ExperimentOptions.NumberRepetitions;
                    trialTable = this.GetTrialTableFromConditions(conditionVars, trialTableOptions);

                case '9DotsCal - Center - 9DotsCal'
                    temp = {[h,0],[-h,0],[0,v],[0,-v],[-h,v],[h,v],[-h,-v],[h,-v]};
                    %random order
                    temp = temp( randsample(length(temp),length(temp)) );
                    %adding center dot
                    this.targetPositions = {temp{1},[0,0]};
                    for i = 2:length(temp)
                        this.targetPositions = {this.targetPositions{:}, temp{i}, [0,0]};
                    end
                    this.targetPositions = {this.targetPositions{:},this.targetPositions{:}};
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
                    trialTableOptions.blocksToRun         = 3;
                    trialTableOptions.blocks                = struct( 'fromCondition', 1, 'toCondition', 15, 'trialsToRun', 15);
                    trialTableOptions.blocks(2)                = struct( 'fromCondition', 16, 'toCondition', 16, 'trialsToRun', 1);
                    trialTableOptions.blocks(3)                = struct( 'fromCondition', 17, 'toCondition', length(targets), 'trialsToRun', 16);

                    trialTableOptions.numberOfTimesRepeatBlockSequence = this.ExperimentOptions.NumberRepetitions;
                    trialTable = this.GetTrialTableFromConditions(conditionVars, trialTableOptions);
                case 'Torsion'
                    this.targetPositions = {[0,0]};
                    targets = 1;
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
                    trialTableOptions.blocks                = struct( 'fromCondition', 1, 'toCondition', 1, 'trialsToRun', 1);

                    trialTableOptions.numberOfTimesRepeatBlockSequence = this.ExperimentOptions.NumberRepetitions;
                    trialTable = this.GetTrialTableFromConditions(conditionVars, trialTableOptions);




            end

        end

        function drawCross(this)
            graph = this.Graph;
            Screen('FillRect', graph.window, this.ExperimentOptions.BackgroundBrightness);

            %-- Draw fixation spot
            [mx, my] = RectCenter(graph.wRect);
            mx = mx + this.newCenter(1);
            my = my + this.newCenter(2);

            this.fixColor = [this.ExperimentOptions.CrossColor_R,this.ExperimentOptions.CrossColor_G,this.ExperimentOptions.CrossColor_B];

            Screen('DrawLine', graph.window,this.fixColor, mx-40, my, mx+40, my,5);
            Screen('DrawLine', graph.window, this.fixColor, mx, my-40, mx, my+40,5);
            Screen('DrawingFinished', graph.window); % Tell PTB that no further drawing commands will follow before Screen('Flip')
            this.Graph.Flip();
        end


        function [trialResult, thisTrialData] = runTrial( this, thisTrialData )

            %for the first time
            if this.iftrial1 == 1

                %save("test",'thisTrialData')

                this.iftrial1 = 0;
                try
                    %-- Find the center of Raster ---------------------------------
                    gamepad = ArumeHardware.GamePad;
                    right = 0;
                    this.newCenter = [this.ExperimentOptions.RasterCenter_x, this.ExperimentOptions.RasterCenter_y];

                    move_newCenter = 1;
                    this.Graph = ArumeCore.Display( );

                    this.Graph.Init( this );

                    while(right == 0)
                        [ ~, ~, right, a, b, x, y] = gamepad.Query;
                        if a == 1, this.newCenter(1) = this.newCenter(1) + move_newCenter;
                        elseif b == 1, this.newCenter(2) = this.newCenter(2) - move_newCenter;
                        elseif x == 1, this.newCenter(2) = this.newCenter(2) + move_newCenter;
                        elseif y == 1,this.newCenter(1) = this.newCenter(1) - move_newCenter;

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
                if strcmp(this.ExperimentOptions.Experiment_Type,'BeaSTCal')
                    % calculating total time of exp based on speed
                    h = this.ExperimentOptions.Calibration_Distance_H;
                    v = this.ExperimentOptions.Calibration_Distance_V;
                    step = this.ExperimentOptions.BeaSTCal_Step_deg;
                    sp = this.ExperimentOptions.BeaSTCal_Speed_degPerSec;
                    % number of steps
                    numStep = floor( v./step);
                    % total length in deg
                    total_deg = h * (numStep+1) + v;
                    % total time for each trial
                    secondsRemaining = sp * total_deg;
                    this.ExperimentOptions.TrialDuration = secondsRemaining;
                elseif strcmp(this.ExperimentOptions.Experiment_Type,'Torsion')
                    h = this.ExperimentOptions.Calibration_Distance_H;
                    v = this.ExperimentOptions.Calibration_Distance_V;
                    sp = this.ExperimentOptions.BeaSTCal_Speed_degPerSec;

                    secondsRemaining    = this.ExperimentOptions.TrialDuration;
              
                else
                    secondsRemaining    = this.ExperimentOptions.TrialDuration;
                end

                thisTrialData.TimeStartLoop = lastFlipTime;

                if ( ~isempty(this.eyeTracker) )
                    thisTrialData.EyeTrackerFrameStartLoop = this.eyeTracker.RecordEvent(sprintf('TRIAL_START_LOOP %d %d %d', thisTrialData.TrialNumber, thisTrialData.Condition, thisTrialData.TargetPosition) );
                end

                if strcmp(this.ExperimentOptions.Experiment_Type,'BeaSTCal')
                    % rec around zero point
                    startPoint = [-h/2, v/2];
                    loc = startPoint;
                    toDown = 0;
                    toLeft = 0;

                    [mx, my] = RectCenter(this.Graph.wRect);
                    mx = this.newCenter(1) + mx;
                    my = this.newCenter(2) + my;
                    secondStartSegment = thisTrialData.TimeStartLoop;

                    while secondsRemaining > 0
                        secondsElapsed      = GetSecs - thisTrialData.TimeStartLoop;
                        secondsRemaining    = this.ExperimentOptions.TrialDuration - secondsElapsed;
                        secondsElapsedInSegment = GetSecs-secondStartSegment;
                        if secondsRemaining > 0
                            % moving straight from left to right
                            if toDown == 0 && toLeft == 0
                                newloc(1) = loc(1) + secondsElapsedInSegment * sp;
                                newloc(2) = loc(2);

                                % if it reached to the end right point
                                if newloc(1) > h/2

                                    loc(1) = h/2;
                                    loc(2) = newloc(2);

                                    % if it reach to the end down
                                    if newloc(2) - step < -v/2
                                        st_v = -v/2;
                                    else
                                        st_v = loc(2) - step;
                                    end
                                    secondStartSegment = GetSecs;

                                    toDown = 1;
                                    toLeft = 1;
                                    continue;
                                end
                            end

                            %moving down
                            if toDown == 1 && newloc(2) ~= -v/2
                                newloc(1) = loc(1);
                                newloc(2) = loc(2) - secondsElapsedInSegment * sp;

                                if newloc(2) < st_v
                                    loc(1) = newloc(1);
                                    loc(2) = st_v;
                                    secondStartSegment = GetSecs;
                                    toDown = 0;
                                    continue;
                                end
                            elseif newloc(2) == -v/2
                                break;
                            end

                            %moving right to left
                            if toDown == 0 && toLeft == 1
                                newloc(1) = loc(1) - secondsElapsedInSegment * sp;
                                newloc(2) = loc(2);

                                % if it reached to the end right point
                                if newloc(1) < -h/2

                                    loc(1) = -h/2;
                                    loc(2) = newloc(2);

                                    % if it reach to the end down
                                    if newloc(2) - step < -v/2
                                        st_v = -v/2;
                                    else
                                        st_v = loc(2) - step;
                                    end
                                    secondStartSegment = GetSecs;
                                    toDown = 1;
                                    toLeft = 0;
                                    continue;
                                end
                            end
                            %loc = newloc;
                            % -----------------------------------------------------------------
                            % --- Drawing of stimulus -----------------------------------------
                            % -----------------------------------------------------------------
                            Screen('FillRect', graph.window, this.ExperimentOptions.BackgroundBrightness);
                            %-- Draw fixation spot

                            targetPix = this.Graph.pxWidth/this.ExperimentOptions.DisplayOptions.ScreenWidth * this.ExperimentOptions.DisplayOptions.ScreenDistance * tand(this.ExperimentOptions.TargetSize);

                            targetLocHPix = this.Graph.pxWidth/this.ExperimentOptions.DisplayOptions.ScreenWidth * this.ExperimentOptions.DisplayOptions.ScreenDistance * tand(newloc(1));
                            targetLocVPix = this.Graph.pxWidth/this.ExperimentOptions.DisplayOptions.ScreenWidth * this.ExperimentOptions.DisplayOptions.ScreenDistance * tand(newloc(2));

                            fixRect = [0 0 targetPix/2 targetPix/2];
                            fixRect = CenterRectOnPointd( fixRect, mx+targetLocHPix/2, my+targetLocVPix/2 );
                            Screen('FillOval', graph.window, this.fixColor, fixRect);

                            %                         fixRect2 = CenterRectOnPointd( fixRect./2, mx+targetHPix/2, my+targetYPix/2 );
                            %                         Screen('FillOval', graph.window, [250,250,250], fixRect2);

                            Screen('DrawingFinished', graph.window); % Tell PTB that no further drawing commands will follow before Screen('Flip')
                        end
                        this.Graph.Flip();
                    end





                elseif strcmp(this.ExperimentOptions.Experiment_Type,'9DotsCal - Center - 9DotsCal')

                    while secondsRemaining > 0

                        secondsElapsed      = GetSecs - thisTrialData.TimeStartLoop;

                        if thisTrialData.BlockNumber == 2
                            secondsRemaining    = this.ExperimentOptions.CenterFixation_Duration - secondsElapsed;
                        else
                            secondsRemaining    = this.ExperimentOptions.TrialDuration - secondsElapsed;
                        end

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






                elseif strcmp(this.ExperimentOptions.Experiment_Type,'Torsion')
                        
                    if ( strcmp(this.ExperimentOptions.Torsion_Direction, 'CCW') )
                        sp = -sp;
                    end
                    % number of random dots
                    numDots = 50;
                    % init position around zero point
                    r = h .* rand(numDots,1); %only using h
                    theta = 2*pi*rand(numDots,1);  
                    %dotsInit = [r r] .* [cos(theta) sin(theta)];
                    
                    [mx, my] = RectCenter(this.Graph.wRect);
                    mx = this.newCenter(1) + mx;
                    my = this.newCenter(2) + my;
                    
                    while secondsRemaining > 0
                        secondsElapsed      = GetSecs - thisTrialData.TimeStartLoop;
                        secondsRemaining    = this.ExperimentOptions.TrialDuration - secondsElapsed;
                        if secondsRemaining > 0
                            dotsLoc = [r r] .* [cos(theta + sp*secondsElapsed) sin(theta + sp*secondsElapsed)];
                            
                            % -----------------------------------------------------------------
                            % --- Drawing of stimulus -----------------------------------------
                            % -----------------------------------------------------------------
                            Screen('FillRect', graph.window, this.ExperimentOptions.BackgroundBrightness);
                            %-- Draw fixation spot

                            targetPix = this.Graph.pxWidth/this.ExperimentOptions.DisplayOptions.ScreenWidth * this.ExperimentOptions.DisplayOptions.ScreenDistance * tand(this.ExperimentOptions.TargetSize);
                            fixRect = [0 0 targetPix/2 targetPix/2];
                            targetLocHPix(1,1) = 0;
                            targetLocVPix(1,1) = 0;
                            for target_i =2:1+length(r)
                                targetLocHPix(1,target_i) = this.Graph.pxWidth/this.ExperimentOptions.DisplayOptions.ScreenWidth * this.ExperimentOptions.DisplayOptions.ScreenDistance * tand(dotsLoc(target_i,1));
                                targetLocVPix(1,target_i) = this.Graph.pxWidth/this.ExperimentOptions.DisplayOptions.ScreenWidth * this.ExperimentOptions.DisplayOptions.ScreenDistance * tand(dotsLoc(target_i,2));

%                                 fixRect(target_i) = CenterRectOnPointd( fixRect, mx+targetLocHPix/2, my+targetLocVPix/2 );
%                                 Screen('FillOval', graph.window, this.fixColor, fixRect);
                            end
                            targetLoc = [mx + targetLocHPix./2; my + targetLocVPix./2];
                            %                         fixRect2 = CenterRectOnPointd( fixRect./2, mx+targetHPix/2, my+targetYPix/2 );
                            %                         Screen('FillOval', graph.window, [250,250,250], fixRect2);
                            Screen('DrawDots', graph.window,targetLoc,targetPix/2,this.fixColor,[0,0],1 )
                            Screen('DrawingFinished', graph.window); % Tell PTB that no further drawing commands will follow before Screen('Flip')
                        end
                        this.Graph.Flip();
                    end








                else

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

                            %                         fixRect2 = CenterRectOnPointd( fixRect./2, mx+targetHPix/2, my+targetYPix/2 );
                            %                         Screen('FillOval', graph.window, [250,250,250], fixRect2);

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




