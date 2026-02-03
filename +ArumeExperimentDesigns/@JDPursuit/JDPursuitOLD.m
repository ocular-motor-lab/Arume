classdef JDPursuit < ArumeExperimentDesigns.EyeTracking
    % DEMO experiment for Arume
    %
    %   1. Copy paste the folder @Demo within +ArumeExperimentDesigns.
    %   2. Rename the folder with the name of the new experiment but keep that @ at the begining!
    %   3. Rename also the file inside to match the name of the folder (without the @ this time).
    %   4. Then change the name of the class inside the folder.
    %
    properties
        % Experiment time properties
        stimOffPursue_seconds = 1.;%1.5;
        stimOnPursue_seconds = 2;

        % Stimulus properties
        stimColor = [0 0 0]; % [255 0 0]; %red
        stimRect_size = [0 0 45 112]; %pixels
        fix_size = 20; %pixels
        LRoffset = 100; %pixels
        allStimShift = 250; %pixels

        % Setting up cow background
        textOrderforTrials = [1 Shuffle([2:19]) 20 1 Shuffle([2:19]) 20 1 Shuffle([2:19]) 20 1 Shuffle([2:19]) 20 1 Shuffle([2:19]) 20 randperm(10) 20 1 1 Shuffle([2:19]) 20 1 Shuffle([2:19]) 20 1 Shuffle([2:19]) 20 1 Shuffle([2:19]) 20 1 Shuffle([2:19]) 20 randperm(10) 20 1 1 Shuffle([2:19]) 20 1 Shuffle([2:19]) 20 1 Shuffle([2:19]) 20 1 Shuffle([2:19]) 20 1 Shuffle([2:19]) 20 randperm(10) 20 1 1 Shuffle([2:19]) 20 1 Shuffle([2:19]) 20 1 Shuffle([2:19]) 20 1 Shuffle([2:19]) 20 1 Shuffle([2:19]) 20 randperm(10) 20 1 1 Shuffle([2:19]) 20 1 Shuffle([2:19]) 20 1 Shuffle([2:19]) 20 1 Shuffle([2:19]) 20]; %made it this way so that the same ones is played in a row

        % Central white circle with black ring
        blackRing_half_length = 120; %increase to make ring bigger, also change 'insideCircle_half_length'
        insideCircle_half_length = ceil(180 / 3); % increase to make white circle bigger
        
        %bools: background
        turnOnRing = 1; %change stimColor
        blankwhite = 0; % 1--ganzfeld, 0--cow

        %tbools: control experiments
        control_Fix_Left = 0;
        control_Fix_noMove = 0;
        
    end

    % ---------------------------------------------------------------------
    % Experiment design methods
    % ---------------------------------------------------------------------
    methods ( Access = protected )

        function dlg = GetOptionsDialog( this, importing )
            dlg = GetOptionsDialog@ArumeExperimentDesigns.EyeTracking(this, importing);

            %% ADD new options
            dlg.Correct_Speed_TEST_deg = { 10 '* (deg/s)' [1 300]};
            dlg.Number_of_Speeds = {7 '* (N (must be odd!))' [1 100] };
            dlg.Rate_lag_or_advance_deg = { 1.2 '* (deg/s)' [1 300]};

            dlg.SpeedStep_deg = { 0.07 '* (deg/s)' [0 300]};

            dlg.NumberOfRepetitions = {2 '* (N)' [1 100] };

            dlg.ScreenWidth_pixels = { 1920 '* (pixels)' [1 5000] };
            dlg.ScreenHeight_pixels = { 1080 '* (pixels)' [1 5000] };

            %% CHANGE DEFAULTS values for existing options

            dlg.UseEyeTracker = 0;
            dlg.Debug.DisplayVariableSelection = 'TrialNumber TrialResult Speed_TEST_pix Speed_middle_pix Speed_ref_pix Response MeanElapsedToFlip StdElapsedToFlip MeanFlipToFlip MaxFlipToFlip NDrops'; % which variables to display every trial in the command line separated by spaces

            dlg.DisplayOptions.ScreenWidth = { 54.5 '* (cm)' [1 3000] };
            dlg.DisplayOptions.ScreenHeight = { 30.5 '* (cm)' [1 3000] };
            dlg.DisplayOptions.ScreenDistance = { 66 '* (cm)' [1 3000] };

            dlg.DisplayOptions.SelectedScreen = 1;

            dlg.HitKeyBeforeTrial = 0;
            dlg.TrialDuration = 10;
            dlg.TrialsBeforeBreak = 150;
            dlg.TrialAbortAction = 'Repeat';
        end

        function trialTable = SetUpTrialTable(this)

            %computing ppd
            theta_W_deg = 2*atand(this.ExperimentOptions.DisplayOptions.ScreenWidth /2 /this.ExperimentOptions.DisplayOptions.ScreenDistance);
            theta_H_deg = 2*atand(this.ExperimentOptions.DisplayOptions.ScreenHeight /2 /this.ExperimentOptions.DisplayOptions.ScreenDistance);

            ppd_W = this.ExperimentOptions.ScreenWidth_pixels/theta_W_deg;
            ppd_H = this.ExperimentOptions.ScreenHeight_pixels/theta_H_deg;

            ppd = (ppd_W + ppd_H)/2;

            this.ExperimentOptions.Correct_Speed_TEST_pix = round(this.ExperimentOptions.Correct_Speed_TEST_deg*ppd);
            this.ExperimentOptions.Rate_lag_or_advance_pix =  round(this.ExperimentOptions.Rate_lag_or_advance_deg*ppd);
            this.ExperimentOptions.SpeedStep_pix = round(this.ExperimentOptions.SpeedStep_deg*ppd);

            t = ArumeCore.TrialTableBuilder();
            t.AddConditionVariable('ppd',ppd);

            center_idx = (this.ExperimentOptions.Number_of_Speeds + 1) / 2;
            t.AddConditionVariable('Speed_TEST_pix',this.ExperimentOptions.Correct_Speed_TEST_pix + ([1:this.ExperimentOptions.Number_of_Speeds] - center_idx) * this.ExperimentOptions.SpeedStep_pix);
            t.AddConditionVariable('StartingFromLeft',[ 0 1]);

            t.AddConditionVariable('Speed_middle_pix',this.ExperimentOptions.Correct_Speed_TEST_pix + this.ExperimentOptions.Rate_lag_or_advance_pix);
            t.AddConditionVariable('Speed_ref_pix',this.ExperimentOptions.Correct_Speed_TEST_pix + this.ExperimentOptions.Rate_lag_or_advance_pix*2);

            trialTable = t.GenerateTrialTable();

            T = table();
            for r = 1: this.ExperimentOptions.NumberOfRepetitions
                T = [T; trialTable];
            end
            numRows = size(T, 1);
            shuffledIndices = randperm(numRows);
            shuffledT = T(shuffledIndices, :);
            trialTable = shuffledT;
        end



        function [trialResult, thisTrialData] = runTrial( this, thisTrialData )

            try

                Enum = ArumeCore.ExperimentDesign.getEnum();
                graph = this.Graph;
                trialResult = Enum.trialResult.CORRECT;
                this.ExperimentOptions.DisplayOptions.PlaySound = 0; %JD add

                lastFlipTime        = GetSecs;
                secondsRemaining    = this.ExperimentOptions.TrialDuration;
                thisTrialData.TimeStartLoop = lastFlipTime;

                %JD added make background
                % if thisTrialData.TrialNumber == 1
                wholeLumRasterCancelJDset = [255 255 255];

                [tex1, tex2, tex3, tex4, tex5, tex6, tex7, tex8, tex9, tex10, tex11, tex12, tex13, tex14, tex15, tex16, tex17, tex18, tex19, tex20, tex99] = ...
                    makeTextures(wholeLumRasterCancelJDset, graph.pxWidth, graph.pxHeight, graph.window);

                %-- Find the center of the screen
                [mx, my] = RectCenter(graph.wRect);

                newTextNumof20 = this.textOrderforTrials(thisTrialData.TrialNumber);
                tex1oF = eval(sprintf('tex%1.f',newTextNumof20));
                [~] = setCowBkgdwithCross(graph, mx, my, this, tex1oF);
                % else
                %     newTextNumof20 = this.textOrderforTrials(thisTrialData.TrialNumber);
                %     tex1oF = eval(sprintf('tex%1.f',newTextNumof20));
                % end

                %Added so that it always places based on where it was at the previous flip, not seconds, so that dropped frames don't judder
                FlipTimes = zeros( 10000 ,2);
                nframe = 0;

                fliptime1 = this.Graph.Flip(this, thisTrialData, 0);
                fliptime = this.Graph.Flip(this, thisTrialData, 0);

                diffflip = fliptime - fliptime1;

                startLoopTime = fliptime;

                if thisTrialData.StartingFromLeft == 0 
                    thisTrialData.Speed_middle_pix = thisTrialData.Speed_middle_pix*-1;
                    thisTrialData.Speed_ref_pix = thisTrialData.Speed_ref_pix*-1;
                    thisTrialData.Speed_TEST_pix = thisTrialData.Speed_TEST_pix*-1;
                end

                while secondsRemaining > 0
                    nframe = nframe+1;

                   % secondsElapsed      = GetSecs - thisTrialData.TimeStartLoop;
                   secondsElapsed      = fliptime + diffflip - startLoopTime;
                   secondsRemaining    = this.ExperimentOptions.TrialDuration - secondsElapsed;

                    FlipTimes(nframe,1) = secondsElapsed;
    
                    % -----------------------------------------------------------------
                    % --- Drawing of stimulus -----------------------------------------
                    % -----------------------------------------------------------------


                    if ( secondsElapsed < this.stimOffPursue_seconds)

                        predicted_travel_distance_pix = ((this.stimOffPursue_seconds+this.stimOnPursue_seconds)*  abs(thisTrialData.Speed_middle_pix));
                        if thisTrialData.StartingFromLeft == 0 %Right
                            this.allStimShift = mx+round(predicted_travel_distance_pix/2);
                            % this.allStimShift = mx+round(predicted_travel_distance_pix); %jd1/29
                        else
                            this.allStimShift = mx-round(predicted_travel_distance_pix/2);
                            % this.allStimShift = mx-round(predicted_travel_distance_pix); %jd1/29
                        end

                        %change background
                        [~] = setCowBkgdwithCross(graph, mx, my, this, tex1oF);

                        if this.control_Fix_Left ==1
                            %LEFT
                            if thisTrialData.StartingFromLeft == 0
                                fixRect = CenterRectOnPointd( this.stimRect_size, secondsElapsed*thisTrialData.Speed_Ref_pix-this.LRoffset +this.allStimShift, my );
                            else
                                fixRect = CenterRectOnPointd( this.stimRect_size, secondsElapsed*thisTrialData.Speed_TEST_pix-this.LRoffset +this.allStimShift, my ); %start left
                            end

                            Screen('FillRect', graph.window,  this.stimColor, fixRect);

                                    % [Adding left--fixation circle in center]----------------------------------------------
                                    fixCirc_size = [0 0 this.fix_size this.fix_size];
                                    fixRect = CenterRectOnPointd( fixCirc_size,   secondsElapsed*thisTrialData.Speed_TEST_pix-this.LRoffset +this.allStimShift, my );
                                    Screen('FillOval', graph.window,  [255 255 255], fixRect);
        
                                    fixCross_size = [0 0 this.fix_size/4 this.fix_size];
                                    fixRect = CenterRectOnPointd( fixCross_size,   secondsElapsed*thisTrialData.Speed_TEST_pix-this.LRoffset +this.allStimShift, my );
                                    Screen('FillOval', graph.window,  [0 0 0], fixRect);
        
                                    fixCross_size = [0 0 this.fix_size this.fix_size/4];
                                    fixRect = CenterRectOnPointd( fixCross_size,   secondsElapsed*thisTrialData.Speed_TEST_pix-this.LRoffset +this.allStimShift, my );
                                    Screen('FillOval', graph.window,  [0 0 0], fixRect);
                                    % ---------------------------------------------- [Adding middle--fixation circle in center]

                            reset_secondsElapsed =secondsElapsed; %exactly amount of time you need to subtract, so that the three are equidistant in the beginning
                            reset_secondsElapsed_middle = secondsElapsed*thisTrialData.Speed_TEST_pix; %how much it traveled during the saccade phase
                        elseif this.control_Fix_noMove == 1
                                    % [Adding left--fixation circle in center]----------------------------------------------
                                    fixCirc_size = [0 0 this.fix_size this.fix_size];
                                    fixRect = CenterRectOnPointd( fixCirc_size,  mx, my);%+ (this.blackRing_half_length-this.insideCircle_half_length)/2+this.insideCircle_half_length );
                                    Screen('FillOval', graph.window,  [255 255 255], fixRect);
        
                                    fixCross_size = [0 0 this.fix_size/4 this.fix_size];
                                    fixRect = CenterRectOnPointd( fixCross_size,  mx, my);%+ (this.blackRing_half_length-this.insideCircle_half_length)/2+this.insideCircle_half_length );
                                    Screen('FillOval', graph.window,  [0 0 0], fixRect);
        
                                    fixCross_size = [0 0 this.fix_size this.fix_size/4];
                                    fixRect = CenterRectOnPointd( fixCross_size, mx, my);%+ (this.blackRing_half_length-this.insideCircle_half_length)/2+this.insideCircle_half_length );
                                    Screen('FillOval', graph.window,  [0 0 0], fixRect);
                                    % ---------------------------------------------- [Adding middle--fixation circle in center]
                        reset_secondsElapsed =secondsElapsed; %exactly amount of time you need to subtract, so that the three are equidistant in the beginning
                        reset_secondsElapsed_middle = secondsElapsed*thisTrialData.Speed_middle_pix; %how much it traveled during the saccade phase
                        
                        predicted_travel_distance_pix = ((this.stimOnPursue_seconds)*  abs(thisTrialData.Speed_middle_pix));
                        predicted_sacdistance_pix = ((this.stimOffPursue_seconds)*  thisTrialData.Speed_middle_pix);
                        if thisTrialData.StartingFromLeft == 0 
                            this.allStimShift = mx+round(predicted_travel_distance_pix/2)-predicted_sacdistance_pix;
                            % this.allStimShift = mx+round(predicted_travel_distance_pix)-predicted_sacdistance_pix -this.fix_size/2;%jd1/29
                        else
                            this.allStimShift = mx-round(predicted_travel_distance_pix/2)-predicted_sacdistance_pix;
                            % this.allStimShift = mx-round(predicted_travel_distance_pix)-predicted_sacdistance_pix +this.fix_size/2;%jd1/29
                        end

                        else %normal experiment
                            %MIDDLE
                            fixRect = CenterRectOnPointd( this.stimRect_size,   secondsElapsed*thisTrialData.Speed_middle_pix +this.allStimShift, my );
                            Screen('FillRect', graph.window,  this.stimColor, fixRect);

                                    % [Adding middle--fixation circle in center]----------------------------------------------
                                    fixCirc_size = [0 0 this.fix_size this.fix_size];
                                    fixRect = CenterRectOnPointd( fixCirc_size,   secondsElapsed*thisTrialData.Speed_middle_pix +this.allStimShift, my );
                                    Screen('FillOval', graph.window,  [255 255 255], fixRect);
        
                                    fixCross_size = [0 0 this.fix_size/4 this.fix_size];
                                    fixRect = CenterRectOnPointd( fixCross_size,   secondsElapsed*thisTrialData.Speed_middle_pix +this.allStimShift, my );
                                    Screen('FillOval', graph.window,  [0 0 0], fixRect);
        
                                    fixCross_size = [0 0 this.fix_size this.fix_size/4];
                                    fixRect = CenterRectOnPointd( fixCross_size,   secondsElapsed*thisTrialData.Speed_middle_pix +this.allStimShift, my );
                                    Screen('FillOval', graph.window,  [0 0 0], fixRect);
                                    % ---------------------------------------------- [Adding middle--fixation circle in center]
                            
                            reset_secondsElapsed =secondsElapsed; %exactly amount of time you need to subtract, so that the three are equidistant in the beginning
                            reset_secondsElapsed_middle = secondsElapsed*thisTrialData.Speed_middle_pix; %how much it traveled during the saccade phase

                        end
                        
                    elseif ( secondsElapsed >= this.stimOffPursue_seconds && secondsElapsed <= (this.stimOffPursue_seconds + this.stimOnPursue_seconds))
                        
                        
                        %change background
                        [~] = setCowBkgdwithCross(graph, mx, my, this, tex1oF);

                        %-- Find the center of the screen
                        [mx, my] = RectCenter(graph.wRect);

                        %LEFT
                        if thisTrialData.StartingFromLeft == 0
                            fixRect = CenterRectOnPointd( this.stimRect_size, (secondsElapsed-reset_secondsElapsed)*  thisTrialData.Speed_ref_pix-this.LRoffset+reset_secondsElapsed_middle +this.allStimShift, my );
                        else
                            fixRect = CenterRectOnPointd( this.stimRect_size, (secondsElapsed-reset_secondsElapsed)*  thisTrialData.Speed_TEST_pix-this.LRoffset+reset_secondsElapsed_middle +this.allStimShift, my );
                        end
                        Screen('FillRect', graph.window,  this.stimColor, fixRect);

                        %MIDDLE
                        fixRect = CenterRectOnPointd( this.stimRect_size, (secondsElapsed-reset_secondsElapsed)*  thisTrialData.Speed_middle_pix+reset_secondsElapsed_middle +this.allStimShift, my );
                        Screen('FillRect', graph.window,  this.stimColor, fixRect);

                                if this.control_Fix_Left ==1
                                    % [Adding left--fixation circle in center]----------------------------------------------
                                    fixCirc_size = [0 0 this.fix_size this.fix_size];
                                    fixRect = CenterRectOnPointd( fixCirc_size, (secondsElapsed-reset_secondsElapsed)*  thisTrialData.Speed_TEST_pix-this.LRoffset+reset_secondsElapsed_middle +this.allStimShift, my );
                                    Screen('FillOval', graph.window,  [255 255 255], fixRect);
        
                                    fixCross_size = [0 0 this.fix_size/4 this.fix_size];
                                    fixRect = CenterRectOnPointd( fixCross_size, (secondsElapsed-reset_secondsElapsed)*  thisTrialData.Speed_TEST_pix-this.LRoffset+reset_secondsElapsed_middle +this.allStimShift, my );
                                    Screen('FillOval', graph.window,  [0 0 0], fixRect);
        
                                    fixCross_size = [0 0 this.fix_size this.fix_size/4];
                                    fixRect = CenterRectOnPointd( fixCross_size, (secondsElapsed-reset_secondsElapsed)*  thisTrialData.Speed_TEST_pix-this.LRoffset+reset_secondsElapsed_middle +this.allStimShift, my );
                                    Screen('FillOval', graph.window,  [0 0 0], fixRect);
                                    % ---------------------------------------------- [Adding middle--fixation circle in center]
                                elseif this.control_Fix_noMove == 1
                                    % [Adding left--fixation circle in center]----------------------------------------------
                                    fixCirc_size = [0 0 this.fix_size this.fix_size];
                                    fixRect = CenterRectOnPointd( fixCirc_size,  mx, my);%+  (this.blackRing_half_length-this.insideCircle_half_length)/2+this.insideCircle_half_length );
                                    Screen('FillOval', graph.window,  [255 255 255], fixRect);
        
                                    fixCross_size = [0 0 this.fix_size/4 this.fix_size];
                                    fixRect = CenterRectOnPointd( fixCross_size,  mx, my);%+ (this.blackRing_half_length-this.insideCircle_half_length)/2+this.insideCircle_half_length );
                                    Screen('FillOval', graph.window,  [0 0 0], fixRect);
        
                                    fixCross_size = [0 0 this.fix_size this.fix_size/4];
                                    fixRect = CenterRectOnPointd( fixCross_size, mx, my);%+  (this.blackRing_half_length-this.insideCircle_half_length)/2+this.insideCircle_half_length );
                                    Screen('FillOval', graph.window,  [0 0 0], fixRect);
                                else
                                    % [Adding middle--fixation circle in center]----------------------------------------------
                                    fixCirc_size = [0 0 this.fix_size this.fix_size];
                                    fixRect = CenterRectOnPointd( fixCirc_size, (secondsElapsed-reset_secondsElapsed)*  thisTrialData.Speed_middle_pix+reset_secondsElapsed_middle +this.allStimShift, my );
                                    Screen('FillOval', graph.window,  [255 255 255], fixRect);
        
                                    fixCross_size = [0 0 this.fix_size/4 this.fix_size];
                                    fixRect = CenterRectOnPointd( fixCross_size, (secondsElapsed-reset_secondsElapsed)*  thisTrialData.Speed_middle_pix+reset_secondsElapsed_middle +this.allStimShift, my );
                                    Screen('FillOval', graph.window,  [0 0 0], fixRect);
        
                                    fixCross_size = [0 0 this.fix_size this.fix_size/4];
                                    fixRect = CenterRectOnPointd( fixCross_size, (secondsElapsed-reset_secondsElapsed)*  thisTrialData.Speed_middle_pix+reset_secondsElapsed_middle +this.allStimShift, my );
                                    Screen('FillOval', graph.window,  [0 0 0], fixRect);
                                    % ---------------------------------------------- [Adding middle--fixation circle in center]
                                end
                        %RIGHT
                        if thisTrialData.StartingFromLeft == 0
                            fixRect = CenterRectOnPointd( this.stimRect_size, (secondsElapsed-reset_secondsElapsed)*  thisTrialData.Speed_TEST_pix+this.LRoffset+reset_secondsElapsed_middle +this.allStimShift, my );
                        else
                            fixRect = CenterRectOnPointd( this.stimRect_size, (secondsElapsed-reset_secondsElapsed)*  thisTrialData.Speed_ref_pix+this.LRoffset+reset_secondsElapsed_middle +this.allStimShift, my );
                        end
                        Screen('FillRect', graph.window,  this.stimColor, fixRect);



                    else
                        
                        %change background
                        [~] = setCowBkgdwithCross(graph, mx, my, this, tex1oF);
                    end


                    % -----------------------------------------------------------------
                    % --- END Drawing of stimulus -------------------------------------
                    % -----------------------------------------------------------------

                    % -----------------------------------------------------------------
                    % -- Flip buffers to refresh screen -------------------------------
                    % -----------------------------------------------------------------
                    fliptime = this.Graph.Flip(this, thisTrialData, secondsRemaining);
                    % -----------------------------------------------------------------
    
                    FlipTimes(nframe,2) = fliptime - startLoopTime;

                    % -----------------------------------------------------------------
                    % --- Collecting responses  ---------------------------------------
                    % -----------------------------------------------------------------

                    if ( secondsElapsed > (this.stimOffPursue_seconds + this.stimOnPursue_seconds) )
                        [keyIsDown, secs, keyCode, deltaSecs] = KbCheck();
                        if ( keyIsDown )
                            keys = find(keyCode);
                            for i=1:length(keys)
                                KbName(keys(i));
                                switch(KbName(keys(i)))
                                    case 'RightArrow'
                                        response = 'R';
                                    case 'LeftArrow'
                                        response = 'L';
                                end
                            end
                        end
                        if (exist('response','var') && ~isempty( response) )
                            thisTrialData.Response = response;
                            thisTrialData.ResponseTime = GetSecs;


                            thisTrialData.MeanElapsedToFlip =  mean((FlipTimes(1:nframe,2)-FlipTimes(1:nframe,1))*1000);
                            thisTrialData.StdElapsedToFlip =    std((FlipTimes(1:nframe,2)-FlipTimes(1:nframe,1))*1000);
                            thisTrialData.MeanFlipToFlip =    mean((FlipTimes(2:nframe,2)-FlipTimes(1:nframe-1,2))*1000);
                            thisTrialData.MaxFlipToFlip =    max((FlipTimes(2:nframe,2)-FlipTimes(1:nframe-1,2))*1000);
                            thisTrialData.NDrops =    sum((FlipTimes(2:nframe,2)-FlipTimes(1:nframe-1,2))*1000 > median((FlipTimes(2:nframe,2)-FlipTimes(1:nframe-1,2))*1000)*1.5);

                            % fprintf('Time from elapsed to flip Mean: %0.1f ms Std: %0.1f ms - Flip to flip mean %0.1f ms max %0.1f ms ndrops= %d \n', ...
                            %     mean((FlipTimes(1:nframe,2)-FlipTimes(1:nframe,1))*1000), ...
                            %     std((FlipTimes(1:nframe,2)-FlipTimes(1:nframe,1))*1000), ...
                            %     mean((FlipTimes(2:nframe,2)-FlipTimes(1:nframe-1,2))*1000), ...
                            %     max((FlipTimes(2:nframe,2)-FlipTimes(1:nframe-1,2))*1000), ...
                            %     sum((FlipTimes(2:nframe,2)-FlipTimes(1:nframe-1,2))*1000 > median((FlipTimes(2:nframe,2)-FlipTimes(1:nframe-1,2))*1000)*1.5));
                            if thisTrialData.TrialNumber == size(this.TrialTable,1)
                                Plot_Psychometric(this)
                            end
                            
                            break;
                        end
                    end
                    % -----------------------------------------------------------------
                    % --- END Collecting responses  -----------------------------------
                    % -----------------------------------------------------------------
                end
            catch ex
                rethrow(ex)
            end

            % Make 20 textures
            function [tex1, tex2, tex3, tex4, tex5, tex6, tex7, tex8, tex9, tex10, tex11, tex12, tex13, tex14, tex15, tex16, tex17, tex18, tex19, tex20, tex99] = ...
                    makeTextures(wholeLumRasterCancelJDset, screen_width, screen_height, w)

                %cd('/Users/Josephine1/Documents/MATLAB/Roorda Lab/Data_and_AnalysisCode/gifMakers_presentations/JDs20noiseBMPS_folder'); % I know this is redundant, but didn't know a better way JD 7/21/23

                folder = fullfile(fileparts(mfilename('fullpath')), 'JDs20noiseBMPS_folder');

                jdnoise1 = imread(fullfile(folder,'JDnoise1.bmp'));
                [jdnoise1] = cropAndMakeBkgLum(jdnoise1,wholeLumRasterCancelJDset, screen_width, screen_height);
                tex1 = Screen('MakeTexture', w, jdnoise1);

                jdnoise2 = imread(fullfile(folder,'JDnoise2.bmp'));
                [jdnoise2] = cropAndMakeBkgLum(jdnoise2,wholeLumRasterCancelJDset, screen_width, screen_height);
                tex2 = Screen('MakeTexture', w, jdnoise2);

                jdnoise3 = imread(fullfile(folder,'JDnoise3.bmp'));
                [jdnoise3] = cropAndMakeBkgLum(jdnoise3,wholeLumRasterCancelJDset, screen_width, screen_height);
                tex3 = Screen('MakeTexture', w, jdnoise3);

                jdnoise4 = imread(fullfile(folder,'JDnoise4.bmp'));
                [jdnoise4] = cropAndMakeBkgLum(jdnoise4,wholeLumRasterCancelJDset, screen_width, screen_height);
                tex4 = Screen('MakeTexture', w, jdnoise4);

                jdnoise5 = imread(fullfile(folder,'JDnoise5.bmp'));
                [jdnoise5] = cropAndMakeBkgLum(jdnoise5,wholeLumRasterCancelJDset, screen_width, screen_height);
                tex5 = Screen('MakeTexture', w, jdnoise5);

                jdnoise6 = imread(fullfile(folder,'JDnoise6.bmp'));
                [jdnoise6] = cropAndMakeBkgLum(jdnoise6,wholeLumRasterCancelJDset, screen_width, screen_height);
                tex6 = Screen('MakeTexture', w, jdnoise6);

                jdnoise7 = imread(fullfile(folder,'JDnoise7.bmp'));
                [jdnoise7] = cropAndMakeBkgLum(jdnoise7,wholeLumRasterCancelJDset, screen_width, screen_height);
                tex7 = Screen('MakeTexture', w, jdnoise7);

                jdnoise8 = imread(fullfile(folder,'JDnoise8.bmp'));
                [jdnoise8] = cropAndMakeBkgLum(jdnoise8,wholeLumRasterCancelJDset, screen_width, screen_height);
                tex8 = Screen('MakeTexture', w, jdnoise8);

                jdnoise9 = imread(fullfile(folder,'JDnoise9.bmp'));
                [jdnoise9] = cropAndMakeBkgLum(jdnoise9,wholeLumRasterCancelJDset, screen_width, screen_height);
                tex9 = Screen('MakeTexture', w, jdnoise9);

                jdnoise10 = imread(fullfile(folder,'JDnoise10.bmp'));
                [jdnoise10] = cropAndMakeBkgLum(jdnoise10,wholeLumRasterCancelJDset, screen_width, screen_height);
                tex10 = Screen('MakeTexture', w, jdnoise10);

                jdnoise11 = imread(fullfile(folder,'JDnoise11.bmp'));
                [jdnoise11] = cropAndMakeBkgLum(jdnoise11,wholeLumRasterCancelJDset, screen_width, screen_height);
                tex11 = Screen('MakeTexture', w, jdnoise11);

                jdnoise12 = imread(fullfile(folder,'JDnoise12.bmp'));
                [jdnoise12] = cropAndMakeBkgLum(jdnoise12,wholeLumRasterCancelJDset, screen_width, screen_height);
                tex12= Screen('MakeTexture', w, jdnoise12);

                jdnoise13 = imread(fullfile(folder,'JDnoise13.bmp'));
                [jdnoise13] = cropAndMakeBkgLum(jdnoise13,wholeLumRasterCancelJDset, screen_width, screen_height);
                tex13 = Screen('MakeTexture', w, jdnoise13);

                jdnoise14 = imread(fullfile(folder,'JDnoise14.bmp'));
                [jdnoise14] = cropAndMakeBkgLum(jdnoise14,wholeLumRasterCancelJDset, screen_width, screen_height);
                tex14 = Screen('MakeTexture', w, jdnoise14);

                jdnoise15 = imread(fullfile(folder,'JDnoise15.bmp'));
                [jdnoise15] = cropAndMakeBkgLum(jdnoise15,wholeLumRasterCancelJDset, screen_width, screen_height);
                tex15 = Screen('MakeTexture', w, jdnoise15);

                jdnoise16 = imread(fullfile(folder,'JDnoise16.bmp'));
                [jdnoise16] = cropAndMakeBkgLum(jdnoise16,wholeLumRasterCancelJDset, screen_width, screen_height);
                tex16 = Screen('MakeTexture', w, jdnoise16);

                jdnoise17 = imread(fullfile(folder,'JDnoise17.bmp'));
                [jdnoise17] = cropAndMakeBkgLum(jdnoise17,wholeLumRasterCancelJDset, screen_width, screen_height);
                tex17 = Screen('MakeTexture', w, jdnoise17);

                jdnoise18 = imread(fullfile(folder,'JDnoise18.bmp'));
                [jdnoise18] = cropAndMakeBkgLum(jdnoise18,wholeLumRasterCancelJDset, screen_width, screen_height);
                tex18 = Screen('MakeTexture', w, jdnoise18);

                jdnoise19 = imread(fullfile(folder,'JDnoise19.bmp'));
                [jdnoise19] = cropAndMakeBkgLum(jdnoise19,wholeLumRasterCancelJDset, screen_width, screen_height);
                tex19 = Screen('MakeTexture', w, jdnoise19);

                jdnoise20 = imread(fullfile(folder,'JDnoise20.bmp'));
                [jdnoise20] = cropAndMakeBkgLum(jdnoise20,wholeLumRasterCancelJDset, screen_width, screen_height);
                tex20 = Screen('MakeTexture', w, jdnoise20);

                whiteblank99 = imread(fullfile(folder,'JDnoise20.bmp'));
                whiteblank99(:,:,:)=255;
                [whiteblank99] = cropAndMakeBkgLum(whiteblank99,wholeLumRasterCancelJDset, screen_width, screen_height);
                tex99 = Screen('MakeTexture', w, whiteblank99);

            end
            % set the raster-cancel-background lum for each texture
            function [noisejpgread_crop] = cropAndMakeBkgLum(noisejpgread,lumJDset, screen_width, screen_height)
                xmin = 260; xmax = 1810;
                ymin = 98;  ymax = 1642;

                width = xmax - xmin;
                height = ymax - ymin;
                rect = [xmin, ymin, width, height];
                noisejpgread_crop = imcrop(noisejpgread, rect);
                % noisejpgread_crop = imcrop(noisejpgread, [240 141  screen_width screen_height]); %this is for jdnoise7323.bmp SCREEN DIMENSIONS, if want crop to best with white on outsidesS 74 not 113
                % % % noisejpgread_crop = imcrop(noisejpgread, [540 241  screen_height screen_height]);
                for frame = 1: size(noisejpgread_crop,1)
                    for c = 1: size(noisejpgread_crop,2)
                        if noisejpgread_crop(frame,c,1) == 255
                            noisejpgread_crop(frame,c,1) = lumJDset(1); %this number is for 4.2e-01 luminance
                            noisejpgread_crop(frame,c,2) = lumJDset(2); %this number is for 4.2e-01 luminance
                            noisejpgread_crop(frame,c,3) = lumJDset(3);%this number is for 4.2e-01 luminance
                        end
                    end
                end

            end

            function [mx] = setCowBkgdwithCross(graph, mx, my, this, tex1oF)

                black = BlackIndex(graph.window);
                white = WhiteIndex(graph.window);

                if this.blankwhite == 1
                    Screen('FillRect', graph.window, white);
                else
                    Screen('FillRect', graph.window, black);
                    destinationRectangle =  []; % so that it fills the whole screen
                    Screen('DrawTexture', graph.window, tex1oF, [],destinationRectangle);

                    if this.turnOnRing == 1
                        blackRing_rgb = [0 0 0];
                        % Drawing the black circle first so that it becomes ring
                        Screen('FillRect',graph.window,blackRing_rgb,... %black oval
                            [1,...
                            my - this.blackRing_half_length,...
                            mx + mx,...
                            my + this.blackRing_half_length]);

                        % Now drawing inside white circle
                        Screen('FillRect', graph.window ...
                            , white,... %these numbers are for 4.2e-01 luminance
                            [1,...
                            my - this.insideCircle_half_length,...
                            mx + mx,...
                            my + this.insideCircle_half_length]); %2deg box
                    end
                end
            end

        end %runTrial end


    end

    % ---------------------------------------------------------------------
    % Plot methods
    % ---------------------------------------------------------------------
    methods ( Access = public )

        function Plot_Psychometric(this)
            folder = fullfile(fileparts(mfilename('fullpath')), 'Palamedes');
            addpath(folder);
            
            trialDataTable = this.Session.currentRun.pastTrialTable;

            %convertings pixels to degrees
            trialDataTable.Speed_TEST_deg = trialDataTable.Speed_TEST_pix./trialDataTable.ppd;
            trialDataTable.Speed_middle_deg = trialDataTable.Speed_middle_pix./trialDataTable.ppd;
            trialDataTable.Speed_ref_deg = trialDataTable.Speed_ref_pix./trialDataTable.ppd;

            %getting rate of gap increasing
            trialDataTable.Rate_test = abs(trialDataTable.Speed_middle_deg - trialDataTable.Speed_TEST_deg);
            trialDataTable.Rate_ref = abs(trialDataTable.Speed_ref_deg - trialDataTable.Speed_middle_deg);

            trialDataTable.ChoseRef(trialDataTable.Speed_ref_pix > 0 & trialDataTable.Response == 'R') = 1;
            trialDataTable.ChoseRef(trialDataTable.Speed_ref_pix < 0 & trialDataTable.Response == 'L') = 1;

            [uniqueVals, ~, idx] = unique(abs(trialDataTable.Speed_TEST_deg));
            % marker_sizes = accumarray(idx,1);

            StimLevels = nan(1,size(uniqueVals,1)); NumPos = nan(1,size(uniqueVals,1)); OutOfNum = nan(1,size(uniqueVals,1));
            for i = 1:size(uniqueVals,1)
                curData = trialDataTable(abs(trialDataTable.Speed_TEST_deg) ==uniqueVals(i), :);

                StimLevels(i) = uniqueVals(i);
                NumPos(i) = sum(curData.ChoseRef);
                OutOfNum(i) = length(curData.ChoseRef);
            end

            PF = @PAL_Logistic;

            paramsFree = [1 1 0 0];
            searchGrid.alpha = mean(StimLevels);
            searchGrid.beta = 3.5;
            searchGrid.gamma = 0.01;
            searchGrid.lambda = 0.01;

            %Perform fit
            [paramsValues] = PAL_PFML_Fit(StimLevels,NumPos,OutOfNum,searchGrid,paramsFree,PF);

            %plotting fit
            ProportionCorrect = NumPos./OutOfNum;
            xEval = [min(StimLevels):max(abs(StimLevels))./1000:max(StimLevels)];
            yEval = PF(paramsValues,xEval);

            PSE = PF(paramsValues, 0.5, 'inv');

            figure('Color','w'); hold on; grid on;
            plot(xEval,yEval,'-','Color',[0 0.7 0], 'LineWidth',4);
            plot(StimLevels,ProportionCorrect, 'ko', 'MarkerFaceColor','k', 'MarkerSize',10);

            B = 400;
            searchGrid.alpha = paramsValues(1) - 0.1:0.001:paramsValues(1)+0.1;
            searchGrid.beta = paramsValues(2);

            [SD] = PAL_PFML_BootstrapParametric(StimLevels,OutOfNum,paramsValues,paramsFree,B,PF,'searchGrid',searchGrid);

            errorbar(PSE,0.5, SD(1), 'horizontal', 'ro','MarkerFaceColor', 'r', 'LineWidth',1.5, 'MarkerSize',15)
            text(PSE+0.5, 0.5, sprintf('PSE = %1.2f \\pm %1.2f ms', PSE, SD(1)));

            %formating graph
            xline(abs(trialDataTable.Speed_middle_deg(1))-trialDataTable.Rate_ref, '--','Color',[0.5 0.5 0.5], 'LineWidth',2);
            ylim([0 1]); yticks([0:0.1:1]);
            ylabel('Reponse == Reference')
            xlabel('Test stimulus speed (deg/s)')
            set(gca,'fontsize', 18)

            % saveas(figure(1), [folder, '/', this.Session.subjectCode, '_', this.Session.sessionCode '_PF.jpg'])
            
            % figure
            % trialDataTable = this.Session.trialDataTable;
            % plot(trialDataTable.Speed_left, double(trialDataTable.Response=='R')+randn(size(trialDataTable.Response))/20,'o');
            % xlabel('Speed')
            % ylabel('Reponse == Right')
        end

    end


end %classdef end