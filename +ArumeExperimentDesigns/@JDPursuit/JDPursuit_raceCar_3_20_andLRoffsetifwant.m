classdef JDPursuit < ArumeExperimentDesigns.EyeTracking
    % Don't forget to turn on sound

    properties
        % Experiment time properties
        stimOffFixate_seconds = 0.5; %seconds
        stimOffPursue_seconds = 1; %seconds
        stimOnPursue_seconds  = 1; %seconds
        stimFadeOut_seconds   = 0.5; %seconds

        % Stimulus properties
        stimColor = [0 0 0]; %[255 255 255]; % [255 0 0]; %red
        stimRect_size =  [0 0 25 30]; % [0 0 45 112]; %pixels %120
        fix_size = 20; %pixels

        allStimShift = 0; %this value here doesn't do anything, it gets updated in the code depending on other parameters

        % Setting up cow background
        textOrderforTrials = [1 Shuffle([2:19]) 20 1 Shuffle([2:19]) 20 1 Shuffle([2:19]) 20 1 Shuffle([2:19]) 20 1 Shuffle([2:19]) 20 randperm(10) 20 1 1 Shuffle([2:19]) 20 1 Shuffle([2:19]) 20 1 Shuffle([2:19]) 20 1 Shuffle([2:19]) 20 1 Shuffle([2:19]) 20 randperm(10) 20 1 1 Shuffle([2:19]) 20 1 Shuffle([2:19]) 20 1 Shuffle([2:19]) 20 1 Shuffle([2:19]) 20 1 Shuffle([2:19]) 20 randperm(10) 20 1 1 Shuffle([2:19]) 20 1 Shuffle([2:19]) 20 1 Shuffle([2:19]) 20 1 Shuffle([2:19]) 20 1 Shuffle([2:19]) 20 randperm(10) 20 1 1 Shuffle([2:19]) 20 1 Shuffle([2:19]) 20 1 Shuffle([2:19]) 20 1 Shuffle([2:19]) 20]; %made it this way so that the same ones is played in a row

        % Central white circle with black ring
        blackRing_half_length = 30; %90; %increase to make ring bigger, also change 'insideCircle_half_length' %120
        insideCircle_half_length = 15;% 40; %ceil(180 / 3); % increase to make white circle bigger
        
        %bools: background
        turnOnRing = 1; % turn on the horizontal bars surrounding stim
        blankblack = 0; % 1--ganzfeld, 0--cow

    end

    % ---------------------------------------------------------------------
    % Experiment design methods
    % ---------------------------------------------------------------------
    methods ( Access = protected )

        function dlg = GetOptionsDialog( this, importing )
            dlg = GetOptionsDialog@ArumeExperimentDesigns.EyeTracking(this, importing);

            %% ADD new options
            dlg.Correct_Speed_TEST_deg = { 4.8 '* (deg/s)' [0 300]}; %4.8 
            dlg.Number_of_Speeds = {9 '* (N (must be odd!))' [1 100] };
            dlg.Rate_lag_or_advance_deg = { 1.2 '* (deg/s)' [0 300]};

            dlg.SpeedStep_deg = { 0.06 '* (deg/s)' [0 300]}; %0.07 %trainingTrials = 2

            dlg.NumberOfRepetitions = {4 '* (N)' [1 100] };

            dlg.LRoffset = {100 '* (100 if grow, 202 if shrink)' [0 500] };
            % dlg.LRoffset_pix_diff = {20 '* (different lengths)' [0 500] };
            dlg.old_control_Fix_Left = 0;
            dlg.control_Fix_noMove = 0;
            dlg.endAtCross = 0;
            dlg.randEye = 0;
            dlg.NOPursueExp = 0;

            dlg.ScreenWidth_pixels = { 1920 '* (pixels)' [1 5000] };
            dlg.ScreenHeight_pixels = { 1080 '* (pixels)' [1 5000] };
            dlg.GrowingOutwards = 1;
            
            % dlg.blankblack = 1;

            %% CHANGE DEFAULTS values for existing options

            dlg.UseEyeTracker = 1;
            dlg.Debug.DisplayVariableSelection = 'TrialNumber TrialResult Speed_TEST_pix Speed_middle_pix Speed_ref_pix Response MeanElapsedToFlip StdElapsedToFlip MeanFlipToFlip MaxFlipToFlip NDrops'; % which variables to display every trial in the command line separated by spaces

            dlg.DisplayOptions.ScreenWidth = { 54.5 '* (cm)' [1 3000] };
            dlg.DisplayOptions.ScreenHeight = { 30.5 '* (cm)' [1 3000] };
            dlg.DisplayOptions.ScreenDistance = { 66 '* (cm)' [1 3000] };
            dlg.DisplayOptions.PlaySound = 1;

            dlg.DisplayOptions.SelectedScreen = 1;

            dlg.HitKeyBeforeTrial = 0;
            dlg.TrialDuration = 10;
            dlg.TrialsBeforeBreak = 150;
            dlg.TrialAbortAction = 'Repeat';
        end

        function trialTable = SetUpTrialTable(this)

            %computing ppd, only horizontal (W)
            theta_W_deg = 2*atand(this.ExperimentOptions.DisplayOptions.ScreenWidth /2 /this.ExperimentOptions.DisplayOptions.ScreenDistance);
            theta_H_deg = 2*atand(this.ExperimentOptions.DisplayOptions.ScreenHeight /2 /this.ExperimentOptions.DisplayOptions.ScreenDistance);

            ppd_W = this.ExperimentOptions.ScreenWidth_pixels/theta_W_deg;
            ppd_H = this.ExperimentOptions.ScreenHeight_pixels/theta_H_deg;

            this.ExperimentOptions.Correct_Speed_TEST_pix = this.ExperimentOptions.Correct_Speed_TEST_deg*ppd_W;
            this.ExperimentOptions.Rate_lag_or_advance_pix =  this.ExperimentOptions.Rate_lag_or_advance_deg*ppd_W;
            this.ExperimentOptions.SpeedStep_pix = this.ExperimentOptions.SpeedStep_deg*ppd_W;

            t = ArumeCore.TrialTableBuilder();
            t.AddConditionVariable('ppd_x',ppd_W);
            t.AddConditionVariable('ppd_y',ppd_H);
            t.AddConditionVariable('LRoffset',this.ExperimentOptions.LRoffset);
            % t.AddConditionVariable('LRoffset_pix_diff', this.ExperimentOptions.LRoffset_pix_diff)
            t.AddConditionVariable('control_Fix_Left', this.ExperimentOptions.old_control_Fix_Left);
            t.AddConditionVariable('control_Fix_noMove', this.ExperimentOptions.control_Fix_noMove);
            t.AddConditionVariable('endAtCross', this.ExperimentOptions.endAtCross);
            t.AddConditionVariable('randEye', this.ExperimentOptions.randEye);
            t.AddConditionVariable('NOPursueExp', this.ExperimentOptions.NOPursueExp);

            center_idx = (this.ExperimentOptions.Number_of_Speeds + 1) / 2;

            %if any are negative, then it will move the wrong way and tests a dif question
            while any(this.ExperimentOptions.Correct_Speed_TEST_pix + ([1:this.ExperimentOptions.Number_of_Speeds] - center_idx) * this.ExperimentOptions.SpeedStep_pix<0)
                this.ExperimentOptions.SpeedStep_pix = this.ExperimentOptions.SpeedStep_pix-1;
            end
            %if any are greater than the ref, then it no longer tests what we're interested in
            while any(this.ExperimentOptions.Correct_Speed_TEST_pix + ([1:this.ExperimentOptions.Number_of_Speeds] - center_idx) * this.ExperimentOptions.SpeedStep_pix>= (this.ExperimentOptions.Correct_Speed_TEST_pix + this.ExperimentOptions.Rate_lag_or_advance_pix*2))
                this.ExperimentOptions.SpeedStep_pix = this.ExperimentOptions.SpeedStep_pix-1;
            end
            this.ExperimentOptions.SpeedStep_deg = this.ExperimentOptions.SpeedStep_pix/ppd_W; %updates if any of the while loops were used

            if this.ExperimentOptions.NOPursueExp == 0
                t.AddConditionVariable('Speed_TEST_pix',this.ExperimentOptions.Correct_Speed_TEST_pix + ([1:this.ExperimentOptions.Number_of_Speeds] - center_idx) * this.ExperimentOptions.SpeedStep_pix);
                t.AddConditionVariable('Speed_middle_pix',this.ExperimentOptions.Correct_Speed_TEST_pix + this.ExperimentOptions.Rate_lag_or_advance_pix);
                t.AddConditionVariable('Speed_ref_pix',this.ExperimentOptions.Correct_Speed_TEST_pix + this.ExperimentOptions.Rate_lag_or_advance_pix*2);
            else
                t.AddConditionVariable('Speed_TEST_pix',-(this.ExperimentOptions.Correct_Speed_TEST_pix + ([1:this.ExperimentOptions.Number_of_Speeds] - center_idx) * this.ExperimentOptions.SpeedStep_pix));
                t.AddConditionVariable('Speed_middle_pix',0);
                t.AddConditionVariable('Speed_ref_pix',this.ExperimentOptions.Correct_Speed_TEST_pix);
            end

            % t.AddConditionVariable('LRoffset',this.ExperimentOptions.LRoffset + ([1:5] - 3) * this.ExperimentOptions.LRoffset_pix_stepsDif);
            t.AddConditionVariable('StartingFromLeft',[ 0 1]); %            t.AddConditionVariable('StartingFromLeft',[ 0 1]);
            t.AddConditionVariable('GrowingOutwards',this.ExperimentOptions.GrowingOutwards); %t.AddConditionVariable('GrowingOutwards',[ 0 1]);
            t.AddConditionVariable('RefIsLead',[ 1]);

            trialTable = t.GenerateTrialTable();

            %making another table to store when we show the front one as TEST----------------------
            t2 = ArumeCore.TrialTableBuilder();
            t2.AddConditionVariable('ppd_x',ppd_W);
            t2.AddConditionVariable('ppd_y',ppd_H);
            t2.AddConditionVariable('LRoffset',this.ExperimentOptions.LRoffset);
            % t2.AddConditionVariable('LRoffset_pix_diff', this.ExperimentOptions.LRoffset_pix_diff)
            t2.AddConditionVariable('control_Fix_Left', this.ExperimentOptions.old_control_Fix_Left);
            t2.AddConditionVariable('control_Fix_noMove', this.ExperimentOptions.control_Fix_noMove);
            t2.AddConditionVariable('endAtCross', this.ExperimentOptions.endAtCross);
            t2.AddConditionVariable('randEye', this.ExperimentOptions.randEye);
            t2.AddConditionVariable('NOPursueExp', this.ExperimentOptions.NOPursueExp);
            
            if this.ExperimentOptions.NOPursueExp == 0
                t2.AddConditionVariable('Speed_TEST_pix',this.ExperimentOptions.Correct_Speed_TEST_pix);
                t2.AddConditionVariable('Speed_middle_pix',this.ExperimentOptions.Correct_Speed_TEST_pix + this.ExperimentOptions.Rate_lag_or_advance_pix);
                t2.AddConditionVariable('Speed_ref_pix',(this.ExperimentOptions.Correct_Speed_TEST_pix + this.ExperimentOptions.Rate_lag_or_advance_pix*2)+ ([1:this.ExperimentOptions.Number_of_Speeds] - center_idx) * this.ExperimentOptions.SpeedStep_pix);
            else
                t2.AddConditionVariable('Speed_TEST_pix',-this.ExperimentOptions.Correct_Speed_TEST_pix);
                t2.AddConditionVariable('Speed_middle_pix',0);
                t2.AddConditionVariable('Speed_ref_pix',this.ExperimentOptions.Correct_Speed_TEST_pix + ([1:this.ExperimentOptions.Number_of_Speeds] - center_idx) * this.ExperimentOptions.SpeedStep_pix);
            end

            % t2.AddConditionVariable('LRoffset',this.ExperimentOptions.LRoffset + ([1:5] - 3) * this.ExperimentOptions.LRoffset_pix_stepsDif);
            t2.AddConditionVariable('StartingFromLeft',[ 0 1]); %            t.AddConditionVariable('StartingFromLeft',[ 0 1]);
            t2.AddConditionVariable('GrowingOutwards',this.ExperimentOptions.GrowingOutwards); %t.AddConditionVariable('GrowingOutwards',[ 0 1]);
            t2.AddConditionVariable('RefIsLead',[ 0]);
            trialTable2 = t2.GenerateTrialTable();
            %----------------------making another table to store when we show the front one as TEST

            T = table();
            for r = 1: this.ExperimentOptions.NumberOfRepetitions
                T = [T; trialTable; trialTable2];
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
                % this.ExperimentOptions.DisplayOptions.PlaySound = 1; %JD add

                if this.blankblack == 1
                    white = BlackIndex(graph.window);
                    black = WhiteIndex(graph.window);
                else
                    black = BlackIndex(graph.window);
                    white = WhiteIndex(graph.window);
                end

                lastFlipTime        = GetSecs;
                secondsRemaining    = this.ExperimentOptions.TrialDuration;

                %JD added make background
                wholeLumRasterCancelJDset = [255 255 255];

                [tex1, tex2, tex3, tex4, tex5, tex6, tex7, tex8, tex9, tex10] = ...
                    makeTextures(wholeLumRasterCancelJDset, graph.pxWidth, graph.pxHeight, graph.window);

                %-- Find the center of the screen
                [mx, my] = RectCenter(graph.wRect);

                newTextNumof20 = this.textOrderforTrials(thisTrialData.TrialNumber);
                if newTextNumof20 > 10
                    newTextNumof20 = newTextNumof20-10;
                end

                tex1oF = eval(sprintf('tex%1.f',newTextNumof20));
                [~] = setCowBkgdwithCross(graph, mx, my, this, tex1oF, white);

                %setting up LRoffset, to be random
                Lag_LRoffset = thisTrialData.LRoffset; %randi(round([thisTrialData.LRoffset-thisTrialData.LRoffset_pix_diff/4 , thisTrialData.LRoffset+3*thisTrialData.LRoffset_pix_diff/4]));

                % % Valid ranges for B
                % lower_range = round(thisTrialData.LRoffset-thisTrialData.LRoffset_pix_diff/4) : (Lag_LRoffset - 15);
                % upper_range = (Lag_LRoffset + 15) : round(thisTrialData.LRoffset+3*thisTrialData.LRoffset_pix_diff/4);
                % 
                % valid_values = [lower_range, upper_range];
                % 
                % if isempty(lower_range) && isempty(upper_range)
                %     Lead_LRoffset = randi(round([thisTrialData.LRoffset-thisTrialData.LRoffset_pix_diff/4 , thisTrialData.LRoffset+3*thisTrialData.LRoffset_pix_diff/4]));
                % else
                %     Lead_LRoffset = valid_values(randi(length(valid_values)));
                % end
                Lead_LRoffset =Lag_LRoffset;
                
                thisTrialData.Lag_LRoffset = Lag_LRoffset;
                thisTrialData.Lead_LRoffset = Lead_LRoffset;

                %Added so that it always places based on where it was at the previous flip, not seconds, so that dropped frames don't judder
                FlipTimes = zeros( 10000 ,2);
                nframe = 0;

                fliptime1 = this.Graph.Flip(this, thisTrialData, 0);
                fliptime = this.Graph.Flip(this, thisTrialData, 0);

                diffflip = fliptime - fliptime1;


                if thisTrialData.StartingFromLeft == 0 
                    thisTrialData.Speed_middle_pix = thisTrialData.Speed_middle_pix*-1;
                    thisTrialData.Speed_ref_pix = thisTrialData.Speed_ref_pix*-1;
                    thisTrialData.Speed_TEST_pix = thisTrialData.Speed_TEST_pix*-1;
                end

                % if ( ~isempty(this.eyeTracker) )
                %     thisTrialData.EyeTrackerFrameStartLoop = this.eyeTracker.RecordEvent("ACTUAL LOOP START");
                % end
                startLoopTime = fliptime;
                thisTrialData.TimeStartLoop = startLoopTime;
                while secondsRemaining > 0
                    nframe = nframe+1;

                   % secondsElapsed      = GetSecs - thisTrialData.TimeStartLoop;
                   secondsElapsed      = fliptime + diffflip - startLoopTime;
                   secondsRemaining    = this.ExperimentOptions.TrialDuration - secondsElapsed;

                    FlipTimes(nframe,1) = secondsElapsed;
    
                    % -----------------------------------------------------------------
                    % --- Drawing of stimulus -----------------------------------------
                    % -----------------------------------------------------------------

                    if ( secondsElapsed < this.stimOffFixate_seconds)
                        
                        %change background
                        [~] = setCowBkgdwithCross(graph, mx, my, this, tex1oF, white);
                           
                        if thisTrialData.control_Fix_noMove == 0 && thisTrialData.NOPursueExp == 0
                            predicted_travel_distance_pix = (this.stimOnPursue_seconds*  abs(thisTrialData.Speed_middle_pix));
                            pursueDist = this.stimOffPursue_seconds*  abs(thisTrialData.Speed_middle_pix);

                            if thisTrialData.StartingFromLeft == 0 %Right
                                this.allStimShift = mx+predicted_travel_distance_pix/2+pursueDist;
                            else
                                this.allStimShift = mx-predicted_travel_distance_pix/2-pursueDist;
                            end

                            %MIDDLE
                            fixRect = CenterRectOnPointd( this.stimRect_size,  this.allStimShift, my );
                            Screen('FillRect', graph.window,  this.stimColor, fixRect);

                                    % [Adding middle--fixation circle in center]----------------------------------------------
                                    fixCirc_size = [0 0 this.fix_size this.fix_size];
                                    fixRect = CenterRectOnPointd( fixCirc_size,  this.allStimShift, my );
                                    Screen('FillOval', graph.window,  black, fixRect);
        
                                    fixCross_size = [0 0 this.fix_size/4 this.fix_size];
                                    fixRect = CenterRectOnPointd( fixCross_size,   this.allStimShift, my );
                                    Screen('FillOval', graph.window,  white, fixRect);
        
                                    fixCross_size = [0 0 this.fix_size this.fix_size/4];
                                    fixRect = CenterRectOnPointd( fixCross_size,   this.allStimShift, my );
                                    Screen('FillOval', graph.window,  white, fixRect);
                                    % ---------------------------------------------- [Adding middle--fixation circle in center]
                        elseif thisTrialData.control_Fix_noMove == 1 || thisTrialData.NOPursueExp == 1
                                   if thisTrialData.randEye == 0
                                       %[Adding left--fixation circle in center]----------------------------------------------
                                       fixCirc_size = [0 0 this.fix_size this.fix_size];
                                       fixRect = CenterRectOnPointd( fixCirc_size,  mx, my);%+ (this.blackRing_half_length-this.insideCircle_half_length)/2+this.insideCircle_half_length );
                                       Screen('FillOval', graph.window,  black, fixRect);

                                       fixCross_size = [0 0 this.fix_size/4 this.fix_size];
                                       fixRect = CenterRectOnPointd( fixCross_size,  mx, my);%+ (this.blackRing_half_length-this.insideCircle_half_length)/2+this.insideCircle_half_length );
                                       Screen('FillOval', graph.window,  white, fixRect);

                                       fixCross_size = [0 0 this.fix_size this.fix_size/4];
                                       fixRect = CenterRectOnPointd( fixCross_size, mx, my);%+ (this.blackRing_half_length-this.insideCircle_half_length)/2+this.insideCircle_half_length );
                                       Screen('FillOval', graph.window,  white, fixRect);
                                       % ---------------------------------------------- [Adding middle--fixation circle in center]
                                   end
                        end

                        first_reset_secondsElapsed =secondsElapsed; %exactly amount of time you need to subtract, so that the three are equidistant in the beginning
                    elseif (secondsElapsed >= this.stimOffFixate_seconds && secondsElapsed < (this.stimOffFixate_seconds+this.stimOffPursue_seconds)) %moving but no stimuli

                        % predicted_travel_distance_pix = ((this.stimOffPursue_seconds+this.stimOnPursue_seconds)*  abs(thisTrialData.Speed_middle_pix)); %whole stim duration is accted%
                        predicted_travel_distance_pix = (this.stimOnPursue_seconds*  abs(thisTrialData.Speed_middle_pix));
                        pursueDist = this.stimOffPursue_seconds*  abs(thisTrialData.Speed_middle_pix);

                        if thisTrialData.StartingFromLeft == 0 %Right
                            this.allStimShift = mx+predicted_travel_distance_pix/2+pursueDist;
                            % this.allStimShift = mx+round(predicted_travel_distance_pix); %jd1/29
                        else
                            this.allStimShift = mx-predicted_travel_distance_pix/2-pursueDist;
                            % this.allStimShift = mx-round(predicted_travel_distance_pix); %jd1/29
                        end

                        %change background
                        [~] = setCowBkgdwithCross(graph, mx, my, this, tex1oF, white);

                        if thisTrialData.control_Fix_Left ==1
                            %LEFT
                            if thisTrialData.StartingFromLeft == 0
                                fixRect = CenterRectOnPointd( this.stimRect_size, (secondsElapsed - first_reset_secondsElapsed)*thisTrialData.Speed_Ref_pix-Lead_LRoffset +this.allStimShift, my );
                            else
                                fixRect = CenterRectOnPointd( this.stimRect_size, (secondsElapsed - first_reset_secondsElapsed)*thisTrialData.Speed_TEST_pix-Lag_LRoffset +this.allStimShift, my ); %start left
                            end

                            Screen('FillRect', graph.window,  this.stimColor, fixRect);

                                    % [Adding left--fixation circle in center]----------------------------------------------
                                    fixCirc_size = [0 0 this.fix_size this.fix_size];
                                    fixRect = CenterRectOnPointd( fixCirc_size,   (secondsElapsed - first_reset_secondsElapsed)*thisTrialData.Speed_TEST_pix-thisTrialData.LRoffset +this.allStimShift, my );
                                    Screen('FillOval', graph.window,  black, fixRect);
        
                                    fixCross_size = [0 0 this.fix_size/4 this.fix_size];
                                    fixRect = CenterRectOnPointd( fixCross_size,   (secondsElapsed - first_reset_secondsElapsed)*thisTrialData.Speed_TEST_pix-thisTrialData.LRoffset +this.allStimShift, my );
                                    Screen('FillOval', graph.window,  white, fixRect);
        
                                    fixCross_size = [0 0 this.fix_size this.fix_size/4];
                                    fixRect = CenterRectOnPointd( fixCross_size,   (secondsElapsed - first_reset_secondsElapsed)*thisTrialData.Speed_TEST_pix-thisTrialData.LRoffset +this.allStimShift, my );
                                    Screen('FillOval', graph.window,  white, fixRect);
                                    % ---------------------------------------------- [Adding middle--fixation circle in center]

                            reset_secondsElapsed =(secondsElapsed - first_reset_secondsElapsed); %exactly amount of time you need to subtract, so that the three are equidistant in the beginning
                            reset_secondsElapsed_middle = (secondsElapsed - first_reset_secondsElapsed)*thisTrialData.Speed_TEST_pix; %how much it traveled during the saccade phase
                        elseif thisTrialData.control_Fix_noMove == 1
                                    if thisTrialData.randEye == 0
                                        % [Adding left--fixation circle in center]----------------------------------------------
                                        fixCirc_size = [0 0 this.fix_size this.fix_size];
                                        fixRect = CenterRectOnPointd( fixCirc_size,  mx, my);%+ (this.blackRing_half_length-this.insideCircle_half_length)/2+this.insideCircle_half_length );
                                        Screen('FillOval', graph.window,  black, fixRect);
        
                                        fixCross_size = [0 0 this.fix_size/4 this.fix_size];
                                        fixRect = CenterRectOnPointd( fixCross_size,  mx, my);%+ (this.blackRing_half_length-this.insideCircle_half_length)/2+this.insideCircle_half_length );
                                        Screen('FillOval', graph.window,  white, fixRect);
        
                                        fixCross_size = [0 0 this.fix_size this.fix_size/4];
                                        fixRect = CenterRectOnPointd( fixCross_size, mx, my);%+ (this.blackRing_half_length-this.insideCircle_half_length)/2+this.insideCircle_half_length );
                                        Screen('FillOval', graph.window,  white, fixRect);
                                        % ---------------------------------------------- [Adding middle--fixation circle in center]
                                    end
                            reset_secondsElapsed =secondsElapsed; %exactly amount of time you need to subtract, so that the three are equidistant in the beginning
                            reset_secondsElapsed_middle = (secondsElapsed - first_reset_secondsElapsed)*thisTrialData.Speed_middle_pix; %how much it traveled during the saccade phase
                            
                            predicted_travel_distance_pix = ((this.stimOnPursue_seconds)*  abs(thisTrialData.Speed_middle_pix));
                            pursueDist = this.stimOffPursue_seconds*  abs(thisTrialData.Speed_middle_pix);
                            if thisTrialData.StartingFromLeft == 0 
                                this.allStimShift = mx+predicted_travel_distance_pix/2+pursueDist;
                                if thisTrialData.endAtCross==1
                                    this.allStimShift = mx+predicted_travel_distance_pix+pursueDist -this.fix_size/2;%jd1/29
                                end
                            else
                                this.allStimShift = mx-predicted_travel_distance_pix/2-pursueDist;
                                 if thisTrialData.endAtCross==1
                                     this.allStimShift = mx-predicted_travel_distance_pix-pursueDist +this.fix_size/2;%jd1/29
                                 end
                            end

                        elseif thisTrialData.NOPursueExp == 1
                             
                            %MIDDLE
                            % [Adding left--fixation circle in center]----------------------------------------------
                            fixCirc_size = [0 0 this.fix_size this.fix_size];
                            fixRect = CenterRectOnPointd( fixCirc_size,  mx, my);%+ (this.blackRing_half_length-this.insideCircle_half_length)/2+this.insideCircle_half_length );
                            Screen('FillOval', graph.window,  black, fixRect);

                            fixCross_size = [0 0 this.fix_size/4 this.fix_size];
                            fixRect = CenterRectOnPointd( fixCross_size,  mx, my);%+ (this.blackRing_half_length-this.insideCircle_half_length)/2+this.insideCircle_half_length );
                            Screen('FillOval', graph.window,  white, fixRect);

                            fixCross_size = [0 0 this.fix_size this.fix_size/4];
                            fixRect = CenterRectOnPointd( fixCross_size, mx, my);%+ (this.blackRing_half_length-this.insideCircle_half_length)/2+this.insideCircle_half_length );
                            Screen('FillOval', graph.window,  white, fixRect);
                            % ---------------------------------------------- [Adding middle--fixation circle in center]

                            reset_secondsElapsed =secondsElapsed ; %exactly amount of time you need to subtract, so that the three are equidistant in the beginning
                            reset_secondsElapsed_middle = (secondsElapsed - first_reset_secondsElapsed)*thisTrialData.Speed_middle_pix; %how much it traveled during the saccade phase

                        else %normal experiment
                            %MIDDLE
                            fixRect = CenterRectOnPointd( this.stimRect_size,   (secondsElapsed - first_reset_secondsElapsed)*thisTrialData.Speed_middle_pix +this.allStimShift, my );
                            Screen('FillRect', graph.window,  this.stimColor, fixRect);

                                    % [Adding middle--fixation circle in center]----------------------------------------------
                                    fixCirc_size = [0 0 this.fix_size this.fix_size];
                                    fixRect = CenterRectOnPointd( fixCirc_size,   (secondsElapsed - first_reset_secondsElapsed)*thisTrialData.Speed_middle_pix +this.allStimShift, my );
                                    Screen('FillOval', graph.window,  black, fixRect);
        
                                    fixCross_size = [0 0 this.fix_size/4 this.fix_size];
                                    fixRect = CenterRectOnPointd( fixCross_size,   (secondsElapsed - first_reset_secondsElapsed)*thisTrialData.Speed_middle_pix +this.allStimShift, my );
                                    Screen('FillOval', graph.window,  white, fixRect);
        
                                    fixCross_size = [0 0 this.fix_size this.fix_size/4];
                                    fixRect = CenterRectOnPointd( fixCross_size,   (secondsElapsed - first_reset_secondsElapsed)*thisTrialData.Speed_middle_pix +this.allStimShift, my );
                                    Screen('FillOval', graph.window,  white, fixRect);
                                    % ---------------------------------------------- [Adding middle--fixation circle in center]
                            
                            reset_secondsElapsed =secondsElapsed ; %exactly amount of time you need to subtract, so that the three are equidistant in the beginning
                            reset_secondsElapsed_middle = (secondsElapsed - first_reset_secondsElapsed)*thisTrialData.Speed_middle_pix; %how much it traveled during the saccade phase

                        end
                    
                    % Turn on Stimuli (white bars)
                    elseif ( secondsElapsed >= (this.stimOffFixate_seconds+this.stimOffPursue_seconds) && secondsElapsed <= (this.stimOffFixate_seconds+this.stimOffPursue_seconds + this.stimOnPursue_seconds))
                        
                        predicted_travel_distance_pix_middle = (this.stimOnPursue_seconds*  abs(thisTrialData.Speed_middle_pix));
                        predicted_travel_distance_pix_REF = (this.stimOnPursue_seconds*  abs(thisTrialData.Speed_ref_pix));
                        predicted_travel_distance_pix_TEST = (this.stimOnPursue_seconds*  abs(thisTrialData.Speed_TEST_pix));

                        predicted_lag = abs(predicted_travel_distance_pix_middle - predicted_travel_distance_pix_TEST);
                        predicted_lead = abs(predicted_travel_distance_pix_REF - predicted_travel_distance_pix_middle);

                        %change background
                        [~] = setCowBkgdwithCross(graph, mx, my, this, tex1oF, white);

                        %-- Find the center of the screen
                        [mx, my] = RectCenter(graph.wRect);

                        %LEFT
                        if thisTrialData.StartingFromLeft == 0
                            if thisTrialData.GrowingOutwards == 1 %normal
                                fixRect = CenterRectOnPointd( this.stimRect_size, (secondsElapsed-reset_secondsElapsed)*  thisTrialData.Speed_ref_pix-Lead_LRoffset+reset_secondsElapsed_middle +this.allStimShift, my );
                            else %white bar shrinks
                                % fixRect = CenterRectOnPointd( this.stimRect_size, (secondsElapsed-reset_secondsElapsed)*  thisTrialData.Speed_middle_pix-thisTrialData.LRoffset+reset_secondsElapsed_middle +this.allStimShift      - predicted_lead* (1-(secondsElapsed-reset_secondsElapsed)/this.stimOnPursue_seconds), my );
                            
                                %making the left one now be Test
                                fixRect = CenterRectOnPointd( this.stimRect_size, (secondsElapsed-reset_secondsElapsed)*  thisTrialData.Speed_middle_pix -Lead_LRoffset+reset_secondsElapsed_middle +this.allStimShift -     predicted_lead* (1-(secondsElapsed-reset_secondsElapsed)/this.stimOnPursue_seconds)+  (predicted_lag-predicted_lead)* (secondsElapsed-reset_secondsElapsed), my ); %this one starts equal and ends unequal
                            end
                        else
                            if thisTrialData.GrowingOutwards == 1 %normal
                                fixRect = CenterRectOnPointd( this.stimRect_size, (secondsElapsed-reset_secondsElapsed)*  thisTrialData.Speed_TEST_pix-Lag_LRoffset+reset_secondsElapsed_middle +this.allStimShift, my );
                            else % white bar shrinks
                                % % % fixRect = CenterRectOnPointd( this.stimRect_size, (secondsElapsed-reset_secondsElapsed)*  thisTrialData.Speed_middle_pix -thisTrialData.LRoffset+reset_secondsElapsed_middle +this.allStimShift -     predicted_lag* (1-(secondsElapsed-reset_secondsElapsed)/this.stimOnPursue_seconds), my ); %this is so that it ends equal,do I want that?
                                % fixRect = CenterRectOnPointd( this.stimRect_size, (secondsElapsed-reset_secondsElapsed)*  thisTrialData.Speed_middle_pix -thisTrialData.LRoffset+reset_secondsElapsed_middle +this.allStimShift -     predicted_lead* (1-(secondsElapsed-reset_secondsElapsed)/this.stimOnPursue_seconds)+  (predicted_lag-predicted_lead)* (secondsElapsed-reset_secondsElapsed), my ); %this one starts equal and ends unequal
                                
                                %making the left one now be Ref
                                fixRect = CenterRectOnPointd( this.stimRect_size, (secondsElapsed-reset_secondsElapsed)*  thisTrialData.Speed_middle_pix-Lag_LRoffset+reset_secondsElapsed_middle +this.allStimShift      - predicted_lead* (1-(secondsElapsed-reset_secondsElapsed)/this.stimOnPursue_seconds), my );
                                
                            end
                        end
                        % fixRect(1) = 0;
                        Screen('FillRect', graph.window,  this.stimColor, fixRect);
                        
                        %MIDDLE
                        fixRect = CenterRectOnPointd( this.stimRect_size, (secondsElapsed-reset_secondsElapsed)*  thisTrialData.Speed_middle_pix+reset_secondsElapsed_middle +this.allStimShift, my );
                        Screen('FillRect', graph.window,  this.stimColor, fixRect);

                                if thisTrialData.control_Fix_Left ==1
                                    % [Adding left--fixation circle in center]----------------------------------------------
                                    fixCirc_size = [0 0 this.fix_size this.fix_size];
                                    fixRect = CenterRectOnPointd( fixCirc_size, (secondsElapsed-reset_secondsElapsed)*  thisTrialData.Speed_TEST_pix-thisTrialData.LRoffset+reset_secondsElapsed_middle +this.allStimShift, my );
                                    Screen('FillOval', graph.window,  black, fixRect);
        
                                    fixCross_size = [0 0 this.fix_size/4 this.fix_size];
                                    fixRect = CenterRectOnPointd( fixCross_size, (secondsElapsed-reset_secondsElapsed)*  thisTrialData.Speed_TEST_pix-thisTrialData.LRoffset+reset_secondsElapsed_middle +this.allStimShift, my );
                                    Screen('FillOval', graph.window,  white, fixRect);
        
                                    fixCross_size = [0 0 this.fix_size this.fix_size/4];
                                    fixRect = CenterRectOnPointd( fixCross_size, (secondsElapsed-reset_secondsElapsed)*  thisTrialData.Speed_TEST_pix-thisTrialData.LRoffset+reset_secondsElapsed_middle +this.allStimShift, my );
                                    Screen('FillOval', graph.window,  white, fixRect);
                                    % ---------------------------------------------- [Adding middle--fixation circle in center]
                                elseif thisTrialData.control_Fix_noMove == 1
                                    if thisTrialData.randEye == 0
                                        % [Adding left--fixation circle in center]----------------------------------------------
                                        fixCirc_size = [0 0 this.fix_size this.fix_size];
                                        fixRect = CenterRectOnPointd( fixCirc_size,  mx, my);%+  (this.blackRing_half_length-this.insideCircle_half_length)/2+this.insideCircle_half_length );
                                        Screen('FillOval', graph.window,  black, fixRect);

                                        fixCross_size = [0 0 this.fix_size/4 this.fix_size];
                                        fixRect = CenterRectOnPointd( fixCross_size,  mx, my);%+ (this.blackRing_half_length-this.insideCircle_half_length)/2+this.insideCircle_half_length );
                                        Screen('FillOval', graph.window,  white, fixRect);

                                        fixCross_size = [0 0 this.fix_size this.fix_size/4];
                                        fixRect = CenterRectOnPointd( fixCross_size, mx, my);%+  (this.blackRing_half_length-this.insideCircle_half_length)/2+this.insideCircle_half_length );
                                        Screen('FillOval', graph.window,  white, fixRect);
                                    end
                                else
                                    % [Adding middle--fixation circle in center]----------------------------------------------
                                    fixCirc_size = [0 0 this.fix_size this.fix_size];
                                    fixRect = CenterRectOnPointd( fixCirc_size, (secondsElapsed-reset_secondsElapsed)*  thisTrialData.Speed_middle_pix+reset_secondsElapsed_middle +this.allStimShift, my );
                                    Screen('FillOval', graph.window,  black, fixRect);
        
                                    fixCross_size = [0 0 this.fix_size/4 this.fix_size];
                                    fixRect = CenterRectOnPointd( fixCross_size, (secondsElapsed-reset_secondsElapsed)*  thisTrialData.Speed_middle_pix+reset_secondsElapsed_middle +this.allStimShift, my );
                                    Screen('FillOval', graph.window,  white, fixRect);
        
                                    fixCross_size = [0 0 this.fix_size this.fix_size/4];
                                    fixRect = CenterRectOnPointd( fixCross_size, (secondsElapsed-reset_secondsElapsed)*  thisTrialData.Speed_middle_pix+reset_secondsElapsed_middle +this.allStimShift, my );
                                    Screen('FillOval', graph.window,  white, fixRect);
                                    % ---------------------------------------------- [Adding middle--fixation circle in center]
                                end
                        %RIGHT
                        if thisTrialData.StartingFromLeft == 0
                            if thisTrialData.GrowingOutwards == 1 %normal
                                fixRect = CenterRectOnPointd( this.stimRect_size, (secondsElapsed-reset_secondsElapsed)*  thisTrialData.Speed_TEST_pix+Lag_LRoffset+reset_secondsElapsed_middle +this.allStimShift, my );
                            else % white bar shrinks
                                % fixRect = CenterRectOnPointd( this.stimRect_size, (secondsElapsed-reset_secondsElapsed)*  thisTrialData.Speed_middle_pix +thisTrialData.LRoffset+reset_secondsElapsed_middle +this.allStimShift +     predicted_lead* (1-(secondsElapsed-reset_secondsElapsed)/this.stimOnPursue_seconds)-  (predicted_lag-predicted_lead)* (secondsElapsed-reset_secondsElapsed), my ); %this one starts equal and ends unequal
                            
                                %making right one now be Ref
                                fixRect = CenterRectOnPointd( this.stimRect_size, (secondsElapsed-reset_secondsElapsed)*  thisTrialData.Speed_middle_pix+Lag_LRoffset+reset_secondsElapsed_middle +this.allStimShift      + predicted_lead* (1-(secondsElapsed-reset_secondsElapsed)/this.stimOnPursue_seconds), my );
                            end
                        else
                            if thisTrialData.GrowingOutwards == 1 %normal
                                fixRect = CenterRectOnPointd( this.stimRect_size, (secondsElapsed-reset_secondsElapsed)*  thisTrialData.Speed_ref_pix+Lead_LRoffset+reset_secondsElapsed_middle +this.allStimShift, my );
                            else % white bar shrinks
                                % fixRect = CenterRectOnPointd( this.stimRect_size, (secondsElapsed-reset_secondsElapsed)*  thisTrialData.Speed_middle_pix+thisTrialData.LRoffset+reset_secondsElapsed_middle +this.allStimShift      + predicted_lead* (1-(secondsElapsed-reset_secondsElapsed)/this.stimOnPursue_seconds), my );
                                
                                %making the right one now be Test
                                fixRect = CenterRectOnPointd( this.stimRect_size, (secondsElapsed-reset_secondsElapsed)*  thisTrialData.Speed_middle_pix +Lead_LRoffset+reset_secondsElapsed_middle +this.allStimShift +     predicted_lead* (1-(secondsElapsed-reset_secondsElapsed)/this.stimOnPursue_seconds)-  (predicted_lag-predicted_lead)* (secondsElapsed-reset_secondsElapsed), my ); %this one starts equal and ends unequal
                            end
                        end
                        % fixRect(3) = mx*2;
                        Screen('FillRect', graph.window,  this.stimColor, fixRect);
        %%%
                                if thisTrialData.control_Fix_noMove == 1 && thisTrialData.randEye == 0
                                    % [Adding left--fixation circle in center]----------------------------------------------
                                    fixCirc_size = [0 0 this.fix_size this.fix_size];
                                    fixRect = CenterRectOnPointd( fixCirc_size,  mx, my);%+  (this.blackRing_half_length-this.insideCircle_half_length)/2+this.insideCircle_half_length );
                                    Screen('FillOval', graph.window,  black, fixRect);

                                    fixCross_size = [0 0 this.fix_size/4 this.fix_size];
                                    fixRect = CenterRectOnPointd( fixCross_size,  mx, my);%+ (this.blackRing_half_length-this.insideCircle_half_length)/2+this.insideCircle_half_length );
                                    Screen('FillOval', graph.window,  white, fixRect);

                                    fixCross_size = [0 0 this.fix_size this.fix_size/4];
                                    fixRect = CenterRectOnPointd( fixCross_size, mx, my);%+  (this.blackRing_half_length-this.insideCircle_half_length)/2+this.insideCircle_half_length );
                                    Screen('FillOval', graph.window,  white, fixRect);
                                end

                        
                                fade_reset_secondsElapsed = secondsElapsed;
                                
                    elseif (secondsElapsed > (this.stimOffFixate_seconds+this.stimOffPursue_seconds + this.stimOnPursue_seconds) && secondsElapsed <= (this.stimOffFixate_seconds+this.stimOffPursue_seconds + this.stimOnPursue_seconds + this.stimFadeOut_seconds))
                        
                        fadeC =  min((secondsElapsed-fade_reset_secondsElapsed) / this.stimFadeOut_seconds, 1);
                        FadeToWhiteColor = [fadeC fadeC fadeC].*255;

                        if this.blankblack == 1
                            FadeToWhiteColor = 255-FadeToWhiteColor;
                        end
                        FadeToWhiteColor = white;
                        %change background
                        [~] = setCowBkgdwithCross(graph, mx, my, this, tex1oF, FadeToWhiteColor);

                        [reset_secondsElapsed] = fadeOut(this, thisTrialData, secondsElapsed, reset_secondsElapsed, fade_reset_secondsElapsed, graph);

                    else % >4s just wait until exp over
                        
                        %change background
                        % [~] = setCowBkgdwithCross(graph, mx, my, this, tex1oF, grey);

                        % % [wo1, wo2, wo3, wo4, wo5, wo6, wo7, wo8, wo9, wo10] = ...
                        % %     makeWOTextures(graph.pxWidth, graph.pxHeight, graph.window);
                        % % 
                        % % %-- Find the center of the screen
                        % % [mx, my] = RectCenter(graph.wRect);
                        % % 
                        % % newTextNumof20 = this.textOrderforTrials(thisTrialData.TrialNumber);
                        % % if newTextNumof20 >10
                        % %     newTextNumof20 = newTextNumof20-10;
                        % % end
                        % % maskTex = eval(sprintf('wo%1.f',newTextNumof20));
                        % % 
                        % % destinationRectangle =  []; % so that it fills the whole screen
                        % % Screen('DrawTexture', graph.window, maskTex, [],destinationRectangle);

                        folder = fullfile(fileparts(mfilename('fullpath')), 'Masks');
                        img = imread(fullfile(folder,'voronoi-v4-1.png'));

                        targetW = graph.pxWidth;
                        targetH = graph.pxHeight;

                        [h, w, ~] = size(img);

                        % Compute scale factor so image fully covers target
                        scale = max(targetW / w, targetH / h);

                        % Resize
                        resized = imresize(img, scale);

                        % Get new size
                        [newH, newW, ~] = size(resized);

                        % Compute center crop coordinates
                        xStart = floor((newW - targetW) / 2) + 1;
                        yStart = floor((newH - targetH) / 2) + 1;

                        cropped = resized(yStart:yStart+targetH-1, ...
                            xStart:xStart+targetW-1, :);

                        maskTex = Screen('MakeTexture', graph.window, cropped);

                        destinationRectangle =  []; % so that it fills the whole screen
                        Screen('DrawTexture', graph.window, maskTex, [],destinationRectangle);
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

                    if ( secondsElapsed > (this.stimOffFixate_seconds+ this.stimOffPursue_seconds + this.stimOnPursue_seconds + this.stimFadeOut_seconds) ) %not including fade out so that you can respond during fade out 
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
                                Plot_Psychometric_2graphs(this)
                                % Plot_EyeData(this)
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

            function [reset_secondsElapsed] = fadeOut(this, thisTrialData, secondsElapsed, reset_secondsElapsed, fade_reset_secondsElapsed, graph)
                %LEFT
                if thisTrialData.StartingFromLeft == 0
                    if thisTrialData.GrowingOutwards == 1 %normal
                        fixRect = CenterRectOnPointd( this.stimRect_size, (secondsElapsed-reset_secondsElapsed)*  thisTrialData.Speed_ref_pix-thisTrialData.Lead_LRoffset+reset_secondsElapsed_middle +this.allStimShift, my );
                    else %white bar shrinks
                        % fixRect = CenterRectOnPointd( this.stimRect_size, (secondsElapsed-reset_secondsElapsed)*  thisTrialData.Speed_middle_pix-thisTrialData.LRoffset+reset_secondsElapsed_middle +this.allStimShift      - predicted_lead* (1-(secondsElapsed-reset_secondsElapsed)/this.stimOnPursue_seconds), my );

                        %making the left one now be Test
                        fixRect = CenterRectOnPointd( this.stimRect_size, (secondsElapsed-reset_secondsElapsed)*  thisTrialData.Speed_middle_pix -thisTrialData.Lead_LRoffset+reset_secondsElapsed_middle +this.allStimShift -     predicted_lead* (1-(secondsElapsed-reset_secondsElapsed)/this.stimOnPursue_seconds)+  (predicted_lag-predicted_lead)* (secondsElapsed-reset_secondsElapsed), my ); %this one starts equal and ends unequal
                    end
                else
                    if thisTrialData.GrowingOutwards == 1 %normal
                        fixRect = CenterRectOnPointd( this.stimRect_size, (secondsElapsed-reset_secondsElapsed)*  thisTrialData.Speed_TEST_pix-thisTrialData.Lag_LRoffset+reset_secondsElapsed_middle +this.allStimShift, my );
                    else % white bar shrinks
                        % % % fixRect = CenterRectOnPointd( this.stimRect_size, (secondsElapsed-reset_secondsElapsed)*  thisTrialData.Speed_middle_pix -thisTrialData.LRoffset+reset_secondsElapsed_middle +this.allStimShift -     predicted_lag* (1-(secondsElapsed-reset_secondsElapsed)/this.stimOnPursue_seconds), my ); %this is so that it ends equal,do I want that?
                        % fixRect = CenterRectOnPointd( this.stimRect_size, (secondsElapsed-reset_secondsElapsed)*  thisTrialData.Speed_middle_pix -thisTrialData.LRoffset+reset_secondsElapsed_middle +this.allStimShift -     predicted_lead* (1-(secondsElapsed-reset_secondsElapsed)/this.stimOnPursue_seconds)+  (predicted_lag-predicted_lead)* (secondsElapsed-reset_secondsElapsed), my ); %this one starts equal and ends unequal

                        %making the left one now be Ref
                        fixRect = CenterRectOnPointd( this.stimRect_size, (secondsElapsed-reset_secondsElapsed)*  thisTrialData.Speed_middle_pix-thisTrialData.Lag_LRoffset+reset_secondsElapsed_middle +this.allStimShift      - predicted_lead* (1-(secondsElapsed-reset_secondsElapsed)/this.stimOnPursue_seconds), my );

                    end
                end
                % fixRect(1) = 0;
                t =   min((secondsElapsed-fade_reset_secondsElapsed) / this.stimFadeOut_seconds, 1);
                FadeToWhiteColor =[t t t].*255; % white; %

                Screen('FillRect', graph.window,  FadeToWhiteColor, fixRect);

                %MIDDLE
                fixRect = CenterRectOnPointd( this.stimRect_size, (secondsElapsed-reset_secondsElapsed)*  thisTrialData.Speed_middle_pix+reset_secondsElapsed_middle +this.allStimShift, my );
                Screen('FillRect', graph.window,  this.stimColor, fixRect);

                if thisTrialData.control_Fix_Left ==1
                            % [Adding left--fixation circle in center]----------------------------------------------
                            fixCirc_size = [0 0 this.fix_size this.fix_size];
                            fixRect = CenterRectOnPointd( fixCirc_size, (secondsElapsed-reset_secondsElapsed)*  thisTrialData.Speed_TEST_pix-thisTrialData.LRoffset+reset_secondsElapsed_middle +this.allStimShift, my );
                            Screen('FillOval', graph.window,  black, fixRect);
        
                            fixCross_size = [0 0 this.fix_size/4 this.fix_size];
                            fixRect = CenterRectOnPointd( fixCross_size, (secondsElapsed-reset_secondsElapsed)*  thisTrialData.Speed_TEST_pix-thisTrialData.LRoffset+reset_secondsElapsed_middle +this.allStimShift, my );
                            Screen('FillOval', graph.window,  white, fixRect);
        
                            fixCross_size = [0 0 this.fix_size this.fix_size/4];
                            fixRect = CenterRectOnPointd( fixCross_size, (secondsElapsed-reset_secondsElapsed)*  thisTrialData.Speed_TEST_pix-thisTrialData.LRoffset+reset_secondsElapsed_middle +this.allStimShift, my );
                            Screen('FillOval', graph.window,  white, fixRect);
                            % ---------------------------------------------- [Adding middle--fixation circle in center]
                elseif thisTrialData.control_Fix_noMove == 1
                            if thisTrialData.randEye == 0
                                % [Adding left--fixation circle in center]----------------------------------------------
                                fixCirc_size = [0 0 this.fix_size this.fix_size];
                                fixRect = CenterRectOnPointd( fixCirc_size,  mx, my);%+  (this.blackRing_half_length-this.insideCircle_half_length)/2+this.insideCircle_half_length );
                                Screen('FillOval', graph.window,  black, fixRect);
        
                                fixCross_size = [0 0 this.fix_size/4 this.fix_size];
                                fixRect = CenterRectOnPointd( fixCross_size,  mx, my);%+ (this.blackRing_half_length-this.insideCircle_half_length)/2+this.insideCircle_half_length );
                                Screen('FillOval', graph.window,  white, fixRect);
        
                                fixCross_size = [0 0 this.fix_size this.fix_size/4];
                                fixRect = CenterRectOnPointd( fixCross_size, mx, my);%+  (this.blackRing_half_length-this.insideCircle_half_length)/2+this.insideCircle_half_length );
                                Screen('FillOval', graph.window,  white, fixRect);
                            end
                else
                            % [Adding middle--fixation circle in center]----------------------------------------------
                            fixCirc_size = [0 0 this.fix_size this.fix_size];
                            fixRect = CenterRectOnPointd( fixCirc_size, (secondsElapsed-reset_secondsElapsed)*  thisTrialData.Speed_middle_pix+reset_secondsElapsed_middle +this.allStimShift, my );
                            Screen('FillOval', graph.window,  black, fixRect);
        
                            fixCross_size = [0 0 this.fix_size/4 this.fix_size];
                            fixRect = CenterRectOnPointd( fixCross_size, (secondsElapsed-reset_secondsElapsed)*  thisTrialData.Speed_middle_pix+reset_secondsElapsed_middle +this.allStimShift, my );
                            Screen('FillOval', graph.window,  white, fixRect);
        
                            fixCross_size = [0 0 this.fix_size this.fix_size/4];
                            fixRect = CenterRectOnPointd( fixCross_size, (secondsElapsed-reset_secondsElapsed)*  thisTrialData.Speed_middle_pix+reset_secondsElapsed_middle +this.allStimShift, my );
                            Screen('FillOval', graph.window,  white, fixRect);
                            % ---------------------------------------------- [Adding middle--fixation circle in center]
                end
                %RIGHT
                if thisTrialData.StartingFromLeft == 0
                    if thisTrialData.GrowingOutwards == 1 %normal
                        fixRect = CenterRectOnPointd( this.stimRect_size, (secondsElapsed-reset_secondsElapsed)*  thisTrialData.Speed_TEST_pix+thisTrialData.Lag_LRoffset+reset_secondsElapsed_middle +this.allStimShift, my );
                    else % white bar shrinks
                        % fixRect = CenterRectOnPointd( this.stimRect_size, (secondsElapsed-reset_secondsElapsed)*  thisTrialData.Speed_middle_pix +thisTrialData.LRoffset+reset_secondsElapsed_middle +this.allStimShift +     predicted_lead* (1-(secondsElapsed-reset_secondsElapsed)/this.stimOnPursue_seconds)-  (predicted_lag-predicted_lead)* (secondsElapsed-reset_secondsElapsed), my ); %this one starts equal and ends unequal

                        %making right one now be Ref
                        fixRect = CenterRectOnPointd( this.stimRect_size, (secondsElapsed-reset_secondsElapsed)*  thisTrialData.Speed_middle_pix+thisTrialData.Lag_LRoffset+reset_secondsElapsed_middle +this.allStimShift      + predicted_lead* (1-(secondsElapsed-reset_secondsElapsed)/this.stimOnPursue_seconds), my );
                    end
                else
                    if thisTrialData.GrowingOutwards == 1 %normal
                        fixRect = CenterRectOnPointd( this.stimRect_size, (secondsElapsed-reset_secondsElapsed)*  thisTrialData.Speed_ref_pix+thisTrialData.Lead_LRoffset+reset_secondsElapsed_middle +this.allStimShift, my );
                    else % white bar shrinks
                        % fixRect = CenterRectOnPointd( this.stimRect_size, (secondsElapsed-reset_secondsElapsed)*  thisTrialData.Speed_middle_pix+thisTrialData.LRoffset+reset_secondsElapsed_middle +this.allStimShift      + predicted_lead* (1-(secondsElapsed-reset_secondsElapsed)/this.stimOnPursue_seconds), my );

                        %making the right one now be Test
                        fixRect = CenterRectOnPointd( this.stimRect_size, (secondsElapsed-reset_secondsElapsed)*  thisTrialData.Speed_middle_pix +thisTrialData.Lead_LRoffset+reset_secondsElapsed_middle +this.allStimShift +     predicted_lead* (1-(secondsElapsed-reset_secondsElapsed)/this.stimOnPursue_seconds)-  (predicted_lag-predicted_lead)* (secondsElapsed-reset_secondsElapsed), my ); %this one starts equal and ends unequal
                    end
                end
                % fixRect(3) = mx*2;
                jdColor = [255 0 0];
                Screen('FillRect', graph.window,  FadeToWhiteColor, fixRect);
                %%%
                        if thisTrialData.control_Fix_noMove == 1 && thisTrialData.randEye == 0
                            % [Adding left--fixation circle in center]----------------------------------------------
                            fixCirc_size = [0 0 this.fix_size this.fix_size];
                            fixRect = CenterRectOnPointd( fixCirc_size,  mx, my);%+  (this.blackRing_half_length-this.insideCircle_half_length)/2+this.insideCircle_half_length );
                            Screen('FillOval', graph.window,  black, fixRect);
        
                            fixCross_size = [0 0 this.fix_size/4 this.fix_size];
                            fixRect = CenterRectOnPointd( fixCross_size,  mx, my);%+ (this.blackRing_half_length-this.insideCircle_half_length)/2+this.insideCircle_half_length );
                            Screen('FillOval', graph.window,  white, fixRect);
        
                            fixCross_size = [0 0 this.fix_size this.fix_size/4];
                            fixRect = CenterRectOnPointd( fixCross_size, mx, my);%+  (this.blackRing_half_length-this.insideCircle_half_length)/2+this.insideCircle_half_length );
                            Screen('FillOval', graph.window,  white, fixRect);
                        end

            end %end fade out

            % Make 20 textures
            function [tex1, tex2, tex3, tex4, tex5, tex6, tex7, tex8, tex9, tex10] = ...
                    makeTextures(wholeLumRasterCancelJDset, screen_width, screen_height, w)

                %cd('/Users/Josephine1/Documents/MATLAB/Roorda Lab/Data_and_AnalysisCode/gifMakers_presentations/JDs20noiseBMPS_folder'); % I know this is redundant, but didn't know a better way JD 7/21/23

                folder = fullfile(fileparts(mfilename('fullpath')), 'JDs20noiseBMPS_folder');

                jdnoise1 = imread(fullfile(folder,'JDnoise1.bmp'));
                [jdnoise1] = cropAndMakeBkgLum(jdnoise1,wholeLumRasterCancelJDset, screen_width, screen_height);
                tex1 = Screen('MakeTexture', w, jdnoise1);

                jdnoise2 = imread(fullfile(folder,'JDnoise2.bmp'));
                [jdnoise2] = cropAndMakeBkgLum(jdnoise2,wholeLumRasterCancelJDset, screen_width, screen_height);
                tex2 = Screen('MakeTexture', w, jdnoise2);

                jdnoise3 = imread(fullfile(folder,'JDnoise20.bmp')); %3 had a dalmation dog face distracting 3/12/26
                [jdnoise3] = cropAndMakeBkgLum(jdnoise3,wholeLumRasterCancelJDset, screen_width, screen_height);
                tex3 = Screen('MakeTexture', w, jdnoise3);

                jdnoise4 = imread(fullfile(folder,'JDnoise4.bmp'));
                [jdnoise4] = cropAndMakeBkgLum(jdnoise4,wholeLumRasterCancelJDset, screen_width, screen_height);
                tex4 = Screen('MakeTexture', w, jdnoise4);

                jdnoise5 = imread(fullfile(folder,'JDnoise5.bmp'));
                [jdnoise5] = cropAndMakeBkgLum(jdnoise5,wholeLumRasterCancelJDset, screen_width, screen_height);
                tex5 = Screen('MakeTexture', w, jdnoise5);

                jdnoise6 = imread(fullfile(folder,'JDnoise11.bmp'));
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

                % jdnoise11 = imread(fullfile(folder,'JDnoise11.bmp'));
                % [jdnoise11] = cropAndMakeBkgLum(jdnoise11,wholeLumRasterCancelJDset, screen_width, screen_height);
                % tex11 = Screen('MakeTexture', w, jdnoise11);
                % 
                % jdnoise12 = imread(fullfile(folder,'JDnoise12.bmp'));
                % [jdnoise12] = cropAndMakeBkgLum(jdnoise12,wholeLumRasterCancelJDset, screen_width, screen_height);
                % tex12= Screen('MakeTexture', w, jdnoise12);
                % 
                % jdnoise13 = imread(fullfile(folder,'JDnoise13.bmp'));
                % [jdnoise13] = cropAndMakeBkgLum(jdnoise13,wholeLumRasterCancelJDset, screen_width, screen_height);
                % tex13 = Screen('MakeTexture', w, jdnoise13);
                % 
                % jdnoise14 = imread(fullfile(folder,'JDnoise14.bmp'));
                % [jdnoise14] = cropAndMakeBkgLum(jdnoise14,wholeLumRasterCancelJDset, screen_width, screen_height);
                % tex14 = Screen('MakeTexture', w, jdnoise14);
                % 
                % jdnoise15 = imread(fullfile(folder,'JDnoise15.bmp'));
                % [jdnoise15] = cropAndMakeBkgLum(jdnoise15,wholeLumRasterCancelJDset, screen_width, screen_height);
                % tex15 = Screen('MakeTexture', w, jdnoise15);
                % 
                % jdnoise16 = imread(fullfile(folder,'JDnoise16.bmp'));
                % [jdnoise16] = cropAndMakeBkgLum(jdnoise16,wholeLumRasterCancelJDset, screen_width, screen_height);
                % tex16 = Screen('MakeTexture', w, jdnoise16);
                % 
                % jdnoise17 = imread(fullfile(folder,'JDnoise17.bmp'));
                % [jdnoise17] = cropAndMakeBkgLum(jdnoise17,wholeLumRasterCancelJDset, screen_width, screen_height);
                % tex17 = Screen('MakeTexture', w, jdnoise17);
                % 
                % jdnoise18 = imread(fullfile(folder,'JDnoise18.bmp'));
                % [jdnoise18] = cropAndMakeBkgLum(jdnoise18,wholeLumRasterCancelJDset, screen_width, screen_height);
                % tex18 = Screen('MakeTexture', w, jdnoise18);
                % 
                % jdnoise19 = imread(fullfile(folder,'JDnoise19.bmp'));
                % [jdnoise19] = cropAndMakeBkgLum(jdnoise19,wholeLumRasterCancelJDset, screen_width, screen_height);
                % tex19 = Screen('MakeTexture', w, jdnoise19);
                % 
                % jdnoise20 = imread(fullfile(folder,'JDnoise20.bmp'));
                % [jdnoise20] = cropAndMakeBkgLum(jdnoise20,wholeLumRasterCancelJDset, screen_width, screen_height);
                % tex20 = Screen('MakeTexture', w, jdnoise20);
                % 
                % whiteblank99 = imread(fullfile(folder,'JDnoise20.bmp'));
                % whiteblank99(:,:,:)=255;
                % [whiteblank99] = cropAndMakeBkgLum(whiteblank99,wholeLumRasterCancelJDset, screen_width, screen_height);
                % tex99 = Screen('MakeTexture', w, whiteblank99);

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
            

            function [mx] = setCowBkgdwithCross(graph, mx, my, this, tex1oF, insideColor)

                blackC = BlackIndex(graph.window);
                whiteC = WhiteIndex(graph.window);

                if this.blankblack == 1
                    Screen('FillRect', graph.window, blackC);
                    
                    Screen('FillRect', graph.window ...
                            , insideColor,... %these numbers are for 4.2e-01 luminance
                            [1,...
                            my - this.insideCircle_half_length,...
                            mx + mx,...
                            my + this.insideCircle_half_length]); %2deg box
                else
                    Screen('FillRect', graph.window, blackC);
                    destinationRectangle =  []; % so that it fills the whole screen
                    Screen('DrawTexture', graph.window, tex1oF, [],destinationRectangle);

                    if this.turnOnRing == 1
                        blackRing_rgb = [1 1 1]; %
                        % Drawing the black circle first so that it becomes ring
                        Screen('FillRect',graph.window,whiteC,... %black oval
                            [1,...
                            my - this.blackRing_half_length,...
                            mx + mx,...
                            my + this.blackRing_half_length]);

                        % Now drawing inside white circle
                        Screen('FillRect', graph.window ...
                            , insideColor,... %these numbers are for 4.2e-01 luminance
                            [1,...
                            my - this.insideCircle_half_length,...
                            mx + mx,...
                            my + this.insideCircle_half_length]); %2deg box
                    end
                end
            end

            % Make 20 textures
            function [tex1, tex2, tex3, tex4, tex5, tex6, tex7, tex8, tex9, tex10] = ...
                    makeWOTextures(screen_width, screen_height, w)
                %cd('/Users/Josephine1/Documents/MATLAB/Roorda Lab/Data_and_AnalysisCode/gifMakers_presentations/JDs20noiseBMPS_folder'); % I know this is redundant, but didn't know a better way JD 7/21/23

                folder = fullfile(fileparts(mfilename('fullpath')), 'Masks');

                mask1 = imread(fullfile(folder,'voronoi-v4-1.png'));
                [mask1] = resizeAndCropWO(mask1, screen_width, screen_height);
                tex1 = Screen('MakeTexture', w, mask1);

                mask2 = imread(fullfile(folder,'voronoi-v4-2.png'));
                [mask2] = resizeAndCropWO(mask2, screen_width, screen_height);
                tex2 = Screen('MakeTexture', w, mask2);

                mask3 = imread(fullfile(folder,'voronoi-v4-3.png'));
                [mask3] = resizeAndCropWO(mask3, screen_width, screen_height);
                tex3 = Screen('MakeTexture', w, mask3);

                mask4 = imread(fullfile(folder,'voronoi-v4-4.png'));
                [mask4] = resizeAndCropWO(mask4, screen_width, screen_height);
                tex4 = Screen('MakeTexture', w, mask4);

                mask5 = imread(fullfile(folder,'voronoi-v4-5.png'));
                [mask5] = resizeAndCropWO(mask5, screen_width, screen_height);
                tex5 = Screen('MakeTexture', w, mask5);

                mask6 = imread(fullfile(folder,'voronoi-v4-6.png'));
                [mask6] = resizeAndCropWO(mask6, screen_width, screen_height);
                tex6 = Screen('MakeTexture', w, mask6);

                mask7 = imread(fullfile(folder,'voronoi-v4-7.png'));
                [mask7] = resizeAndCropWO(mask7, screen_width, screen_height);
                tex7 = Screen('MakeTexture', w, mask7);

                mask8 = imread(fullfile(folder,'voronoi-v4-8.png'));
                [mask8] = resizeAndCropWO(mask8, screen_width, screen_height);
                tex8 = Screen('MakeTexture', w, mask8);

                mask9 = imread(fullfile(folder,'voronoi-v4-9.png'));
                [mask9] = resizeAndCropWO(mask9, screen_width, screen_height);
                tex9 = Screen('MakeTexture', w, mask9);

                mask10 = imread(fullfile(folder,'voronoi-v4-10.png'));
                [mask10] = resizeAndCropWO(mask10, screen_width, screen_height);
                tex10 = Screen('MakeTexture', w, mask10);
            end

            function [cropped] = resizeAndCropWO(img, screen_width, screen_height)
                % targetW = graph.pxWidth;
                % targetH = graph.pxHeight;

                [h, w, ~] = size(img);

                % Compute scale factor so image fully covers target
                scale = max(screen_width / w, screen_height / h);

                % Resize
                resized = imresize(img, scale);

                % Get new size
                [newH, newW, ~] = size(resized);

                % Compute center crop coordinates
                xStart = floor((newW - screen_width) / 2) + 1;
                yStart = floor((newH - screen_height) / 2) + 1;

                cropped = resized(yStart:yStart+screen_height-1, ...
                    xStart:xStart+screen_width-1, :);

                % maskTex = Screen('MakeTexture', graph.window, cropped);
            end

        end %runTrial end


    end

    % ---------------------------------------------------------------------
    % Plot methods
    % ---------------------------------------------------------------------
    methods ( Access = public )

        function Plot_EyeData(this)

            figure

            hold
            for i=1:height(this.TrialTable)
                WhichTrial = i;

                samplesDataTable = this.Session.samplesDataTable;
                trialDataTable = this.Session.trialDataTable;

                samplesInTrial = find(samplesDataTable.TrialNumber==WhichTrial);

                % time point in openiris time and trial psyhchtoolbox time that is common
                % is samplesDataTable.Time(samplesInTrial(1)) and
                % trialDataTable(WhichTrial,:).TimeTrialStart (NOTE there is still some
                % latency betwen open iris and psychtoolbox that we don't measusre order of
                % magnitude probably 10 to 100ms) but will depend on how busy the open iris
                % computer is.

                t = samplesDataTable.Time(samplesInTrial) - samplesDataTable.Time(samplesInTrial(1)) - (trialDataTable(WhichTrial,:).TimeStartLoop-trialDataTable(WhichTrial,:).TimeTrialStart);

                plot(t+ samplesDataTable.Time(samplesInTrial(1)) ,samplesDataTable.LeftX(samplesInTrial));
                line([0.5 0.5], [-15 15])
                line([1.5 1.5], [-15 15])
                line([3.5 3.5], [-15 15])

                timeWindow = find(t>1.8 & t<3.3);
                mdl = fitlm(t(timeWindow), samplesDataTable.LeftX(samplesInTrial(timeWindow)));
                
                plot(t(timeWindow)+ samplesDataTable.Time(samplesInTrial(1)) , mdl.predict(t(timeWindow)),'r' )

                % plot(t, x0 + (t-t0)*v)
                plot(t(timeWindow) + samplesDataTable.Time(samplesInTrial(1)) , samplesDataTable.LeftX(samplesInTrial(timeWindow(1))) + (t(timeWindow)-t(timeWindow(1)))*trialDataTable.Speed_middle_pix(WhichTrial)/trialDataTable.ppd_x(WhichTrial),'g')


                plot(t+ samplesDataTable.Time(samplesInTrial(1)) , samplesDataTable.QuickPhase(samplesInTrial)*1000)
            end
        end

        function Plot_Psychometric(this)
            folder = fullfile(fileparts(mfilename('fullpath')), 'Palamedes');
            addpath(folder);

            trialDataTable = this.Session.currentRun.pastTrialTable;

            %convertings pixels to degrees
            trialDataTable.Speed_TEST_deg = trialDataTable.Speed_TEST_pix./trialDataTable.ppd_x;
            trialDataTable.Speed_middle_deg = trialDataTable.Speed_middle_pix./trialDataTable.ppd_x;
            trialDataTable.Speed_ref_deg = trialDataTable.Speed_ref_pix./trialDataTable.ppd_x;

            %getting rate of gap increasing
            trialDataTable.Rate_test = abs(trialDataTable.Speed_middle_deg - trialDataTable.Speed_TEST_deg);
            trialDataTable.Rate_ref = abs(trialDataTable.Speed_ref_deg - trialDataTable.Speed_middle_deg);

            numConditions = unique(trialDataTable.RefIsLead);

            T = table();
            for nc = 1: size(numConditions,1)
                curT = trialDataTable(trialDataTable.RefIsLead == numConditions(nc),:);

                if numConditions(nc) == 1
                    curT.ChoseRef(curT.StartingFromLeft ==1 & curT.Response == 'L' & curT.GrowingOutwards == 1) = 1;
                    curT.ChoseRef(curT.StartingFromLeft ==0 & curT.Response == 'R' & curT.GrowingOutwards == 1) = 1;

                    %I added for table shrinking...
                    curT.ChoseRef(curT.StartingFromLeft ==1 & curT.Response == 'R' & curT.GrowingOutwards == 0) = 1;
                    curT.ChoseRef(curT.StartingFromLeft ==0 & curT.Response == 'L' & curT.GrowingOutwards == 0) = 1;

                    curT.RateLag = abs(curT.Speed_middle_deg) - abs(curT.Speed_TEST_deg);
                    curT.RateLead = abs(curT.Speed_ref_deg) - abs(curT.Speed_middle_deg);

                    curT.RateRefFinal = curT.RateLead - curT.RateLead(1);
                    curT.RateTestFinal = curT.RateLag - curT.RateLead(1);
                    curT.RateTestFinal = curT.RateTestFinal *-1; %to flip vertical axis

                    T = [curT ;T];

                else
                    curT.ChoseRef(curT.StartingFromLeft ==1 & curT.Response == 'L' & curT.GrowingOutwards == 1) = 1;
                    curT.ChoseRef(curT.StartingFromLeft ==0 & curT.Response == 'R' & curT.GrowingOutwards == 1) = 1;

                    %I added for table shrinking...
                    curT.ChoseRef(curT.StartingFromLeft ==1 & curT.Response == 'R' & curT.GrowingOutwards == 0) = 1;
                    curT.ChoseRef(curT.StartingFromLeft ==0 & curT.Response == 'L' & curT.GrowingOutwards == 0) = 1;

                    curT.RateLag = abs(curT.Speed_middle_deg) - abs(curT.Speed_TEST_deg);
                    curT.RateLead = abs(curT.Speed_ref_deg) - abs(curT.Speed_middle_deg);

                    curT.RateRefFinal = curT.RateLag - curT.RateLag(1);
                    curT.RateTestFinal = curT.RateLead - curT.RateLag(1);

                    T = [curT ;T];
                end
            end

            curT = T;
            curT.RateTestFinal =  round(curT.RateTestFinal,4);
            [uniqueVals, ~, idx] = unique(curT.RateTestFinal);
            % marker_sizes = accumarray(idx,1);

            StimLevels = nan(1,size(uniqueVals,1)); NumPos = nan(1,size(uniqueVals,1)); OutOfNum = nan(1,size(uniqueVals,1));
            for i = 1:size(uniqueVals,1)
                curData = curT(curT.RateTestFinal==uniqueVals(i), :);

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
            if curT.GrowingOutwards(1) == 0
                [paramsValues] = PAL_PFML_Fit(StimLevels,NumPos,OutOfNum,searchGrid,paramsFree,PF);
            else
                NumPos_flipped = OutOfNum - NumPos;
                [paramsValues] = PAL_PFML_Fit(StimLevels,NumPos_flipped,OutOfNum,searchGrid,paramsFree,PF);
            end

            %plotting fit
            ProportionCorrect = NumPos./OutOfNum;
            xEval = [min(StimLevels):max(abs(StimLevels))./1000:max(StimLevels)];
            yEval = PF(paramsValues,xEval);
            if curT.GrowingOutwards(1) == 1
                yEval = 1 - yEval;
            end

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
            % if trialDataTable.NOPursueExp ==0
            %     if numConditions(nc) == 1
            %         xline(abs(curT.Speed_middle_deg(1))-curT.Rate_ref(1), '--','Color',[0.5 0.5 0.5], 'LineWidth',2);
            %     else
            %         xline(abs(curT.Speed_middle_deg(1))+curT.Rate_test(1), '--','Color',[0.5 0.5 0.5], 'LineWidth',2);
            %     end
            % end
            xline(T.RateRefFinal, '--','Color',[0.5 0.5 0.5], 'LineWidth',2);

            ylim([0 1]); yticks([0:0.1:1]);

            % if numConditions(nc) == 1
            %     xlabel('Lagging stimulus speed (deg/s)')
            %     ylabel('Reponse == Leading Stimulus')
            % else
                xlabel('Test stimulus rate of lag or lead (deg/s)')
                ylabel('Reponse == Reference')
            % end

            set(gca,'fontsize', 18)
            title(strrep(this.Session.shortName,'_',' '))
            saveFolder = this.Session.folder;
            saveas(gcf, [saveFolder, '/', this.Session.subjectCode, '_', this.Session.sessionCode '_PF_','RefIsLead_',num2str(numConditions(nc))','_indiv.jpg'])


        end

        function Plot_Psychometric_2graphs(this)
            folder = fullfile(fileparts(mfilename('fullpath')), 'Palamedes');
            addpath(folder);

            trialDataTable = this.Session.currentRun.pastTrialTable;

            %convertings pixels to degrees
            trialDataTable.Speed_TEST_deg = trialDataTable.Speed_TEST_pix./trialDataTable.ppd_x;
            trialDataTable.Speed_middle_deg = trialDataTable.Speed_middle_pix./trialDataTable.ppd_x;
            trialDataTable.Speed_ref_deg = trialDataTable.Speed_ref_pix./trialDataTable.ppd_x;

            %getting rate of gap increasing
            trialDataTable.Rate_test = abs(trialDataTable.Speed_middle_deg - trialDataTable.Speed_TEST_deg);
            trialDataTable.Rate_ref = abs(trialDataTable.Speed_ref_deg - trialDataTable.Speed_middle_deg);

            numConditions = unique(trialDataTable.RefIsLead);


            for nc = 1: size(numConditions,1)
                curT = trialDataTable(trialDataTable.RefIsLead == numConditions(nc),:);

                if numConditions(nc) == 1
                    curT.ChoseRef(curT.StartingFromLeft ==1 & curT.Response == 'R' & curT.GrowingOutwards == 1) = 1;
                    curT.ChoseRef(curT.StartingFromLeft ==0 & curT.Response == 'L' & curT.GrowingOutwards == 1) = 1;

                    %I added for table shrinking...
                    curT.ChoseRef(curT.StartingFromLeft ==1 & curT.Response == 'L' & curT.GrowingOutwards == 0) = 1;
                    curT.ChoseRef(curT.StartingFromLeft ==0 & curT.Response == 'R' & curT.GrowingOutwards == 0) = 1;

                    [uniqueVals, ~, idx] = unique(abs(curT.Speed_TEST_deg));
                    % marker_sizes = accumarray(idx,1);

                    StimLevels = nan(1,size(uniqueVals,1)); NumPos = nan(1,size(uniqueVals,1)); OutOfNum = nan(1,size(uniqueVals,1));
                    for i = 1:size(uniqueVals,1)
                        curData = curT(abs(curT.Speed_TEST_deg) ==uniqueVals(i), :);

                        StimLevels(i) = uniqueVals(i);
                        NumPos(i) = sum(curData.ChoseRef);
                        OutOfNum(i) = length(curData.ChoseRef);
                    end
                else
                    curT.ChoseRef(curT.StartingFromLeft ==1 & curT.Response == 'L' & curT.GrowingOutwards == 1) = 1;
                    curT.ChoseRef(curT.StartingFromLeft ==0 & curT.Response == 'R' & curT.GrowingOutwards == 1) = 1;

                    %I added for table shrinking...
                    curT.ChoseRef(curT.StartingFromLeft ==1 & curT.Response == 'R' & curT.GrowingOutwards == 0) = 1;
                    curT.ChoseRef(curT.StartingFromLeft ==0 & curT.Response == 'L' & curT.GrowingOutwards == 0) = 1;

                    [uniqueVals, ~, idx] = unique(abs(curT.Speed_ref_deg));
                    % marker_sizes = accumarray(idx,1);

                    StimLevels = nan(1,size(uniqueVals,1)); NumPos = nan(1,size(uniqueVals,1)); OutOfNum = nan(1,size(uniqueVals,1));
                    for i = 1:size(uniqueVals,1)
                        curData = curT(abs(curT.Speed_ref_deg) ==uniqueVals(i), :);

                        StimLevels(i) = uniqueVals(i);
                        NumPos(i) = sum(curData.ChoseRef);
                        OutOfNum(i) = length(curData.ChoseRef);
                    end
                end

                PF = @PAL_Logistic;

                paramsFree = [1 1 0 0];
                searchGrid.alpha = mean(StimLevels);
                searchGrid.beta = 3.5;
                searchGrid.gamma = 0.01;
                searchGrid.lambda = 0.01;


                %Perform fit
                if (numConditions(nc) == 1 && curT.GrowingOutwards(1) == 1) || (numConditions(nc) == 0 && curT.GrowingOutwards(1) == 0)
                    [paramsValues] = PAL_PFML_Fit(StimLevels,NumPos,OutOfNum,searchGrid,paramsFree,PF);
                else
                    NumPos_flipped = OutOfNum - NumPos;
                    [paramsValues] = PAL_PFML_Fit(StimLevels,NumPos_flipped,OutOfNum,searchGrid,paramsFree,PF);
                end

                %plotting fit
                ProportionCorrect = NumPos./OutOfNum;
                xEval = [min(StimLevels):max(abs(StimLevels))./1000:max(StimLevels)];
                yEval = PF(paramsValues,xEval);
                if (numConditions(nc) == 0 && curT.GrowingOutwards(1) == 1) || (numConditions(nc) == 1 && curT.GrowingOutwards(1) == 0)
                    yEval = 1 - yEval;
                end

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
                if trialDataTable.NOPursueExp ==0
                    if numConditions(nc) == 1
                        xline(abs(curT.Speed_middle_deg(1))-curT.Rate_ref(1), '--','Color',[0.5 0.5 0.5], 'LineWidth',2);
                    else
                        xline(abs(curT.Speed_middle_deg(1))+curT.Rate_test(1), '--','Color',[0.5 0.5 0.5], 'LineWidth',2);
                    end
                end
                ylim([0 1]); yticks([0:0.1:1]);

                if numConditions(nc) == 1
                    xlabel('Lagging stimulus speed (deg/s)')
                    ylabel('Reponse == Leading Stimulus')
                else
                    xlabel('Leading stimulus speed (deg/s)')
                    ylabel('Reponse == Lagging Stimulus')
                end

                set(gca,'fontsize', 18)
                title(strrep(this.Session.shortName,'_',' '))
                saveFolder = this.Session.folder;
                saveas(gcf, [saveFolder, '/', this.Session.subjectCode, '_', this.Session.sessionCode '_PF_','RefIsLead_',num2str(numConditions(nc))','.jpg'])
            end

            % figure
            % trialDataTable = this.Session.trialDataTable;
            % plot(trialDataTable.Speed_left, double(trialDataTable.Response=='R')+randn(size(trialDataTable.Response))/20,'o');
            % xlabel('Speed')
            % ylabel('Reponse == Right')
        end

    end


end %classdef end