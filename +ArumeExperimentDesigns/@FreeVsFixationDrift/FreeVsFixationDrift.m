classdef FreeVsFixationDrift < ArumeExperimentDesigns.EyeTracking
    %OPTOKINETICTORSION Summary of this class goes here
    %   Detailed explanation goes here

    properties
        fixRad = 20;
        fixColor = [255 0 0];
        targetPositions =[];
        stimTexture = [];
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
            dlg.StimulusContrast0to100 = {60 '* (%)' [0 100] };
            dlg.StimSizeDeg = {15 '* (deg)' [1 100] };

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

            % %% trial table for fixation targets
            % t = ArumeCore.TrialTableBuilder();
            % 
            % t.AddConditionVariable("TargetPosition", { ...
            %     [0 0], [0 10], [10 0], [10 10], [-10 10], [-10 -10], [10 -10] [-10 0], [0 -10] ...
            %     [0 2], [0 4], [0 6], [0 8], [0 -2], [0 -4], [0 -6], [0 -8], ...
            %     [2 0], [4 0], [6 0], [8 0], [-2 0], [-4 0], [-6 0], [-8 0], ...
            %     [2 2], [4 4], [6 6], [8 8], [-2 2], [-4 4], [-6 6], [-8 8], ...
            %     [2 -2], [4 -4], [6 -6], [8 -8], [-2 -2], [-4 -4], [-6 -6], [-8 -8] ...
            %     });
            % 
            % t.AddBlock(1:height(t.ConditionTable), nReps);
            % 
            % trialSequence = 'Random';
            % blockSequence = 'Sequential';
            % blockSequenceRepeatitions = 1;
            % abortAction = 'Repeat';
            % trialsPerSession = 1000;  % You can use a large number if not splitting
            % 
            % trialTable = t.GenerateTrialTable(trialSequence, blockSequence, blockSequenceRepeatitions, abortAction, trialsPerSession);
            
            %% trial table for freeviewing images
            t = ArumeCore.TrialTableBuilder();
            nReps = this.ExperimentOptions.NumberRepetitions;

            t.AddConditionVariable("Image", {'01', '02', '03', '04', '05', '06', '07', '08', '09', '10'}); %% currently allows for ten images, numbered 01-10. add more if needed

            t.AddBlock(1:height(t.ConditionTable), nReps); %% same number of repetitions for the fixation targets and the images

            trialSequence = 'Random';
            blockSequence = 'Sequential';
            blockSequenceRepeatitions = 1;
            abortAction = 'Repeat';
            trialsPerSession = 1000;
            trialTable = t.GenerateTrialTable(trialSequence, blockSequence, blockSequenceRepeatitions, abortAction, trialsPerSession);
        end

        function [trialResult, thisTrialData] = runPreTrial( this, thisTrialData )
            Enum = ArumeCore.ExperimentDesign.getEnum();
            trialResult = Enum.trialResult.CORRECT;
            
            % JORGE AT THE MEETING
            %experimentFolder = fileparts(mfilename('fullpath'));
            %imageFile = fullfile(experimentFolder,[thisTrialData.Image '.jpg']);
            % END JORGE
            test = string(thisTrialData.Image);
            imageFile = fullfile(fileparts(mfilename('fullpath')),[test + ".jpeg"]);
            
            I = imread(imageFile);
                
            monitorWidthPix     = this.Graph.wRect(3);
            monitorWidthCm      = this.ExperimentOptions.DisplayOptions.ScreenWidth;
            monitorDistanceCm   = this.ExperimentOptions.DisplayOptions.ScreenDistance;
            stimSizeDeg         = this.ExperimentOptions.StimSizeDeg;

            % we will asume that pixels are square
            monitorWidthDeg     = 2*atand(monitorWidthCm/monitorDistanceCm/2);
            % asuming linearity (not completely true for very large displays
            %             pixelsPerDeg        = monitorWidthPix/monitorWidthDeg;
            %             stimSizePix         = pixelsPerDeg * stimSizeDeg;

            % non linear aproximation
            stimSizeCm  = 2*tand(stimSizeDeg/2)*monitorDistanceCm
            %stimSizePix = stimSizeCm/monitorWidthCm*monitorWidthPix;
            stimSizePix = (monitorWidthPix/monitorWidthCm)*stimSizeCm

            Isquare = uint8(double(I(:,(size(I,2) - size(I,1))/2+(1:(size(I,1))),:,:))*this.ExperimentOptions.StimulusContrast0to100/100);
            Isquare = imresize(Isquare, [stimSizePix stimSizePix], 'bilinear');
            this.stimTexture = Screen('MakeTexture', this.Graph.window, Isquare);
            
        end

        function [trialResult, thisTrialData] = runTrial( this, thisTrialData )
           
            Enum = ArumeCore.ExperimentDesign.getEnum();
            graph = this.Graph;
            
            trialDuration = this.ExperimentOptions.TrialDuration;
            
            %-- add here the trial code
            Screen('FillRect', graph.window, 128);
            
            
            lastFlipTime                        = Screen('Flip', graph.window);
            secondsRemaining                    = trialDuration;
            thisTrialData.TimeStartLoop         = lastFlipTime;
            if ( ~isempty(this.eyeTracker) )
                thisTrialData.EyeTrackerFrameStartLoop = this.eyeTracker.RecordEvent(sprintf('TRIAL_START_LOOP %d %d', thisTrialData.TrialNumber, thisTrialData.Condition) );
            end
            while secondsRemaining > 0
                
                secondsElapsed      = GetSecs - thisTrialData.TimeStartLoop;
                secondsRemaining    = trialDuration - secondsElapsed;
                
                % -----------------------------------------------------------------
                % --- Drawing of stimulus -----------------------------------------
                % -----------------------------------------------------------------


                %-- Find the center of the screen
                [mx, my] = RectCenter(graph.wRect);
                fixRect = [0 0 10 10];
                fixRect = CenterRectOnPointd( fixRect, mx, my );

                Screen('DrawTexture', this.Graph.window, this.stimTexture, [],[]);

                this.Graph.Flip(this, thisTrialData, secondsRemaining);
                % -----------------------------------------------------------------
                % --- END Drawing of stimulus -------------------------------------
                % -----------------------------------------------------------------
                
                
            end
            
            trialResult = Enum.trialResult.CORRECT;

        end
            
    end
end