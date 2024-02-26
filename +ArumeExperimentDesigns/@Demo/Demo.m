classdef Demo < ArumeExperimentDesigns.EyeTracking
    % DEMO experiment for Arume 
    %
    %   1. Copy paste the folder @Demo within +ArumeExperimentDesigns.
    %   2. Rename the folder with the name of the new experiment but keep that @ at the begining!
    %   3. Rename also the file inside to match the name of the folder (without the @ this time).
    %   4. Then change the name of the class inside the folder.
    %     
    properties
        fixRad = 20;
        fixColor = [255 0 0];
    end
    
    % ---------------------------------------------------------------------
    % Experiment design methods
    % ---------------------------------------------------------------------
    methods ( Access = protected )
        function dlg = GetOptionsDialog( this, importing )
            dlg = GetOptionsDialog@ArumeExperimentDesigns.EyeTracking(this, importing);
            
            %% ADD new options
            dlg.Max_Speed = { 30 '* (deg/s)' [0 100] };
            dlg.Number_of_Speeds = {3 '* (N)' [1 100] };
            
            dlg.Number_of_Dots = { 2000 '* (deg/s)' [10 10000] };
            dlg.Max_Radius = { 40 '* (deg)' [1 100] };
            dlg.Min_Radius = { 1 '* (deg)' [0 100] };

            dlg.Min_Dot_Diam = {0.1  '* (deg)' [0.01 100] };
            dlg.Max_Dot_Diam = {0.4  '* (deg)' [0.01 100] };
            dlg.Number_of_Dot_Sizes = {5 '* (N)' [1 100] };
            
            dlg.NumberOfRepetitions = {8 '* (N)' [1 100] };
            
            dlg.Do_Blank = { {'0','{1}'} };
            
            dlg.TargetSize = 0.5;
            
            dlg.BackgroundBrightness = 0;
            
            %% CHANGE DEFAULTS values for existing options

            dlg.UseEyeTracker = 0;
            dlg.Debug.DisplayVariableSelection = 'TrialNumber TrialResult Speed Stimulus'; % which variables to display every trial in the command line separated by spaces

            dlg.DisplayOptions.ScreenWidth = { 121 '* (cm)' [1 3000] };
            dlg.DisplayOptions.ScreenHeight = { 68 '* (cm)' [1 3000] };
            dlg.DisplayOptions.ScreenDistance = { 60 '* (cm)' [1 3000] };
            
            dlg.HitKeyBeforeTrial = 0;
            dlg.TrialDuration = 10;
            dlg.TrialsBeforeBreak = 15;
            dlg.TrialAbortAction = 'Repeat';
        end
        
        function trialTable = SetUpTrialTable(this)
            
            t = ArumeCore.TrialTableBuilder();

            t.AddConditionVariable('Speed',this.ExperimentOptions.Max_Speed/this.ExperimentOptions.Number_of_Speeds * [0:this.ExperimentOptions.Number_of_Speeds]);
            t.AddConditionVariable('Direction',{'CW' 'CCW'})
            if ( this.ExperimentOptions.Do_Blank )
                t.AddConditionVariable('Stimulus',{'Blank' 'Dots'});
            else
                t.AddConditionVariable('Stimulus',{'Dots'});
            end

            trialTable = t.GenerateTrialTable();
        end
        

        function [trialResult, thisTrialData] = runTrial( this, thisTrialData )

            try


                Enum = ArumeCore.ExperimentDesign.getEnum();
                graph = this.Graph;
                trialResult = Enum.trialResult.CORRECT;


                lastFlipTime        = GetSecs;
                secondsRemaining    = this.ExperimentOptions.TrialDuration;
                thisTrialData.TimeStartLoop = lastFlipTime;


                while secondsRemaining > 0

                    secondsElapsed      = GetSecs - thisTrialData.TimeStartLoop;
                    secondsRemaining    = this.ExperimentOptions.TrialDuration - secondsElapsed;


                    % -----------------------------------------------------------------
                    % --- Drawing of stimulus -----------------------------------------
                    % -----------------------------------------------------------------


                    %-- Find the center of the screen
                    [mx, my] = RectCenter(graph.wRect);

                    fixRect = [0 0 10 10];
                    fixRect = CenterRectOnPointd( fixRect, mx, my );
                    Screen('FillOval', graph.window,  this.fixColor, fixRect);


                    % -----------------------------------------------------------------
                    % --- END Drawing of stimulus -------------------------------------
                    % -----------------------------------------------------------------
                    
                    % -----------------------------------------------------------------
                    % -- Flip buffers to refresh screen -------------------------------
                    % -----------------------------------------------------------------
                    this.Graph.Flip(this, thisTrialData, secondsRemaining);
                    % -----------------------------------------------------------------


                    % -----------------------------------------------------------------
                    % --- Collecting responses  ---------------------------------------
                    % -----------------------------------------------------------------

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
                    if ( ~isempty( response) )
                        thisTrialData.Response = response;
                        thisTrialData.ResponseTime = GetSecs;

                        break;
                    end
                    % -----------------------------------------------------------------
                    % --- END Collecting responses  -----------------------------------
                    % -----------------------------------------------------------------


                end
            catch ex
                rethrow(ex)
            end
            
        end        
    end
    
end