classdef OptokineticTorsion < ArumeExperimentDesigns.EyeTracking
    %OPTOKINETICTORSION Summary of this class goes here
    %   Detailed explanation goes here
    
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
            
            
            dlg.Display.ScreenWidth = { 121 '* (cm)' [1 3000] };
            dlg.Display.ScreenHeight = { 68 '* (cm)' [1 3000] };
            dlg.Display.ScreenDistance = { 60 '* (cm)' [1 3000] };
            
            dlg.Trial_Duration =  { 30 '* (s)' [1 100] };
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
            
            %% CHANGE DEFAULTS 
            dlg.UseEyeTracker = 0;
            dlg.Debug.DisplayVariableSelection = 'TrialNumber TrialResult Speed Stimulus'; % which variables to display every trial in the command line separated by spaces
            
                dlg.HitKeyBeforeTrial = 0;
                dlg.TrialDuration = 10;
                dlg.TrialsBeforeBreak = 15;
                dlg.TrialAbortAction = 'Repeat';
        end
        
        function trialTable = SetUpTrialTable(this)
            
            %-- condition variables ---------------------------------------
            i= 0;
            
            i = i+1;
            conditionVars(i).name   = 'Speed';
            conditionVars(i).values = this.ExperimentOptions.Max_Speed/this.ExperimentOptions.Number_of_Speeds * [0:this.ExperimentOptions.Number_of_Speeds];
            
            i = i+1;
            conditionVars(i).name   = 'Direction';
            conditionVars(i).values = {'CW' 'CCW'};
            
            i = i+1;
            conditionVars(i).name   = 'Stimulus';
            if ( this.ExperimentOptions.Do_Blank ) 
                conditionVars(i).values = {'Blank' 'Dots'};
            else
                conditionVars(i).values = {'Dots'};
            end
            
            trialTableOptions = this.GetDefaultTrialTableOptions();
            trialTableOptions.trialSequence = 'Random';
            trialTableOptions.trialAbortAction = 'Delay';
            trialTableOptions.trialsPerSession = 1000;
            trialTableOptions.numberOfTimesRepeatBlockSequence = this.ExperimentOptions.NumberOfRepetitions;
            trialTable = this.GetTrialTableFromConditions(conditionVars, trialTableOptions);
        end
        
        
        function [trialResult, thisTrialData] = runTrial( this, thisTrialData )
            
            try
                
                Enum = ArumeCore.ExperimentDesign.getEnum();
                graph = this.Graph;
                trialResult = Enum.trialResult.CORRECT;
                
                
                lastFlipTime        = GetSecs;
                secondsRemaining    = this.ExperimentOptions.TrialDuration;
                thisTrialData.TimeStartLoop = lastFlipTime;
                
                % prepare dots
                
                mon_width   = this.Graph.pxWidth/this.ExperimentOptions.DisplayOptions.ScreenWidth;   % horizontal dimension of viewable screen (cm)
                v_dist      = this.ExperimentOptions.DisplayOptions.ScreenDistance ;   % viewing distance (cm)
                
                dot_speed   = thisTrialData.Speed;    % dot speed (deg/sec)
                ndots       = this.ExperimentOptions.Number_of_Dots; % number of dots
                
                max_d       = this.ExperimentOptions.Max_Radius;   % maximum radius of  annulus (degrees)
                min_d       = this.ExperimentOptions.Min_Radius;    % minumum
                
                differentsizes = this.ExperimentOptions.Number_of_Dot_Sizes; % Use different sizes for each point if >= 1. Use one common size if == 0.
                dot_w       = this.ExperimentOptions.Min_Dot_Diam;  % width of dot (deg)
                                
                ppd = pi * (this.Graph.wRect(3)-this.Graph.wRect(1)) / atan(mon_width/v_dist/2) / 360;    % pixels per degree
                s = dot_w * ppd;                                        % dot size (pixels)
                
                rmax = max_d * ppd;	% maximum radius of annulus (pixels from center)
                rmin = min_d * ppd; % minimum
                r = rmin + (rmax-rmin) * sqrt(rand(ndots,1));	% r
                t = 2*pi*rand(ndots,1);                         % theta polar coordinate
                
                cs = [cos(t), sin(t)];
                xy = [r r] .* cs;   % dot positions in Cartesian coordinates (pixels from center)
                xymatrix = transpose(xy);
                 
                mdir = 1;    % motion direction (in or out) for each dot
                if ( thisTrialData.Direction == 'CCW' )
                    mdir = -1;
                end
                dt = 2*pi/360/this.Graph.frameRate * thisTrialData.Speed*mdir;                       % change in theta per frame (radians)
                % Create a vector with different point sizes for each single dot, if
                % requested:
                if (differentsizes>0)
                    s=(1+rand(1, ndots)*(differentsizes-1))*s;
                end
                [center(1), center(2)] = RectCenter(this.Graph.wRect);
                
                % end prepare dots
                
                if ( ~isempty(this.eyeTracker) )
                    thisTrialData.EyeTrackerFrameStartLoop = this.eyeTracker.RecordEvent(sprintf('TRIAL_START_LOOP %d %d', thisTrialData.TrialNumber, thisTrialData.Condition) );
                end
                
                while secondsRemaining > 0
                    
                    secondsElapsed      = GetSecs - thisTrialData.TimeStartLoop;
                    secondsRemaining    = this.ExperimentOptions.TrialDuration - secondsElapsed;
                    
                    
                    % -----------------------------------------------------------------
                    % --- Drawing of stimulus -----------------------------------------
                    % -----------------------------------------------------------------
                    switch(thisTrialData.Stimulus)
                        case 'Dots'
                             s( s>63) = 63;
                            Screen('DrawDots', graph.window, xymatrix, s, WhiteIndex(graph.window), center,1);  % change 1 to 0 to draw square dots
                    end
                    
                    t = t + dt;                         % update theta
                    xy = [r r] .* [cos(t), sin(t)];     % compute new positions
                    xymatrix = transpose(xy);
                    
                    
                    
                    %-- Draw fixation spot
                    [mx, my] = RectCenter(this.Graph.wRect);
                    
                    targetPix = this.Graph.pxWidth/this.ExperimentOptions.DisplayOptions.ScreenWidth * this.ExperimentOptions.DisplayOptions.ScreenDistance * tan(this.ExperimentOptions.TargetSize/180*pi);
                    fixRect = [0 0 targetPix targetPix];
                    fixRect = CenterRectOnPointd( fixRect, mx, my );
                    Screen('FillOval', graph.window, this.fixColor, fixRect);
                    
                    Screen('DrawingFinished', graph.window); % Tell PTB that no further drawing commands will follow before Screen('Flip')
                    
                    
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
    
end