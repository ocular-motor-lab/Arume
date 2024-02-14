classdef Stereoacuity_MethodOfConstantStimuli < ArumeExperimentDesigns.EyeTracking
    % DEMO experiment for Arume
    %
    %   1. Copy paste the folder @Demo within +ArumeExperimentDesigns.
    %   2. Rename the folder with the name of the new experiment but keep that @ at the begining!
    %   3. Rename also the file inside to match the name of the folder (without the @ this time).
    %   4. Then change the name of the class inside the folder.
    %
    properties
        stimTextureLeft = [];
        stimTextureRight = [];
        targetColor = [150 150 150];
    end
    
    % ---------------------------------------------------------------------
    % Experiment design methods
    % ---------------------------------------------------------------------
    methods ( Access = protected )
        function dlg = GetOptionsDialog( this, importing )
            dlg = GetOptionsDialog@ArumeExperimentDesigns.EyeTracking(this, importing);
            
            %% ADD new options
            dlg.Number_of_Dots = { 750 '* (num)' [10 10000] }; %750
            dlg.Size_of_Dots = { 1 '* (pix)' [1 100] };
            dlg.MaxStimDeg = {3.5 '* (deg)' [1 100] }; %3.5
            dlg.MinStimDeg = {2 '* (deg)' [1 100] };
            dlg.FixationSpotSize = { 0.25 '* (diameter_in_deg)' [0 5] };
            dlg.TimeStimOn = { 0.2 '* (sec)' [0 60] }; 
            dlg.InitFixDuration = { 1 '* (sec)' [0 60] };
            
            dlg.NumberOfRepetitions = {20 '* (N)' [1 200] }; 
            dlg.BackgroundBrightness = 0;
            
            %% CHANGE DEFAULTS values for existing options
            
            dlg.UseEyeTracker = 0;
            dlg.Debug.DisplayVariableSelection = 'TrialNumber TrialResult RotateDots DisparityArcMin GuessedCorrectly'; % which variables to display every trial in the command line separated by spaces
            
            dlg.DisplayOptions.ScreenWidth = { 59.5 '* (cm)' [1 3000] };
            dlg.DisplayOptions.ScreenHeight = { 34 '* (cm)' [1 3000] };
            dlg.DisplayOptions.ScreenDistance = { 57 '* (cm)' [1 3000] };
            dlg.DisplayOptions.StereoMode = { 4 '* (mode)' [0 9] }; 
            dlg.DisplayOptions.SelectedScreen = { 1 '* (screen)' [0 5] };
            
            dlg.HitKeyBeforeTrial = 0;
            dlg.TrialDuration = 90;
            dlg.TrialsBeforeBreak = 150;
            dlg.TrialAbortAction = 'Repeat';
        end
        
        function trialTable = SetUpTrialTable(this)
            
            %-- condition variables ---------------------------------------
            i= 0;
             
            i = i+1;
            conditionVars(i).name   = 'Disparities';
            conditionVars(i).values = [0.1:0.2:0.9];
            %conditionVars(i).values = ones(size([0.1:0.2:0.9]))*30;
            
            i = i+1;
            conditionVars(i).name   = 'RotateDots';
            conditionVars(i).values = [0 5 10 30]; %[0 10 45];
            
            i = i+1;
            conditionVars(i).name   = 'SignDisparity';
            conditionVars(i).values = [-1 1]; 
            
            trialTableOptions = this.GetDefaultTrialTableOptions();
            trialTableOptions.trialSequence = 'Random';
            trialTableOptions.trialAbortAction = 'Repeat'; % Repeat, Delay, Drop
            trialTableOptions.trialsPerSession = 1000;
            trialTableOptions.numberOfTimesRepeatBlockSequence = this.ExperimentOptions.NumberOfRepetitions;
            trialTable = this.GetTrialTableFromConditions(conditionVars, trialTableOptions);
        end
        
        
        function [trialResult, thisTrialData] = runTrial( this, thisTrialData )
            
            try
                
                Enum = ArumeCore.ExperimentDesign.getEnum();
                graph = this.Graph;
                trialResult = Enum.trialResult.CORRECT;
                Screen('FillRect', graph.window, 0); % not sure if needed
                ShowCursor();
                
                % Screen and monitor settings
                screenWidth_pix = this.Graph.wRect(3); % screen width of one of the two eyes in pix
                moniterWidth_deg = (atan2d(this.ExperimentOptions.DisplayOptions.ScreenWidth/2, this.ExperimentOptions.DisplayOptions.ScreenDistance)) * 2;
                pixPerDeg = (screenWidth_pix*2) / moniterWidth_deg;
                                
                % Get fixation spot size in pix
                fixSizePix = pixPerDeg * this.ExperimentOptions.FixationSpotSize;
                
                % How big should the dots stimulus be in pix
                dots = zeros(3, this.ExperimentOptions.Number_of_Dots);
                stimWindow_pix = this.ExperimentOptions.MaxStimDeg * pixPerDeg;
                stimRingMin_pix = this.ExperimentOptions.MinStimDeg * pixPerDeg;
                
                % Initialize disparity settings:
                thisTrialData.DisparityArcMin = thisTrialData.Disparities * thisTrialData.SignDisparity;
                disparity_deg = thisTrialData.DisparityArcMin/60;
                disparityNeeded_pix = pixPerDeg*disparity_deg;

                % Initialize dots
                dots(1, :) = stimWindow_pix*rand(1, this.ExperimentOptions.Number_of_Dots) - (stimWindow_pix/2); % SR x coords
                dots(2, :) = stimWindow_pix*rand(1, this.ExperimentOptions.Number_of_Dots) - (stimWindow_pix/2); % SR y coords
                
                % Make the dot stimulus circular 
                distFromCenter = sqrt((dots(1,:)).^2 + (dots(2,:)).^2);
                while isempty(distFromCenter(distFromCenter>stimWindow_pix/2 | distFromCenter<stimRingMin_pix/2)) == 0 % while there are dots that are outside of the desired circle
                    idxs=find(distFromCenter>stimWindow_pix/2 | distFromCenter<stimRingMin_pix/2);
                    dots(1, idxs) = stimWindow_pix*rand(1, length(idxs)) - (stimWindow_pix/2); % resample those dots
                    dots(2, idxs) = stimWindow_pix*rand(1, length(idxs)) - (stimWindow_pix/2); 
                    distFromCenter = sqrt((dots(1,:)).^2 + (dots(2,:)).^2);
                end
                dots(3, :) = (ones(size(dots,2),1)')*disparityNeeded_pix; % how much the dots will shift by in pixels
                
                % Right and left shifted dots 
                leftStimDots = [dots(1,:)+(dots(3,:)/2); dots(2,:)]; %dots(1:2, :) + [dots(3, :)/2; zeros(1, numDots)]; 
                rightStimDots = [dots(1,:)-(dots(3,:)/2); dots(2,:)]; 
                
                % Rotating the dots 
                leftDistFromCenter = sqrt((leftStimDots(1,:)).^2 + (leftStimDots(2,:)).^2); %dist is measured from the ORIGIN (0,0) which is where the fixation dot is so we're rotating around fixation dot
                leftThetaDeg = atan2d(leftStimDots(2,:),leftStimDots(1,:));
                leftPolarPtX = cosd(leftThetaDeg + thisTrialData.RotateDots) .* leftDistFromCenter;
                leftPolarPtY = sind(leftThetaDeg + thisTrialData.RotateDots) .* leftDistFromCenter;
                rightDistFromCenter = sqrt((rightStimDots(1,:)).^2 + (rightStimDots(2,:)).^2); %where leftStimDots(1,:) is the x coord and leftStimDots(2,:) is the y coord
                rightThetaDeg = atan2d(rightStimDots(2,:),rightStimDots(1,:));
                rightPolarPtX = cosd(rightThetaDeg + thisTrialData.RotateDots) .* rightDistFromCenter;
                rightPolarPtY = sind(rightThetaDeg + thisTrialData.RotateDots) .* rightDistFromCenter;
                leftStimDots = [leftPolarPtX;leftPolarPtY];
                rightStimDots = [rightPolarPtX;rightPolarPtY];
                
                % What the response should be
                if thisTrialData.DisparityArcMin > 0
                    thisTrialData.CorrectResponse = 'F';
                elseif thisTrialData.DisparityArcMin < 0
                    thisTrialData.CorrectResponse = 'B';
                end
                
                % For the while loop trial start
                lastFlipTime                        = Screen('Flip', graph.window);
                secondsRemaining                    = this.ExperimentOptions.TrialDuration;
                thisTrialData.TimeStartLoop         = lastFlipTime;
                
                response = []; %initialize this
                
                while secondsRemaining > 0
                    
                    secondsElapsed      = GetSecs - thisTrialData.TimeStartLoop;
                    secondsRemaining    = this.ExperimentOptions.TrialDuration - secondsElapsed;
                    
                    
                    % -----------------------------------------------------------------
                    % --- Drawing of stimulus -----------------------------------------
                    % -----------------------------------------------------------------
                    
                    
                    % If it's during the time when the stimulus (dots) is on, then show the stimulus plus the fixation dot
                    if ( secondsElapsed > this.ExperimentOptions.InitFixDuration && secondsElapsed < this.ExperimentOptions.InitFixDuration + this.ExperimentOptions.TimeStimOn) % then show dots + fixation dot
                        % Select left-eye image buffer for drawing:
                        Screen('SelectStereoDrawBuffer', this.Graph.window, 0);
                        
                        % Draw left stim:
                        Screen('DrawDots', this.Graph.window, leftStimDots, this.ExperimentOptions.Size_of_Dots, [], this.Graph.wRect(3:4)/2, 1); % the 1 at the end means dot type where 1 2 or 3 is circular
                        Screen('DrawDots', this.Graph.window, [0;0], fixSizePix, this.targetColor, this.Graph.wRect(3:4)/2, 1); % fixation spot
                        Screen('FrameRect', this.Graph.window, [1 0 0], [], 5);
                        
                        % Select right-eye image buffer for drawing:
                        Screen('SelectStereoDrawBuffer', this.Graph.window, 1);
                        
                        % Draw right stim:
                        Screen('DrawDots', this.Graph.window, rightStimDots, this.ExperimentOptions.Size_of_Dots, [], this.Graph.wRect(3:4)/2, 1);
                        Screen('DrawDots', this.Graph.window, [0;0], fixSizePix, this.targetColor, this.Graph.wRect(3:4)/2, 1); % fixation spot
                        Screen('FrameRect', this.Graph.window, [0 1 0], [], 5);
                        
                    end
                    
                    % If it's the initial fixation time, or if it's after the stimulus time has expired, then just show the fixation dot
                    if ( secondsElapsed <= this.ExperimentOptions.InitFixDuration ||  secondsElapsed > this.ExperimentOptions.InitFixDuration + this.ExperimentOptions.TimeStimOn) 
                        % Select left-eye image buffer for drawing:
                        Screen('SelectStereoDrawBuffer', this.Graph.window, 0);
                        
                        % Draw left stim:
                        Screen('DrawDots', this.Graph.window, [0;0], fixSizePix, this.targetColor, this.Graph.wRect(3:4)/2, 1); % fixation spot
                        Screen('FrameRect', this.Graph.window, [1 0 0], [], 5);
                        
                        % Select right-eye image buffer for drawing:
                        Screen('SelectStereoDrawBuffer', this.Graph.window, 1);
                        
                        % Draw right stim:
                        Screen('DrawDots', this.Graph.window, [0;0], fixSizePix, this.targetColor, this.Graph.wRect(3:4)/2, 1); % fixation spot
                        Screen('FrameRect', this.Graph.window, [0 1 0], [], 5);
                        
                    end
                                        
                    
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
                    
                    % Only collect responses after the stim has been displayed
                    if secondsElapsed > this.ExperimentOptions.InitFixDuration
                        [keyIsDown, secs, keyCode, deltaSecs] = KbCheck();
                        if ( keyIsDown )
                            keys = find(keyCode);
                            for i=1:length(keys)
                                KbName(keys(i));
                                switch(KbName(keys(i)))
                                    case 'RightArrow'
                                        response = 'F';
                                    case 'LeftArrow'
                                        response = 'B';
                                    case 'DownArrow'
                                        response = 'F';
                                    case 'UpArrow'
                                        response = 'B';
                                end
                            end
                            if ( ~isempty(response) ) % if there is a response, break this trial and start the next
                                thisTrialData.Response = response;
                                thisTrialData.ResponseTime = GetSecs;
                                break;
                            end
                        end
                    end
                end
                
                if ( isempty(response) )
                    thisTrialData.Response = 'N';
                    thisTrialData.ResponseTime = GetSecs;
                    trialResult = Enum.trialResult.ABORT;
                else
                    trialResult = Enum.trialResult.CORRECT;
                end
                
                
            catch ex
                rethrow(ex)
            end
            
        end
        
       function [trialResult, thisTrialData] = runPostTrial( this, thisTrialData )
           % Record if the subject guessed correctly or not 
           if thisTrialData.Response == thisTrialData.CorrectResponse
               thisTrialData.GuessedCorrectly = 1;
           elseif thisTrialData.Response ~= thisTrialData.CorrectResponse
               thisTrialData.GuessedCorrectly = 0;
           end
           
           % Move this forward
           trialResult = thisTrialData.TrialResult;
       end
            
             
        
    end

    

    methods
        
        function [out] = Plot_ConstantStim(this)
            %%
            trialTable = this.Session.trialDataTable;
            
            
            trialTable.RespondedFront = ones(size(trialTable,1),1);
            trialTable.RespondedFront(trialTable.Response == 'B') = trialTable.RespondedFront(trialTable.Response == 'B') *0;
            RotateDotsCond = unique(trialTable.RotateDots);
            here=1;
            
                figure
                
                for arotation = 1:length(RotateDotsCond)
                    idxs = find(trialTable.RotateDots == RotateDotsCond(arotation));
                    temp=sortrows(array2table([trialTable.DisparityArcMin(idxs) trialTable.RespondedFront(idxs)],'VariableNames',{'DisparityArcMin','RespondedFront'}));
                    
                    % Fit model
                    modelspec = 'RespondedFront ~ DisparityArcMin';
                    mdl = fitglm(temp(:,{'RespondedFront', 'DisparityArcMin'}), modelspec, 'Distribution', 'binomial');
                    a=[-1:0.01:0, 0.01:0.01:1];
                    p = predict(mdl,a');
                    
                    % Get average response for that disparity
                    temp.meanedResp=zeros(height(temp),1);
                    [uniqueDisparities,~] = unique(temp.DisparityArcMin);
                    for i = 1:length(uniqueDisparities)
                        idxs=find(temp.DisparityArcMin == uniqueDisparities(i));
                        temp.meanedResp(idxs) = mean(temp.RespondedFront(idxs));
                    end
                    
                    % Get the numbers
                    alpha = mdl.Coefficients.Estimate(2);
                    beta=mdl.Coefficients.Estimate(1);
                    p1=0.5;
                    x1 = (log(p1/(1-p1))-beta)/alpha; %PSE % log(p/1-p) = ax+b where we know p and a and b and are trying to get x.
                    p2=.25;
                    x2 = (log(p2/(1-p2))-beta)/alpha;
                    p3=.75;
                    x3 = (log(p3/(1-p3))-beta)/alpha;
                    slope = (p3-p2) / (x3-x2);
                    threshold=(x3-x2)/2;
                  
                    % Plot
                    subplot(1,(length(RotateDotsCond)),arotation)
                    plot(a,p) % plot prediction
                    hold on;
                    plot(temp.DisparityArcMin,temp.meanedResp,'o'); hold on
                    ylim([0 1])
                    xlabel('Disparity (arcmin)')
                    ylabel('Proportion Front')
                    %text(min(xlim)+0.05, max(ylim)-0.1, sprintf('PSE: %.2f', x1), 'Horiz','left', 'Vert','bottom')
                    text(min(xlim)+0.05, max(ylim)-0.1, sprintf('Slope: %.2f', slope), 'Horiz','left', 'Vert','bottom')
                    text(min(xlim)+0.05, max(ylim)-0.15, sprintf('Threshold(arcmin): %.2f', threshold), 'Horiz','left', 'Vert','bottom')
                    title(sprintf('Rotation: %s',string(RotateDotsCond(arotation))))
                    
                    here=here+1;
                end
            



        end
    end
end