classdef OpticFlow_DualTask < ArumeExperimentDesigns.EyeTracking

    % OpticFlow_DualTaskexperiment for Arume 

    properties
        camera
        uicomponents
        audio
        shapes
    end
    
    % ---------------------------------------------------------------------
    % Experiment design methods
    % ---------------------------------------------------------------------
    methods ( Access = protected )
        function dlg = GetOptionsDialog( this, importing )
            dlg = GetOptionsDialog@ArumeExperimentDesigns.EyeTracking(this, importing);

            %% ADD new options
            dlg.TaskType = {{'Visual Search','Heading','Both','All interleaved'}};
            dlg.UniformityOptions = {{'Only Uniform' '{Only non-uniform}' 'both'}};
            dlg.NumShapes = {2000, 'Number of Shapes',[1000,100000]};
            dlg.ShapeSz = {12, 'Shape size (circle diam in px)',[1,100]}; % px - preferable odd
            dlg.HFOV = {66.35, 'horizontal FoV (degrees)',[1,10000]}; % based on view distance of 1.3000 meters
            dlg.FCP = {60, 'far clipping plane',[1,1000]}; % in meters
            dlg.ObserverHeight = {1.23, 'Observer Height',[0.1,10]}; % Center of the screen in meters
            dlg.ShapeSizeCue = {{'0','{1}'}, 'Size Cue?'};
            dlg.NumberTargets = {20, 'Number of Target Shapes'};

            % condition parameters
            % dlg.useeyelink = { {'{0}','1'}, 'Use Eyelink'};
            dlg.HeadingChanges = {[-15, -12.5, -10, -7.5, -5, -2.5, 2.5, 5, 7.5, 10, 12.5, 15], 'Heading Deltas (degrees)',[-100,100]}; % in meters per second walking = 1.31, jogging = 3.25, running = 5.76, driving = 13.41
            dlg.HeadingChangeDuration = {2, 'Heading change duration',[0,10]};
            dlg.Smoothing = {{'Gaussian','Linear','None'}};
            dlg.WalkSpeed = {3.25, 'Locomotion Speed (m/sec)',[0,100]}; % in meters per second walking = 1.31, jogging = 3.25, running = 5.76, driving = 13.41
            dlg.ShapeLifetime = {2, 'Shape Lifetime (secs)',[0,20]}; % in secs
            dlg.NumberTrials = {1, 'Number of Trials Per Condition',[1,10000]};
            dlg.AuditoryFeedback = {{'0','{1}'}, 'Auditory Feedback?'};

            %% CHANGE DEFAULTS values for existing options
            dlg.UseEyeTracker = { {'0','{1}' }};
            dlg.EyeTracker      = { {'OpenIris' 'Fove' '{Eyelink}' 'Mouse sim'} };
            dlg.Debug.DisplayVariableSelection = 'TrialNumber TrialResult Speed Stimulus'; % which variables to display every trial in the command line separated by spaces

            dlg.DisplayOptions.ScreenWidth = { 121 '* (cm)' [1 3000]};
            dlg.DisplayOptions.ScreenHeight = { 68 '* (cm)' [1 3000]};
            dlg.DisplayOptions.ScreenDistance = { 60 '* (cm)' [1 3000]};
            
            dlg.HitKeyBeforeTrial = { {'{0}','1'} };
            dlg.TrialDuration = {8, 'Trial duration',[1,100]}; % in secs
            dlg.TrialsBeforeBreak = 10; %             dlg.numberblocks = {5, 'Number of Blocks', [1,1000]};
            dlg.TrialsBeforeCalibration = 10; %             dlg.numberblocks = {5, 'Number of Blocks', [1,1000]};
            dlg.TrialAbortAction = 'Repeat';
        end
        
        function trialTable = SetUpTrialTable(this)
            
            %-- condition variables ---------------------------------------
            t = ArumeCore.TrialTableBuilder();

            switch(this.ExperimentOptions.UniformityOptions)
                case 'Only Uniform'
                    t.AddConditionVariable( 'Density', {'Uniform'})
                case 'Only non-uniform'
                    t.AddConditionVariable( 'Density', {'NonUniform'})
                case 'both'
                    t.AddConditionVariable( 'Density',{'Uniform' 'NonUniform'})
            end

            switch (this.ExperimentOptions.TaskType) 
                case 'Visual Search'
                    t.AddConditionVariable('Task',{'Visual Search'});
                case 'Heading'
                    t.AddConditionVariable('Task',{'Heading'});
                case 'Both'
                    t.AddConditionVariable('Task',{'Both'});
                case 'All interleaved'
                    t.AddConditionVariable('Task',{'Visual Search','Heading','Both'});
            end

            t.AddConditionVariable('WalkingSpeed',this.ExperimentOptions.WalkSpeed);
            t.AddConditionVariable('HeadingChange',this.ExperimentOptions.HeadingChanges);

            trialSequence = 'Random';
            blockSequence =  'Sequential';
            blockSequenceRepetitions = this.ExperimentOptions.NumberTrials; % basically how many repetitions of each unique condition. block wording confusing
            abortAction = 'Delay';
            trialsPerSession = this.ExperimentOptions.NumberTrials*height(t.ConditionTable); % full experiment
            trialTable = t.GenerateTrialTable(trialSequence, blockSequence, blockSequenceRepetitions, abortAction, trialsPerSession);

            this.ExperimentOptions.nTrialsTotal = height(trialTable);

            % set up delta-onset time
            mintonset = 1; % heading change cannot happen immediately
            mintoffset = 1; % and the trial cannot end immediately after heading change is completed
            trange = this.ExperimentOptions.TrialDuration-this.ExperimentOptions.HeadingChangeDuration-mintoffset-mintonset;
            trialTable.HeadingChangeOnsetTime = rand(height(trialTable),1).*trange + mintonset;

            % set up target type
            targetTypes = categorical({'red squares' 'red circles' 'green squares' 'green circles'});
            trialTable.SearchTarget = targetTypes(randi(length(targetTypes),this.ExperimentOptions.nTrialsTotal,1))';

            % set up target present/absent. Arume hates logicals so use doubles
            trialTable.TargetPresent = double(rand(this.ExperimentOptions.nTrialsTotal,1)>.5);

            % set up target present/absent
            responseOrder = categorical({'HeadingFirst' 'SearchFirst'});
            trialTable.ResponseOrder = responseOrder(randi(length(responseOrder),this.ExperimentOptions.nTrialsTotal,1))';

            % and trial number
            trialTable.Trial = (1:this.ExperimentOptions.nTrialsTotal)';

            % and block derivative (block change)
            trialTable.BlockChange = double([diff(trialTable.BlockNumber);0]);

        end
        
        % run initialization before the first trial is run
        % Use this function to initialize things that need to be
        % initialized before running but don't need to be initialized for
        % every single trial
        function shouldContinue = initBeforeRunning( this )

            % create camera object with rendering properties
            this = setUpShapeAppearanceAndCameraProjectionMatrix(this);

            % UI screen coordinates and shapes
            this = setUpUIVariables(this);

            % give feedback in the form of a sound effect?
            if this.ExperimentOptions.AuditoryFeedback
                this = getAudioFeedbackFiles(this);
            end
            shouldContinue = 1;
        end
        
        % runPreTrial
        % use this to prepare things before the trial starts
        function [trialResult, thisTrialData] = runPreTrial(this, thisTrialData )

            % initialize the camera to the starting position every trial
            this.camera.pos = [0, this.ExperimentOptions.ObserverHeight, 0, 0];

            % precompute the x-z position for every frame in the trial
            this = createLocomotionTrajectory(this,thisTrialData);

            % create array of randomly positioned distractor locations
            this = initializeShapePlacement(this,thisTrialData);

            % create a smaller array of targets
            this = initializeTargetShapes(this,thisTrialData);

            % and show pre-trial fixation
            showFixationCross(this, thisTrialData)

            Enum = ArumeCore.ExperimentDesign.getEnum();
            trialResult = Enum.trialResult.CORRECT;
        end
        
        

        function [trialResult, thisTrialData] = runTrial( this, thisTrialData )

            try

                Enum = ArumeCore.ExperimentDesign.getEnum();
                trialResult = Enum.trialResult.CORRECT;
                lastFlipTime        = GetSecs;
                thisTrialData.TimeStartLoop = lastFlipTime;

                % present optic flow stimulus. No response recorded during
                % movie presentation
                [this,thisTrialData] = presentStimulus(this,thisTrialData);

                % gdata = this.eyeTracker.GetCurrentData();
                
            catch ex
                rethrow(ex)
            end
            
        end        


        % runPostTrial
        function [trialResult, thisTrialData] = runPostTrial(this, thisTrialData)

            % request response from observer
            switch thisTrialData.Task
                case 'Visual Search'
                    [this, thisTrialData]  = getVisualSearchResponse(this, thisTrialData);
        
                case 'Heading'
                    [this, thisTrialData]  = getHeadingResponse(this, thisTrialData);
        
                case 'Both'
                    % randomly choose which one goes first though. 
                    switch thisTrialData.ResponseOrder
                        case 'SearchFirst'
                            [this, thisTrialData]  = getVisualSearchResponse(this, thisTrialData);
                            [this, thisTrialData]  = getHeadingResponse(this, thisTrialData);

                        case 'HeadingFirst'
                            [this, thisTrialData]  = getHeadingResponse(this, thisTrialData);
                            [this, thisTrialData]  = getVisualSearchResponse(this, thisTrialData); %#ok<*ASGLU>
                    end
            end

            Enum = ArumeCore.ExperimentDesign.getEnum();
            trialResult = Enum.trialResult.CORRECT;

            % % % % % % % % % % % % % % % % % % % % endOfTrialSequence(this,thisTrialData);
        end
        
        % run cleaning up after the session is completed or interrupted
        function cleanAfterRunning(this)

            % close audio buffers
            if this.ExperimentOptions.AuditoryFeedback
                if ( ~isempty(this.audio))
                    try
                        PsychPortAudio('Close', this.audio.pahandlecorrect);
                        PsychPortAudio('Close', this.audio.pahandleincorrect);
                    catch
                    end
                end
            end

        end
        
    end


    % methods
    %     function Plot_MattTesting(this)
    %         samples = this.Session.samplesDataTable;
    %         trials = this.Session.trialDataTable;
    % 
    %         figure
    %     end
    % end
    
end