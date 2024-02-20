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
            dlg.TaskType = {{'Visual Search','Heading','Both','{Single Tasks Interleaved}','All interleaved'}};
            % dlg.BlockTasks = {{'0','{1}'},'Block together tasks?'};
            dlg.NumberTaskBlocks = {3, 'How many task blocks?',[1,10000]}; % based on view distance of 1.3000 meters

            dlg.UniformityOptions = {{'Only Uniform' '{Only non-uniform}' 'both'}};
            dlg.NumShapes = {2000, 'Number of Shapes',[1000,100000]};
            dlg.ShapeSz = {12, 'Shape size (circle diam in px)',[1,100]}; % px - preferable odd
            dlg.HFOV = {66.35, 'horizontal FoV (degrees)',[1,10000]}; % based on view distance of 1.3000 meters
            dlg.FCP = {60, 'far clipping plane',[1,1000]}; % in meters
            dlg.ObserverHeight = {1.23, 'Observer Height',[0.1,10]}; % Center of the screen in meters
            dlg.ShapeSizeCue = {{'0','{1}'}, 'Size Cue?'};
            dlg.NumberTargets = {15, 'Number of Target Shapes'};

            % condition parameters
            dlg.HeadingChanges = {[-15, -12.5, -10, -7.5, -5, -2.5, 2.5, 5, 7.5, 10, 12.5, 15], 'Heading Deltas (degrees)',[-100,100]}; % in meters per second walking = 1.31, jogging = 3.25, running = 5.76, driving = 13.41
            dlg.HeadingChangeDuration = {2, 'Heading change duration',[0,10]};
            dlg.Smoothing = {{'Gaussian','Linear','None'}};
            dlg.WalkSpeed = {[3.25,5.76], 'Locomotion Speed (m/sec)',[0,100]}; % in meters per second walking = 1.31, jogging = 3.25, running = 5.76, driving = 13.41
            dlg.ShapeLifetime = {1, 'Shape Lifetime (secs)',[0,20]}; % in secs
            dlg.NumberTrials = {3, 'Number of Trials Per Condition',[1,10000]};
            dlg.AuditoryFeedback = {{'0','{1}'}, 'Auditory Feedback?'};

            %% CHANGE DEFAULTS values for existing options
            dlg.UseEyeTracker = { {'0','{1}' }};
            dlg.EyeTracker      = { {'OpenIris' 'Fove' '{Eyelink}' 'Mouse sim'} };
            dlg.Debug.DisplayVariableSelection = 'TrialNumber TrialResult Speed Stimulus'; % which variables to display every trial in the command line separated by spaces

            % overwrite default display parameters
            dlg.DisplayOptions.ScreenWidth = { 170 '* (cm)' [1 3000]};
            dlg.DisplayOptions.ScreenHeight = { 96 '* (cm)' [1 3000]};
            dlg.DisplayOptions.ScreenDistance = { 130 '* (cm)' [1 3000]};
            dlg.DisplayOptions.PlaySound = {{'{0}','1'},'Display Sound?'};
            
            dlg.HitKeyBeforeTrial = { {'{0}','1'} };
            dlg.TrialDuration = {8, 'Trial duration',[1,100]}; % in secs
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
                case 'Single Tasks Interleaved'
                    t.AddConditionVariable('Task',{'Visual Search','Heading'});
                case 'All interleaved'
                    t.AddConditionVariable('Task',{'Visual Search','Heading','Both'});
            end

            t.AddConditionVariable('WalkingSpeed',this.ExperimentOptions.WalkSpeed);
            t.AddConditionVariable('HeadingChange',this.ExperimentOptions.HeadingChanges);

            % Specify the blocking structure of the task
            nt =  this.ExperimentOptions.NumberTrials;
            nb = this.ExperimentOptions.NumberTaskBlocks;

            % we want to always have the specified number of trials, even
            % if they can't be evenly split into blocks. Using the below
            % code, we split the blocks up as evenly as possible (last
            % n blocks have 1 fewer trial
            ntpb = zeros(1,nb);
            ntpb = ntpb + fix(nt/nb);
            ntpb(1:mod(nt,nb)) = ntpb(1:mod(nt,nb)) + 1;
            
            switch (this.ExperimentOptions.TaskType) 
                case 'Visual Search'
                    for i = 1:this.ExperimentOptions.NumberTaskBlocks
                        t.AddBlock(find(t.ConditionTable.Task=='Visual Search'), ntpb(i));
                    end
                case 'Heading'
                    for i = 1:this.ExperimentOptions.NumberTaskBlocks
                        t.AddBlock(find(t.ConditionTable.Task=='Heading'), ntpb(i));
                    end
                case 'Both'
                    for i = 1:this.ExperimentOptions.NumberTaskBlocks
                        t.AddBlock(find(t.ConditionTable.Task=='Both'), ntpb(i));
                    end
                case 'Single Tasks Interleaved'
                    % we don't really want complete randomization, because
                    % we want the task to change regularly
                    tasks = categorical({'Visual Search','Heading'});
                    taskorder = tasks(randperm(length(tasks)));
                    for i = 1:this.ExperimentOptions.NumberTaskBlocks
                        for j = 1:length(tasks)
                            t.AddBlock(find(t.ConditionTable.Task==taskorder(j)), ntpb(i));
                        end
                    end
                case 'All interleaved'
                    tasks = categorical({'Visual Search','Heading','Both'});
                    taskorder = tasks(randperm(length(tasks)));
                    for i = 1:this.ExperimentOptions.NumberTaskBlocks
                        for j = 1:length(tasks)
                            t.AddBlock(find(t.ConditionTable.Task==taskorder(j)), ntpb(i));
                        end
                    end
            end

            trialSequence = 'Random';
            blockSequence =  'Sequential';
            blockSequenceRepetitions = 1; % we handle the blocking repetitions manually above
            abortAction = 'Delay';
            trialsPerSession = realmax; % full experiment
            trialTable = t.GenerateTrialTable(trialSequence, blockSequence, blockSequenceRepetitions, abortAction, trialsPerSession);

            this.ExperimentOptions.nTrialsTotal = height(trialTable);
            fprintf('\n\n THERE ARE %i TRIALS/ROWS IN THE CONSTRUCTED TABLE \n\n',this.ExperimentOptions.nTrialsTotal)

            % set up delta-onset time
            mintonset = 1; % heading change cannot happen immediately
            mintoffset = 1; % and the trial cannot end immediately after heading change is completed
            trange = this.ExperimentOptions.TrialDuration-this.ExperimentOptions.HeadingChangeDuration-mintoffset-mintonset;
            trialTable.HeadingChangeOnsetTime = rand(height(trialTable),1).*trange + mintonset;

            % set up target type. all equally probable
            targetTypes = categorical({'red squares' 'red circles' 'green squares' 'green circles'});
            idxs = repmat(1:length(targetTypes),1,ceil(this.ExperimentOptions.nTrialsTotal/length(targetTypes)));
            idxs = idxs(randperm(length(idxs)));
            idxs = idxs(1:this.ExperimentOptions.nTrialsTotal);
            trialTable.SearchTarget = targetTypes(idxs)';

            % set up target present/absent. Arume hates logicals so use doubles
            idxs = repmat([0,1],1,ceil(this.ExperimentOptions.nTrialsTotal/2));
            idxs = idxs(randperm(length(idxs)));
            idxs = idxs(1:this.ExperimentOptions.nTrialsTotal);
            trialTable.TargetPresent = idxs';

            % set up target present/absent
            responseOrder = categorical({'HeadingFirst' 'SearchFirst'});
            idxs = repmat(1:length(responseOrder),1,ceil(this.ExperimentOptions.nTrialsTotal/length(responseOrder)));
            idxs = idxs(randperm(length(idxs)));
            idxs = idxs(1:this.ExperimentOptions.nTrialsTotal);
            trialTable.ResponseOrder = responseOrder(idxs)';

            % and trial number
            trialTable.Trial = (1:this.ExperimentOptions.nTrialsTotal)';

            % and block derivative (block change)
            trialTable.BlockChange = double([diff(trialTable.BlockNumber);0]);

            % make sure that the break is always half way through a
            % task-block
            midthroughblock = length(t.Blocks.ConditionsIncluded{1})*t.Blocks.TrialsPerCondition(1)/2;

            if floor(midthroughblock) ~= midthroughblock
                warning('Breaks are not perfectly aligned with blocks')
                midthroughblock = ceil(midthroughblock);
            end

            this.ExperimentOptions.TrialsBeforeBreak = midthroughblock; %             dlg.numberblocks = {5, 'Number of Blocks', [1,1000]};
            this.ExperimentOptions.TrialsBeforeCalibration = midthroughblock; %             dlg.numberblocks = {5, 'Number of Blocks', [1,1000]};
            
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

            % show block number (Arume doesn't do this by default)
            endOfTrialSequence(this,thisTrialData);
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