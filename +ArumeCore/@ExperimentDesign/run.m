function run(this)

% --------------------------------------------------------------------
% -- EXPERIMENT LOOP -------------------------------------------------
% --------------------------------------------------------------------
% This is the main method that controls the flow of the
% experiment. It functions as a finite state matchine

Enum = ArumeCore.ExperimentDesign.getEnum();


% --------------------------------------------------------------------
% possible states of the loop
% --------------------------------------------------------------------
INITIALIZNG_HARDWARE = 0;
INITIALIZNG_EXPERIMENT = 1;
IDLE = 2;
CALIBRATING = 3;
STARTING_RECORDING = 4;
RUNNING_TRIALS = 5;
FINILIZING_EXPERIMENT = 6;
SESSIONFINISHED = 7;
BREAK = 8;
DOWNLOADING_DATA = 9;
FINALIZING_HARDWARE = 10;
TRY_FINALIZING_AFTER_ERROR = 11;
% --------------------------------------------------------------------
% end possible states
% --------------------------------------------------------------------

% initialize variables
state = INITIALIZNG_HARDWARE;
trialsSinceBreak = 0;
trialsSinceCalibration = 0;
calibrationCounter = [];
lastError = [];
while(1)
    try
        switch( state )
            case INITIALIZNG_HARDWARE

                % initialize psychtoolbox
                this.Graph = ArumeCore.PTB(this.ExperimentOptions.Debug.DebugMode, this.ExperimentOptions.DisplayOptions);

                shouldContinue = 1;

                % initialize eye tracker
                if ( this.ExperimentOptions.UseEyeTracker )

                    switch(this.ExperimentOptions.EyeTracker)
                        case 'OpenIris'
                            this.eyeTracker = ArumeHardware.VOG();
                        case 'Eyelink'
                            this.eyeTracker = ArumeHardware.EyeTrackerEyelink();
                            this.eyeTracker.experimentOptions = this.ExperimentOptions;
                    end

                    shouldContinue = this.eyeTracker.Connect(this.Graph);

                    if ( shouldContinue )
                        this.eyeTracker.SetSessionName(this.Session.name);
                    else
                        this.eyeTracker = [];
                    end
                end

                if ( shouldContinue )
                    state = INITIALIZNG_EXPERIMENT;
                else
                    state = FINILIZING_EXPERIMENT;
                end

            case INITIALIZNG_EXPERIMENT
                shouldContinue = this.initBeforeRunning();

                if ( shouldContinue )
                    state = CALIBRATING;
                else
                    state = FINILIZING_EXPERIMENT;
                end

            case IDLE
                result = this.Graph.DlgSelect( ...
                    'What do you want to do next:', ...
                    { 'n' 'c' 'q' 'e'}, ...
                    { 'Continue with next trial' 'Calibrate' 'Quit' 'Show last error'} , [],[]);

                switch( result )
                    case 'n'
                        state = RUNNING_TRIALS;
                    case 'c'
                        state = CALIBRATING;
                    case {'q' 0}
                        dlgResult = this.Graph.DlgYesNo( 'Are you sure you want to exit?',[],[],20,20);
                        if( dlgResult )
                            state = FINILIZING_EXPERIMENT;
                        end
                    case 'e'
                        if ( ~isempty(lastError) )
                            disp(lastError.getReport)
                        else
                            disp('NO ERRROS');
                        end
                end

            case CALIBRATING

                calibrationSuccessful = 1;
                if ( ~isempty(this.eyeTracker))
                    this.eyeTracker.StopRecording();
                    calibrationSuccessful =  this.eyeTracker.Calibrate();
                end

                if ( calibrationSuccessful)
                    trialsSinceCalibration = 0;
                    if ( ~isempty(this.Session.currentRun.pastTrialTable))
                        calibrationCounter = max(this.Session.currentRun.pastTrialTable.CalibrationCounter) + 1;
                    else
                        calibrationCounter = 1;
                    end
                    state = STARTING_RECORDING;
                else
                    result = this.Graph.DlgSelect( ...
                        'Calibration was not successful what do you want to do?', ...
                        { 'n' 'c' 'q'}, ...
                        { 'Continue with next trial anyway' 'Try to calibrate again' 'Quit experiment'} , [],[]);

                    switch( result )
                        case 'n'
                            state = STARTING_RECORDING;
                        case 'c'
                            state = CALIBRATING;
                        case 'q'
                            dlgResult = this.Graph.DlgYesNo( 'Are you sure you want to exit?',[],[],20,20);
                            if( dlgResult )
                                state = FINILIZING_EXPERIMENT;
                            else
                                state = IDLE;
                            end
                    end
                end

            case STARTING_RECORDING

                if ( ~isempty(this.eyeTracker))
                    this.eyeTracker.StartRecording();
                end

                state = RUNNING_TRIALS;

            case BREAK
                this.Graph.DlgHitKey( 'Time for a break! hit a key to continue',[],[]);
                trialsSinceBreak = 0;
                state = IDLE;

            case RUNNING_TRIALS
                % force to hit a key to continue if the
                % previous trial was an abort or if the
                % experiment is set to ask for hit key before
                % every trial
                if ( (~isempty(this.Session.currentRun.pastTrialTable) && this.Session.currentRun.pastTrialTable.TrialResult(end) == Enum.trialResult.ABORT) ...
                        || this.ExperimentOptions.HitKeyBeforeTrial )
                    dlgResult = this.Graph.DlgHitKey( 'Hit a key to continue',[],[]);
                    if ( ~dlgResult )
                        state = IDLE;
                        continue;
                    end
                end

                try
                    commandwindow;

                    nCorrectTrials = 0;
                    if ( any(strcmp(this.Session.currentRun.pastTrialTable.Properties.VariableNames,'TrialResult')))
                        nCorrectTrials = sum(this.Session.currentRun.pastTrialTable.TrialResult == 'CORRECT');
                    end

                    %-- find which condition to run and the variable values for that condition
                    thisTrialData = table();
                    thisTrialData.TrialNumber  = nCorrectTrials+1;
                    thisTrialData.DateTimeTrialStart = string(datetime);
                    thisTrialData.CalibrationCounter = calibrationCounter;
                    thisTrialData = [thisTrialData this.Session.currentRun.futureTrialTable(1,:)];

                    fprintf('\nARUME :: TRIAL %d START (%d TOTAL) ...\n', nCorrectTrials+1, height(this.Session.currentRun.originalFutureTrialTable));

                    %------------------------------------------------------------
                    % -- PRE TRIAL ----------------------------------------------
                    %------------------------------------------------------------
                    thisTrialData.TimePreTrialStart = GetSecs;

                    [trialResult, thisTrialData] = this.runPreTrial( thisTrialData );
                    thisTrialData.TrialResult = trialResult;
                    thisTrialData.TimePreTrialStop = GetSecs;

                    %------------------------------------------------------------
                    % -- TRIAL --------------------------------------------------
                    %------------------------------------------------------------
                    if ( trialResult == Enum.trialResult.CORRECT )

                        thisTrialData.TimeTrialStart = GetSecs;

                        if ( ~isempty(this.eyeTracker))

                            if ( ~this.eyeTracker.IsRecording )
                                ME = MException('ArumeHardware.VOG:NotRecording', 'The eye tracker is not recording.');
                                throw(ME);
                            end

                            [framenumber, eyetrackertime] = this.eyeTracker.RecordEvent( ...
                                sprintf('TRIAL_START [trial=%d, condition=%d]', ...
                                thisTrialData.TrialNumber, thisTrialData.Condition) );

                            thisTrialData.EyeTrackerFrameNumberTrialStart = framenumber;
                            thisTrialData.EyeTrackerTimeTrialStart = eyetrackertime;

                            % Keep track of how many eye tracking files this session is
                            % split in and mark this trial with the correct file number
                            % from the linked files list
                            % TODO: a bit ugly. It would be
                            % good to clean this up.
                            if ( isempty( this.Session.currentRun.LinkedFiles) )
                                thisTrialData.FileNumber = 1;
                            else
                                if ( ischar(this.Session.currentRun.LinkedFiles.vogDataFile) )
                                    thisTrialData.FileNumber = 2;
                                else
                                    thisTrialData.FileNumber = length(this.Session.currentRun.LinkedFiles.vogDataFile)+1;
                                end
                            end
                        end

                        if ( ~isempty(this.Graph) )
                            this.Graph.ResetFlipTimes();
                        end
                        [trialResult, thisTrialData] = this.runTrial( thisTrialData );
                        thisTrialData.TrialResult = trialResult;
                        if ( ~isempty(this.Graph) )
                            thisTrialData.NumFlips = this.Graph.NumFlips;
                            thisTrialData.NumSlowFlips = this.Graph.NumSlowFlips;
                            thisTrialData.NumSuperSlowFlips = this.Graph.NumSuperSlowFlips;
                        end

                        thisTrialData.TimeTrialStop = GetSecs;

                        if (~isempty(this.eyeTracker) )
                            [framenumber, eyetrackertime] = this.eyeTracker.RecordEvent(....
                                sprintf('TRIAL_STOP [trial=%d, condition=%d]', ...
                                thisTrialData.TrialNumber, thisTrialData.Condition) );
                            thisTrialData.EyeTrackerFrameNumberTrialStop = framenumber;
                            thisTrialData.EyeTrackerTimeTrialStop = eyetrackertime;

                            if ( ~this.eyeTracker.IsRecording )
                                ME = MException('ArumeHardware.VOG:NotRecording', 'The eye tracker is not recording.');
                                throw(ME);
                            end
                        end
                    end

                    %------------------------------------------------------------
                    % -- POST TRIAL ---------------------------------------------
                    %------------------------------------------------------------
                    if ( trialResult == Enum.trialResult.CORRECT )

                        if ( this.ExperimentOptions.DisplayOptions.PlaySound)
                            this.PlaySound(thisTrialData.TrialResult);
                        end

                        thisTrialData.TimePostTrialStart = GetSecs;

                        [trialResult, thisTrialData] = this.runPostTrial( thisTrialData );
                        thisTrialData.TrialResult = trialResult;

                        thisTrialData.TimePostTrialStop = GetSecs;
                    end


                catch lastError
                    if ( streq(lastError.identifier, 'PSYCORTEX:USERQUIT' ) )
                        thisTrialData.TrialResult = Enum.trialResult.QUIT;
                        lastError = [];
                    else
                        thisTrialData.TrialResult = Enum.trialResult.ERROR;
                        thisTrialData.ErrorMessage = string(lastError.message);
                        % display error

                        beep
                        cprintf('red', '\n')
                        cprintf('red', '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n')
                        cprintf('red', '!!!!!!!!!!!!! ARUME ERROR: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n')
                        cprintf('red', '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n')
                        % disp(err.getReport);
                        disp(lastError.message);
                        cprintf('red', '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n')
                        cprintf('red', '!!!!!!!!!!!!! END ARUME ERROR: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n')
                        cprintf('red', '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n')
                        cprintf('red', '\n')
                    end
                end


                % -- Update past trial table
                this.Session.currentRun.AddPastTrialData(thisTrialData);

                if ( this.ExperimentOptions.DisplayOptions.ShowTrialTable)
                    % -- Display trial Table for last 20 trials
                    data = this.Session.currentRun.pastTrialTable;
                    varSelection = intersect(strsplit(this.ExperimentOptions.Debug.DisplayVariableSelection,' '),data.Properties.VariableNames,'stable');
                    if ( ~this.ExperimentOptions.Debug.DebugMode )
                        disp(data(max(1,end-20):end,varSelection));
                    else
                        disp(data);
                    end
                end


                if ( thisTrialData.TrialResult == Enum.trialResult.CORRECT )
                    %-- remove the condition that has just run from the future conditions list
                    this.Session.currentRun.futureTrialTable(1,:) = [];

                    if ( mod( nCorrectTrials, this.ExperimentOptions.TrialsBeforeSaving ) == 0 )
                        %-- save to disk temporary data
                        this.Session.save();
                    end

                    trialsSinceBreak = trialsSinceBreak + 1;
                else
                    %-- what to do in case of abort
                    switch(this.TrialTable.Properties.UserData.trialTableOptions.trialAbortAction) % TODO: save the trial abort action somehwere else
                        case 'Repeat'
                            % do nothing
                        case 'Delay'
                            % randomly get one of the future conditions in the current block
                            % and switch it with the next
                            currentblock = this.Session.currentRun.futureTrialTable.BlockNumber(1);
                            currentblockSeqNumber = this.Session.currentRun.futureTrialTable.BlockSequenceNumber(1);
                            futureConditionsInCurrentBlock = this.Session.currentRun.futureTrialTable(this.Session.currentRun.futureTrialTable.BlockNumber==currentblock & this.Session.currentRun.futureTrialTable.BlockSequenceNumber==currentblockSeqNumber,:);

                            newPosition = ceil(rand(1)*(height(futureConditionsInCurrentBlock)-1))+1;
                            c = futureConditionsInCurrentBlock(1,:);
                            futureConditionsInCurrentBlock(1,:) = futureConditionsInCurrentBlock(newPosition,:);
                            futureConditionsInCurrentBlock(newPosition,:) = c;
                            this.Session.currentRun.futureTrialTable(this.Session.currentRun.futureTrialTable.BlockNumber==currentblock & this.Session.currentRun.futureTrialTable.BlockSequenceNumber==currentblockSeqNumber,:) = futureConditionsInCurrentBlock;
                        case 'Drop'
                            %-- remove the condition that has just run from the future conditions list
                            this.Session.currentRun.futureTrialTable(1,:) = [];
                    end
                end

                %-- handle errors
                switch ( thisTrialData.TrialResult )
                    case Enum.trialResult.ERROR
                        state = IDLE;
                        continue;
                    case Enum.trialResult.QUIT
                        state = IDLE;
                        continue;
                end

                % -- Experiment or session finished ?
                if ( trialsSinceCalibration >= this.ExperimentOptions.TrialsBeforeCalibration )
                    state = CALIBRATING;
                end
                if ( trialsSinceBreak >= this.ExperimentOptions.TrialsBeforeBreak )
                    state = BREAK;
                end
                if ( ~isempty(this.Session.currentRun.futureTrialTable) && ~isempty(this.Session.currentRun.pastTrialTable) )
                    if ( this.Session.currentRun.pastTrialTable.Session(end) ~= this.Session.currentRun.futureTrialTable.Session(1) )
                        state = SESSIONFINISHED;
                    end
                end
                if ( isempty(this.Session.currentRun.futureTrialTable) )
                    state = FINILIZING_EXPERIMENT;
                end

            case SESSIONFINISHED
                cprintf('blue', '---------------------------------------------------------\n')
                cprintf('blue', '---------------------------------------------------------\n')
                cprintf('*blue', 'Session part finished! closing down and saving data ...\n');
                cprintf('blue', '---------------------------------------------------------\n')
                cprintf('blue', '---------------------------------------------------------\n')

                state = DOWNLOADING_DATA;

            case FINILIZING_EXPERIMENT
                cprintf('blue', '---------------------------------------------------------\n')
                cprintf('blue', '---------------------------------------------------------\n')
                cprintf('*blue', 'Experimental session finished! closing down and saving data ...\n');
                cprintf('blue', '---------------------------------------------------------\n')
                cprintf('blue', '---------------------------------------------------------\n')

                state = DOWNLOADING_DATA;

            case DOWNLOADING_DATA

                if (~isempty(this.eyeTracker))
                    this.eyeTracker.StopRecording();

                    disp('Downloading eye tracking files...');
                    now = datetime;
                    now.Format = 'uuuuMMdd_HHmmss';
                    files = this.eyeTracker.DownloadFile( this.Session.dataPath, strcat(this.Session.name, "_", string(now), ".edf"));

                    if (~isempty( files) )
                        switch(this.ExperimentOptions.EyeTracker)
                            case 'OpenIris'
                                disp(files{1});
                                disp(files{2});
                                if (length(files) > 2 )
                                    disp(files{3});
                                end
                                disp('Finished downloading');

                                this.Session.addFile('vogDataFile', files{1});
                                this.Session.addFile('vogCalibrationFile', files{2});
                                if (length(files) > 2 )
                                    this.Session.addFile('vogEventsFile', files{3});
                                end
                            case 'Eyelink'
                                this.Session.addExistingFile('vogDataFile', files{1});
                        end
                    else
                        disp('No eye tracking files downloaded!');
                    end
                end

                state = FINALIZING_HARDWARE;

            case FINALIZING_HARDWARE

                this.cleanAfterRunning();

                if ( ~isempty(this.eyeTracker))
                    this.eyeTracker.Disconnect();
                end

                this.Graph.Clear();

                this.Graph = [];
                disp('ARUME:: Done closing display and connections!');

                break; % finish loop

            case TRY_FINALIZING_AFTER_ERROR
                % This state is to try to finalize whatever we
                % can but to not make them depend on each other
                % so much. Specially we want to make sure
                % psytoolbox closes and frees the screen

                try
                    this.cleanAfterRunning();
                catch
                end

                try
                    if ( ~isempty(this.eyeTracker))
                        this.eyeTracker.Disconnect();
                    end
                catch
                end

                try
                    this.Graph.Clear();
                catch
                end


        end
    catch lastError
        beep
        cprintf('*red', '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n')
        cprintf('*red', '!!!!!!!!!!!!! ARUME ERROR: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n')
        cprintf('*red', '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n')
        disp(lastError.getReport);
        disp(lastError.message);
        cprintf('*red', '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n')
        cprintf('*red', '!!!!!!!!!!!!! END ARUME ERROR: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n')
        cprintf('*red', '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n')

        if ( state == FINALIZING_HARDWARE )
            state = TRY_FINALIZING_AFTER_ERROR;
        else
            state = FINALIZING_HARDWARE;
        end
    end
end
% --------------------------------------------------------------------
%% -------------------- END EXPERIMENT LOOP ---------------------------
% --------------------------------------------------------------------
end
