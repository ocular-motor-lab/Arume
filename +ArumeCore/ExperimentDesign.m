classdef ExperimentDesign < handle
    %EXPERIMENTDESIGN Base class for all experiment designs (paradigms).
    % All experiment designs must inherit from this class and must override
    % some of the methods.
    %
    % A experiment design contains the main trail flow but also a lot of
    % options regarding configuration of the experiment, randomization,
    % etc.
    
    properties( SetAccess = protected)
        Name                % Name of the experiment design

        Session = [];       % The session that is currently running this experiment design
        ExperimentOptions   = [];  % Options of this specific experiment design
        TrialTable          = [];

        Graph               = [];   % Display handle (usually psychtoolbox).
        eyeTracker          = [];   % Eye tracker handle
    end
        
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % THESE ARE THE METHODS THAT SHOULD BE IMPLEMENTED BY NEW EXPERIMENT
    % DESIGNS
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access=protected) % Methods to be overriden in experiments
        
        % Gets the options that be set in the UI when creating a new
        % session of this experiment (in structdlg format)
        % Some common options will be added
        function dlg = GetOptionsDialog( this, importing )
            dlg = [];
            
            if ( ~importing)
                dlg.Debug.DebugMode = { {'{0}','1'} };
                dlg.Debug.DisplayVariableSelection = 'TrialNumber TrialResult'; % which variables to display every trial in the command line separated by spaces
            end

            dlg.DisplayOptions = ArumeCore.PTB.GetDisplayOptions();

            dlg.HitKeyBeforeTrial = { {'{0}','1'} };
            dlg.TrialDuration = 10;
            dlg.TrialsBeforeBreak = 1000;
            dlg.TrialsBeforeCalibration = 1000;
            dlg.TrialsBeforeSaving = 1;

            dlg.UseEyeTracker   = { {'0' '{1}'} };
            dlg.EyeTracker      = { {'{OpenIris}' 'Fove' 'Eyelink' 'Mouse sim'} };

            if ( exist('importing','var') && importing )
                dlg.DataFiles = { {['uigetfile(''' fullfile(pwd,'*.txt') ''',''MultiSelect'', ''on'')']} };
                dlg.EventFiles = { {['uigetfile(''' fullfile(pwd,'*.txt') ''',''MultiSelect'', ''on'')']} };
                dlg.CalibrationFiles = { {['uigetfile(''' fullfile(pwd,'*.cal') ''',''MultiSelect'', ''on'')']} };
            end
        end
        
        % Set up the trial table when a new session is created
        function trialTable = SetUpTrialTable( this )
            t = ArumeCore.TrialTableBuilder();
            trialTable = t.GenerateTrialTable();
        end
        
        % run initialization before the first trial is run
        % Use this function to initialize things that need to be
        % initialized before running but don't need to be initialized for
        % every single trial
        function shouldContinue = initBeforeRunning( this )
            shouldContinue = 1;
        end
        
        % runPreTrial
        % use this to prepare things before the trial starts
        function [trialResult, thisTrialData] = runPreTrial(this, thisTrialData )
            Enum = ArumeCore.ExperimentDesign.getEnum();
            trialResult = Enum.trialResult.CORRECT;
        end
        
        % runTrial
        function [trialResult, thisTrialData] = runTrial( this, thisTrialData)
            Enum = ArumeCore.ExperimentDesign.getEnum();
            trialResult = Enum.trialResult.CORRECT;
        end
        
        % runPostTrial
        function [trialResult, thisTrialData] = runPostTrial(this, thisTrialData)
            Enum = ArumeCore.ExperimentDesign.getEnum();
            trialResult = Enum.trialResult.CORRECT;
        end
        
        % run cleaning up after the session is completed or interrupted
        function cleanAfterRunning(this)
        end
        
    end
        
    % --------------------------------------------------------------------
    % Analysis methods --------------------------------------------------
    % --------------------------------------------------------------------
    
    methods ( Access = public )
        
        function optionsDlg = GetAnalysisOptionsDialog(this)
            optionsDlg = VOGAnalysis.GetParameterOptions();
        end

        function [analysisResults, samplesDataTable, trialDataTable, sessionTable]  = RunDataAnalyses(this, analysisResults, samplesDataTable, trialDataTable, sessionTable, options)
        end
        
        function ImportSession( this )
            newRun = ArumeCore.ExperimentRun.SetUpNewRun( this );
            vars = newRun.futureTrialTable;
            vars.TrialResult = categorical(cellstr('CORRECT'));
            vars.TrialNumber = 1;
            vars.FileNumber = 1;
            newRun.AddPastTrialData(vars);
            newRun.futureTrialTable(:,:) = [];
            this.Session.importCurrentRun(newRun);
            
            
            dataFiles = this.ExperimentOptions.DataFiles;
            eventFiles = this.ExperimentOptions.EventFiles;
            calibrationFiles = this.ExperimentOptions.CalibrationFiles;
            if ( ~iscell(dataFiles) )
                dataFiles = {dataFiles};
            end
            if ( ~iscell(eventFiles) )
                eventFiles = {eventFiles};
            end
            if ( ~iscell(calibrationFiles) )
                calibrationFiles = {calibrationFiles};
            end
            
            for i=1:length(dataFiles)
                if (exist(dataFiles{i},'file') )
                    this.Session.addFile('vogDataFile', dataFiles{i});
                end
            end
            for i=1:length(eventFiles)
                if (exist(eventFiles{i},'file') )
                    this.Session.addFile('vogEventsFile', eventFiles{i});
                end
            end
            for i=1:length(calibrationFiles)
                if (exist(calibrationFiles{i},'file') )
                    this.Session.addFile('vogCalibrationFile', calibrationFiles{i});
                end
            end
        end
    end

    methods (Access = private)

        function  [samplesDataTable, trialDataTable, sessionTable] = prepareTablesForAnalysis( this, options)
            Enum = ArumeCore.ExperimentDesign.getEnum();
            
            if ( isempty(  this.Session.currentRun ) )
                return;
            end
            
            if ( ~options.Prepare_For_Analysis_And_Plots )
                return;
            end

            %% 0) Create the basic trial data table (without custom experiment stuff)
            trialDataTable = this.Session.currentRun.pastTrialTable;

            % remove errors and aborts for analysis
            if (~isempty(trialDataTable))
                % Trial attempt is just a continuos unique number for
                % each past trial.
                trialDataTable.TrialAttempt = (1:height(trialDataTable))';

                % just in case for old data. TrialResult used to be
                % numeric. Now it is categorical but the categories
                % match the old numbers+1;
                if ( ~iscategorical(trialDataTable.TrialResult) )
                    trialDataTable.TrialResult = Enum.trialResult.PossibleResults(trialDataTable.TrialResult+1);
                end
                
                KEEP_ONLY_CORRECT_TRIALS = 1;
                if ( KEEP_ONLY_CORRECT_TRIALS )

                    % in old files TrialNumber counted all trials not just
                    % correct trials. So we fix it for code down the line
                    % it could also be missing
                    if ( ~any(strcmp(trialDataTable.Properties.VariableNames,'TrialNumber')) || ...
                            sum(trialDataTable.TrialResult == Enum.trialResult.CORRECT) < max(trialDataTable.TrialNumber) )
                        % rebuild trial number as a counter of past correct
                        % trials plus one
                        trialDataTable.TrialNumber = cumsum([1;trialDataTable.TrialResult(1:end-1) == Enum.trialResult.CORRECT]);
                    end

                    % keep only correct trials from now on
                    % TODO: rethink this. Depending on how the experiment
                    % is programmed it may be interesting to look at the
                    % aborts.
                    trialDataTable(trialDataTable.TrialResult ~= Enum.trialResult.CORRECT ,:) = [];
                end

                % merge the columns in trials with the ones already
                % present in the trialDataTable.
                % It is only necessary to rerun this stage zero if
                % this.trialDataTable is not empty because there may be
                % changes on the code. Otherwise we could change it to
                % get here only if trialDataTable is empty.
                if ( ~isempty(this.Session.trialDataTable) )
                    rightVariables = setdiff(this.Session.trialDataTable.Properties.VariableNames, trialDataTable.Properties.VariableNames);
                    trialDataTable =  outerjoin(trialDataTable, this.Session.trialDataTable, 'Keys', 'TrialNumber', 'MergeKeys',true, 'RightVariables', rightVariables );
                end
            end

            %% 1) Prepare the sample data table
            samplesDataTable = this.Session.samplesDataTable;

            if ( isempty(samplesDataTable) )
                % In most cases this will just be from EyeTracking
                % experiment but there could be others that have a
                % different way to load sample data.
                [samplesDataTable, cleanedData, calibratedData, rawData] = this.PrepareSamplesDataTableEyeTracking(options);
                % TODO: I don't like this here. It should be moved
                % to session. But may have memory problems at some
                % point
                this.Session.WriteVariableIfNotEmpty(rawData,'rawDataTable');
                this.Session.WriteVariableIfNotEmpty(cleanedData,'cleanedData');
                this.Session.WriteVariableIfNotEmpty(calibratedData,'calibratedData');
            end
            cprintf('blue', '++ ARUME::Done with samplesDataTable.\n');

            %% 2) Prepare the trial data table
            [trialDataTable,samplesDataTable] = this.PrepareTrialDataTableEyeTracking(trialDataTable, samplesDataTable,  options);
            cprintf('blue', '++ ARUME::Done with trialDataTable.\n');

            %% 3) Prepare session data table
            newSessionDataTable = this.Session.GetBasicSessionDataTable();
            newSessionDataTable = this.PrepareSessionDataTableEyeTracking(newSessionDataTable, options);
            newSessionDataTable.LastAnalysisDateTime = char(string(datetime('now')));

            options = FlattenStructure(options); % eliminate strcuts with the struct so it can be made into a row of a table
            opts = fieldnames(options);
            s = this.GetExperimentOptionsDialog(1);
            for i=1:length(opts)
                if ( isempty(options.(opts{i})))
                    newSessionDataTable.(['AnalysisOption_' opts{i}]) = {''};
                elseif ( ~ischar( options.(opts{i})) && numel(options.(opts{i})) <= 1)
                    newSessionDataTable.(['AnalysisOption_' opts{i}]) = options.(opts{i});
                elseif (isfield( s, opts{i}) && iscell(s.(opts{i})) && iscell(s.(opts{i}){1}) && length(s.(opts{i}){1}) >1)
                    newSessionDataTable.(['AnalysisOption_' opts{i}]) = categorical(cellstr(options.(opts{i})));
                elseif (~ischar(options.(opts{i})) && numel(options.(opts{i})) > 1 )
                    newSessionDataTable.(['AnalysisOption_' opts{i}]) = {options.(opts{i})};
                else
                    newSessionDataTable.(['AnalysisOption_' opts{i}]) = string(options.(opts{i}));
                end
            end

            sessionTable = newSessionDataTable;
        end
        
        function [samplesDataTable, cleanedData, calibratedData, rawData] = PrepareSamplesDataTableEyeTracking(this, options)
            
            if ( isfield(this.ExperimentOptions,'EyeTracker') )
                eyeTrackerType = this.ExperimentOptions.EyeTracker;
            else
                eyeTrackerType = 'OpenIris';
            end
            
            switch(eyeTrackerType)
                case 'OpenIris'
                    calibrationsForEachTrial = [];
                    
                    % if this session is not a calibration
                    if ( ~strcmp(this.Session.experimentDesign.Name, 'Calibration') )

                        % TODO: FIND A BETTER WAY TO GET ALL THE RELATED
                        % SESSIONS
                        arume = Arume('nogui');
                        calibrationSessions = arume.currentProject.findSessionBySubjectAndExperiment(this.Session.subjectCode, 'Calibration');
                        calibrationTables = {};
                        calibrationCRTables = {};
                        calibrationTimes = NaT(0);
                        calibrationNames = {};
                        for i=1:length(calibrationSessions)
               
                            if ( ~isfield( calibrationSessions(i).analysisResults, 'calibrationTable') || isempty(calibrationSessions(i).analysisResults.calibrationTable))
                                % analyze the calibration sessions just in case
                                % they have not before
                                analysisOptions = arume.getAnalysisOptionsDefault( calibrationSessions(i) );
                                calibrationSessions(i).runAnalysis(analysisOptions);
                            end

                            if ( isfield( calibrationSessions(i).analysisResults, 'calibrationTable') )
                                calibrationTables{i} = calibrationSessions(i).analysisResults.calibrationTable;
                                calibrationCRTables{i} = calibrationSessions(i).analysisResults.calibrationTableCR;
                            else
                                calibrationTables{i} = table();
                                calibrationCRTables{i} = table();
                            end
                            calibrationTimes(i) = datetime(calibrationSessions(i).currentRun.pastTrialTable.DateTimeTrialStart{end});
                            calibrationNames{i} =  calibrationSessions(i).name;
                        end

                        calibrations = table(string(calibrationNames'), calibrationTables', calibrationCRTables', calibrationTimes','VariableNames',{'SessionName','CalibrationTable','CalibrationCRTable','DateTime'});
                        calibrations = sortrows(calibrations,'DateTime');

                        % loop through trials to find the relavant calibration
                        calibrationsForEachTrial = nan(height(this.Session.currentRun.pastTrialTable),1);
                        for i=1:height(this.Session.currentRun.pastTrialTable)
                            if ( isempty( calibrations))
                                continue;
                            end
                            trialStartTime = datetime(this.Session.currentRun.pastTrialTable.DateTimeTrialStart(i));

                            pastClosestCalibration = find((trialStartTime - calibrations.DateTime)>0,1, 'last');
                            if ( isempty( pastClosestCalibration))
                                continue;
                            end

                            if (i==1)
                                if ( (trialStartTime-calibrations.DateTime(pastClosestCalibration)) < minutes(5) )
                                    calibrationsForEachTrial(i) = pastClosestCalibration;
                                end
                            else
                                previousTrialCalibration = calibrationsForEachTrial(i-1);
                                % TODO consider the case when you take a
                                % break and forget to do a calibration
                                % before restarting. Right now we will keep
                                % the calibration from the previous trial
                                
                                if( previousTrialCalibration == pastClosestCalibration)
                                    % this is the case for a following trial
                                    % after a calibration
                                    calibrationsForEachTrial(i) = previousTrialCalibration;
                                else
                                    % this is the case for trial following a
                                    % break when one or more calibrations where
                                    % performed
                                    calibrationsForEachTrial(i) = pastClosestCalibration;
                                end
                            end
                        end

                        % if we did not find any calibration for the trials
                        % we behave as if there were no calibrations
                        if ( all(isnan(calibrationsForEachTrial)))
                            calibrationsForEachTrial = [];
                        end
                    end


                    samplesDataTable = table();
                    cleanedData = table();
                    calibratedData = table();
                    rawData = table();
                    
                    if ( ~isprop(this.Session.currentRun, 'LinkedFiles' ) || ~isfield(this.Session.currentRun.LinkedFiles,  'vogDataFile') )
                        return;
                    end
                    
                    dataFiles = this.Session.currentRun.LinkedFiles.vogDataFile;
                    calibrationFiles = this.Session.currentRun.LinkedFiles.vogCalibrationFile;
                    
                    if (~iscell(dataFiles) )
                        dataFiles = {dataFiles};
                    end
                    
                    if (~iscell(calibrationFiles) )
                        calibrationFiles = {calibrationFiles};
                    end
                    
                    if (length(calibrationFiles) == 1)
                        calibrationFiles = repmat(calibrationFiles,size(dataFiles));
                    elseif length(calibrationFiles) ~= length(dataFiles)
                        error('ERROR preparing sample data set: The session should have the same number of calibration files as data files or 1 calibration file');
                    end
                    




                    if ( isempty(calibrationsForEachTrial) )
                        [samplesDataTable, cleanedData, calibratedData, rawData] = VOGAnalysis.LoadCleanAndResampleData(this.Session.dataPath, dataFiles, calibrationFiles, options);
                    else
                        % Dealing with trials without a calibration. Add an
                        % empty calibration to the table and change the NaN
                        % indices for the index of the empty table
                        calibrations.CalibrationTable{end+1} = table();
                        calibrationsForEachTrial(isnan(calibrationsForEachTrial)) = height(calibrations);
                        
                        calibrationTables = [calibrations(calibrationsForEachTrial,:) this.Session.currentRun.pastTrialTable(:,'FileNumber')];
                        [~,idx] = unique(calibrationTables.SessionName);
                        calibrfationTablesPerFile = calibrationTables(idx,:);
                        [samplesDataTable, cleanedData, calibratedData, rawData] = VOGAnalysis.LoadCleanAndResampleDataArumeMultiCalibration(this.Session.dataPath, dataFiles, calibrationFiles, calibrfationTablesPerFile , options);
                    end

                case 'Fove'
                    
                    samplesDataTable = table();
                    cleanedData = table();
                    calibratedData = table();
                    rawData = table();
                    
                    if ( ~isprop(this.Session.currentRun, 'LinkedFiles' ) || ~isfield(this.Session.currentRun.LinkedFiles,  'foveDataFile') )
                        return;
                    end
                    
                    dataFiles = this.Session.currentRun.LinkedFiles.foveDataFile;
                    if (~iscell(dataFiles) )
                        dataFiles = {dataFiles};
                    end
                   
                    [samplesDataTable, cleanedData, rawData] = VOGAnalysis.LoadCleanAndResampleDataFOVE(this.Session.dataPath, dataFiles, options);
                case 'Eyelink'
                    
                    samplesDataTable = table();
                    cleanedData = table();
                    calibratedData = table();
                    rawData = table();
                    
                    if ( ~isprop(this.Session.currentRun, 'LinkedFiles' ) || ~isfield(this.Session.currentRun.LinkedFiles,  'vogDataFile') )
                        return;
                    end
                    
                    dataFiles = this.Session.currentRun.LinkedFiles.vogDataFile;
                    if (~iscell(dataFiles) )
                        dataFiles = {dataFiles};
                    end
                   
                    [samplesDataTable, cleanedData, rawData] = VOGAnalysis.LoadCleanAndResampleDataEyelink(this.Session.dataPath, dataFiles, options);

            end
        end
        
        function [trialDataTable samplesDataTable] = PrepareTrialDataTableEyeTracking( this, trialDataTable, samplesDataTable, options)
            
            if ( isfield(this.ExperimentOptions,'EyeTracker') )
                eyeTrackerType = this.ExperimentOptions.EyeTracker;
            else
                eyeTrackerType = 'OpenIris';
            end
            
            if ( ~isempty( samplesDataTable ) )
                
                
                switch(eyeTrackerType)
                    case 'OpenIris'
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % DEAL WITH OLD FILES THAT HAVE SOME MISSING FIELDS
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                        % Old versions did not have a file number column. This code
                        % recovers one from the events file.
                        if ( ~any(strcmp(trialDataTable.Properties.VariableNames,'FileNumber')) )
                            % Find the file number that corresponds with each trial.
                            if ( isfield(this.Session.currentRun.LinkedFiles, 'vogEventsFile'))
                                ev = VOGAnalysis.ReadEventFiles(this.Session.dataPath, this.Session.currentRun.LinkedFiles.vogEventsFile);
                                trialDataTable.FileNumber = ev.FileNumber(this.Session.currentRun.pastTrialTable.TrialResult=='CORRECT');
                            else
                                % if there are no event files assume all the data
                                % comes from one single file
                                trialDataTable.FileNumber = ones(size(trialDataTable.TrialNumber));
                            end
                        end
                        
                        % Old versions did not have a timestamp from the eye
                        % tracker to line up trials precisely. This reads the old
                        % event file to try to line them up.
                        if ( ~any(strcmp(trialDataTable.Properties.VariableNames,'EyeTrackerFrameNumberTrialStart')) )
                            if ( height(trialDataTable) == 1 )
                                trialDataTable.EyeTrackerFrameNumberTrialStart = samplesDataTable.RawFrameNumber(1);
                                trialDataTable.EyeTrackerFrameNumberTrialStop = samplesDataTable.RawFrameNumber(end);
                            else
                                events = readtable(fullfile(this.Session.folder, this.Session.currentRun.LinkedFiles.vogEventsFile),'Delimiter',' ');
                                % get frame number from event table
                                % get trial duration from trialDataTable
                                % calculate frame number off of those two things
                                events = events(this.Session.currentRun.pastTrialTable.TrialResult=='CORRECT',:);
                                trialDataTable.EyeTrackerFrameNumberTrialStart = events.Var2 - samplesDataTable.LeftCameraRawFrameNumber(1)+1;
                                if ( min(trialDataTable.EyeTrackerFrameNumberTrialStart) < 0 )
                                    % crappy fix for files recorded around july-aug
                                    % 2018. The data files and the event files have
                                    % different frame numbers so they cannot be
                                    % lined up exactly.
                                    daysForFirstTrialStart = datenum(events.Var1{1},'yyyy-mm-dd-HH:MM:SS');
                                    a = regexp(this.Session.currentRun.LinkedFiles.vogEventsFile,'.+PostProc\-(?<date>.+)\-events\.txt', 'names');
                                    daysForFileOpening = datenum(a.date,'yyyymmmdd-HHMMSS');
                                    secondsFromFileOpeningToFirstTrial = (daysForFirstTrialStart-daysForFileOpening)*24*60*60;
                                    frameNumberFirstTrialStart = min( secondsFromFileOpeningToFirstTrial*100,  max(samplesDataTable.FrameNumber) - (trialDataTable.TimeTrialStop(end)-trialDataTable.TimeTrialStart(1))*100);
                                    trialDataTable.EyeTrackerFrameNumberTrialStart = (trialDataTable.TimeTrialStart-trialDataTable.TimeTrialStart(1))*100 + frameNumberFirstTrialStart;
                                end
                                trialDataTable.EyeTrackerFrameNumberTrialStop = (trialDataTable.TimeTrialStop - trialDataTable.TimeTrialStart)*100 + trialDataTable.EyeTrackerFrameNumberTrialStart;
                            end
                        end
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % END OF DEALING WITH OLD FILES
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                        
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % LINE UP TRIALS AND SAMPLES DATA
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                        % Find the samples that mark the begining and ends of trials
                        trialDataTable.SampleStartTrial = nan(size(trialDataTable.TrialNumber));
                        trialDataTable.SampleStopTrial = nan(size(trialDataTable.TrialNumber));
                        if ( ~any(strcmp(samplesDataTable.Properties.VariableNames, 'FileNumber')))
                            samplesDataTable.FileNumber = ones(size(samplesDataTable.Time));
                        end
                        for i=1:height(trialDataTable)
                            trialDataTable.SampleStartTrial(i) = find(samplesDataTable.FileNumber' == trialDataTable.FileNumber(i) & samplesDataTable.RawFrameNumber'>=trialDataTable.EyeTrackerFrameNumberTrialStart(i),1,'first');
                            trialDataTable.SampleStopTrial(i) = find(samplesDataTable.FileNumber' == trialDataTable.FileNumber(i) & samplesDataTable.RawFrameNumber'<=trialDataTable.EyeTrackerFrameNumberTrialStop(i),1,'last');
                        end
                    case 'Fove'
                        % Find the samples that mark the begining and ends of trials
                        trialDataTable.SampleStartTrial = nan(size(trialDataTable.TrialNumber));
                        trialDataTable.SampleStopTrial = nan(size(trialDataTable.TrialNumber));
                        if ( ~any(strcmp(samplesDataTable.Properties.VariableNames, 'FileNumber')))
                            samplesDataTable.FileNumber = ones(size(samplesDataTable.Time));
                        end
                        for i=1:height(trialDataTable)
                            trialDataTable.SampleStartTrial(i) = find(samplesDataTable.FileNumber' == trialDataTable.FileNumber(i) & samplesDataTable.RawTime'>=trialDataTable.TrialStartTime(i),1,'first');
                            trialDataTable.SampleStopTrial(i) = find(samplesDataTable.FileNumber' == trialDataTable.FileNumber(i) & samplesDataTable.RawTime'<=trialDataTable.TrialEndTime(i),1,'last');
                        end

                        % Add average head position data to the trialtable
                        trialDataTable.HeadYaw = nan(height(trialDataTable),1);
                        trialDataTable.HeadPitch = nan(height(trialDataTable),1);
                        trialDataTable.HeadRoll = nan(height(trialDataTable),1);
                        for i=1:height(trialDataTable)
                            idx = trialDataTable.SampleStartTrial(i):trialDataTable.SampleStopTrial(i);

                            trialDataTable.HeadYaw(i) = mean(samplesDataTable.HeadYaw(idx),1,"omitnan");
                            trialDataTable.HeadPitch(i) = mean(samplesDataTable.HeadPitch(idx),1,"omitnan");
                            trialDataTable.HeadRoll(i) = mean(samplesDataTable.HeadRoll(idx),1,"omitnan");
                        end
                    case 'Eyelink'

                        % Eyelink does not provide frame numbers in real
                        % time, only timestamps. The frame numbers are
                        % created artificially after reading the edf but
                        % are still useful so we will regenate them here
                        % in the trial table
                        for i=1:height(trialDataTable)
                            i1 = find(samplesDataTable.RawTime>=trialDataTable.EyeTrackerTimeTrialStart(i),1,'first');
                            i2 = find(samplesDataTable.RawTime<=trialDataTable.EyeTrackerTimeTrialStop(i),1,'last');
                            trialDataTable.EyeTrackerFrameNumberTrialStart(i) = samplesDataTable.RawFrameNumber(i1);
                            trialDataTable.EyeTrackerFrameNumberTrialStop(i) = samplesDataTable.RawFrameNumber(i2); 
                        end

                        % Find the samples that mark the begining and ends of trials
                        trialDataTable.SampleStartTrial = nan(size(trialDataTable.TrialNumber));
                        trialDataTable.SampleStopTrial = nan(size(trialDataTable.TrialNumber));
                        if ( ~any(strcmp(samplesDataTable.Properties.VariableNames, 'FileNumber')))
                            samplesDataTable.FileNumber = ones(size(samplesDataTable.Time));
                        end
                        for i=1:height(trialDataTable)
                            trialDataTable.SampleStartTrial(i) = find(samplesDataTable.FileNumber' == trialDataTable.FileNumber(i) & samplesDataTable.RawFrameNumber'>=trialDataTable.EyeTrackerFrameNumberTrialStart(i),1,'first');
                            trialDataTable.SampleStopTrial(i) = find(samplesDataTable.FileNumber' == trialDataTable.FileNumber(i) & samplesDataTable.RawFrameNumber'<=trialDataTable.EyeTrackerFrameNumberTrialStop(i),1,'last');
                        end

                end                
                
                % Build a column for the samples with the trial number
                samplesDataTable.TrialNumber = nan(size(samplesDataTable.FrameNumber));
                for i=1:height(trialDataTable)
                    idx = trialDataTable.SampleStartTrial(i):trialDataTable.SampleStopTrial(i);
                    samplesDataTable.TrialNumber(idx) = trialDataTable.TrialNumber(i);
                end
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % CALCULATE AVERAGE EYE MOVEMENT STATS FOR EACH TRIAL
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                LRdataVars = {'X' 'Y' 'T'};
                for i=1:length(LRdataVars)
                    if ( ~any(strcmp(samplesDataTable.Properties.VariableNames,['Left' LRdataVars{i}])) )
                        samplesDataTable.(['Left' LRdataVars{i}]) = nan(size(samplesDataTable.(['Right' LRdataVars{i}])));
                    end
                    if ( ~any(strcmp(samplesDataTable.Properties.VariableNames,['Right' LRdataVars{i}])) )
                        samplesDataTable.(['Right' LRdataVars{i}]) = nan(size(samplesDataTable.(['Left' LRdataVars{i}])));
                    end
                    
                    % average both eyes
                    samplesDataTable.(LRdataVars{i}) = mean(samplesDataTable{:,{['Left' LRdataVars{i}],['Right' LRdataVars{i}]}},2,'omitnan');
                end
                
                dataVars = { 'X' 'Y' 'T' 'LeftX' 'LeftY' 'LeftT' 'RightX' 'RightY' 'RightT'};
                
                if ( any(strcmp(samplesDataTable.Properties.VariableNames,'LeftBadData')) && any(strcmp(samplesDataTable.Properties.VariableNames,'RightBadData')) )
                    samplesDataTable.GoodData = ~samplesDataTable.LeftBadData & ~samplesDataTable.RightBadData;
                elseif ( any(strcmp(samplesDataTable.Properties.VariableNames,'LeftBadData')) )
                    samplesDataTable.GoodData = ~samplesDataTable.LeftBadData ;
                elseif ( any(strcmp(samplesDataTable.Properties.VariableNames,'RightBadData')) )
                    samplesDataTable.GoodData =  ~samplesDataTable.RightBadData;
                end
                isInTrial = ~isnan(samplesDataTable.TrialNumber);
                
                stats = grpstats(...
                    samplesDataTable(isInTrial,:), ...     % Selected rows of data
                    'TrialNumber', ...                  % GROUP VARIABLE
                    {'median' 'mean' 'std'}, ...        % Stats to calculate
                    'DataVars', dataVars );             % Vars to do stats on
                leftRightBadData = intersect({'LeftBadData', 'RightBadData'},samplesDataTable.Properties.VariableNames);
                stats2 = grpstats(...
                    samplesDataTable(isInTrial,:), ...     % Selected rows of data
                    'TrialNumber', ...                  % GROUP VARIABLE
                    {'sum'}, ...             % Vars to do stats on
                    'DataVars', {'GoodData', leftRightBadData{:}});	% Vars to do stats on
                stats2.Properties.VariableNames{'GroupCount'} = 'count_GoodSamples';
                
                samplerate = samplesDataTable.Properties.UserData.sampleRate;
                varsBeforeJoin = trialDataTable.Properties.VariableNames;
                trialDataTable = outerjoin(trialDataTable, stats,'Keys','TrialNumber','MergeKeys',true, 'LeftVariables', setdiff(trialDataTable.Properties.VariableNames, stats.Properties.VariableNames) );
                trialDataTable = outerjoin(trialDataTable, stats2,'Keys','TrialNumber','MergeKeys',true, 'LeftVariables', setdiff(trialDataTable.Properties.VariableNames, stats2.Properties.VariableNames) );
                
                trialDataTable.TotalGoodSec = trialDataTable.sum_GoodData/samplerate;
                trialDataTable.TotalSec = trialDataTable.count_GoodSamples/samplerate;
                
                % keep a list of the variables added so it is easier to do
                % states across trials
                trialDataTable.Properties.UserData.EyeTrackingPrepareTrialDataTableVariables = setdiff(trialDataTable.Properties.VariableNames, varsBeforeJoin);
            end
        end
        
        function sessionDataTable = PrepareSessionDataTableEyeTracking(this, sessionDataTable, options)
            return
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % CALCULATE AVERAGE EYE MOVEMENT ACROSS TRIALS
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            trialDataTable = this.Session.trialDataTable;
            varsToGroup = trialDataTable.Properties.UserData.EyeTrackingPrepareTrialDataTableVariables;
            
            % IT IS TOO SLOW NOW! 11/23/2021 CANNOT DO ALL THIS WHEN THERE
            % ARE MANY CONDITIONS OR MANY LEVELS
            return;

            % Also calculate average across conditions for each of the
            % values of each condition variables. 
            % Also create a variable that containes all the combinations of
            % all conditions and avarege across those
            ConditionVarsNames = {};
            condition = [];
            for i=1:length(this.Session.experimentDesign.ConditionVars) %TODO CHANGE
                if ( numel(this.Session.experimentDesign.ConditionVars(i).values)>1)
                    ConditionVarsNames{end+1} = this.Session.experimentDesign.ConditionVars(i).name;
                    if (isempty(condition) )
                        condition = string(trialDataTable{:,ConditionVarsNames(i)});
                    else
                        condition = strcat(condition,'_', string(trialDataTable{:,ConditionVarsNames(i)}));
                    end
                end
            end
            if ( ~isempty(ConditionVarsNames))
                trialDataTable.Condition = condition;
                ConditionVarsNames = horzcat({'Condition'}, ConditionVarsNames);
                for i=1:length(ConditionVarsNames)
                    st = grpstats(trialDataTable, ...     % Selected rows of data
                        ConditionVarsNames{i}, ...                  % GROUP VARIABLE
                        {'mean' 'std'}, ...        % Stats to calculate
                        'DataVars', varsToGroup);             % Vars to do stats on
                    b = unstack(st,setdiff(st.Properties.VariableNames, ConditionVarsNames{i}),ConditionVarsNames{i});
                    b.Properties.RowNames = {};
                    sessionDataTable(:,intersect(b.Properties.VariableNames, sessionDataTable.Properties.VariableNames)) = [];
                    sessionDataTable = horzcat(sessionDataTable, b);
                end
            end
        end

        function [analysisResults, samplesDataTable, trialDataTable, sessionDataTable]  = RunDataAnalysesEyeTracking(this, analysisResults, samplesDataTable, trialDataTable, sessionDataTable, options)
            
            updateTrialsAndSessionTables = false;
            
            if ( isfield(options,'Detect_Quik_and_Slow_Phases') && options.Detect_Quik_and_Slow_Phases )
                samplesDataTable = VOGAnalysis.DetectQuickPhases(samplesDataTable, options);
                samplesDataTable = VOGAnalysis.DetectSlowPhases(samplesDataTable, options);
                updateTrialsAndSessionTables = true;
            end
            
            if ( isfield(options,'Mark_Data') && options.Mark_Data )
                samplesDataTable = VOGDataExplorer.MarkData(samplesDataTable);
                updateTrialsAndSessionTables = true;
            end
            
            if ( 0 )
                T = samplesDataTable.Properties.UserData.sampleRate;
                analysisResults.SPV.Time = samplesDataTable.Time(1:T:(end-T/2));
                fields = {'LeftX', 'LeftY' 'LeftT' 'RightX' 'RightY' 'RightT'};
                
                t = samplesDataTable.Time;
                                
                %
                % calculate binocular spv
                %
                LRdataVars = {'X' 'Y' 'T'};
                for j =1:length(LRdataVars)
                    
                    [vleft, xleft] = VOGAnalysis.GetSPV_SimpleQP(t, samplesDataTable.(['Left' LRdataVars{j}]), samplesDataTable.QuickPhase );
                    [vright, xright] = VOGAnalysis.GetSPV_SimpleQP(t, samplesDataTable.(['Right' LRdataVars{j}]), samplesDataTable.QuickPhase );
                    
                    vmed = nanmedfilt(mean([vleft, vright],2),T,1/2,'omitnan');
                    samplesDataTable.(['SPV' LRdataVars{j}]) = vmed;
                end
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % CALCULATE AVERAGE SAMPLE PROPERTIES ACROSS TRIALS
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                
                dataVars = { 'SPVX' 'SPVY' 'SPVT'};
                isInTrial = ~isnan(samplesDataTable.TrialNumber);
                
                stats = grpstats(...
                    samplesDataTable(isInTrial,:), ...     % Selected rows of data
                    'TrialNumber', ...                  % GROUP VARIABLE
                    {'median' 'mean' 'std'}, ...        % Stats to calculate
                    'DataVars', dataVars );             % Vars to do stats on
                
                trialDataTable = outerjoin(trialDataTable, stats,'Keys','TrialNumber','MergeKeys',true, 'LeftVariables', setdiff(trialDataTable.Properties.VariableNames, stats.Properties.VariableNames) );
            end
            
            if (updateTrialsAndSessionTables)
                if ( width(this.TrialTable) > 5 ) % TODO: a bit ugly! 
                    % experiments that are imported behave a bit different
                    % than experiments that are run original with arume. 
                    ConditionVarsNames = this.TrialTable.Properties.VariableNames(6:end);
                else
                    ConditionVarsNames = this.Session.currentRun.pastTrialTable.Properties.VariableNames(6:end); 
                end
                condition = [];
                for i=1:length(ConditionVarsNames)
                    conditionVarLevels = categories(categorical(this.Session.currentRun.pastTrialTable{:,ConditionVarsNames{i}}));
                    if ( numel(conditionVarLevels)>1)
                        if (isempty(condition) )
                            condition = string(trialDataTable{:,ConditionVarsNames(i)});
                        else
                            condition = strcat(condition,'_', string(trialDataTable{:,ConditionVarsNames(i)}));
                        end
                    end
                end
                
                % add columns to quick and slow phases for trial number and
                % also the values of each condition variable that
                % corresponds with that trial
                [qp, sp] = VOGAnalysis.GetQuickAndSlowPhaseTable(samplesDataTable);
                qpDataVars = qp.Properties.VariableNames;
                spDataVars = sp.Properties.VariableNames;
                
                warning('off','MATLAB:table:RowsAddedNewVars');
                qp.TrialNumber = samplesDataTable.TrialNumber(qp.StartIndex);
                sp.TrialNumber = samplesDataTable.TrialNumber(sp.StartIndex);
                for i=1:numel(ConditionVarsNames)
                    if ( iscategorical(trialDataTable{qp.TrialNumber(~isnan(qp.TrialNumber)),ConditionVarsNames{i}}))
                        qp{~isnan(qp.TrialNumber),ConditionVarsNames{i}} =  categorical(trialDataTable{qp.TrialNumber(~isnan(qp.TrialNumber)),ConditionVarsNames{i}});
                    elseif iscellstr(trialDataTable{qp.TrialNumber(~isnan(qp.TrialNumber)),ConditionVarsNames{i}})
                        qp{~isnan(qp.TrialNumber),ConditionVarsNames{i}} =  categorical(trialDataTable{qp.TrialNumber(~isnan(qp.TrialNumber)),ConditionVarsNames{i}});
                    else
                        qp{:,ConditionVarsNames{i}} = nan(height(qp),1);
                        qp{~isnan(qp.TrialNumber),ConditionVarsNames{i}} =  trialDataTable{qp.TrialNumber(~isnan(qp.TrialNumber)),ConditionVarsNames{i}};
                    end
                    sp{~isnan(sp.TrialNumber),ConditionVarsNames{i}} =  categorical(trialDataTable{sp.TrialNumber(~isnan(sp.TrialNumber)),ConditionVarsNames{i}});
                end
                warning('on','MATLAB:table:RowsAddedNewVars');
                
                analysisResults.QuickPhases = qp;
                analysisResults.SlowPhases = sp;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % CALCULATE AVERAGE QP AND SP PROPERTIES FOR EACH TRIAL
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Add average properties of QP and SP to the trial table
                goodQP = qp(~isnan(qp.TrialNumber),:);
                stats = grpstats(goodQP, ...     % Selected rows of data
                    'TrialNumber', ...                  % GROUP VARIABLE
                    {'median' 'mean' 'std'}, ...        % Stats to calculate
                    'DataVars', qpDataVars);             % Vars to do stats on
                stats.Properties.VariableNames(2:end) = strcat('QP_', stats.Properties.VariableNames(2:end));
                trialDataTable = outerjoin(trialDataTable, stats,'Keys','TrialNumber','MergeKeys',true, 'LeftVariables', setdiff(trialDataTable.Properties.VariableNames, stats.Properties.VariableNames) );
                
                goodSP = sp(~isnan(sp.TrialNumber),:);
                stats = grpstats(goodSP, ...     % Selected rows of data
                    'TrialNumber', ...                  % GROUP VARIABLE
                    {'median' 'mean' 'std'}, ...        % Stats to calculate
                    'DataVars', spDataVars);             % Vars to do stats on
                stats.Properties.VariableNames(2:end) = strcat('SP_', stats.Properties.VariableNames(2:end));
                trialDataTable = outerjoin(trialDataTable, stats,'Keys','TrialNumber','MergeKeys',true, 'LeftVariables', setdiff(trialDataTable.Properties.VariableNames, stats.Properties.VariableNames) );
                
                trialDataTable.QP_Rate = trialDataTable.QP_GroupCount ./ trialDataTable.TotalGoodSec;
                
                
                

                % IT IS TOO SLOW NOW! 11/23/2021 CANNOT DO ALL THIS WHEN THERE
                % ARE MANY CONDITIONS OR MANY LEVELS
                return;

                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % CALCULATE AVERAGE QP AND SP PROPERTIES ACROSS TRIALS
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if ( ~isempty(condition) )
                    varsToGroup = strcat('QP_mean_', qpDataVars(~contains(qpDataVars,'Left') & ~contains(qpDataVars,'Right')));
                    varsToGroup = horzcat({'QP_Rate' 'median_SPVX','median_SPVY','median_SPVT'}, varsToGroup);
                    varsToGroup = intersect(varsToGroup, trialDataTable.Properties.VariableNames);

                    trialDataTable.Condition = condition; % combination of condition variables
                    ConditionVarsNames = horzcat({'Condition'}, ConditionVarsNames);
                    for i=1:length(ConditionVarsNames)
                        stats = grpstats(trialDataTable, ...     % Selected rows of data
                            ConditionVarsNames{i}, ...                  % GROUP VARIABLE
                            {'mean' 'std'}, ...        % Stats to calculate
                            'DataVars', varsToGroup);             % Vars to do stats on
                        conditionStats = unstack(stats,setdiff(stats.Properties.VariableNames, ConditionVarsNames{i}),ConditionVarsNames{i});
                        conditionStats.Properties.RowNames = {};
                        sessionDataTable(:,intersect(conditionStats.Properties.VariableNames, sessionDataTable.Properties.VariableNames)) = [];
                        sessionDataTable = horzcat(sessionDataTable, conditionStats);
                    end
                end
                
                cprintf('blue','++ EyeTracking :: Done with eye tracking analsyis\n');
            end
        end
                
        
        function markData(this, options)
            options.Prepare_For_Analysis_And_Plots = 0;
            options.Preclear_Samples_Table = 0;
            options.Preclear_Trial_Table = 0;
            options.Preclear_Session_Table = 0;
            options.Mark_Data = 1;
            
            this.Session.runAnalysis(options);
        end
    end

    methods(Access=public,Sealed=true) % 

        function [analysisResults, samplesDataTable, trialDataTable, sessionTable] = RunExperimentAnalysis(this, options)
            [samplesDataTable, trialDataTable, sessionTable]  = this.prepareTablesForAnalysis(options);
            analysisResults  = struct();
            [analysisResults, samplesDataTable, trialDataTable, sessionTable]  = this.RunDataAnalyses(analysisResults, samplesDataTable, trialDataTable, sessionTable, options);
        end
        
        function this = ExperimentDesign()
            className = class(this);
            this.Name = className(find(className=='.',1, 'last')+1:end);
        end

        function trialTable = GetTrialTable(this)
            trialTable = this.TrialTable;
        end
        
        function trialTableOptions = GetDefaultTrialTableOptions(this)
            % OBSOLTE METHOD use TrialTableBuilder

            % Trial sequence and blocking
            trialTableOptions = [];
            trialTableOptions.trialSequence       = 'Sequential';	% Sequential, Random, Random with repetition, ...
            trialTableOptions.trialAbortAction    = 'Repeat';     % Repeat, Delay, Drop
            trialTableOptions.trialsPerSession    = 1;
            trialTableOptions.trialsBeforeBreak   = 1;
            
            trialTableOptions.blockSequence       = 'Sequential';	% Sequential, Random, Random with repetition, ...numberOfTimesRepeatBlockSequence = 1;
            trialTableOptions.blocksToRun         = 1;
            trialTableOptions.blocks              = [];
            trialTableOptions.numberOfTimesRepeatBlockSequence = 1;
        end
            
        function trialTable = GetTrialTableFromConditions(this, conditionVars, trialTableOptions)
            % OBSOLETE!!! USE TrialTableBuilderInstead
            %
            % Create the matrix with all the possible combinations of
            % condition variables. Each combination is a condition
            % total number of conditions is the product of the number of
            % values of each condition variable
            nConditions = 1;
            for iVar = 1:length(conditionVars)
                nConditions = nConditions * length(conditionVars(iVar).values);
            end
            
            conditionMatrix = [];
            
            %-- recursion to create the condition matrix
            % for each variable, we repeat the previous matrix as many
            % times as values the current variable has and in each
            % repetition we add a new column with one of the values of the
            % current variable
            % example: var1 = {a b} var2 = {e f g}
            % step 1: matrix = [ a ;
            %                    b ];
            % step 2: matrix = [ a e ;
            %                    b e ;
            %                    a f ;
            %                    b f ;
            %                    a g ;
            %                    b g ];
            for iVar = 1:length(conditionVars)
                nValues(iVar) = length(conditionVars(iVar).values);
                conditionMatrix = [ repmat(conditionMatrix,nValues(iVar),1)  ceil((1:prod(nValues))/prod(nValues(1:end-1)))' ];
            end
            
            % if the blocks are empty add one that includes all the
            % conditions
            if ( isempty( trialTableOptions.blocks) )
                trialTableOptions.blocks = struct( 'fromCondition', 1, 'toCondition', size(conditionMatrix,1), 'trialsToRun', size(conditionMatrix,1)  );
                trialTableOptions.blocksToRun = 1;
            end
        
        
            blockSeqWithRepeats = [];
            for iRepeatBlockSequence = 1:trialTableOptions.numberOfTimesRepeatBlockSequence
        
                % generate the sequence of blocks, a total of
                % parameters.blocksToRun blocks will be run
                nBlocks = length(trialTableOptions.blocks);
                blockSeq = [];
                switch(trialTableOptions.blockSequence)
                    case 'Sequential'
                        blockSeq = mod( (1:trialTableOptions.blocksToRun)-1,  nBlocks ) + 1;
                    case 'Random'
                        [~, theBlocks] = sort( rand(1,trialTableOptions.blocksToRun) ); % get a random shuffle of 1 ... blocks to run
                        blockSeq = mod( theBlocks-1,  nBlocks ) + 1; % limit the random sequence to 1 ... nBlocks
                    case 'Random with repetition'
                        blockSeq = ceil( rand(1,trialTableOptions.blocksToRun) * nBlocks ); % just get random block numbers
                    case 'Manual'
                        blockSeq = [];
                        
                        while length(blockSeq) ~= trialTableOptions.blocksToRun
                            S.Block_Sequence = [1:trialTableOptions.blocksToRun];
                            S = StructDlg( S, ['Block Sequence'], [],  CorrGui.get_default_dlg_pos() );
                            blockSeq =  S.Block_Sequence;
                        end
                        %                     if length(parameters.manualBlockSequence) == parameters.blocksToRun;
                        %                         %                         blockSequence = parameters.manualBlockSequence;
                        %
                        %                     else
                        %                         disp(['Error with the manual block sequence. Please fix.']);
                        %                     end
                end
                blockSeq = [blockSeq;ones(size(blockSeq))*iRepeatBlockSequence];
                blockSeqWithRepeats = [blockSeqWithRepeats blockSeq];
            end
            blockSeq = blockSeqWithRepeats;
            
            futureConditions = [];
            for iblock=1:size(blockSeq,2)
                i = blockSeq(1,iblock);
                possibleConditions = trialTableOptions.blocks(i).fromCondition : trialTableOptions.blocks(i).toCondition; % the possible conditions to select from in this block
                nConditions = length(possibleConditions);
                nTrials = trialTableOptions.blocks(i).trialsToRun;
                
                switch( trialTableOptions.trialSequence )
                    case 'Sequential'
                        trialSeq = possibleConditions( mod( (1:nTrials)-1,  nConditions ) + 1 );
                    case 'Random'
                        [~, conditions] = sort( rand(1,nTrials) ); % get a random shuffle of 1 ... nTrials
                        conditionIndexes = mod( conditions-1,  nConditions ) + 1; % limit the random sequence to 1 ... nConditions
                        trialSeq = possibleConditions( conditionIndexes ); % limit the random sequence to fromCondition ... toCondition for this block
                    case 'Random with repetition' 
                        trialSeq = possibleConditions( ceil( rand(1,nTrials) * nConditions ) ); % nTrialss numbers between 1 and nConditions
                end
                futureConditions = cat(1,futureConditions, [trialSeq' ones(size(trialSeq'))*iblock  ones(size(trialSeq'))*i ones(size(trialSeq'))*blockSeq(2,iblock)] );
            end
            
            newTrialTable = table();
            newTrialTable.Condition = futureConditions(:,1);
            newTrialTable.BlockNumber = futureConditions(:,2);
            newTrialTable.BlockSequenceNumber = futureConditions(:,3);
            newTrialTable.BlockSequenceRepeat = futureConditions(:,4);
            newTrialTable.Session = ceil((1:height(newTrialTable))/min(height(newTrialTable), trialTableOptions.trialsPerSession))';
            
            variableTable = table();
            for i=1:height(newTrialTable)
                variables = [];
                for iVar=1:length(conditionVars)
                    varName = conditionVars(iVar).name;
                    varValues = conditionVars(iVar).values;
                    if iscell( varValues )
                        variables.(varName) = categorical(varValues(conditionMatrix(newTrialTable.Condition(i),iVar)));
                    else
                        variables.(varName) = varValues(conditionMatrix(newTrialTable.Condition(i),iVar));
                    end
                end
            
                variableTable = cat(1,variableTable,struct2table(variables,'AsArray',true));
            end
            
            trialTable = [newTrialTable variableTable];
            
            trialTable.Properties.UserData.conditionVars = conditionVars;
            trialTable.Properties.UserData.trialTableOptions = trialTableOptions;
        end
        
        function options = GetDefaultOptions(this)
            options = StructDlg(this.GetAnalysisOptionsDialog(),'',[],[],'off');
        end
        
        function UpdateExperimentOptions(this, newOptions)
            this.ExperimentOptions = newOptions;
        end
        
        function [dataTable, idx, selectedFilters] = FilterTableByConditionVariable(this, dataTable, Select_Conditions, columns, columnNames, DataToIncludeFilter1, DataToIncludeFilter2, datafilter)
            
            % if there are no filters just get everything
            if ( ~exist('DataToIncludeFilter1','var') )
                DataToIncludeFilter1 = 'All';
            end
            if ( ~exist('DataToIncludeFilter2','var') )
                DataToIncludeFilter2 = 'All';
            end

            % TODO: not great right now. But get all the columns from the
            % trial table that have condition variables. 
            if(size(this.Session.experimentDesign.TrialTable,2)>6) % if not imported
                ConditionVars = this.Session.experimentDesign.TrialTable.Properties.VariableNames(6:end);
            else % if imported (TODO: not very good)
                ConditionVars = this.Session.currentRun.pastTrialTable.Properties.VariableNames(6:end);
            end

            % get all the possible values of the condition variables. But
            % only if they have less than 10 possible values. Otherwise it
            % gets too cumbersome
            if (ischar(dataTable) )
                switch(dataTable)
                    case 'get_filters'
                        Select_Conditions = struct();
                        Select_Conditions.All = { {'0', '{1}'}};
                        for i=1:length(ConditionVars)
                            name = ConditionVars{i};
                            if ( isnumeric(this.Session.currentRun.pastTrialTable{:,ConditionVars{i}}))
                                values = unique(this.Session.currentRun.pastTrialTable{:,ConditionVars{i}});
                            else
                                values = categories(categorical(this.Session.currentRun.pastTrialTable{:,ConditionVars{i}}));
                            end
                            if ( numel(values) < 10 )
                                for j=1:numel(values)
                                    if (~ismissing(values(j)))
                                        Select_Conditions.(strcat(name, '_', genvarname(string(values(j))))) = { {'{0}', '1'}};
                                    end
                                end
                            end
                        end
                        dataTable = Select_Conditions;
                        return;
                end
            end
            
            % Add the All condition variable column to allow that filter to
            % work like all others and not need special code later
            dataTable.All = ones(height(dataTable),1);
            
            % Get the actual variable name and the value of that variable
            % for the group filter
            Select_ConditionsFilters = struct();
            Select_ConditionsFilters.All.VarName = 'All';
            Select_ConditionsFilters.All.VarValue = 1;
            for i=1:length(ConditionVars)
                name = ConditionVars{i};
                if ( isnumeric(this.Session.currentRun.pastTrialTable{:,ConditionVars{i}}))
                    values = unique(this.Session.currentRun.pastTrialTable{:,ConditionVars{i}});
                else
                    values = categories(categorical(this.Session.currentRun.pastTrialTable{:,ConditionVars{i}}));
                end
                if ( numel(values) < 10 )
                    for j=1:numel(values)
                        if (~ismissing(values(j)))
                            Select_ConditionsFilters.(strcat(name, '_', genvarname(string(values(j))))).VarName = name;
                            Select_ConditionsFilters.(strcat(name, '_', genvarname(string(values(j))))).VarValue = values(j);
                        end
                    end
                end
            end

            % Get the actual variable name and the value of that variable
            % for the data to include filters
            switch(DataToIncludeFilter1)
                case 'All'
                    DataToIncludeName1 = 'All';
                    DataToIncludeValue1 = 1;
                otherwise
                    DataToIncludeNameValue1 = strsplit(DataToIncludeFilter1,'_');
                    DataToIncludeName1 = DataToIncludeNameValue1{1};
                    DataToIncludeValue1 = DataToIncludeNameValue1{2};
            end
            switch(DataToIncludeFilter2)
                case 'All'
                    DataToIncludeName2 = 'All';
                    DataToIncludeValue2 = 1;
                otherwise
                    DataToIncludeNameValue2 = strsplit(DataToIncludeFilter2,'_');
                    DataToIncludeName2 = DataToIncludeNameValue2{1};
                    DataToIncludeValue2 = DataToIncludeNameValue2{2};
            end
            
            % find which filters have been selected
            selectedFilters = {};
            
            filters = fieldnames(Select_Conditions);
            for i=1:length(filters)
                if ( Select_Conditions.(filters{i}) )
                    selectedFilters{end+1} = filters{i};
                end
            end
            
            % create a table with one row per filter and in each row the
            % filter name and the filtered table
            sessionDataTable = dataTable;
            dataTable = table();
            idx = table();
            for i=1:length(selectedFilters)
                idxf = sessionDataTable.(Select_ConditionsFilters.(selectedFilters{i}).VarName) == Select_ConditionsFilters.(selectedFilters{i}).VarValue;
                idxf = idxf & sessionDataTable.(DataToIncludeName1) == DataToIncludeValue1 & sessionDataTable.(DataToIncludeName2) == DataToIncludeValue2;
                dataTable{i, {'Data' 'Condition' 'Idx'}} = {sessionDataTable(idxf,:), selectedFilters{i}, find(idxf)};
            end
            
            % get the specific components we want and if more than one
            % replicate the rows as many time as needed but with the data
            % from the particular component
            if ( exist('columns','var'))
                dataTableTemp = table();
                for j=1:height(dataTable)
                    props = dataTable(j,:);
                    
                    props = repmat(props, numel(columns), 1);
                    props.Component = columnNames';
                    for iComp=1:numel(columns)
                        props{iComp,'Data'} = {props.Data{iComp}.(columns{iComp})};
                    end
                    dataTableTemp = vertcat(dataTableTemp, props);
                end
                dataTable = dataTableTemp;
                
                
                dataTable.Component = categorical(cellstr(dataTable.Component));
            end

            % apply data filters. For each row get only the elements of the
            % table that meet a particular condition
            if ( exist('datafilter','var'))
                for j=1:height(dataTable)
                    data = sessionDataTable(dataTable.Idx{j},:);
                    rowsToKeep = eval(datafilter);
                    dataTable.Data{j}  = dataTable.Data{j}(rowsToKeep,:);
                    dataTable.Idx{j}  = dataTable.Idx{j}(rowsToKeep,:);
                end
            end

            dataTable.Condition = categorical(cellstr(dataTable.Condition));
        end
    end
    
    % --------------------------------------------------------------------
    %% PUBLIC and sealed METHODS ------------------------------------------
    % --------------------------------------------------------------------
    % to be called from gui or command line
    % --------------------------------------------------------------------
    methods(Sealed = true) % MAIN public methods
        
        %
        % Options to set at runtime, this options will appear as a dialog
        % when creating a new session. If one experiment inherits from another
        % one it is a good idea to first call GetExperimentDesignOptions from
        % the parent class to get the options and then add new ones.
        %
        % This options may also appear when importing a session and it is
        % possible that the parameters that want to be displaied in that
        % case are different
        function dlg = GetExperimentOptionsDialog( this, importing )
            if ( ~exist( 'importing', 'var') )
                importing = 0;
            end
            dlg = this.GetOptionsDialog(importing);
        end
        
        function init(this, session, options, importing)
            this.Session = session;
            if ( ~exist( 'importing', 'var') )
                importing = 0;
            end
            
            %-- Check if all the options are there, if not add the default
            % values. This is important to mantain past compatibility if
            % options are added in the future.
            % TODO deal with nested structures
            optionsDlg = this.GetOptionsDialog(importing);
            if ( ~isempty( optionsDlg ) && (isempty(options) || ~isempty(setdiff(fieldnames(optionsDlg), fieldnames(options)))) )
                
                defaultOptions = StructDlg(optionsDlg,'',[],[],'off');
                
                fields = fieldnames(defaultOptions);
                for i=1:length(fields)
                    if ( ~isfield(options, fields{i}))
                        options.(fields{i}) = defaultOptions.(fields{i});
                    end
                end
            end
            
            this.ExperimentOptions = options;
            
            newTrialTable = this.SetUpTrialTable();
            
            % TODO: Check trialTable
            
            this.TrialTable = newTrialTable;
        end
        
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
                                { 'n' 'c' 'q'}, ...
                                { 'Continue with next trial' 'Calibrate' 'Quit'} , [],[]);
                            
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
                                
                            
                            catch err
                                if ( streq(err.identifier, 'PSYCORTEX:USERQUIT' ) )
                                    thisTrialData.TrialResult = Enum.trialResult.QUIT;
                                else
                                    thisTrialData.TrialResult = Enum.trialResult.ERROR;
                                    thisTrialData.ErrorMessage = string(err.message);
                                    % display error
                                    
                                    beep
                                    cprintf('red', '\n')
                                    cprintf('red', '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n')
                                    cprintf('red', '!!!!!!!!!!!!! ARUME ERROR: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n')
                                    cprintf('red', '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n')
                                    disp(err.getReport);
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
        
        function abortExperiment(this)
            throw(MException('PSYCORTEX:USERQUIT', ''));
        end
    end
        
    methods (Access=private) % Private methods
        
        function PlaySound(this,trialResult)
            
            Enum = ArumeCore.ExperimentDesign.getEnum();
            
            %% make a sound for the end of the trial
            fs = 8000;
            T = 0.1; % 2 seconds duration
            t = 0:(1/fs):T;
            if ( trialResult == Enum.trialResult.CORRECT )
                f = 500;
            else
                f = 250;
            end
            y = sin(2*pi*f*t);
            sound(y, fs);
        end
    end
    
    methods ( Access = public ) % Eye tracking Plot methods
        
        function Plot_VOG_RawData(this)
            switch(this.ExperimentOptions.EyeTracker) 
                case 'OpenIris'
                    data = this.Session.rawDataTable;
                    VOGAnalysis.PlotRawTraces(data);
                case 'Fove'
                    data = this.Session.rawDataTable;
                    VOGAnalysis.PlotRawTraces(data,'Fove');
            end
        end
        
        function Plot_VOG_Data_Explorer(this)
            VOGDataExplorer.Open(this.Session.samplesDataTable);
        end
            
        function [out, options] = Plot_VOG_Traces(this, options)
            
            out = [];
            if ( nargin == 1 )
                [~, options] = this. Plot_VOG_Traces('get_defaults');
            end
            
            if ( ischar(options) )
                command = options;
                switch( command)
                    case 'get_options'
                        options = VOGAnalysis.PlotTraces('get_options');
                        options.Show_Trials = { {'0','{1}'} };
                        return;
                    case 'get_defaults'
                        optionsDlg = VOGAnalysis.PlotTraces('get_options');
                        options = StructDlg(optionsDlg,'',[],[],'off');
                        return
                end
            end
            
            % Start plot
            
            data = this.Session.samplesDataTable;
            
            h = VOGAnalysis.PlotTraces(data, options);
            
            % Add vertical lines for trial begin and end
            if ( options.Show_Trials )
                for i=1:length(h)
                    axes(h(i));
                    yl = get(h(i),'ylim');
                    
                    trialStarts = nan(size(data.Time));
                    trialStarts(this.Session.trialDataTable.SampleStartTrial-1) = yl(1);
                    trialStarts(this.Session.trialDataTable.SampleStartTrial) = yl(2);
                    
                    trialStops = nan(size(data.Time));
                    trialStops(this.Session.trialDataTable.SampleStopTrial-1) = yl(1);
                    trialStops(this.Session.trialDataTable.SampleStopTrial) = yl(2);
                    
                    plot(data.Time, trialStarts,'color',[0.8 0.8 0.8]);
                    plot(data.Time, trialStops,'--','color',[0.8 0.8 0.8]);
                end
            end
        end
        
        function Plot_VOG_DebugCleaning(this)
            
            rawdata = this.Session.rawDataTable;
            data = this.Session.samplesDataTable;
            
            VOGAnalysis.PlotCleanAndResampledData(rawdata,data);
        end
        
        function Plot_VOG_SaccadeTraces(this)
            data = this.Session.samplesDataTable;
            VOGAnalysis.PlotQuickPhaseDebug(data)
        end
                
        function [out, options] = PlotAggregate_VOG_QuickPhase_MainSequence(this, sessions, options)
            
            out = [];
            if ( nargin == 1 )
                options = this.Plot_VOG_MainSequence('get_defaults');
            end
            
            if ( ischar(sessions) )
                command = sessions;
                switch( command)
                    case 'get_options'
                        options = VOGAnalysis.PlotMainsequence('get_options');
                        options.Component = { '{XY}|X|Y|T|All|X and Y' };
                        filterNames = fieldnames(this.FilterTableByConditionVariable('get_filters'));
                        options.DataToInclude = {filterNames};
                        options.Select_Trial_Conditions = this.FilterTableByConditionVariable('get_filters');
                        options.Figures_Axes_Lines_Order = {{...
                            '{Sessions-Conditions-Components}' 'Sessions-Components-Conditions' ...
                            'Conditions-Sessions-Components' 'Conditions-Components-Sessions' ...
                            'Components-Sessions-Conditions'  'Components-Conditions-Sessions'}};
                        return;
                    case 'get_defaults'
                        optionsDlg = VOGAnalysis.PlotMainsequence('get_options');
                        options = StructDlg(optionsDlg,'',[],[],'off');
                        return
                end
            end

            % Pick the variables that corresponds with the options for the
            % components. This can be single or multiple. 
            componentNames = {options.Component};
                switch(options.Component)
                    case 'XY'
                        Xcomponents = {'Amplitude'};
                        Ycomponents = {'PeakSpeed'};
                    case 'X'
                        Xcomponents = {'X_Amplitude'};
                        Ycomponents = {'X_PeakSpeed'};
                    case 'Y'
                        Xcomponents = {'Y_Amplitude'};
                        Ycomponents = {'Y_PeakSpeed'};
                    case 'T'
                        Xcomponents = {'T_Amplitude'};
                        Ycomponents = {'T_PeakSpeed'};
                    case 'All'
                        Xcomponents = {'X_Amplitude' 'Y_Amplitude' 'T_Amplitude'};
                        Ycomponents = {'X_PeakSpeed' 'Y_PeakSpeed' 'T_PeakSpeed'};
                        componentNames  = {'Horizotal', 'Vertical', 'Torsional'};
                    case 'X and Y'
                        Xcomponents = {'X_Amplitude' 'Y_Amplitude'};
                        Ycomponents = {'X_PeakSpeed' 'Y_PeakSpeed'};
                        componentNames  = {'Horizotal', 'Vertical'};
                end

            Xallprops = table();
            for i=1:length(sessions)
                [sessionProps, ~] = this.FilterTableByConditionVariable(sessions(i).analysisResults.QuickPhases, options.Select_Trial_Conditions, Xcomponents, componentNames, options.DataToInclude);
                sessionProps.Session = categorical(cellstr(repmat(sessions(i).shortName,height(sessionProps),1)));
                Xallprops = vertcat(Xallprops, sessionProps);
            end
            
            Yallprops = table();
            for i=1:length(sessions)
                [sessionProps, ~] = this.FilterTableByConditionVariable(sessions(i).analysisResults.QuickPhases, options.Select_Trial_Conditions, Ycomponents, componentNames, options.DataToInclude);
                sessionProps.Session = categorical(cellstr(repmat(sessions(i).shortName,height(sessionProps),1)));
                Yallprops = vertcat(Yallprops, sessionProps);
            end
                        
            
            
%             out = VOGAnalysis.PlotMainsequence(options, xdata, ydata );
%             if (~isempty(legendText) )
%                 legend(out.forLegend, legendText,'box','off');
%             end
% %             title(['Main sequence - ', strrep(filters, '_', ' ')]);

                        
            nplot1 = [1 1 1 2 2 2 2 2 3 2 3 3 3 3 4 4 4 4 4 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5];
            nplot2 = [1 2 3 2 3 3 4 4 3 5 4 4 5 5 4 5 5 5 5 5 5 5 5 5 6 6 6 6 6 6 7 7 7 7 7 7];
            
            
            
            
            
            
            
                        
            COLUMNS = {'Session','Condition', 'Component'};
            switch(options.Figures_Axes_Lines_Order)
                case 'Sessions-Conditions-Components'
                    COLUMNS = {'Session','Condition', 'Component'};
                case 'Sessions-Components-Conditions'
                    COLUMNS = {'Session', 'Component','Condition'};
                case 'Conditions-Sessions-Components'
                    COLUMNS = {'Condition', 'Session', 'Component'};
                case 'Conditions-Components-Sessions'
                    COLUMNS = {'Condition', 'Component', 'Session'};
                case 'Components-Sessions-Conditions'
                    COLUMNS = {'Component', 'Session','Condition'};
                case 'Components-Conditions-Sessions'
                    COLUMNS = {'Component','Condition', 'Session'};
            end

            ELEMENTS = {};
            for i=1:3
                ELEMENTS{i} = unique(Xallprops.(COLUMNS{i}),'stable');
            end
            for i=1:length(ELEMENTS{1})
                figure('name',string(ELEMENTS{1}(i)))
                for j=1:length(ELEMENTS{2})
                    ax = subplot(nplot1(length(ELEMENTS{2})),nplot2(length(ELEMENTS{2})),j,'nextplot','add');
                    xdata = Xallprops(Xallprops.(COLUMNS{1})==ELEMENTS{1}(i) & Xallprops.(COLUMNS{2})==ELEMENTS{2}(j),:);
                    ydata = Yallprops(Yallprops.(COLUMNS{1})==ELEMENTS{1}(i) & Yallprops.(COLUMNS{2})==ELEMENTS{2}(j),:);                    
                    title = char(strcat('Main seq. - ', strrep(string(ELEMENTS{2}(j)), '_', ' ')));
                    out = VOGAnalysis.PlotMainsequence(ax, xdata.Data{1}, ydata.Data{1}, title);
                end
            end
            legend(out.forLegend, strrep(string(ELEMENTS{3}),'_', ' '),'box','off');

        end
          


        function [out, options] = PlotAggregate_VOG_QuickPhase_Distribution(this, sessions, options)
            
            out = [];
            if ( nargin == 1 )
                % if passing no parameters get the default options
                options = this.PlotAggregate_VOG_QuickPhase_Distribution('get_defaults');
            end
            
            if ( ischar(sessions) )
                command = sessions;
                switch( command)
                    case 'get_options'
                        options = VOGAnalysis.PlotHistogram('get_options');
                        options.Feature =  {'{Amplitude}|PeakSpeed|Displacement|Direction'};
                        options.Component = { '{XY}|X|Y|T|All|X and Y' };
                        filterNames = fieldnames(this.FilterTableByConditionVariable('get_filters'));
                        options.DataToInclude = {filterNames};
                        options.Select_Trial_Conditions = this.FilterTableByConditionVariable('get_filters');
                        options.Figures_Axes_Lines_Order = {{...
                            '{Sessions-Conditions-Components}' 'Sessions-Components-Conditions' ...
                            'Conditions-Sessions-Components' 'Conditions-Components-Sessions' ...
                            'Components-Sessions-Conditions'  'Components-Conditions-Sessions'}};
                        return;
                    case 'get_defaults'
                        optionsDlg = VOGAnalysis.PlotAggregate_VOG_QuickPhase_Distribution('get_options');
                        options = StructDlg(optionsDlg,'',[],[],'off');
                        return
                end
            end
            
            % Get the right units for the labels depending on the feature
            % thatis being plotted
            switch(options.Feature)
                case 'Amplitude'
                    units = 'Deg';
                case 'PeakSpeed'
                    units = 'Deg/s';
                case 'Displacement' 
                    units = 'Deg';
                case 'Direction'
                    units = 'Deg';
            end

            % Pick the variables that corresponds with the options for the
            % components. This can be single or multiple. 
            switch(options.Component)
                case 'XY'
                    components = {options.Feature};
                    componentNames = {'Polar'};
                case 'X'
                    components = {['X_' options.Feature]};
                    componentNames = {'Horizontal'};
                case 'Y'
                    components = {['Y' options.Feature]};
                    componentNames = {'Vertical'};
                case 'T'
                    components = {['T' options.Feature]};
                    componentNames = {'Torsion'};
                case 'All'
                    components = {['X_' options.Feature], ['Y_' options.Feature], ['T_' options.Feature]};
                    componentNames = {'Horizontal', 'Vertical', 'Torsion'};
                case 'X and Y'
                    components = {['X_' options.Feature], ['Y_' options.Feature]};
                    componentNames = {'Horizontal', 'Vertical'};
            end
            
            allprops = table();
            for i=1:length(sessions)
                [sessionProps, ~] = this.FilterTableByConditionVariable(sessions(i).analysisResults.QuickPhases, options.Select_Trial_Conditions, components, componentNames, options.DataToInclude);
                sessionProps.Session = categorical(cellstr(repmat(sessions(i).shortName,height(sessionProps),1)));
                allprops = vertcat(allprops, sessionProps);
            end
                        
            nplot1 = [1 1 1 2 2 2 2 2 3 2 3 3 3 3 4 4 4 4 4 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5];
            nplot2 = [1 2 3 2 3 3 4 4 3 5 4 4 5 5 4 5 5 5 5 5 5 5 5 5 6 6 6 6 6 6 7 7 7 7 7 7];
            

            COLUMNS = {'Session','Condition', 'Component'};
            switch(options.Figures_Axes_Lines_Order)
                case 'Sessions-Conditions-Components'
                    COLUMNS = {'Session','Condition', 'Component'};
                case 'Sessions-Components-Conditions'
                    COLUMNS = {'Session', 'Component','Condition'};
                case 'Conditions-Sessions-Components'
                    COLUMNS = {'Condition', 'Session', 'Component'};
                case 'Conditions-Components-Sessions'
                    COLUMNS = {'Condition', 'Component', 'Session'};
                case 'Components-Sessions-Conditions'
                    COLUMNS = {'Component', 'Session','Condition'};
                case 'Components-Conditions-Sessions'
                    COLUMNS = {'Component','Condition', 'Session'};
            end

            ELEMENTS = {};
            for i=1:3
                ELEMENTS{i} = unique(allprops.(COLUMNS{i}),'stable');
            end
            for i=1:length(ELEMENTS{1})
                figure('name',string(ELEMENTS{1}(i)))
                for j=1:length(ELEMENTS{2})
                    ax = subplot(nplot1(length(ELEMENTS{2})),nplot2(length(ELEMENTS{2})),j);
                    xdata = allprops(allprops.(COLUMNS{1})==ELEMENTS{1}(i) & allprops.(COLUMNS{2})==ELEMENTS{2}(j),:);
                    title = strcat(options.Feature, ' distribution - ', strrep(string(ELEMENTS{2}(j)), '_', ' '));
                    xlab = [options.Feature '(' units ')'];
                    out = VOGAnalysis.PlotHistogram(ax, options, xdata.Data, [], title, xlab );
                end
            end
            legend(out.forLegend, strrep(string(ELEMENTS{3}),'_', ' '),'box','off');

        end
        
        
        function [out, options] = PlotAggregate_VOG_QuickPhase_Polar_Distribution(this, sessions, options)
            
            out = [];
            if ( nargin == 1 )
                % if passing no parameters get the default options
                options = this.PlotAggregate_VOG_QuickPhase_Distribution('get_defaults');
            end
            
            if ( ischar(sessions) )
                command = sessions;
                switch( command)
                    case 'get_options'
                        options = struct();%VOGAnalysis.PlotHistogram('get_options');
                        options.Number_of_bins = 36;
                        options.Feature =  {'{Direction}|NOTIMPLEMENTED'};
                        options.Component = { '{Both eyes combined}|Left eye|Right eye|Left and right eye' };
                        options.Normalize = { '{No}|By area|By max' };
                        options.Average = { {'{0}','1'} };
                        filterNames = fieldnames(this.FilterTableByConditionVariable('get_filters'));
                        options.Trials_to_Include_Filter1 = {filterNames};
                        options.Trials_to_Include_Filter2 = {filterNames};
                        options.Trial_Groups_To_Plot = this.FilterTableByConditionVariable('get_filters');
                        options.Figures_Axes_Lines_Order = {{...
                            '{Sessions-Groups-Components}' 'Sessions-Components-Groups' ...
                            'Groups-Sessions-Components' 'Groups-Components-Sessions' ...
                            'Components-Sessions-Groups'  'Components-Groups-Sessions'}};
                        options.QuickPhase_Filter = {{...
                            '{data.Amplitude > 0}' ...
                            'data.Amplitude <= 1' 'data.Amplitude <= 3' 'data.Amplitude <= 5' 'data.Amplitude < 10' ...
                            'data.Amplitude > 1' 'data.Amplitude > 3' 'data.Amplitude > 5' 'data.Amplitude > 10'
                            }};
                        return;
                    case 'get_defaults'
                        optionsDlg = VOGAnalysis.PlotAggregate_VOG_QuickPhase_Distribution('get_options');
                        options = StructDlg(optionsDlg,'',[],[],'off');
                        return
                end
            end
            
            % Get the right units for the labels depending on the feature
            % thatis being plotted
            switch(options.Feature)
                case 'Direction'
                    units = 'Deg';
            end

            % Pick the variables that corresponds with the options for the
            % components. This can be single or multiple. 
            switch(options.Component)
                case 'Both eyes combined'
                    options.Component = { '{}|||' };
                    components = {options.Feature};
                    componentNames = {'Both eyes'};
                case 'Left eye'
                    components = {['Left_' options.Feature]};
                    componentNames = {'Left eye'};
                case 'Right eye'
                    components = {['Right_' options.Feature]};
                    componentNames = {'Right eye'};
                case 'Left and right eye'
                    components = {['Left_' options.Feature], ['Right_' options.Feature]};
                    componentNames = {'Left eye', 'Right eye'};
            end
            
            allprops = table();
            for i=1:length(sessions)
                [sessionProps, ~] = this.FilterTableByConditionVariable(sessions(i).analysisResults.QuickPhases, options.Trial_Groups_To_Plot, components, componentNames, options.Trials_to_Include_Filter1, options.Trials_to_Include_Filter2, options.QuickPhase_Filter);
                sessionProps.Session = categorical(cellstr(repmat(sessions(i).shortName,height(sessionProps),1)));
                allprops = vertcat(allprops, sessionProps);
            end
                        
            nplot1 = [1 1 1 2 2 2 2 2 3 2 3 3 3 3 4 4 4 4 4 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5];
            nplot2 = [1 2 3 2 3 3 4 4 3 5 4 4 5 5 4 5 5 5 5 5 5 5 5 5 6 6 6 6 6 6 7 7 7 7 7 7];
            

            COLUMNS = {'Session','Condition', 'Component'};
            switch(options.Figures_Axes_Lines_Order)
                case 'Sessions-Groups-Components'
                    COLUMNS = {'Session','Condition', 'Component'};
                case 'Sessions-Components-Groups'
                    COLUMNS = {'Session', 'Component','Condition'};
                case 'Groups-Sessions-Components'
                    COLUMNS = {'Condition', 'Session', 'Component'};
                case 'Groups-Components-Sessions'
                    COLUMNS = {'Condition', 'Component', 'Session'};
                case 'Components-Sessions-Groups'
                    COLUMNS = {'Component', 'Session','Condition'};
                case 'Components-Groups-Sessions'
                    COLUMNS = {'Component','Condition', 'Session'};
            end

            if ( ~options.Average)
                ELEMENTS = {};
                for i=1:3
                    ELEMENTS{i} = unique(allprops.(COLUMNS{i}),'stable');
                end
                hforLegend = [];
                for i=1:length(ELEMENTS{1})
                    figure('name',string(ELEMENTS{1}(i)),'color','w')
                    for j=1:length(ELEMENTS{2})
                        ax1 = subplot(nplot1(length(ELEMENTS{2})),nplot2(length(ELEMENTS{2})),j);
                        ax = polaraxes('Units',ax1.Units,'Position',ax1.Position, 'nextplot','add'); % https://www.mathworks.com/matlabcentral/answers/443441-can-i-plot-multiple-polar-histograms-together
                        delete(ax1);
                        xdata = allprops(allprops.(COLUMNS{1})==ELEMENTS{1}(i) & allprops.(COLUMNS{2})==ELEMENTS{2}(j),:);
                        tit = strcat(options.Feature, ' distribution - ', strrep(string(ELEMENTS{2}(j)), '_', ' '));
                        xlab = [options.Feature '(' units ')'];
                        for k=1:length(xdata.Data)
                            angles = rad2deg(xdata.Data{k});
                            binsize = 360/options.Number_of_bins;

                            binedges = (-180-binsize/2):binsize:(180-binsize/2); % shift the bins to have one bin centered in zero
                            bincenters = (binedges(1:end-1) + binedges(2:end))/2;

                            angles(angles>max(binedges)) = -360+angles(angles>max(binedges)); % so the circular binning works with our bin shift

                            h = histcounts(angles, binedges);
                            switch(options.Normalize)
                                case 'By area'
                                    h = h/sum(h);
                                case 'By max'
                                    h = h/max(h);
                            end

                            radBinsCenter = deg2rad(bincenters);

                            h = polarplot(ax, radBinsCenter([1:end 1]), h([1:end 1]));

                            hforLegend(k) = h;

                        end

                        title(tit);
                    end
                end
                ll = legend(hforLegend, strrep(string(ELEMENTS{3}),'_', ' '),'box','off');
                set(ll,'Location','northeast');
            else
                ELEMENTS = {};
                for i=1:3
                    ELEMENTS{i} = unique(allprops.(COLUMNS{i}),'stable');
                end
                hforLegend = [];
                figure('color','w')
                for i=1:length(ELEMENTS{1})
                    tit = ['Average ' options.Feature ' distribution - ' strrep(char(ELEMENTS{1}(i)), '_', ' ')];
                    ax1 = subplot(nplot1(length(ELEMENTS{1})),nplot2(length(ELEMENTS{1})),i);
                    ax = polaraxes('Units',ax1.Units,'Position',ax1.Position, 'nextplot','add'); % https://www.mathworks.com/matlabcentral/answers/443441-can-i-plot-multiple-polar-histograms-together
                    delete(ax1);
                    for j=1:length(ELEMENTS{2})
                        xdata = allprops(allprops.(COLUMNS{1})==ELEMENTS{1}(i) & allprops.(COLUMNS{2})==ELEMENTS{2}(j),:);
                        allH = [];
                        for k=1:length(xdata.Data)
                            angles = rad2deg(xdata.Data{k});
                            binsize = 360/options.Number_of_bins;

                            binedges = (-180-binsize/2):binsize:(180-binsize/2); % shift the bins to have one bin centered in zero
                            bincenters = (binedges(1:end-1) + binedges(2:end))/2;

                            angles(angles>max(binedges)) = -360+angles(angles>max(binedges)); % so the circular binning works with our bin shift

                            h = histcounts(angles, binedges);
                            switch(options.Normalize)
                                case 'By area'
                                    h = h/sum(h);
                                case 'By max'
                                    h = h/max(h);
                            end

                            radBinsCenter = deg2rad(bincenters);

                            allH = vertcat(allH, h);
                        end

                        havg = mean(allH, 1);
                        h = polarplot(ax, radBinsCenter([1:end 1]), havg([1:end 1]),'linewidth',2);

                        xlab = [options.Feature '(' units ')'];
                        hforLegend(j) = h;
                        title(tit);
                    end
                end
                ll = legend(hforLegend, strrep(string(ELEMENTS{2}),'_', ' '),'box','off');
                set(ll,'Location','northeast');
            end

        end
        
        
        function [out, options] = PlotAggregate(this, sessions, options)
            
        end
    end
    
    methods ( Static = true ) % Static methods
        
        function options = GetDefaultExperimentOptions(experiment, importing)
            if ( ~exist('importing','var') )
                importing = 0;
            end
            experiment = ArumeCore.ExperimentDesign.Create(experiment);
            optionsDlg = experiment.GetExperimentOptionsDialog( importing);
            options = StructDlg(optionsDlg,'',[],[],'off');
        end
            
        function experimentList = GetExperimentList()
            experimentList = {};
            
            % It is important to not navigate the expPackage object in the
            % loop. It is vry slow. Keep the objects in their own variables
            % and used them like that.
            classList = meta.package.fromName('ArumeExperimentDesigns').ClassList;
            for i=1:length(classList)
                c = classList(i);
                if (~c.Abstract)
                    experimentList{end+1} = strrep( c.Name, 'ArumeExperimentDesigns.','');
                end
            end
        end
        
        function experiment = Create(experimentName)
            
            if ( exist( ['ArumeExperimentDesigns.' experimentName],  'class') )
                % Create the experiment design object
                experiment = ArumeExperimentDesigns.(experimentName)();
            elseif ( exist( ['AlconExperimentDesigns.' experimentName],  'class') )
                % Create the experiment design object
                experiment = AlconExperimentDesigns.(experimentName)();
            else
                % Create the experiment design object
                experiment = ArumeCore.ExperimentDesign();
            end
        end
        
        function Enum = getEnum()
            % -- possible trial results
            Enum.trialResult.CORRECT = categorical(cellstr('CORRECT')); % Trial finished correctly
            Enum.trialResult.ABORT = categorical(cellstr('ABORT'));   % Trial not finished, wrong key pressed, subject did not fixate, etc
            Enum.trialResult.ERROR = categorical(cellstr('ERROR'));   % Error during the trial
            Enum.trialResult.QUIT = categorical(cellstr('QUIT'));    % Escape was pressed during the trial
            Enum.trialResult.SOFTABORT = categorical(cellstr('SOFTABORT')); % Like an abort but does not go to hitkey to continue
            Enum.trialResult.PossibleResults = [...
                Enum.trialResult.CORRECT ...
                Enum.trialResult.ABORT ...
                Enum.trialResult.ERROR ...
                Enum.trialResult.QUIT ...
                Enum.trialResult.SOFTABORT]';
                
            
            % -- useful key codes
            try
                
                KbName('UnifyKeyNames');
                Enum.keys.SPACE     = KbName('space');
                Enum.keys.ESCAPE    = KbName('ESCAPE');
                Enum.keys.RETURN    = KbName('return');
                % Enum.keys.BACKSPACE = KbName('backspace');
                %
                %             Enum.keys.TAB       = KbName('tab');
                %             Enum.keys.SHIFT     = KbName('shift');
                %             Enum.keys.CONTROL   = KbName('control');
                %             Enum.keys.ALT       = KbName('alt');
                %             Enum.keys.END       = KbName('end');
                %             Enum.keys.HOME      = KbName('home');
                
                Enum.keys.LEFT      = KbName('LeftArrow');
                Enum.keys.UP        = KbName('UpArrow');
                Enum.keys.RIGHT     = KbName('RightArrow');
                Enum.keys.DOWN      = KbName('DownArrow');
            catch
            end
            
        end
        
    end
end

