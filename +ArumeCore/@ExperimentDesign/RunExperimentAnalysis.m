function [analysisResults, samplesDataTable, trialDataTable, sessionTable] = RunExperimentAnalysis(this, options)
Enum = ArumeCore.ExperimentDesign.getEnum();

if ( isempty(  this.Session.currentRun ) )
    return;
end

if ( options.Prepare_For_Analysis_And_Plots )

    %% 1) Create the basic trial data table (without custom experiment stuff)
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

    %% 2) Prepare the sample data table

    [samplesDataTable, cleanedData, calibratedData, rawData] = this.EyeTrackingLoadSamplesData(options);
    this.Session.WriteVariableIfNotEmpty(rawData,'rawDataTable');
    this.Session.WriteVariableIfNotEmpty(cleanedData,'cleanedData');
    this.Session.WriteVariableIfNotEmpty(calibratedData,'calibratedData');

    cprintf('blue', '++ ARUME::Done with samplesDataTable.\n');

    %% 3) Prepare the trial data table
    if ( ~isempty(samplesDataTable) )
        [sampleStartTrial, sampleStopTrial, trialDataTable, samplesDataTable] = this.EyeTrackingSyncTrialsAndSamples(trialDataTable, samplesDataTable,  options);
        trialDataTable.SampleStartTrial = sampleStartTrial;
        trialDataTable.SampleStopTrial = sampleStopTrial;
        % Build a column for the samples with the trial number
        samplesDataTable.TrialNumber = nan(size(samplesDataTable.FrameNumber));
        for i=1:height(trialDataTable)
            idx = trialDataTable.SampleStartTrial(i):trialDataTable.SampleStopTrial(i);
            samplesDataTable.TrialNumber(idx) = trialDataTable.TrialNumber(i);
        end
    end

    [trialDataTable] = this.EyeTrackingGetTrialStats(trialDataTable, samplesDataTable,  options);
    cprintf('blue', '++ ARUME::Done with trialDataTable.\n');

    %% 4) Prepare session data table
    newSessionDataTable = this.Session.GetBasicSessionDataTable();
    newSessionDataTable = this.EyeTrackingGetSessionStats(newSessionDataTable, options);
    newSessionDataTable.LastAnalysisDateTime = char(string(datetime('now')));

    optionsf = FlattenStructure(options); % eliminate strcuts with the struct so it can be made into a row of a table
    opts = fieldnames(optionsf);
    s = this.GetExperimentOptionsDialog(1);
    for i=1:length(opts)
        varName = ['AnalysisOption_' opts{i}];

        % If too long, keep the END of the name
        if length(varName) > namelengthmax
            warning('Variable name too long. Truncating end of: %s', varName);
            % Keep the last part that fits within MATLAB's limit
            varName = varName(end-namelengthmax+1:end);
        end

        % Assign to table
        if isempty(optionsf.(opts{i}))
            newSessionDataTable.(varName) = {''};
        elseif ~ischar(optionsf.(opts{i})) && numel(optionsf.(opts{i})) <= 1
            newSessionDataTable.(varName) = optionsf.(opts{i});
        elseif isfield(s, opts{i}) && iscell(s.(opts{i})) && iscell(s.(opts{i}){1}) && length(s.(opts{i}){1}) > 1
            newSessionDataTable.(varName) = categorical(cellstr(optionsf.(opts{i})));
        elseif ~ischar(optionsf.(opts{i})) && numel(optionsf.(opts{i})) > 1
            newSessionDataTable.(varName) = {optionsf.(opts{i})};
        else
            newSessionDataTable.(varName) = string(optionsf.(opts{i}));
        end
    end

    sessionTable = newSessionDataTable;

    analysisResults  = struct();
end

%% 5) Run analysis for the experiment
[analysisResults, samplesDataTable, trialDataTable, sessionTable]  = this.RunDataAnalyses(analysisResults, samplesDataTable, trialDataTable, sessionTable, options);

[analysisResults, samplesDataTable, trialDataTable, sessionTable]  = this.RunDataAnalysesEyeTracking(analysisResults, samplesDataTable, trialDataTable, sessionTable, options);
            


%% 6) Save data to disk
this.Session.WriteVariableIfNotEmpty(samplesDataTable,'samplesDataTable');
this.Session.WriteVariableIfNotEmpty(trialDataTable,'trialDataTable');
this.Session.WriteVariableIfNotEmpty(sessionTable,'sessionDataTable');

% save the fields of AnalysisResults into separate variables
if ( isstruct(analysisResults))
    fields=fieldnames(analysisResults);
    for i=1:length(fields)
        field = fields{i};
        this.Session.WriteVariableIfNotEmpty(analysisResults.(field),['AnalysisResults_' field]);
    end
else
    this.Session.WriteVariableIfNotEmpty(analysisResults,'AnalysisResults');
end

end