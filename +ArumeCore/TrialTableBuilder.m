classdef TrialTableBuilder < handle
    %TRIALTABLEBUILDER Class to help build trial tables for psychophysical
    % experiments based on condition variables.
    %
    % See ArumeHardware.TrialTableBuilder.Demo for an example of how to run it
    %
    % Example:
    %
    % t = ArumeCore.TrialTableBuilder();
    %
    % t.AddConditionVariable( 'Distance', ["near" "far"]);
    % t.AddConditionVariable( 'Direction', {'Left' 'Right' });
    % t.AddConditionVariable( 'Tilt', [-30 0 30]);
    %
    % Add three blocks. One with all the upright trials, one with the rest,
    % and another one with upright trials. Running only one repeatition of
    % each upright trial and 3 repeatitions of the other trials,
    % t.AddBlock(find(t.ConditionTable.Tilt==0), 1);
    % t.AddBlock(find(t.ConditionTable.Tilt~=0), 3);
    % t.AddBlock(find(t.ConditionTable.Tilt==0), 1);
    %
    % trialSequence = 'Random';
    % blockSequence =  'Sequential';
    % blockSequenceRepeatitions = 2;
    % abortAction = 'Repeat';
    % trialsPerSession = 10;
    % trialTable = t.GenerateTrialTable(trialSequence, blockSequence, blockSequenceRepeatitions, abortAction,trialsPerSession);

    properties (SetAccess=private)

        % Table of conditions with columns: Condition, BlockNumber,
        % BlockSequenceNumber, BlockSequenceRepeat, Session
        ConditionTable

        % Table of condition variables with columns: Name (string), Values
        % (cell)
        ConditionVariables = table( 'Size', [0 2],  ...
            'VariableNames', {'Name' 'Values'}, ...
            'VariableTypes', {'string', 'cell'});

        % Table of blocks with columns: ConditionsIncluded (cell),
        % TrialsPerCondition (double)end
        Blocks = table( 'Size', [0 2],  ...
            'VariableNames', {'ConditionsIncluded' 'TrialsPerCondition'}, ...
            'VariableTypes', {'cell', 'double'});
    end

    methods % GETTER METHODS
        function ConditionTable = get.ConditionTable(this)
            ConditionTable = ArumeCore.TrialTableBuilder.MakeConditionTable(this.ConditionVariables);
        end
    end

    methods
        function AddConditionVariable(this, name, values)
            %
            % Adds a new condition variable. A condition variable is a
            % variable that have one among a set of different values during
            % a experiment. It will be use to create all the possible
            % conditions as combinations of the values of the all the
            % condition variables.
            %
            %   AddConditionVariable( name, values)
            %
            %   Params
            %
            %       name: name of the variable
            %       values: possible values of the variable, usually an
            %           array of numbers or an array of strings. But it can
            %           also be a 1 Dimensional cell containing anything
            %
            %
            %   Example:
            %
            %   t = ArumeCore.TrialTableBuilder();
            %   t.AddConditionVariable( "Direction", ["Left" "Right"])
            %   t.AddConditionVariable( "Distance", {'near' 'far'})
            %   t.AddConditionVariable( "Position", [-20 0 20])
            %   t.AddConditionVariable( "TargetPosition", {[0 0] [0 10] [10 0] [10 10] [-10 10] [-10 -10] [10 -10]})
            %
            %
            this.ConditionVariables(end+1, ["Name" "Values"]) = table(string(name),{values});
        end

        function AddBlock(this, conditionsIncluded, trialsPerCondition)
            % Adds a new block of conditions to help with flexibility on
            % the randomization and structure of the experiment.
            %
            %   AddBlock(this, conditionsIncluded, trialsPerCondition)
            %
            %   Params
            %
            %       conditionsIncluded: array of condition numbers to be
            %           included in this block
            %       trialsPerCondition: how many times to run a trial for
            %           each condition in this block.
            %
            %
            %   Example:
            %
            %   t = ArumeCore.TrialTableBuilder();
            %   t.AddConditionVariable( "Distance", ["near" "far"])
            %   t.AddConditionVariable( "Position", [-30 0 30])
            %
            %
            % Add one block that will run only trials for task one
            % condition but 3 repeatitions of all conditions
            %
            % firstConditionNear = find(t.ConditionTable.Distance=='near',1);
            % lastConditionNear = find(t.ConditionTable.Distance=='near',1,'last');
            % t.AddBlock(firstConditionNear:lastConditionNear, 3);
            %
            % % Add another block with the rest of the conditions
            % % that will run only once per trial
            % t.AddBlock((lastConditionNear+1):height(t.ConditionTable), 1);
            %
            % % Add another block with only -30 tilt conditions
            % t.AddBlock(find(t.ConditionTable.Position==-30), 3);
            %
            %
            if ( iscolumn(conditionsIncluded))
                conditionsIncluded = conditionsIncluded';
            end
            this.Blocks(end+1, ["ConditionsIncluded" "TrialsPerCondition"]) = table({conditionsIncluded}, trialsPerCondition);
        end

        function trialTable = GenerateTrialTable(this, trialSequence, blockSequence, numberOfTimesRepeatBlockSequence, trialAbortAction, trialsPerSession)
            %
            % Generates the trial table after configuring the condition
            % variables and the blocks
            %
            %   GenerateTrialTable(this, trialSequence, blockSequence, numberOfTimesRepeatBlockSequence, trialAbortAction, trialsPerSession)
            %
            % Params
            %       trialSequence: Sequencing (randomization) of trials within a
            %           block. Possible values:
            %               'Sequential': Conditions run in sequential
            %                   order
            %               'Random': Ensures balanced number of trials per
            %                   condition
            %               'Random with repetition': Just a random
            %                   selection of the conditions
            %       blockSequence: Sequencing (randomization) of blocks.
            %           Possible values:
            %               'Sequential': Conditions run in sequential
            %                   order
            %               'Random': Shuffles the order of the blocks
            %       numberOfTimesRepeatBlockSequence: number of times the
            %           entire block sequence is repeated
            %       trialAbortAction: what to do when a trial is aborted.
            %           Possible values:
            %               'Repeat': Repeat the trial again.
            %               'Delay': Repeat the trial sometime in the
            %                   future but within the same block.
            %               'Drop': Ignore the trial and move on with the
            %                   nextone
            %       trialsPerSession: how many trials to run in each
            %           session part of an experiment. In arume, the experiment
            %           will automatically exit when all the trials of a
            %           session are run.
            %
            % Example:
            %
            %
            % trialSequence = 'Random';
            % blockSequence =  'Sequential';
            % blockSequenceRepeatitions = 2;
            % abortAction = 'Repeat';
            % trialsPerSession = 10;
            % trialTable = t.GenerateTrialTable(trialSequence, blockSequence, blockSequenceRepeatitions, abortAction,trialsPerSession);



            conditionTable = this.ConditionTable;
            if ( isempty(conditionTable))
                vars = this.ConditionVariables;
                vars(end+1, ["Name" "Values"]) = table("Var1",{1});
                conditionTable = ArumeCore.TrialTableBuilder.MakeConditionTable(vars);
            end
            blocks = this.Blocks;

            % if the blocks are empty add one that includes all the
            % conditions
            if ( isempty( blocks) )
                blocks(end+1, ["ConditionsIncluded" "TrialsPerCondition"]) = table({1:height(conditionTable)}, 1);
            end


            if (~exist('trialSequence','var'))
                trialSequence = 'Sequential';
            end
            if (~exist('blockSequence','var'))
                blockSequence = 'Sequential';
            end

            if (~exist('numberOfTimesRepeatBlockSequence','var'))
                numberOfTimesRepeatBlockSequence = 1;
            end
            if (~exist('trialAbortAction','var'))
                trialAbortAction = 'Repeat';
            end
            if (~exist('trialsPerSession','var'))
                trialsPerSession = 10000;
            end

            trialTableOptions.trialSequence = trialSequence;
            trialTableOptions.blockSequence = blockSequence;
            trialTableOptions.numberOfTimesRepeatBlockSequence = numberOfTimesRepeatBlockSequence;
            trialTableOptions.trialsPerSession = 10000;
            trialTableOptions.trialAbortAction = trialAbortAction;     % Repeat, Delay, Drop
            trialTableOptions.trialsPerSession = trialsPerSession;

            trialTable = ArumeCore.TrialTableBuilder.MakeTrialTable(conditionTable, blocks, trialTableOptions);
        end
    end

    methods(Static)
        function conditionTable = MakeConditionTable(conditionVariables)
            % Create the matrix with all the possible combinations of
            % condition variables. Each combination is a condition
            % total number of conditions is the product of the number of
            % values of each condition variable
            nConditions = 1;
            for iVar = 1:height(conditionVariables)
                nConditions = nConditions * length(conditionVariables{iVar, 'Values'}{1});
            end

            conditionTable = table();

            %-- recursion to create the condition matrix
            % for each variable, we repeat the previous matrix as many
            % times as values the current variable has and in each
            % repetition we add a new column with one of the values of the
            % current variable
            % we start with the last variable so the table has the big
            % blocks first and the small blocks later
            %
            % example: var1 = {a b} var2 = {e f g}
            % step 1: matrix = [ e ;
            %                    f ;
            %                    g ];
            % step 2: matrix = [ a e ;
            %                    a f ;
            %                    a g ;
            %                    b e ;
            %                    b f ;
            %                    b g ];
            nValues = ones(height(conditionVariables),1);
            for iVar = 1:height(conditionVariables)
                % index of the variables starting from the last
                iVarRev = height(conditionVariables) - iVar + 1;
                values = conditionVariables{iVarRev, 'Values'}{1};
                if (isstring(values) || iscellstr(values))
                    % make string and cell char variables to be
                    % categoricals
                    values = categorical(values);
                end
                nValues(iVar) = length(values);
                newVariable = table( values(ceil((1:prod(nValues(1:iVar)))/prod(nValues(1:iVar-1))))', 'VariableNames', conditionVariables.Name(iVarRev));
                conditionTable = horzcat( newVariable, repmat(conditionTable, nValues(iVar),1) );
            end

            % Add the condition number to the table
            conditionNumber = (1:height(conditionTable))';
            conditionTable = [table(conditionNumber) conditionTable];

            conditionTable.Properties.UserData.ConditionVariables = conditionVariables;
        end

        function trialTable = MakeTrialTable(conditionTable, blocks, trialTableOptions)

            % first create the block sequence
            blockSeqWithRepeats = [];
            for iRepeatBlockSequence = 1:trialTableOptions.numberOfTimesRepeatBlockSequence

                % generate the sequence of blocks
                nBlocks = height(blocks);
                blockSeq = [];
                switch(trialTableOptions.blockSequence)
                    case 'Sequential'
                        blockSeq = 1:nBlocks;
                    case 'Random'
                        [~, theBlocks] = sort( rand(1,nBlocks) ); % get a random shuffle of 1 ... nBlocks
                        blockSeq = mod( theBlocks-1,  nBlocks ) + 1; % limit the random sequence to 1 ... nBlocks
                end
                blockSeq = [blockSeq;ones(size(blockSeq))*iRepeatBlockSequence];
                blockSeqWithRepeats = [blockSeqWithRepeats blockSeq];
            end
            blockSeq = blockSeqWithRepeats;

            % then create the sequence of conditions going block by block
            futureConditions = [];
            for iblock=1:size(blockSeq,2)
                i = blockSeq(1,iblock);
                possibleConditions = blocks.ConditionsIncluded{i}; % the possible conditions to select from in this block
                nConditions = length(possibleConditions);
                nTrials = floor(nConditions * blocks.TrialsPerCondition(i));

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

            % finally use the sequence of conditions to build the trial
            % table with the values of all the conditions

            trialTable = table();
            trialTable.Condition = futureConditions(:,1);
            trialTable.BlockNumber = futureConditions(:,2);
            trialTable.BlockSequenceNumber = futureConditions(:,3);
            trialTable.BlockSequenceRepeat = futureConditions(:,4);
            trialTable.Session = ceil((1:height(trialTable))/min(height(trialTable), trialTableOptions.trialsPerSession))';

            trialTable = horzcat(trialTable, conditionTable(futureConditions(:,1),2:end));

            trialTable.Properties.UserData.ConditionVariables = conditionTable.Properties.UserData.ConditionVariables;
            trialTable.Properties.UserData.Blocks = blocks;
            trialTable.Properties.UserData.ConditionTable = conditionTable;
            trialTable.Properties.UserData.trialTableOptions = trialTableOptions;

            trialTable.Properties.VariableDescriptions("Condition") = "Condition number corresponding with the row in the condition table with the particular values of the condition variables.";
            trialTable.Properties.VariableDescriptions("BlockNumber") = "Block number, just counting how many blocks are in the trial table starting at 1 and going up sequentially.";
            trialTable.Properties.VariableDescriptions("BlockSequenceNumber") = "Block sequence number corresponding with the row number in the blocks table used in Trial Table Builder.";
            trialTable.Properties.VariableDescriptions("BlockSequenceRepeat") = "Number of block sequence repeat this trial belongs to.";
            trialTable.Properties.VariableDescriptions("Session") = "Session number.";
        end

        function trialTable = Demo()
            %%

            t = ArumeCore.TrialTableBuilder();

            t.AddConditionVariable( 'Distance', ["near" "far"]);
            t.AddConditionVariable( 'Direction', {'Left' 'Right' });
            t.AddConditionVariable( 'Tilt', [-30 0 30]);

            %  >> t.ConditionTable
            %
            % ans =
            %   12×4 table
            %
            %     ConditionNumber    Distance    Direction    Tilt
            %     _______________    ________    _________    ____
            %
            %            1             near        Left       -30
            %            2             near        Left         0
            %            3             near        Left        30
            %            4             near        Right      -30
            %            5             near        Right        0
            %            6             near        Right       30
            %            7             far         Left       -30
            %            8             far         Left         0
            %            9             far         Left        30
            %           10             far         Right      -30
            %           11             far         Right        0
            %           12             far         Right       30


            % Add one block that will run only trials for task one
            % condition but 3 repeatitions of all conditions
            firstConditionNear = find(t.ConditionTable.Distance=='near',1);
            lastConditionNear = find(t.ConditionTable.Distance=='near',1,'last');
            t.AddBlock(firstConditionNear:lastConditionNear, 3);

            % Add another block with the rest of the conditions
            % that will run only once per trial
            t.AddBlock((lastConditionNear+1):height(t.ConditionTable), 1);

            % Add another block with only -30 tilt conditions
            t.AddBlock(find(t.ConditionTable.Tilt==-30), 3);


            trialSequence = 'Random';
            blockSequence =  'Random';
            blockSequenceRepeatitions = 2;
            abortAction = 'Repeat';
            trialsPerSession = 10;
            trialTable = t.GenerateTrialTable(trialSequence, blockSequence, blockSequenceRepeatitions, abortAction,trialsPerSession);

            % Add a random variable to the trial table that is not
            % randomized as a condition variable
            trialTable.Delay = randn(height(trialTable),1)*10+200;

            % trialTable =
            %
            %   72×9 table
            %
            %     Condition    BlockNumber    BlockSequenceNumber    BlockSequenceRepeat    Session    Distance    Direction    Tilt    Delay
            %     _________    ___________    ___________________    ___________________    _______    ________    _________    ____    ______
            %
            %         8             1                  2                      1                1         far         Left         0     180.98
            %         7             1                  2                      1                1         far         Left       -30     223.74
            %         9             1                  2                      1                1         far         Left        30     197.67
            %        10             1                  2                      1                1         far         Right      -30     204.04
            %        11             1                  2                      1                1         far         Right        0     211.92
            %        12             1                  2                      1                1         far         Right       30     183.15
            %         1             2                  1                      1                1         near        Left       -30     204.13
            %         3             2                  1                      1                1         near        Left        30     205.02
            %         4             2                  1                      1                1         near        Right      -30     200.83
            %         2             2                  1                      1                1         near        Left         0     201.58
            %         1             2                  1                      1                2         near        Left       -30     194.72
            %         4             2                  1                      1                2         near        Right      -30     207.23
            %         1             2                  1                      1                2         near        Left       -30      191.5
            %
            %         :             :                  :                      :                :          :            :         :        :
            %
            %         1             5                  3                      2                6         near        Left       -30     183.84
            %         7             5                  3                      2                7         far         Left       -30     200.19
            %         1             5                  3                      2                7         near        Left       -30     195.74
            %        10             5                  3                      2                7         far         Right      -30     179.66
            %        10             5                  3                      2                7         far         Right      -30     186.71
            %         1             5                  3                      2                7         near        Left       -30      196.8
            %         4             5                  3                      2                7         near        Right      -30     208.25
            %        10             6                  2                      2                7         far         Right      -30     197.71
            %        12             6                  2                      2                7         far         Right       30     202.62
            %         9             6                  2                      2                7         far         Left        30     191.02
            %         7             6                  2                      2                7         far         Left       -30     178.43
            %        11             6                  2                      2                8         far         Right        0     200.94
            %         8             6                  2                      2                8         far         Left         0     190.49
            %
            % 	Display all 72 rows.

        end
    end
end