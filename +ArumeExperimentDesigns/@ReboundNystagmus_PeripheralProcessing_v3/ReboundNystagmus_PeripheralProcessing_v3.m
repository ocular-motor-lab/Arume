classdef ReboundNystagmus_PeripheralProcessing_v3 < ArumeExperimentDesigns.EyeTracking
    % ReboundNystagmus_PeripheralProcessing experiment for Arume
    %
    % 11/13/2022 - Added a break screen and rest state variable
    % 3/3/2023   - Refactored the landmark task code as a finite-state
    % machine.
    % 3/6/2023   - Add exp design implementation.
    % 5/11/2023  - 
    %  
    %
    % Coded by Terence Tyson, 11/1/2022
    properties
        fixRad   = 20;
        fixColor = [255 0 0];
    end

    % ---------------------------------------------------------------------
    % Experiment design methods
    % ---------------------------------------------------------------------
    methods ( Access = protected )
        function dlg = GetOptionsDialog( this, importing )
            dlg = GetOptionsDialog@ArumeExperimentDesigns.EyeTracking(this, importing);

            %% ADD new options
            %dlg.Max_Speed        = { 30 '* (deg/s)' [0 100] };
            %dlg.Number_of_Speeds = {3 '* (N)' [1 100] };

            %dlg.Number_of_Dots = { 2000 '* (deg/s)' [10 10000] };
            %dlg.Max_Radius     = { 40 '* (deg)' [1 100] };
            %dlg.Min_Radius     = { 1 '* (deg)' [0 100] };

            %dlg.Min_Dot_Diam        = {0.1  '* (deg)' [0.01 100] };
            %dlg.Max_Dot_Diam        = {0.4  '* (deg)' [0.01 100] };
            %dlg.Number_of_Dot_Sizes = {5 '* (N)' [1 100] };

            dlg.NumberOfRepetitions = {8 '* (N)' [1 100] };

            dlg.Do_Blank = { {'0','{1}'} };

            dlg.TargetSize = 0.5;

            dlg.BackgroundBrightness = 0;

            %% CHANGE DEFAULTS values for existing options
 
            dlg.UseEyeTracker = 0;
            dlg.Debug.DisplayVariableSelection = 'TrialNumber TrialResult controlORexperimental lineDistance  Response gazeHoldingSide distanceAdjSide'; % which variables to display every trial in the command line separated by spaces

%             dlg.DisplayOptions.ScreenWidth    = { 121 '* (cm)' [1 3000] };
%             dlg.DisplayOptions.ScreenHeight   = { 68 '* (cm)' [1 3000] };
%             dlg.DisplayOptions.ScreenDistance = { 60 '* (cm)' [1 3000] };
            dlg.DisplayOptions.SelectedScreen = 1; % This should be the large monitor in the lab

            dlg.HitKeyBeforeTrial         = 1;
            %dlg.TrialsBeforeBreak         = 20;
            dlg.TrialsBeforeBreak         = 400;
            dlg.TrialAbortAction          = 'Repeat';
            %dlg.TrialDurSequence          = [5 35 35 36 36 36.2];
            dlg.TrialDurSequence          = [0 5 5 35 35 36 36 36.55]; % 550 ms here is about ~xx ms
            %dlg.TrialDurSequence_Short    = [0 5 5 10 10 11 11 11.55];
            dlg.TrialDurSequence_Short    = [0 1.5 1.5 6.5 6.5 7.5 7.5 8.05];
            %dlg.TrialDurSequence_Control  = [0 5 5 5.2];
            dlg.TrialDurSequence_Control  = [0 5 5 5.55];
            dlg.TrialDuration             = 50;
            dlg.RestDuration              = 0;
            dlg.flickerRate               = [];
            dlg.refreshRate               = 60; % 60 Hz monitor refresh rate
            %dlg.LeftORRight_Gaze  = 'L';


            % additional parameters on the setup
            whichDisplay                    = 'Lab'; % lab or office?

            switch whichDisplay
                case 'Office'
                    dlg.viewingDistance     = 33.4;   % in cm
                    dlg.targetEccentricity  = 40;     % in deg
                    %dlg.targetHorDisplace   = 28;     % in cm
                    dlg.targetHorDisplace   = ...
                        dlg.viewingDistance*tan(dlg.targetEccentricity*(pi/180)); % in cm

                    dlg.screenWidth         = 59.5;   % in cm
                    dlg.screenHeight        = 33.5;   % in cm
                    dlg.screenResolutionHor = 3840;   % in pixels
                    dlg.screenResolutionVer = 2160;   % in pixels
                    dlg.lineFixed           = 20;     % in deg
                    dlg.textScreenFontSize  = 34;     % font size
                case 'Lab'
                    dlg.viewingDistance     = 83;     % in cm
                    %dlg.viewingDistance     = 86;     % in cm
                    dlg.targetEccentricity  = 40;     % in deg
                    dlg.targetHorDisplace   = ...
                        dlg.viewingDistance*tan(dlg.targetEccentricity*(pi/180)); % in cm
                    %dlg.targetHorDisplace   = 69.6453;% in cm
                    dlg.screenWidth         = 144.78; % in cm
                    dlg.screenHeight        = 82.55;  % in cm
                    dlg.screenResolutionHor = 3840;   % in pixels
                    dlg.screenResolutionVer = 2160;   % in pixels
                    dlg.lineFixed           = 20;     % in deg
                    dlg.textScreenFontSize  = 34;     % font size
             end


            dlg.debugFlag           = NaN;
        end

        function trialTable = SetUpTrialTable(this)

            %-- condition variables ---------------------------------------
            i = 0;

            % experimental variable 1: line distance (11 conditions)
            i = i+1;
            conditionVars(i).name   = 'lineDistance';

            % subject ID: LT1
%             conditionVars(i).values = {'18.75','19','19.25','19.5','19.75',...
%                  '20','20.25','20.5','20.75','21','21.25'}; % LT1, control1 (medium)
%             conditionVars(i).values = {'17.5','18','18.5','19','19.5',...
%                  '20','20.5','21','21.5','22','22.5'}; % LT1, control2 (easy)
%             conditionVars(i).values = {'19','19.2','19.4','19.6','19.8',...
%                  '20','20.2','20.4','20.6','20.8','21'}; % LT1, control3
                 %(difficult)
%             conditionVars(i).values = {'18.75','19','19.25','19.5','19.75',...
%                   '20','20.25','20.5','20.75','21','21.25'}; % LT1, experiment1 (medium)
%
%           second experiment with wider range (easy)
             conditionVars(i).values = {'17.5','18','18.5','19','19.5',...
                  '20','20.5','21','21.5','22','22.5'}; % LT1, experiment2 (easy)


            % subject ID: LT2
%            conditionVars(i).values = {'15','16','17','18','19',...
%                '20','21','22','23','24','25'}; % LT2, control_1
%            conditionVars(i).values = {'12.5','14','15.5','17','18.5',...
%                '20','21.5','23','24.5','26','27.5'}; % LT2, control_2
%            conditionVars(i).values = {'16.25','17','17.75','18.5','19.25',...
%                '20','20.75','21.5','22.25','23','23.75'}; % LT2, control_3
%             conditionVars(i).values = {'15','16','17','18','19',...
%                             '20','21','22','23','24','25'}; % LT2, experiment_1 (easy)

            % subject ID: LT3
%             conditionVars(i).values = {'18.75','19','19.25','19.5','19.75',...
%                  '20','20.25','20.5','20.75','21','21.25'}; % LT3, control1 (medium)
%             conditionVars(i).values = {'17.5','18','18.5','19','19.5',...
%                  '20','20.5','21','21.5','22','22.5'}; % LT3, control2 (easy)
%             conditionVars(i).values = {'19','19.2','19.4','19.6','19.8',...
%                  '20','20.2','20.4','20.6','20.8','21'}; % LT3, control3
                % (difficult)
%             conditionVars(i).values = {'18.75','19','19.25','19.5','19.75',...
%                   '20','20.25','20.5','20.75','21','21.25'}; % LT3, experiment1 (medium)

                        % subject ID: LT4
%             conditionVars(i).values = {'17.5','18','18.5','19','19.5',...
%                  '20','20.5','21','21.5','22','22.5'}; % LT4, control2 (easy)
%             conditionVars(i).values = {'18.75','19','19.25','19.5','19.75',...
%                  '20','20.25','20.5','20.75','21','21.25'}; % LT4, control1 (medium)
%             conditionVars(i).values = {'19','19.2','19.4','19.6','19.8',...
%                  '20','20.2','20.4','20.6','20.8','21'}; % LT4, control3
                % (difficult)
%                 conditionVars(i).values = {'17.5','18','18.5','19','19.5',...
%                   '20','20.5','21','21.5','22','22.5'}; % LT4, control2 (easy)    

                        % subject ID: LT8
%             conditionVars(i).values = {'19','19.2','19.4','19.6','19.8',...
%                  '20','20.2','20.4','20.6','20.8','21'}; % LT8, control1
                % (difficult)
%             conditionVars(i).values = {'17.5','18','18.5','19','19.5',...
%                  '20','20.5','21','21.5','22','22.5'}; % LT8, control2 (easy)
%             conditionVars(i).values = {'18.75','19','19.25','19.5','19.75',...
%                  '20','20.25','20.5','20.75','21','21.25'}; % LT8, control3 (medium)
%              conditionVars(i).values = {'17.5','18','18.5','19','19.5',...
%                   '20','20.5','21','21.5','22','22.5'}; % LT8, experiment1 (easy)

                        % subject ID: LT6
%             conditionVars(i).values = {'18.75','19','19.25','19.5','19.75',...
%                  '20','20.25','20.5','20.75','21','21.25'}; % LT6, control1 (medium)
%             conditionVars(i).values = {'19','19.2','19.4','19.6','19.8',...
%                  '20','20.2','20.4','20.6','20.8','21'}; % LT6, control2
                % (difficult)
%             conditionVars(i).values = {'17.5','18','18.5','19','19.5',...
%                  '20','20.5','21','21.5','22','22.5'}; % LT6, control3 (easy)
%               conditionVars(i).values = {'17.5','18','18.5','19','19.5',...
%                  '20','20.5','21','21.5','22','22.5'}; % LT6, experiment1 (easy)

                        % subject ID: LT7
%             conditionVars(i).values = {'19','19.2','19.4','19.6','19.8',...
%                  '20','20.2','20.4','20.6','20.8','21'}; % LT7, control1
                %--(difficult)                  
%             conditionVars(i).values = {'17.5','18','18.5','19','19.5',...
%                  '20','20.5','21','21.5','22','22.5'}; % LT7, control2 (easy)
%             conditionVars(i).values = {'18.75','19','19.25','19.5','19.75',...
%                  '20','20.25','20.5','20.75','21','21.25'}; % LT7, control3 (medium)
%              conditionVars(i).values = {'17.5','18','18.5','19','19.5',...
%                   '20','20.5','21','21.5','22','22.5'}; % LT7, experiment1 (easy)

                        % subject ID: LT5
%             conditionVars(i).values = {'19','19.2','19.4','19.6','19.8',...
%                  '20','20.2','20.4','20.6','20.8','21'}; % LT5, control3
%                 (difficult)            
%             conditionVars(i).values = {'18.75','19','19.25','19.5','19.75',...
%                  '20','20.25','20.5','20.75','21','21.25'}; % LT5, control1 (medium)            
%             conditionVars(i).values = {'17.5','18','18.5','19','19.5',...
%                  '20','20.5','21','21.5','22','22.5'}; % LT5, control2 (easy)
%              conditionVars(i).values = {'17.5','18','18.5','19','19.5',...
%                   '20','20.5','21','21.5','22','22.5'}; % LT5, experiment1 (easy)          
%             

                       % subject ID: LT10
%               conditionVars(i).values = {'19','19.2','19.4','19.6','19.8',...
%                   '20','20.2','20.4','20.6','20.8','21'}; % LT10, control1
                % (difficult)
%             conditionVars(i).values = {'17.5','18','18.5','19','19.5',...
%                  '20','20.5','21','21.5','22','22.5'}; % LT10, control2 (easy)
%             conditionVars(i).values = {'18.75','19','19.25','19.5','19.75',...
%                  '20','20.25','20.5','20.75','21','21.25'}; % LT10, control3 (medium)
%              conditionVars(i).values = {'17.5','18','18.5','19','19.5',...
%                   '20','20.5','21','21.5','22','22.5'}; % LT10, experiment1 (easy)

                     %subject ID: LT9
%             conditionVars(i).values = {'17.5','18','18.5','19','19.5',...
%                  '20','20.5','21','21.5','22','22.5'}; % LT9, control1 
                % (easy)
%             conditionVars(i).values = {'18.75','19','19.25','19.5','19.75',...
%                  '20','20.25','20.5','20.75','21','21.25'}; % LT9, control2 
                % (medium)
%             conditionVars(i).values = {'19','19.2','19.4','19.6','19.8',...
%                  '20','20.2','20.4','20.6','20.8','21'}; % LT9, control3
%                 % (difficult)
%               conditionVars(i).values = {'17.5','18','18.5','19','19.5',...
%                  '20','20.5','21','21.5','22','22.5'}; % LT11, control1 
% %                 % (easy)

%                     %subject ID: LT11
%             conditionVars(i).values = {'19','19.2','19.4','19.6','19.8',...
%                  '20','20.2','20.4','20.6','20.8','21'}; % LT11, control1
                % (difficult)              
%             conditionVars(i).values = {'17.5','18','18.5','19','19.5',...
%                  '20','20.5','21','21.5','22','22.5'}; % LT11, control2 
                % (easy)
%             conditionVars(i).values = {'18.75','19','19.25','19.5','19.75',...
%                  '20','20.25','20.5','20.75','21','21.25'}; % LT11, control3 
                % (medium)

                    % subject ID: LT12
%             conditionVars(i).values = {'17.5','18','18.5','19','19.5',...
%                  '20','20.5','21','21.5','22','22.5'}; % LT12, control1 (easy)
%             conditionVars(i).values = {'19','19.2','19.4','19.6','19.8',...
%                  '20','20.2','20.4','20.6','20.8','21'}; % LT12  , control2
%                 (difficult)
%             conditionVars(i).values = {'18.75','19','19.25','19.5','19.75',...
%                  '20','20.25','20.5','20.75','21','21.25'}; % LT12, control3 (medium)
%              conditionVars(i).values = {'17.5','18','18.5','19','19.5',...
%                   '20','20.5','21','21.5','22','22.5'}; % LT12, experiment1 (easy)

                    % subject ID: LT13
%              conditionVars(i).values = {'17.5','18','18.5','19','19.5',...
%                   '20','20.5','21','21.5','22','22.5'}; % LT12, control1 (easy)
%             conditionVars(i).values = {'19','19.2','19.4','19.6','19.8',...
%                  '20','20.2','20.4','20.6','20.8','21'}; % LT12  , control2
                %(difficult)
%              conditionVars(i).values = {'18.75','19','19.25','19.5','19.75',...
%                   '20','20.25','20.5','20.75','21','21.25'}; % LT12, control3 (medium)
%              conditionVars(i).values = {'18.75','19','19.25','19.5','19.75',...
%                   '20','20.25','20.5','20.75','21','21.25'}; % LT12, experiment1 (medium)


                              
            % experimental variable 2: control or experimental (2




            % conditions)
            i = i + 1;
            conditionVars(i).name   = 'controlORexperimental';
            conditionVars(i).values = {'control','experimental'};
            %conditionVars(i).values = {'experimental'};
            %conditionVars(i).values = {'control'};

            % experimental variable 3: leftward or rightward gaze holding
            % (2 conditions)
            i = i + 1;
            conditionVars(i).name   = 'gazeHoldingSide';
            conditionVars(i).values   = {'leftwardGazeHolding','rightwardGazeHolding'};

            % trial table options
            trialTableOptions                   = this.GetDefaultTrialTableOptions();
            trialTableOptions.trialSequence     = 'Random';
            trialTableOptions.trialAbortAction  = 'Delay';
            trialTableOptions.trialsPerSession  = 1000;
            trialTableOptions.numberOfTimesRepeatBlockSequence = this.ExperimentOptions.NumberOfRepetitions;
            %trialTableOptions.experimentType    = 'controlOnly'; % fullExperiment or controlOnly
             trialTableOptions.experimentType    = 'fullExperiment'; % fullExperiment or controlOnly
            trialTable = this.GetTrialTableFromConditions(conditionVars, trialTableOptions);

            % organize this table to get the sequencing correct
            % trialTable
            
            % now randomize left right condition for a given line distance
%             leftRightOrder{1} = [1 0 1 0 1 0 1 0 1 0 1];
%             leftRightOrder{2} = [1 0 1 0 1 0 1 0 1 0 1];
%             leftRightOrder{1} = leftRightOrder{1}(randperm(length(leftRightOrder{1})));
%             leftRightOrder{2} = leftRightOrder{2}(randperm(length(leftRightOrder{2})));  

            %leftRightAdjOrder = [1 0 1 0 1 0 1 0 1 0 1];
            leftRightAdjOrder{1} = [1 0 1 0 1 0 1 0 1 0 1];
            leftRightAdjOrder{2} = [1 0 1 0 1 0 0 0 1 0 1];

            % add eccentric gaze holding duration column            
            trialTable.gazeHoldingDuration = ones(height(trialTable),1);
%            uniqueDist                     = unique(trialTable.lineDistance); % unique distances          
            
            % add a column for the side that the distance was adjusted
            trialTable.distanceAdjSide     = zeros(height(trialTable),1);

            % new code here (4/17/2023)
            nReps = 8; 
            uniqueDist = unique(trialTable.lineDistance); % unique distances
            gazeHoldingCategories = {'leftwardGazeHolding','rightwardGazeHolding'};
            newTable = [];
            for i = 1:nReps
                % gaze side permutation
                gazeSideInd_Permuted = randperm(2);

                % loop through each gaze side for a given repetition
                for j = 1:length(gazeSideInd_Permuted)
                    % exp/control condition permutation
                    expCondInd_Permuted = randperm(2);

                    % loop through each exp condition for a given gaze side
                    for k = 1:length(expCondInd_Permuted)
                        % line distance permutation
                        %lineDistance_Permuted = randperm(length(uniqueDist));                      

                        % if exp condition, insert 30 sec trial dur at head
                        % and 5 sec trial dur for the proceeding 10 trials
                        if expCondInd_Permuted(k) == 1
                            % find data for given block number, gaze side,
                            % and is experimental condition
                            tempTable = trialTable(trialTable.BlockNumber == i & ...
                            trialTable.controlORexperimental == ...
                            categorical({'experimental'}) & ...
                            trialTable.gazeHoldingSide == categorical(cellstr(gazeHoldingCategories{gazeSideInd_Permuted(j)})),:);
                            tempTable.gazeHoldingDuration(1)     = 30;
                            tempTable.gazeHoldingDuration(2:end) = 5;
                        % if control condition, insert 1 sec trial dur for
                        % the 11 trials
                        elseif expCondInd_Permuted(k) == 2
                            % find data for given block number, gaze side,
                            % and is experimental condition
                            tempTable = trialTable(trialTable.BlockNumber == i & ...
                            trialTable.controlORexperimental == ...
                            categorical({'control'}) & ...
                            trialTable.gazeHoldingSide == categorical(cellstr(gazeHoldingCategories{gazeSideInd_Permuted(j)})),:);
                            tempTable.gazeHoldingDuration(:)     = 1;

                        end
                        % randomly assign the side the distance is manipulated,
                        % permutation
%                             permSide                        = randperm(nReps);
%                             distanceAdjustmentSide_Permuted = leftRightAdjOrder(permSide);
%                             %for L = 1:length(distanceAdjustmentSide_Permuted)
%                             tempTable.distanceAdjSide = distanceAdjustmentSide_Permuted;
                        %end

                        % save data to new table
                        newTable  = [newTable;tempTable];
                        tempTable = [];

                    end
                end
            end
            % add a column that randomizes which side to adjust the line
            % distances
%             i = 1;
%             j = 11;
%             k = 11; % increment size
%             iter = 1;
% 
%             stayLoop = 1;
%             permBinary = [1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0]';
%             permBinary = permBinary(randperm(32));
%             while j <= 352
%                 permAdjSide = leftRightAdjOrder{permBinary(iter)+1}(randperm(11))';
%                 while stayLoop
%                     tempTable2      = newTable.distanceAdjSide;
%                     tempTable2(i:j) = permAdjSide;
%                     if sum(tempTable2) <= 176
%                         newTable.distanceAdjSide(i:j) = permAdjSide;
%                         i    = j + 1;
%                         j    = j + k;
%                         iter = iter + 1;
%                         break
%                     elseif sum(tempTable2) > 176
%                         randInd = randperm(11,1);
%                         permAdjSide(randInd) = 0;
%                     end
%                 end
%             end

            % initialize the binary array representing the sides in which
            % the line distance is varied; this will be used to create the
            % random permutation of this array.
            lineSideVary = [0 1 0 1];

            % define variables to iterate through
            lineDist        = unique(newTable.lineDistance);
            expControlCond  = unique(newTable.controlORexperimental); 
            gazeHoldingSide = unique(newTable.gazeHoldingSide); 
            for i = 1:length(lineDist) % for a given line distance
                for j = 1:length(expControlCond) % for a give condition (control or experiment)
                    for k = 1:length(gazeHoldingSide) % for a given gaze-holding side
                        % generate random permutation
                        permInd        = [randperm(length(lineSideVary)) randperm(length(lineSideVary))];
                        lineSide4Table = [lineSideVary(permInd(1:4))';...
                            lineSideVary(permInd(5:end))'];
                        % redefine the given rows with the new
                        % distanceAdjSide
                        newTable.distanceAdjSide(newTable.lineDistance == lineDist(i) & ...
                            newTable.controlORexperimental == expControlCond(j) & ...
                            newTable.gazeHoldingSide == gazeHoldingSide(k),:) = lineSide4Table; 
                    end
                end
            end


            trialTable = newTable;

            % If control-only experiment, then get control trials of the
            % trial table.
            switch trialTableOptions.experimentType
                case 'controlOnly'
                    trialTable = trialTable(trialTable.controlORexperimental == ...
                        categorical({'control'}),:);
            end

        end
  


        function [trialResult, thisTrialData] = runTrial( this, thisTrialData )
            try
                Enum = ArumeCore.ExperimentDesign.getEnum();
                graph = this.Graph;
                trialResult = Enum.trialResult.CORRECT;

                lastFlipTime        = GetSecs;
                secondsRemaining    = this.ExperimentOptions.TrialDuration;
                thisTrialData.TimeStartLoop = lastFlipTime;

                % initialize flicker stim start time
                timeStamps     = nan(this.ExperimentOptions.refreshRate*...
                    this.ExperimentOptions.TrialDuration,1);
                timeInd        = 1;
                
                % trial level flags
                waitTime4Button = 0;
                unlimitedTime   = 0; % for 2AFC response
                restState       = 0; % rest state
%                 beepState       = 1;
%                 expState        = 1;

                % initialize experimentalState for finite state machine
                % implementation
                experimentalState = 1;
                targetFlashState  = 1;

                while secondsRemaining > 0
                    secondsNow          = GetSecs;
                    secondsElapsed      = secondsNow - thisTrialData.TimeStartLoop;
                    if ~unlimitedTime
                        secondsRemaining    = this.ExperimentOptions.TrialDuration - secondsElapsed;
                        timeStamps(timeInd) = secondsNow;
                        if timeInd > 1
                            secondsperIteration = secondsperIteration + (secondsNow - timeStamps(timeInd-1));
                            waitTime4Button     = waitTime4Button - (secondsNow - timeStamps(timeInd-1));
                        else
                            secondsperIteration = 0;
                        end
                        timeInd             = timeInd + 1;
                    end

                    % -----------------------------------------------------------------
                    % --- Drawing of stimulus -----------------------------------------
                    % -----------------------------------------------------------------

                    %-- Find the center of the screen
                    [mx, my] = RectCenter(graph.wRect);

                    fixRect = [0 0 10 10];
                    fixRect = CenterRectOnPointd( fixRect, mx, my );
                    %fixRect =  [440   190   450   300];
                    %fixRect =  [500   290   510   300]';

                    % 
                    bisectionFlag    = 0;
                    bisectionFlag2   = 0;
                    eccentricFixFlag = 0;
                    
                    switch thisTrialData.controlORexperimental
                        case 'experimental'
                            switch experimentalState
                                % State 1: Trial start, Center fixation (0
                                % deg), 5 second initial fixations for the long
                                % gaze holding trials and 1 second initial
                                % fixations for the short gaze holding
                                % trials
                                case 1
                                    if secondsRemaining <= this.ExperimentOptions.TrialDuration-this.ExperimentOptions.TrialDurSequence(1) && ...
                                    secondsRemaining > this.ExperimentOptions.TrialDuration-this.ExperimentOptions.TrialDurSequence(2) && ...
                                    thisTrialData.gazeHoldingDuration == 30
                                        Beeper(500,0.4,0.15); % start trial beep
                                        experimentalState = 2;
                                        targetFlashState  = 1; % target flashing
                                    elseif secondsRemaining <= this.ExperimentOptions.TrialDuration-this.ExperimentOptions.TrialDurSequence_Short(1) && ...
                                    secondsRemaining > this.ExperimentOptions.TrialDuration-this.ExperimentOptions.TrialDurSequence_Short(2) && ...
                                    thisTrialData.gazeHoldingDuration == 5
                                        Beeper(500,0.4,0.15); % start trial beep
                                        experimentalState = 2;
                                        targetFlashState  = 1; % target flashing
                                    end
                                % State 2: Eccentric Gaze Holding (40 deg)
                                case 2
                                    % if gaze holding is long (30 seconds)
                                    if secondsRemaining <= this.ExperimentOptions.TrialDuration-this.ExperimentOptions.TrialDurSequence(3) && ...
                                    secondsRemaining > this.ExperimentOptions.TrialDuration-this.ExperimentOptions.TrialDurSequence(4) && ...
                                    thisTrialData.gazeHoldingDuration == 30
                                        experimentalState = 3;
                                        targetFlashState  = 2; % target continuous
                                        Beeper(500,0.4,0.15);  % initial beep to indicate state change
                                    % if gaze holding is short (5 seconds)
                                    elseif secondsRemaining <= this.ExperimentOptions.TrialDuration-this.ExperimentOptions.TrialDurSequence_Short(3) && ...
                                    secondsRemaining > this.ExperimentOptions.TrialDuration-this.ExperimentOptions.TrialDurSequence_Short(4) && ...
                                    thisTrialData.gazeHoldingDuration == 5
                                        experimentalState = 3;
                                        targetFlashState  = 2; % target continuous
                                        Beeper(500,0.4,0.15);  % initial beep to indicate state change
                                    end
                                % State 3: 1 sec interval to move gaze position from
                                % eccentric target back to central target   
                                case 3                                  
                                    if secondsRemaining <= this.ExperimentOptions.TrialDuration-this.ExperimentOptions.TrialDurSequence(5) && ...
                                        secondsRemaining > this.ExperimentOptions.TrialDuration-this.ExperimentOptions.TrialDurSequence(6) && ...
                                        thisTrialData.gazeHoldingDuration == 30
                                        experimentalState = 4;
                                        targetFlashState  = 5; % Target is on for 500 ms and off for 500 ms
                                        Beeper(500,0.4,0.15);
                                    elseif secondsRemaining <= this.ExperimentOptions.TrialDuration-this.ExperimentOptions.TrialDurSequence_Short(5) && ...
                                        secondsRemaining > this.ExperimentOptions.TrialDuration-this.ExperimentOptions.TrialDurSequence_Short(6) && ...
                                        thisTrialData.gazeHoldingDuration == 5
                                        experimentalState = 4;
                                        targetFlashState  = 5; % Target is on for 500 ms and off for 500 ms
                                        Beeper(500,0.4,0.15);
                                    end     
                                % State 4: Primary Position, bisection task   
                                case 4
                                    if secondsRemaining <= this.ExperimentOptions.TrialDuration-this.ExperimentOptions.TrialDurSequence(7) && ...
                                    secondsRemaining > this.ExperimentOptions.TrialDuration-this.ExperimentOptions.TrialDurSequence(8) && ...
                                        thisTrialData.gazeHoldingDuration == 30
                                        Beeper(1000,0.4,0.15);
                                        experimentalState = 5;
                                        targetFlashState  = 3; % Bisection stimulus
                                        %GetSecs
%                                         [fixRect,bisectionFlag] = biSection_Stimulus(fixRect,...
%                                         thisTrialData.lineDistance(end),...
%                                         this.ExperimentOptions);
                                    elseif secondsRemaining <= this.ExperimentOptions.TrialDuration-this.ExperimentOptions.TrialDurSequence_Short(7) && ...
                                    secondsRemaining > this.ExperimentOptions.TrialDuration-this.ExperimentOptions.TrialDurSequence_Short(8) && ...
                                        thisTrialData.gazeHoldingDuration == 5
                                        Beeper(1000,0.4,0.15);
                                        experimentalState = 5;
                                        targetFlashState  = 3; % Bisection stimulus
                                    end
 
                                % State 5: Fixation at central position (0 deg) 
                                % unlimited time
                                case 5
                                    if secondsRemaining < this.ExperimentOptions.TrialDuration-this.ExperimentOptions.TrialDurSequence(8) && ...
                                            thisTrialData.gazeHoldingDuration == 30
                                        %experimentalState = 6;
                                        targetFlashState  = 4; % blank screen
                                        unlimitedTime     = 1; % unlimited time to respond
                                        bisectionFlag2    = 1; % collect response for task
                                    elseif secondsRemaining < this.ExperimentOptions.TrialDuration-this.ExperimentOptions.TrialDurSequence_Short(8) && ...
                                            thisTrialData.gazeHoldingDuration == 5
                                        targetFlashState  = 4; % blank screen
                                        unlimitedTime     = 1; % unlimited time to respond
                                        bisectionFlag2    = 1; % collect response for task
                                    else
                                        %GetSecs
                                    end
%                                 % State 6: End trial
%                                 % REST (eyes closed)
%                                 case 6
%                                     response = 'NA';
%                                     thisTrialData.Response = response;
%                                     thisTrialData.ResponseTime = GetSecs;
%                                     break 
                                
                            end
                        case 'control'
                            switch experimentalState
                                % State 1: 5 sec gaze in primary position
                                case 1
                                    if secondsRemaining <= this.ExperimentOptions.TrialDuration-this.ExperimentOptions.TrialDurSequence_Control(1) && ...
                                            secondsRemaining > this.ExperimentOptions.TrialDuration-this.ExperimentOptions.TrialDurSequence_Control(2)                                     
                                        Beeper(500,0.4,0.15); % start trial beep
                                        experimentalState = 2;
                                        targetFlashState  = 1; % target flashing
                                    end
                                % State 2: Primary Position, bisection task
                                case 2
                                    if secondsRemaining <= this.ExperimentOptions.TrialDuration-this.ExperimentOptions.TrialDurSequence_Control(3) && ...
                                        secondsRemaining > this.ExperimentOptions.TrialDuration-this.ExperimentOptions.TrialDurSequence_Control(4)
                                        Beeper(1000,0.4,0.15);
                                        experimentalState = 3;
                                        targetFlashState  = 3; % bisection stimulus
                                        GetSecs
                            
                                    end
                                % State 3: Primary position, question  
                                case 3
                                    if secondsRemaining < this.ExperimentOptions.TrialDuration-this.ExperimentOptions.TrialDurSequence_Control(4)
                                        %experimentalState = 4;
                                        targetFlashState  = 4; % blank screen
                                        unlimitedTime     = 1; % unlimited time to respond
                                        bisectionFlag2    = 1; % collect response for task
                                    else
                                        %GetSecs
                                        1;
                                    end
                            end

                    end

                    % TARGET FINITE-STATE MACHINE % -----------------------
                    % Only draw for 16.7 ms, then off for 983.3 ms
                    switch targetFlashState
                        case 1 % Target State 1: Flashing
                            if secondsperIteration <= 0.0167 && secondsperIteration >= 0
                                % Draw fixation position at central position
                                % (0 deg, listing secondary position, disconjugate)
                                Screen('FillOval', graph.window,  this.fixColor, fixRect);
                            elseif mod(secondsperIteration,0.0167) < (0.0167)*(0.75) && ...
                                    secondsperIteration >= 0 && ...
                                    secondsperIteration <= 0.0167*2
                                % If secondsperIteration is delayed by half
                                % a frame or less, then reset the counter
                                % by subtracting off the delay
                                secondsperIteration = secondsperIteration - ...
                                    mod(secondsperIteration,0.0167);
                                Screen('FillOval', graph.window,  this.fixColor, fixRect);
                            end
                            if secondsperIteration > 0.0167 %&& expState < 3
                                secondsperIteration = -0.98330;
                            end
                        case 2 % Target State 2: Continuous
                            % Draw fixation position in eccentric position
                            % (40 deg, listing secondary position, conjugate)
                            % calculate the fixation eccentricity in pixels
                            eccentricFix_Pixels = this.ExperimentOptions.targetEccentricity*(this.ExperimentOptions.targetHorDisplace/...
                                this.ExperimentOptions.targetEccentricity)*...
                                (this.ExperimentOptions.screenResolutionHor/this.ExperimentOptions.screenWidth);
                            % check if debugger mode off or on
                            switch this.ExperimentOptions.Debug.DebugMode
                                case 0 % debugger mode off
                                    switch thisTrialData.gazeHoldingSide
                                        case 'rightwardGazeHolding'
                                            fixRect = [fixRect(1)+eccentricFix_Pixels ...
                                                fixRect(2) ...
                                                fixRect(3)+eccentricFix_Pixels ...
                                                fixRect(4)]';
                                        case 'leftwardGazeHolding'
                                            fixRect = [fixRect(1)-eccentricFix_Pixels ...
                                                fixRect(2) ...
                                                fixRect(3)-eccentricFix_Pixels ...
                                                fixRect(4)]';
                                    end
                                otherwise
                                    % Scaled to the debugger mode window (scaling factor
                                    % in denominator (i.e., 7.322)
                                    fixRect= [fixRect(1)+round(436.5/7.322) fixRect(2) ...
                                        fixRect(3)+round(436.5/7.322) fixRect(4)]';
                            end  
                            Screen('FillOval', graph.window,  this.fixColor, fixRect);
                        case 3 % Target State 3: Bisection Stimulus
                            % Bisection stuff goes here
                            % create bisection stimulus
                            [fixRect,bisectionFlag] = biSection_Stimulus(fixRect,...
                                thisTrialData.lineDistance(end),...
                                this.ExperimentOptions,...
                                thisTrialData.gazeHoldingSide,...
                                thisTrialData.distanceAdjSide);
                            % Draw bisection stimulus
                            Screen('FillRect', graph.window,  this.fixColor, fixRect);
                        case 4 % Target State 4: Blank Screen
                            Screen('TextSize', graph.window, 36);
                            Screen('TextFont', graph.window, 'Courier');                          
                        case 5 % Target State 5: First part of half-second target state
                            Screen('FillOval', graph.window,  this.fixColor, fixRect);
                            secondsperIteration = -0.500; 
                            targetFlashState = 6;
                        case 6 % Target State 6: Second part of half-second target state
                            if secondsperIteration <= 0.5 && secondsperIteration >= 0
                                1;
                            else
                                Screen('FillOval', graph.window,  this.fixColor, fixRect);
                            end
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
                    %if bisectionFlag || bisectionFlag2
                    if bisectionFlag2
                        % time between iterations of KbCheck loop
                        [keyIsDown, secs, keyCode, deltaSecs] = KbCheck();
                        response = '';
                        if ( keyIsDown ) && waitTime4Button <= 0
                            keys = find(keyCode);
                            for i=1:length(keys)
                                KbName(keys(i));
                                switch(KbName(keys(i)))
                                    case 'RightArrow'
                                        response = 'R';
                                        % if right key then, have kbcheck
                                        % wait for a full second before
                                        % another check.
                                        waitTime4Button = 1;
                                    case 'LeftArrow'
                                        response = 'L';
                                        waitTime4Button = 1;
                                end
                            end

                        end
                        if ( ~isempty( response) )
                            thisTrialData.Response     = response;
                            thisTrialData.ResponseTime = GetSecs;

                            % Subject rests for the remainder of the trial
                            while secondsRemaining > 0
                                switch thisTrialData.controlORexperimental
                                    case 'experimental'
                                        this.ExperimentOptions.RestDuration = 0;
                                    case 'control'
                                        this.ExperimentOptions.RestDuration = 0;
                                end
                                if ~restState
                                    TimeStartLoop_Rest = GetSecs;
                                    restState          = 1;
                                end
                                secondsNow        = GetSecs;
                                secondsElapsed    = secondsNow - TimeStartLoop_Rest;
                                secondsRemaining  = this.ExperimentOptions.RestDuration - secondsElapsed;
                                this.Graph.Flip(this, thisTrialData, secondsRemaining);
                            end

                            break
                        end
%                             thisTrialData.Response = '';
%                             thisTrialData.ResponseTime = nan;

                    end
                    % -----------------------------------------------------------------
                    % --- END Collecting responses  -----------------------------------
                    % -----------------------------------------------------------------
                end
            catch ex
                rethrow(ex)
            end

            % -------------------------------------------------------------
            % --- NESTED HELPER FUNCTIONS ---------------------------------
            % -------------------------------------------------------------
            function [fixRect,bisectionFlag] = biSection_Stimulus(fixRect,...
                    lineDistance,ExperimentOptions, gazeHoldingSide, distanceAdjSide)
                % pixel 2 degree
                fixRect_Temp  = fixRect;
                fixRect_Temp2 = fixRect;

                % calculate eccentricity in pixels
                lineDistance = double(string(lineDistance))*(ExperimentOptions.targetHorDisplace/...
                    ExperimentOptions.targetEccentricity)*...
                    (ExperimentOptions.screenResolutionHor/ExperimentOptions.screenWidth);
                lineFixed      = ExperimentOptions.lineFixed*(ExperimentOptions.targetHorDisplace/...
                    ExperimentOptions.targetEccentricity)*...
                    (ExperimentOptions.screenResolutionHor/ExperimentOptions.screenWidth);

                % middle line
                fixRect_Temp  = fixRect;
                fixRect       = [fixRect_Temp(1) ;fixRect_Temp(2)-50 ;fixRect_Temp(3) ;fixRect_Temp(4)+50];

                % compute the displacements of the vertical bars of the
                % bisection stimulus probe
                switch  gazeHoldingSide
                    % left line (varied), right line (fixed); rightward
                    % gaze-holding trials
                    case 'rightwardGazeHolding'
                        switch ExperimentOptions.Debug.DebugMode
                            case 0   % debugger mode off
                                switch distanceAdjSide
                                    case 0 % leftward adjustment
                                        fixRect       = [fixRect [fixRect_Temp(1)-lineDistance;...
                                            fixRect_Temp(2)-50; ...
                                            fixRect_Temp(3)-lineDistance ;...
                                            fixRect_Temp(4)+50]];
                                    case 1
                                        fixRect       = [fixRect [fixRect_Temp(1)+lineDistance;...
                                            fixRect_Temp(2)-50; ...
                                            fixRect_Temp(3)+lineDistance ;...
                                            fixRect_Temp(4)+50]];
                                end
                                
                            otherwise % debugger mode
                                fixRect       = [fixRect [fixRect_Temp(1)-lineDistance;...
                                    fixRect_Temp(2)-50; ...
                                    fixRect_Temp(3)-lineDistance;...
                                    fixRect_Temp(4)+50]];
                        end
                        switch distanceAdjSide
                            case 0 % leftward adjustment
                                % right line (fixed)
                                fixRect       = [fixRect [fixRect_Temp(1)+lineFixed; fixRect_Temp(2)-50; ...
                                fixRect_Temp(3)+lineFixed; fixRect_Temp(4)+50]];
                            case 1 % rightward adjustment
                                % left line (fixed)
                                fixRect       = [fixRect [fixRect_Temp(1)-lineFixed; fixRect_Temp(2)-50; ...
                                fixRect_Temp(3)-lineFixed; fixRect_Temp(4)+50]];
                        end
                    % right line (varied), left line (fixed); leftward
                    % gaze-holding trials
                    case 'leftwardGazeHolding'
                        switch ExperimentOptions.Debug.DebugMode
                            case 0   % debugger mode off
                                switch distanceAdjSide
                                    case 0 % leftward adjustment
                                        fixRect       = [fixRect [fixRect_Temp(1)-lineDistance;...
                                            fixRect_Temp(2)-50; ...
                                            fixRect_Temp(3)-lineDistance ;...
                                            fixRect_Temp(4)+50]];
                                    case 1 % rightward adjustment
                                        fixRect       = [fixRect [fixRect_Temp(1)+lineDistance;...
                                            fixRect_Temp(2)-50; ...
                                            fixRect_Temp(3)+lineDistance ;...
                                            fixRect_Temp(4)+50]];
                                end
                            otherwise % debugger mode
                                fixRect       = [fixRect [fixRect_Temp(1)+lineDistance;...
                                    fixRect_Temp(2)-50; ...
                                    fixRect_Temp(3)+lineDistance;...
                                    fixRect_Temp(4)+50]];
                        end
                        switch distanceAdjSide
                            case 0 % leftward adjustment
                                % right line (fixed)
                                fixRect       = [fixRect [fixRect_Temp(1)+lineFixed; fixRect_Temp(2)-50; ...
                                fixRect_Temp(3)+lineFixed; fixRect_Temp(4)+50]];
                            case 1 % rightward adjustment
                                % left line (fixed)
                                fixRect       = [fixRect [fixRect_Temp(1)-lineFixed; fixRect_Temp(2)-50; ...
                                fixRect_Temp(3)-lineFixed; fixRect_Temp(4)+50]];
                        end
                        
                end

%                 % right line (fixed)
%                 fixRect       = [fixRect [fixRect_Temp(1)+lineFixed; fixRect_Temp(2)-50; ...
%                     fixRect_Temp(3)+lineFixed; fixRect_Temp(4)+50]];
                % change bisection flag to be true
                bisectionFlag = 1;
            end
        end
    end
end
