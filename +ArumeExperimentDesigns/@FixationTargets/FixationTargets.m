classdef FixationTargets < ArumeExperimentDesigns.EyeTracking
    %OPTOKINETICTORSION Summary of this class goes here
    %   Detailed explanation goes here

    properties
        fixRad = 20;
        fixColor = [255 0 0];
        targetPositions =[];
    end

    % ---------------------------------------------------------------------
    % Experiment design methods
    % ---------------------------------------------------------------------
    methods ( Access = protected )
        function dlg = GetOptionsDialog( this, importing )
            if( ~exist( 'importing', 'var' ) )
                importing = 0;
            end
            dlg = GetOptionsDialog@ArumeExperimentDesigns.EyeTracking(this, importing);
            dlg.Debug.DisplayVariableSelection = 'TrialNumber TrialResult TargetPosition'; % which variables to display every trial in the command line separated by spaces

            dlg.DisplayOptions.ScreenWidth = { 55 '* (cm)' [1 3000] };
            dlg.DisplayOptions.ScreenHeight = { 31 '* (cm)' [1 3000] };
            dlg.DisplayOptions.ScreenDistance = { 67 '* (cm)' [1 3000] };

            dlg.TrialDuration =  { 10 '* (s)' [1 100] };
            dlg.NumberRepetitions = 10;

            dlg.TargetSize = 1;
            dlg.Calibration_Type = { {'Center dot' '5 dots' '{9 dots}' '13 dots' '17 dots'} };
            dlg.Calibration_Distance_H = { 10 '* (deg)' [1 3000] };
            dlg.Calibration_Distance_V = { 10 '* (deg)' [1 3000] };

            dlg.BackgroundBrightness = 0;
        end


        function trialTable = SetUpTrialTable(this)
            Screen('Preference', 'SkipSyncTests', 0);                  % DO NOT skip sync tests (set to 1 or 2 for debugging only)
            Screen('Preference', 'VisualDebugLevel', 1);               % Minimize startup splash
            Screen('Preference', 'SuppressAllWarnings', 0);            % Show all warnings
            Screen('Preference', 'ConserveVRAM', 4096);

            h = this.ExperimentOptions.Calibration_Distance_H;
            v = this.ExperimentOptions.Calibration_Distance_V;
            temp = 1.5;

            switch(this.ExperimentOptions.Calibration_Type)
                case 'Center dot'
                    this.targetPositions = {[0,0]};
                case '5 dots'
                    this.targetPositions = {[0,0],[h,v],[h,-v],[-h,v],[-h,-v]};
                case '9 dots'
                    this.targetPositions = {[0,0],[h,0],[-h,0],[0,v],[0,-v],[h,v],[h,-v],[-h,v],[-h,-v]};
                case '13 dots'
                    this.targetPositions = {[0,0],[h,0],[-h,0],[0,v],[0,-v],[h,v],[h,-v],[-h,v],[-h,-v],...
                        [h/temp,v/temp],[h/temp,-v/temp],[-h/temp,v/temp],[-h/temp,-v/temp]};
                case '17 dots'
                    this.targetPositions = {[0,0],[h,0],[-h,0],[0,v],[0,-v],[h,v],[h,-v],[-h,v],[-h,-v],...
                        [h/temp,0],[-h/temp,0],[0,v/temp],[0,-v/temp],[h/temp,v/temp],[h/temp,-v/temp],[-h/temp,v/temp],[-h/temp,-v/temp]};
            end
            t = ArumeCore.TrialTableBuilder();

            t.AddConditionVariable("TargetPosition", { ...
                [0 0], [0 10], [10 0], [10 10], [-10 10], [-10 -10], [10 -10] [-10 0], [0 -10] ...
                [0 2], [0 4], [0 6], [0 8], [0 -2], [0 -4], [0 -6], [0 -8], ...
                [2 0], [4 0], [6 0], [8 0], [-2 0], [-4 0], [-6 0], [-8 0], ...
                [2 2], [4 4], [6 6], [8 8], [-2 2], [-4 4], [-6 6], [-8 8], ...
                [2 -2], [4 -4], [6 -6], [8 -8], [-2 -2], [-4 -4], [-6 -6], [-8 -8] ...
                });

            % Add all conditions to a single block with N repetitions
            nReps = this.ExperimentOptions.NumberRepetitions;
            t.AddBlock(1:height(t.ConditionTable), nReps);

            trialSequence = 'Random';
            blockSequence = 'Sequential';
            blockSequenceRepeatitions = 1;
            abortAction = 'Repeat';
            trialsPerSession = 1000;  % You can use a large number if not splitting

            trialTable = t.GenerateTrialTable(trialSequence, blockSequence, blockSequenceRepeatitions, abortAction, trialsPerSession);
        end


        function [trialResult, thisTrialData] = runTrial(this, thisTrialData)
            try
                Enum = ArumeCore.ExperimentDesign.getEnum();
                graph = this.Graph;
                trialResult = Enum.trialResult.CORRECT;

                lastFlipTime        = GetSecs;
                secondsRemaining    = this.ExperimentOptions.TrialDuration;
                thisTrialData.TimeStartLoop = lastFlipTime;

                if (~isempty(this.eyeTracker))
                    thisTrialData.EyeTrackerFrameStartLoop = this.eyeTracker.RecordEvent( ...
                        sprintf('TRIAL_START_LOOP %d %d [%d %d]', ...
                        thisTrialData.TrialNumber, ...
                        thisTrialData.Condition, ...
                        thisTrialData.TargetPosition(1), ...
                        thisTrialData.TargetPosition(2)) );
                end

                while secondsRemaining > -5
                    secondsElapsed      = GetSecs - thisTrialData.TimeStartLoop;
                    secondsRemaining    = this.ExperimentOptions.TrialDuration - secondsElapsed;

                    if secondsRemaining > 0
                        %-- Draw fixation spot as a white cross (+)
                        Screen('FillRect', graph.window, this.ExperimentOptions.BackgroundBrightness);

                        [mx, my] = RectCenter(this.Graph.wRect);

                        % Get the stimulus position in degrees from the trial table
                        targetDeg = thisTrialData.TargetPosition{1};
                        dx = targetDeg(1);
                        dy = targetDeg(2);

                        % === LOG POSITION INTO TRIAL DATA ===
                        thisTrialData.StimulusPosition_X = targetDeg(1);
                        thisTrialData.StimulusPosition_Y = targetDeg(2);

                        % Convert visual degrees to pixels
                        pixelsPerDegree = this.Graph.pxWidth / this.ExperimentOptions.DisplayOptions.ScreenWidth * ...
                            this.ExperimentOptions.DisplayOptions.ScreenDistance;

                        targetHPix = pixelsPerDegree * tand(dx);
                        targetYPix = pixelsPerDegree * tand(dy);

                        % Fixation location
                        fixX = mx + targetHPix / 2;
                        fixY = my + targetYPix / 2;

                        % Fixation cross parameters
                        targetSizeDeg = this.ExperimentOptions.TargetSize;
                        crossLength = pixelsPerDegree * tand(targetSizeDeg); % in pixels
                        crossThickness = 2;
                        crossColor = [255, 255, 255];

                        % Define cross lines centered on fixX, fixY
                        crossCoords = [ ...
                            -crossLength/2, 0; ...
                            crossLength/2, 0; ...
                            0, -crossLength/2; ...
                            0,  crossLength/2 ...
                            ]';

                        Screen('DrawLines', graph.window, crossCoords, crossThickness, crossColor, [fixX, fixY], 2);

                        Screen('DrawingFinished', graph.window);
                    else
                        Screen('FillRect', graph.window, this.ExperimentOptions.BackgroundBrightness);
                        Screen('DrawingFinished', graph.window);
                    end

                    % Flip screen buffer
                    this.Graph.Flip();
                end

            catch ex
                rethrow(ex)
            end
        end

    end

    methods ( Access = public )
    end

    % ---------------------------------------------------------------------
    % Plot methods
    % ---------------------------------------------------------------------
    methods ( Access = public )
        function [out, options] = PlotAggregate_MicroCorrelations(this, sessions, options)


            out = [];
            if ( nargin == 1 )
                options = this.PlotAggregate_MicroCorrelations('get_defaults');
            end

            if ( ischar(sessions) )
                command = sessions;
                switch( command)
                    case 'get_options'
                        options =[];
                        return
                    case 'get_defaults'
                        optionsDlg = VOGAnalysis.PlotAggregate_MicroCorrelations('get_options');
                        options = StructDlg(optionsDlg,'',[],[],'off');
                        return
                end
            end


            data = [];
            % for each subj, ignoring 000
            for subj = 1:length(sessions)
                %load data
                AnalysisResults_QuickPhases = sessions(subj).analysisResults.QuickPhases;

                % pick the center positions
                center_indx = AnalysisResults_QuickPhases.Left_X_MeanPosition <= 2 &...
                    AnalysisResults_QuickPhases.Left_X_MeanPosition >= -2 &...
                    AnalysisResults_QuickPhases.Left_Y_MeanPosition <= 2 &...
                    AnalysisResults_QuickPhases.Left_Y_MeanPosition >= -2;

                d = AnalysisResults_QuickPhases(center_indx,:);

                % create the table
                data(subj).X_Vergence = d.Left_X_Displacement - d.Right_X_Displacement;
                data(subj).Y_Vergence = d.Left_Y_Displacement - d.Right_Y_Displacement;
                data(subj).T_Vergence = d.Left_T_Displacement - d.Right_T_Displacement;

                data(subj).X_Version = (d.Left_X_Displacement + d.Right_X_Displacement)./2;
                data(subj).Y_Version = (d.Left_Y_Displacement + d.Right_Y_Displacement)./2;
                data(subj).T_Version = (d.Left_T_Displacement + d.Right_T_Displacement)./2;

                data(subj).Left_X = d.Left_X_Displacement;
                data(subj).Left_Y = d.Left_Y_Displacement;
                data(subj).Left_T = d.Left_T_Displacement;

                data(subj).Right_X = d.Right_X_Displacement;
                data(subj).Right_Y = d.Right_Y_Displacement;
                data(subj).Right_T = d.Right_T_Displacement;


            end


            %% create the corr table

            % combine all data
            dataAll.X_Vergence = [];
            dataAll.Y_Vergence = [];
            dataAll.T_Vergence = [];

            dataAll.X_Version = [];
            dataAll.Y_Version = [];
            dataAll.T_Version = [];

            dataAll.Left_X = [];
            dataAll.Left_Y = [];
            dataAll.Left_T = [];

            dataAll.Right_X = [];
            dataAll.Right_Y = [];
            dataAll.Right_T = [];

            for subj = 1:length(data)
                dataAll.X_Vergence = [dataAll.X_Vergence; data(subj).X_Vergence];
                dataAll.Y_Vergence = [dataAll.Y_Vergence; data(subj).Y_Vergence];
                dataAll.T_Vergence = [dataAll.T_Vergence; data(subj).T_Vergence];

                dataAll.Y_Version = [dataAll.Y_Version; data(subj).Y_Version];
                dataAll.X_Version = [dataAll.X_Version; data(subj).X_Version];
                dataAll.T_Version = [dataAll.T_Version; data(subj).T_Version];

                dataAll.Left_X = [dataAll.Left_X; data(subj).Left_X];
                dataAll.Left_Y = [dataAll.Left_Y; data(subj).Left_Y];
                dataAll.Left_T = [dataAll.Left_T; data(subj).Left_T];

                dataAll.Right_X = [dataAll.Right_X; data(subj).Right_X];
                dataAll.Right_Y = [dataAll.Right_Y; data(subj).Right_Y];
                dataAll.Right_T = [dataAll.Right_T; data(subj).Right_T];
            end

            dataAll_VV = [...
                dataAll.X_Version,...
                dataAll.Y_Version,...
                dataAll.X_Vergence,...
                dataAll.Y_Vergence,...
                dataAll.T_Vergence,...
                dataAll.T_Version...
                ];

            dataAll_LR = [...
                dataAll.Left_X,...
                dataAll.Left_Y,...
                dataAll.Left_T,...
                dataAll.Right_X,...
                dataAll.Right_Y,...
                dataAll.Right_T...
                ];

            [R,P] = corr(dataAll_VV,'rows','complete');

            x = [-1:0.5:1];
            c = 0;
            figure,
            for i = 1:6
                for j = 1:6
                    c = c + 1;
                    h(i,j) = subplot(6,6,c);
                    if i <= j, continue; end

                    %indx = (dataAll_(:,j) < 1 & dataAll_(:,j) > -1) & (dataAll_(:,i) < 1 & dataAll_(:,i) > -1);

                    fitlm_{i,j} = fitlm(dataAll_VV(:,j),dataAll_VV(:,i),'RobustOpts','on');
                    %fitlm_{i,j}.Coefficients = robustfit(dataAll_VV(:,j),dataAll_VV(:,i));

                    plot(x,x*fitlm_{i,j}.Coefficients{2,1} + fitlm_{i,j}.Coefficients{1,1},'r--','LineWidth',2),
                    %plot(fitlm_{i,j}),legend off
                    hold on
                    plot(dataAll_VV(:,j),dataAll_VV(:,i),'b.');
                    xlim([-1,1])
                    ylim([-1,1])
                    title(strcat("R = ",num2str(round(R(c),2)),", P = ",num2str(P(c)) ,", Slope = ",num2str( fitlm_{i,j}.Coefficients{2,1} ) ))


                end
            end

            linkaxes(h(:))

            n = [...
                "X Version",...
                "Y Version",...
                "X Vergence",...
                "Y Vergence",...
                "T Vergence",...
                "T Version"...
                ];

            c=1;
            for i = 1:6:36
                subplot(6,6,i)
                ylabel(n(c))
                c=c+1;
            end


            for i = 31:36
                subplot(6,6,i)
                xlabel(n(i-30))
            end

            sgtitle("Center Fixational Micro Saccades Positional Displacement (degree)")


            figure,
            subplot(1,2,1)
            heatmap(R)
            title("R")

            subplot(1,2,2)
            heatmap(P)
            title("P")
            sgtitle("Center Fixational Micro Saccades Correlation parameters for Vergence and Version")

            figure,
            bar(mean(dataAll_VV,'omitnan'))
            hold on
            errorbar(mean(dataAll_VV,'omitnan'),std(dataAll_VV,'omitnan'),'o')
            xticklabels(n)
            ylabel('Amplitude (degree)')


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%% LEFT RIGHT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % [R,P] = corr(dataAll_LR,'rows','complete');
            %
            %             x = [-1:0.5:1];
            %             c = 0;
            %             figure,
            %             for i = 1:6
            %                 for j = 1:6
            %                     c = c + 1;
            %                     h2(i,j) = subplot(6,6,c);
            %                     if i >= j, continue; end
            %
            %                     %indx = (dataAll_(:,j) < 1 & dataAll_(:,j) > -1) & (dataAll_(:,i) < 1 & dataAll_(:,i) > -1);
            %
            %                     fitlm_{i,j} = fitlm(dataAll_LR(:,j),dataAll_LR(:,i));
            %
            %                     plot(x,x*fitlm_{i,j}.Coefficients{2,1} + fitlm_{i,j}.Coefficients{1,1},'r--','LineWidth',2),
            %                     %plot(fitlm_{i,j}),legend off
            %                     hold on
            %                     plot(dataAll_LR(:,j),dataAll_LR(:,i),'b.');
            %                     xlim([-1,1])
            %                     ylim([-1,1])
            %                     title(strcat("R = ",num2str(round(R(c),2)), ", P = ",num2str(P(c)) ,", Slope = ",num2str( fitlm_{i,j}.Coefficients{2,1} ) ))
            %
            %
            %                 end
            %             end
            %
            %             linkaxes(h2(:))
            %
            %             n = [...
            %                 "Left X",...
            %                 "Left Y",...
            %                 "Left T",...
            %                 "Right X",...
            %                 "Right Y",...
            %                 "Right T"...
            %                 ];
            %
            %             c=1;
            %             for i = 1:6:36
            %                 subplot(6,6,i)
            %                 ylabel(n(c))
            %                 c=c+1;
            %             end
            %
            %
            %             for i = 31:36
            %                 subplot(6,6,i)
            %                 xlabel(n(i-30))
            %             end
            %
            %             sgtitle("Center Fixational Micro Saccades Positional Displacement (degree)")
            %
            %
            %             figure,
            %             subplot(1,2,1)
            %             heatmap(R)
            %             title("R")
            %
            %             subplot(1,2,2)
            %             heatmap(P)
            %             title("P")
            %             sgtitle("Center Fixational Micro Saccades Correlation parameters for Left and Right")
            %




        end

        function Plot_Listings(this)
            figure
            s = this.Session.samplesDataTable;
            plot(s.LeftX, s.LeftT,'o')
            hold
            plot(s.RightX, s.RightT,'o')
            set(gca,'xlim',[-20 20],'ylim',[-20 20])
            xlabel('Horizontal')
            ylabel('Torsion')
        end

         function [out, options] = Plot_XY(this)

            s = this.Session.samplesDataTable;

            figure
            plot(s.LeftX, s.LeftY,'o','MarkerSize',3);
            hold
            plot(s.RightX, s.RightY,'o','MarkerSize',3);

            set(gca,'xlim',[-40 40], 'ylim', [-30 30])
         end
         
         function [out, options] = Plot_PositionVelocity(this)

             p = this.Session.analysisResults.SlowPhases.X_MeanPosition;
             pp = round(p/2.5)*2.5;
             [means,pred,grp] = grpstats(this.Session.analysisResults.SlowPhases.X_MeanVelocity,pp, ["mean","predci","gname"],"Alpha",0.1);
             figure
             subplot(2,2,1)
             plot(str2double(grp), means,'o')
             set(gca,'xlim',[-15 15], 'ylim', [-3 3])
             xlabel('Horizontal position (deg)')
             ylabel('Horizontal drift velocity (deg/s)')

             p = this.Session.analysisResults.SlowPhases.Y_MeanPosition;
             pp = round(p/2.5)*2.5;
             [means,pred,grp] = grpstats(this.Session.analysisResults.SlowPhases.Y_MeanVelocity,pp, ["mean","predci","gname"],"Alpha",0.1);
             
             subplot(2,2,2)
             plot(str2double(grp), means,'o')
             xlabel('Vertical position (deg)')
             ylabel('Vertical drift velocity (deg/s)')
             set(gca,'xlim',[-15 15], 'ylim', [-3 3])



             p = this.Session.analysisResults.SlowPhases.Left_X_MeanPosition - this.Session.analysisResults.SlowPhases.Right_X_MeanPosition;
             v = this.Session.analysisResults.SlowPhases.Left_X_MeanVelocity - this.Session.analysisResults.SlowPhases.Right_X_MeanVelocity;
             pp = round(p*2.5)/2.5;
             [means,pred,grp] = grpstats(v,pp, ["mean","predci","gname"],"Alpha",0.1);
             
             subplot(2,2,3)
             plot(str2double(grp), means,'o')
             set(gca,'xlim',[-5 5], 'ylim', [-3 3])
             xlabel('Vergence position (deg)')
             ylabel('Vergence drift velocity (deg/s)')

             % p = this.Session.analysisResults.SlowPhases.Y_MeanPosition;
             % pp = round(p/2.5)*2.5;
             % [means,pred,grp] = grpstats(this.Session.analysisResults.SlowPhases.Y_MeanVelocity,pp, ["mean","predci","gname"],"Alpha",0.1);
             % 
             % subplot(2,2,4)
             % plot(str2double(grp), means,'o')

             % set(gca,'xlim',[-15 15], 'ylim', [-3 3])
         end
         
    end
end




