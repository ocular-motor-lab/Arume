classdef ArumeGui < matlab.apps.AppBase
    %ARUMEGUI Summary of app class goes here
    %   Detailed explanation goes here

    properties (Access = public)

        % main controller
        arumeController

        % figure handle
        figureHandle        matlab.ui.Figure

        % control handles
        sessionTree
        infoBox
        sessionTable
        trialTable
        commentsTextBox

        % panel handles
        leftPanel
        tabSessions
        rightPanel
        tabSessionInfo
        tabSessionTable
        tabTrialTable

        % Menu items
        menuProject
        menuProjectNewProject
        menuProjectLoadProject
        menuProjectSaveProjectBackup
        menuProjectLoadProjectBackup
        menuProjectLoadRecentProject
        menuProjectCloseProject

        menuSession
        menuSessionNewSession
        menuSessionImportSession
        menuSessionEditSettings
        menuSessionRename
        menuSessionDelete
        menuSessionCopy
        menuSessionSendDataToWorkspace
        menuSessionNewSessionAggregate

        menuRun
        menuRunStartSession
        menuRunResumeSession
        menuRunRestartSession
        menuResumeSessionFrom

        menuAnalyze
        menuAnalyzeRunAnalyses
        menuAnalyzeMarkData 
        menuAnalyzeClearAll


        menuPlot

        menuTools
        menuToolsBiteBarGui
        menuToolsOpenProjectFolderInExplorer
    end

    %% Constructor
    methods
        function app = ArumeGui( arumeController )

            app.arumeController = arumeController;

            app.InitUI();
        end

        function InitUI(app )

            screenSize = get(groot,'ScreenSize');
            screenWidth = screenSize(3);
            screenHeight = screenSize(4);
            w = screenWidth*0.5;
            h = screenHeight*0.5;
            left = screenWidth/2-w/2;
            bottom = screenHeight/2-h/2;


            %  Construct the figure
            app.figureHandle = uifigure();
            app.figureHandle.Name              = 'Arume';
            app.figureHandle.Position          =  [left bottom w h];
            app.figureHandle.Tag               = 'Arume';
            app.figureHandle.CloseRequestFcn   = @app.figureCloseRequest;
            app.figureHandle.AutoResizeChildren = 'off';
            app.figureHandle.SizeChangedFcn    = @app.figureResizeFcn;
            app.figureHandle.UserData          = app;

            app.InitMenu();

            %  Construct panels

            app.leftPanel = uitabgroup(app.figureHandle);

            app.rightPanel = uitabgroup(app.figureHandle);
            app.figureResizeFcn();

            app.tabSessions = uitab(app.leftPanel);
            app.tabSessions.Title = 'Sessions';

            app.tabSessionInfo = uitab(app.rightPanel);
            app.tabSessionInfo.Title = 'Session info';

            %             app.tabSessionTable = uitab(app.rightPanel);
            %             app.tabSessionTable.Title = 'Session table';
            %             app.rightPanel.SelectionChangedFcn = @app.tabRightPanelCallBack;
            %
            %             app.tabTrialTable = uitab(app.rightPanel);
            %             app.tabTrialTable.Title = 'Trial table';


            %  Construct the components
            app.sessionTree = uitree( app.tabSessions);
            app.sessionTree.Position =  [1 1 app.tabSessions.Position(3)-3 app.tabSessions.Position(4)-35];
            app.sessionTree.FontName = 'consolas';
            app.sessionTree.Multiselect = 'on';
            app.sessionTree.SelectionChangedFcn = @app.sessionTreeCallBack;

            app.infoBox = uitextarea(app.tabSessionInfo);
            app.infoBox.FontName = 'consolas';
            app.infoBox.HorizontalAlignment = 'Left';
            app.infoBox.Editable = 'off';
            app.infoBox.Value = '';
            app.infoBox.Position = [1 app.tabSessionInfo.Position(4)/5+2 app.tabSessionInfo.Position(3)-3 app.tabSessionInfo.Position(4)*4/5-35];
            app.infoBox.BackgroundColor = 'w';

            %             app.sessionTable = uitable(app.tabSessionTable);
            %             app.sessionTable.Position = [1 1 app.tabSessionInfo.Position(3)-3 app.tabSessionInfo.Position(4)-35];
            %
            %             app.trialTable = uitable(app.tabTrialTable);
            %             app.trialTable.Position = [1 1 app.tabSessionInfo.Position(3)-3 app.tabSessionInfo.Position(4)-35];


            app.commentsTextBox = 	uitextarea(app.tabSessionInfo);
            app.commentsTextBox.FontName = 'consolas';
            app.commentsTextBox.HorizontalAlignment = 'Left';
            app.commentsTextBox.Value = 'Session notes:';
            app.commentsTextBox.Editable = 'on';
            app.commentsTextBox.BackgroundColor = [1 1 0.8];
            app.commentsTextBox.ValueChangedFcn = @app.commentsTextBoxCallBack;
            app.commentsTextBox.Position = [1 1 app.tabSessionInfo.Position(3)-3 app.tabSessionInfo.Position(4)*1/5];

            % This is to avoid a close all closing the GUI
            set(app.figureHandle, 'handlevisibility', 'off');

            app.UpdateGui();


            % Register the app with App Designer
            registerApp(app, app.figureHandle)

            if nargout == 0
                clear app
            end
        end

        function InitMenu(app)

            SIMPLEMODE = app.arumeController.SIMPLEMODE;

            % menu
            set(app.figureHandle,'MenuBar','none');

            app.menuProject = uimenu(app.figureHandle);
            app.menuProject.Text = 'Project';
            app.menuProject.Callback = @app.menuProjectCallback;

            app.menuProjectNewProject = uimenu(app.menuProject);
            app.menuProjectNewProject.Text = 'New project ...';
            app.menuProjectNewProject.Callback = @app.newProject;

            app.menuProjectLoadProject = uimenu(app.menuProject);
            app.menuProjectLoadProject.Text = 'Load project ...';
            app.menuProjectLoadProject.Callback = @app.loadProject;

            app.menuProjectLoadRecentProject = uimenu(app.menuProject);
            app.menuProjectLoadRecentProject.Text = 'Load recent project';

            app.menuProjectCloseProject = uimenu(app.menuProject);
            app.menuProjectCloseProject.Text = 'Close project';
            app.menuProjectCloseProject.Callback =  @app.closeProject;

            app.menuProjectSaveProjectBackup = uimenu(app.menuProject);
            app.menuProjectSaveProjectBackup.Text = 'Backup project ...';
            app.menuProjectSaveProjectBackup.Separator = 'on';
            app.menuProjectSaveProjectBackup.Callback = @app.saveProjectBackup;

            app.menuProjectLoadProjectBackup = uimenu(app.menuProject);
            app.menuProjectLoadProjectBackup.Text = 'Restore project backup ...';
            app.menuProjectLoadProjectBackup.Callback = @app.loadProjectBackup;


            app.menuSession = uimenu(app.figureHandle);
            app.menuSession.Text = 'Session';

            app.menuSessionNewSession = uimenu(app.menuSession);
            app.menuSessionNewSession.Text = 'New session';
            app.menuSessionNewSession.Callback = @app.newSession;

            app.menuSessionImportSession = uimenu(app.menuSession);
            app.menuSessionImportSession.Text = 'Import session';
            app.menuSessionImportSession.Callback =  @app.importSession;

            app.menuSessionCopy = uimenu(app.menuSession);
            app.menuSessionCopy.Label = 'Copy sessions (no data)...';
            app.menuSessionCopy.Callback = @app.CopySessions;
            app.menuSessionCopy.Separator = 'on';

            app.menuSessionDelete = uimenu(app.menuSession);
            app.menuSessionDelete.Label = 'Delete sessions ...';
            app.menuSessionDelete.Callback = @app.DeleteSessions;

            app.menuSessionRename = uimenu(app.menuSession);
            app.menuSessionRename.Label = 'Rename sessions ...';
            app.menuSessionRename.Callback = @app.RenameSessions;

            app.menuSessionEditSettings = uimenu(app.menuSession);
            app.menuSessionEditSettings.Label = 'Edit settings ...';
            app.menuSessionEditSettings.Callback = @app.EditSessionSettings;

            app.menuSessionSendDataToWorkspace = uimenu(app.menuSession);
            app.menuSessionSendDataToWorkspace.Label = 'Send data to workspace ...';
            app.menuSessionSendDataToWorkspace.Callback = @app.SendDataToWorkspace;

            app.menuSessionNewSessionAggregate = uimenu(app.menuSession);
            app.menuSessionNewSessionAggregate.Label = 'New session aggreage ...';
            app.menuSessionNewSessionAggregate.Callback = @app.NewSessionAggregate;
            app.menuSessionNewSessionAggregate.Separator = 'on';


            app.menuRun = uimenu(app.figureHandle);
            app.menuRun.Text = 'Run';

            app.menuRunStartSession = uimenu(app.menuRun);
            app.menuRunStartSession.Text = 'Start session...';
            app.menuRunStartSession.Callback = @app.startSession;

            app.menuRunResumeSession = uimenu(app.menuRun);
            app.menuRunResumeSession.Text = 'Resume session';
            app.menuRunResumeSession.Callback = @app.resumeSession;

            app.menuRunRestartSession = uimenu(app.menuRun);
            app.menuRunRestartSession.Text = 'Restart session';
            app.menuRunRestartSession.Callback = @app.restartSession;

            app.menuResumeSessionFrom = uimenu(app.menuRun);
            app.menuResumeSessionFrom.Text = 'Resume session from ...';
            app.menuResumeSessionFrom.Separator = 'on' ;


            app.menuAnalyze = uimenu(app.figureHandle);
            app.menuAnalyze.Text = 'Analyze';

            app.menuAnalyzeRunAnalyses = uimenu(app.menuAnalyze);
            app.menuAnalyzeRunAnalyses.Text = 'Run data analyses ...';
            app.menuAnalyzeRunAnalyses.Callback = @app.RunDataAnalyses;

            app.menuAnalyzeMarkData = uimenu(app.menuAnalyze);
            app.menuAnalyzeMarkData.Text = 'Mark data ...';
            app.menuAnalyzeMarkData.Separator = 'on';
            app.menuAnalyzeMarkData.Callback = @app.MarkData;


            app.menuAnalyzeClearAll = uimenu(app.menuAnalyze);
            app.menuAnalyzeClearAll.Text = 'Clear all data';
            app.menuAnalyzeClearAll.Separator = 'on';
            app.menuAnalyzeClearAll.Callback = @app.ClearAllPrepareAndAnalyses;



            app.menuPlot = uimenu(app.figureHandle);
            app.menuPlot.Text = 'Plot';

            app.menuTools = uimenu(app.figureHandle);
            app.menuTools.Text = 'Tools';

            app.menuToolsBiteBarGui = uimenu(app.menuTools);
            app.menuToolsBiteBarGui.Text = 'Bite bar GUI';
            app.menuToolsBiteBarGui.Callback = @BitebarGUI;


            app.menuToolsOpenProjectFolderInExplorer = uimenu(app.menuTools);
            app.menuToolsOpenProjectFolderInExplorer.Text = 'Open project folder in explorer...';
            app.menuToolsOpenProjectFolderInExplorer.Callback = @app.OpenProjectFolderInExplorer;


            if ( SIMPLEMODE )
                app.menuProjectSaveProjectBackup.Visible = 'off';
                app.menuProjectLoadProjectBackup.Visible = 'off';

                app.menuSessionImportSession.Visible = 'off';
                app.menuSessionSendDataToWorkspace.Visible = 'off';
                app.menuSessionNewSessionAggregate.Visible = 'off';

                app.menuResumeSessionFrom.Visible = 'off';

                app.menuAnalyzeMarkData.Visible = 'off';
                app.menuAnalyzeClearAll.Visible = 'off';
                app.menuToolsBiteBarGui.Visible = 'off';
            end

        end

        function UpdateGui( app, fastOption )
            if ( ~exist('fastOption','var') )
                fastOption = 0;
            end

            if ( isempty( app.arumeController ))
                return;
            end

            % update top box info
            if ( ~isempty( app.arumeController.currentProject ) )
                app.figureHandle.Name = sprintf('Arume - Project: %s', app.arumeController.currentProject.path);
            else
                app.figureHandle.Name = sprintf('Arume');
            end

            if ( ~fastOption )
                app.updateSessionTree();
            end

            % update info box
            app.updateInfoBox();

            % update comments text box
            if ( ~isempty( app.arumeController.currentSession ) )
                app.commentsTextBox.Enable = 'on';
                app.commentsTextBox.Value = cellstr(app.arumeController.currentSession.comment);
            else
                app.commentsTextBox.Enable = 'off';
                app.commentsTextBox.Value = '';
            end

            % update menu

            % top level menus
            if ( ~isempty( app.arumeController.currentSession ) )
                app.menuAnalyze.Enable = 'on';
                app.menuPlot.Enable = 'on';
                if ( length( app.arumeController.selectedSessions )==1 )
                    app.menuRun.Enable = 'on';
                    app.menuAnalyzeMarkData.Enable = 'on';
                else
                    app.menuRun.Enable = 'off';
                    app.menuAnalyzeMarkData.Enable = 'off';
                end
            else
                app.menuRun.Enable = 'off';
                app.menuAnalyze.Enable = 'off';
                app.menuPlot.Enable = 'off';
            end


            % sub menus

            if ( ~isempty( app.arumeController.currentProject ) )
                set(app.menuProjectCloseProject, 'Enable', 'on');
                set(app.menuProjectSaveProjectBackup, 'Enable', 'on');
                set(app.menuSession, 'Enable', 'on');
                set(app.menuSessionNewSession, 'Enable', 'on');

            else
                set(app.menuProjectCloseProject, 'Enable', 'off');
                set(app.menuProjectSaveProjectBackup, 'Enable', 'off');
                set(app.menuSession, 'Enable', 'off');
                set(app.menuSessionNewSession, 'Enable', 'off');
            end

            if ( ~isempty( app.arumeController.currentSession ) && numel(app.arumeController.selectedSessions) == 1 )

                set(app.menuRunStartSession, 'Enable', 'off');
                set(app.menuRunRestartSession, 'Enable', 'off');
                set(app.menuRunResumeSession, 'Enable', 'off');

                if ( ~app.arumeController.currentSession.isStarted && ~app.arumeController.currentSession.isFinished )
                    set(app.menuRunStartSession, 'Enable', 'on');
                end
                if ( app.arumeController.currentSession.isStarted )
                    set(app.menuRunRestartSession, 'Enable', 'on');
                end
                if ( app.arumeController.currentSession.isStarted && ~app.arumeController.currentSession.isFinished )
                    set(app.menuRunResumeSession, 'Enable', 'on');
                end


                %if ( app.arumeController.currentSession.isFinished )
                set(app.menuAnalyzeMarkData, 'Enable', 'on');
                %else
                %    set(app.menuAnalyzeMarkData, 'Enable', 'off');
                %end

                set(app.menuSessionRename, 'Enable', 'on');
                set(app.menuSessionDelete, 'Enable', 'On');

            elseif ( ~isempty( app.arumeController.currentSession ) && numel(app.arumeController.selectedSessions) >= 1 )
                set(app.menuSessionDelete, 'Enable', 'on');
            else

                set(app.menuRun, 'Enable', 'off');

                set(app.menuSessionDelete, 'Enable', 'off');
                set(app.menuAnalyzeMarkData, 'Enable', 'off');
            end

            % Update past runs

            delete(get(app.menuResumeSessionFrom,'children'));
            set(app.menuResumeSessionFrom, 'Enable', 'off');

            session = app.arumeController.currentSession;
            if (~isempty(session) )
                for i=length(session.pastRuns):-1:1
                    if ( ~isempty( session.pastRuns(i).pastTrialTable ) && ...
                            any(strcmp(session.pastRuns(i).pastTrialTable.Properties.VariableNames, 'TrialNumber')) && ...
                            any(strcmp(session.pastRuns(i).pastTrialTable.Properties.VariableNames, 'DateTimeTrialStart') ))
                        label = [];
                        try
                            label = sprintf('[%d] Trial %d interrupted on %s', i, session.pastRuns(i).pastTrialTable{end,'TrialNumber'},session.pastRuns(i).pastTrialTable{end,'DateTimeTrialStart'}{1});
                        end
                        if ( ~isempty(label) )
                            uimenu(app.menuResumeSessionFrom, ...
                                'Label'     , label, ...
                                'Callback'  , @app.resumeSessionFrom);
                            set(app.menuResumeSessionFrom, 'Enable', 'on');
                        end
                    end
                end
            end

            switch(app.rightPanel.SelectedTab.Title)
                case 'Session table'
                    if ( ~isempty( app.arumeController.currentProject ) )
                        app.sessionTable.Data = app.arumeController.currentProject.sessionsTable;
                    end
                case 'Trial table'
                    if ( ~isempty( app.arumeController.currentSession) && ~isempty(app.arumeController.currentSession.currentRun))
                        app.trialTable.Data = app.arumeController.currentSession.currentRun.pastTrialTable;
                    end
            end

            % update plot menu
            delete(get(app.menuPlot,'children'));

            [plotsList, plotsListAgg]  = app.arumeController.GetPlotList();

            if ( app.arumeController.SIMPLEMODE )
                plotsList = {'VOG_Traces'};
                plotsListAgg = {};
            end

            for i=1:length(plotsList)
                uimenu(app.menuPlot, ...
                    'Label'     , plotsList{i}, ...
                    'Callback'  , @app.Plot);
            end

            for i=1:length(plotsListAgg)
                if ( i==1)
                    separator = 'on';
                else
                    separator = 'off';
                end
                uimenu(app.menuPlot, ...
                    'Label'     , plotsListAgg{i}, ...
                    'Separator' , separator, ...
                    'Callback'  , @app.Plot);
            end
            % end update plot menu
        end
    end


    %%  Callbacks
    methods

        % Figure

        function figureCloseRequest( app, ~, ~ )
            if ( app.closeProject() )
                delete(app.figureHandle)
                Arume('clear');
            end
        end

        function figureResizeFcn( app, ~, ~ )
            figurePosition = app.figureHandle.Position;

            %             app.figureHandle.Position = [app.figureHandle.Position(1:3) max(app.figureHandle.Position(4),600)];
            w = figurePosition(3);  % figure width
            h = figurePosition(4);  % figure height
            h = h;

            m = 2;      % margin between panels
            lw = 300;   % left panel width

            app.leftPanel.Position = [1 1 lw h-2];
            app.rightPanel.Position = [lw+3 1 (w-lw-4) h-2];
        end

        % Project

        function menuProjectCallback( app, ~, ~ )

            % Clean up and refill the recent projects menu

            delete(get(app.menuProjectLoadRecentProject,'children'));

            recentProjects = app.arumeController.configuration.recentProjects;
            recentProjects = recentProjects(~contains(recentProjects','.aruprj'));

            for i=1:min(length(recentProjects),10)
                uimenu(app.menuProjectLoadRecentProject, ...
                    'Text'     , recentProjects{i}, ...
                    'Callback'  , @app.loadProject);
            end
        end

        function newProject(app, ~, ~ )
            if ( app.closeProject() )

                P.Path = app.arumeController.configuration.defaultDataFolder;
                P.Name = 'ProjectName';

                while(1)

                    sDlg.Path = { {['uigetdir(''' P.Path ''')']} };
                    sDlg.Name = P.Name;

                    P = StructDlg(sDlg, 'New project');
                    if ( isempty( P ) )
                        return
                    end

                    if ( ~exist( P.Path, 'dir' ) )
                        uiwait(msgbox('The folder selected does not exist', 'Error', 'Modal'));
                        continue;
                    end

                    if ( ~ArumeCore.Project.IsValidProjectName(P.Name) )
                        uiwait(msgbox('The project name is not valid (no spaces or special signs)', 'Error', 'Modal'));
                        continue;
                    end

                    if ( exist( fullfile(P.Path, P.Name), 'dir') )
                        uiwait(msgbox('There is already a project with that name in that folder.', 'Error', 'Modal'));
                        continue;
                    end

                    break;
                end

                if ( ~isempty( app.arumeController.currentProject ) )
                    app.arumeController.currentProject.save();
                end

                app.arumeController.newProject( P.Path, P.Name);
                app.UpdateGui();
            end
        end

        function loadProject(app, source, ~ )

            if ( app.closeProject() )

                if ( app.menuProjectLoadProject == source )
                    pathname = uigetdir(app.arumeController.configuration.defaultDataFolder, 'Pick a project folder');
                else % load a recent project
                    pathname = get(source,'Label');
                end

                if ( isempty(pathname) || (isscalar(pathname) && (~pathname)) || ~exist(pathname,'dir')  )
                    msgbox('File does not exist');
                    return
                end

                h = waitbar(0,'Please wait..');

                waitbar(1/3)
                if ( ~isempty(app.arumeController.currentProject) )
                    app.arumeController.currentProject.save();
                    app.arumeController.closeProject();
                end

                waitbar(2/3)
                app.arumeController.loadProject(pathname);
                waitbar(3/3)
                app.UpdateGui();
                close(h)
            end

        end

        function loadProjectBackup(app, ~, ~ )

            if ( app.closeProject() )

                [filename, pathname] = uigetfile({'*.zip;*.aruprj', 'Arume backup files (*.zip, *.aruprj'}, 'Pick a project backup');
                if ( isempty(filename) )
                    return
                end

                backupFile = fullfile(pathname, filename);
                newParentPath = uigetdir(app.arumeController.configuration.defaultDataFolder, 'Pick the parent folder for the restored project');
                if ( ~newParentPath  )
                    return
                end

                if ( ~isempty(app.arumeController.currentProject) )
                    app.arumeController.currentProject.save();
                end

                h=waitbar(0,'Please wait..');
                waitbar(1/2)
                app.arumeController.loadProjectBackup(backupFile, newParentPath);
                waitbar(2/2)
                app.UpdateGui();
                close(h)
            end

        end

        function saveProjectBackup(app, ~, ~ )

            file = fullfile(app.arumeController.configuration.defaultDataFolder, [app.arumeController.currentProject.name '-backup-'  datestr(now,'yyyy-mm-dd') '.zip']);

            [filename, pathname] = uiputfile(file, 'Pick a project backup');
            if ( isempty(filename) || ( isscalar(filename) && filename == 0) )
                return
            end

            backupFile = fullfile(pathname, filename);

            h=waitbar(0,'Please wait..');
            waitbar(1/2)
            app.arumeController.saveProjectBackup(backupFile);
            waitbar(2/2)
            app.UpdateGui();
            close(h)

        end

        function closed = closeProject(app, ~, ~ )
            closed = true;
            if (~isempty( app.arumeController.currentProject) )

                choice = questdlg('Do you want to close the current project?', ...
                    'Closing', ...
                    'Yes','Backup before closing','No','No');

                switch choice
                    case 'Backup first'
                        app.saveProjectBackup();
                    case 'No'
                        closed = false;
                        return;
                end

                app.arumeController.closeProject();
                app.UpdateGui();
            end
        end

        % Session

        function sessionTreeCallBack( app, ~, eventdata )

            % first if there are subject nodes selected select all the
            % sessions
            shouldSelect = [];
            for i =1:length( eventdata.SelectedNodes)
                node = eventdata.SelectedNodes(i);
                if ( isempty(node.NodeData) )
                    shouldSelect = cat(1,shouldSelect, node.Children);
                end
            end
            app.sessionTree.SelectedNodes = unique(cat(1,app.sessionTree.SelectedNodes,shouldSelect));


            nodes = app.sessionTree.SelectedNodes;
            sessionListBoxCurrentValue = nan(length(nodes),1);
            for i =1:length( nodes)
                node = nodes(i);
                if ( ~isempty(node.NodeData) )
                    [~,j] = app.arumeController.currentProject.findSessionByIDNumber( node.NodeData );
                    sessionListBoxCurrentValue(i) = j;
                end
            end
            sessionListBoxCurrentValue(isnan(sessionListBoxCurrentValue)) = [];

            if ( sessionListBoxCurrentValue > 0 )
                selectedSessions = app.arumeController.currentProject.sessions(sort(sessionListBoxCurrentValue));
                app.arumeController.selectedSessions = selectedSessions;
            else
                app.arumeController.selectedSessions = [];
            end

            app.UpdateGui(1);
        end

        function newSession( app, ~, ~ )

            possibleExperiments = sort(ArumeCore.ExperimentDesign.GetExperimentList());
            if (~isempty(app.arumeController.currentProject.sessions) )
                lastExperiment = app.arumeController.currentProject.sessions(end).experimentDesign.Name;
            else
                lastExperiment = possibleExperiments{1};
            end

            session.Experiment = lastExperiment;
            session.Subject_Code = '000';
            session.Session_Code = 'Z';

            while(1)
                sessionDlg.Experiment = {possibleExperiments};
                sessionDlg.Experiment{1}{strcmp(possibleExperiments,lastExperiment)} = ['{'  lastExperiment '}'];

                sessionDlg.Subject_Code = session.Subject_Code;
                sessionDlg.Session_Code = session.Session_Code;

                session = StructDlg(sessionDlg, 'New Session');
                if ( isempty( session ) )
                    return
                end

                if ( ~ArumeCore.Session.IsValidSubjectCode(session.Subject_Code) )
                    uiwait(msgbox('The subject code is not valid', 'Error', 'Modal'));
                    continue;
                end
                if ( ~ArumeCore.Session.IsValidSessionCode(session.Session_Code) )
                    uiwait(msgbox('The session code is not valid', 'Error', 'Modal'));
                    continue;
                end

                % Check if session already exists
                if ( isempty(app.arumeController.currentProject.findSession( session.Subject_Code, session.Session_Code)))
                    break;
                else
                    uiwait(msgbox('There is already a session with app name/code', 'Error', 'Modal'));
                end
            end


            % Show the dialog for experiment options if necessary
            experiment = ArumeCore.ExperimentDesign.Create(session.Experiment);
            optionsDlg = experiment.GetExperimentOptionsDialog( );
            if ( ~isempty( optionsDlg) )
                options = StructDlg(optionsDlg, 'Edit experiment options');
                if ( ~isempty( options ) )
                    app.arumeController.newSession( session.Experiment, session.Subject_Code, session.Session_Code, options );
                end
            else
                options = [];
                app.arumeController.newSession( session.Experiment, session.Subject_Code, session.Session_Code, options );
            end

            app.UpdateGui();
        end

        function importSession( app, ~, ~ )
            possibleExperiments = sort(ArumeCore.ExperimentDesign.GetExperimentList());

            if (~isempty(app.arumeController.currentProject.sessions) )
                lastExperiment = app.arumeController.currentProject.sessions(end).experimentDesign.Name;
            else
                lastExperiment = possibleExperiments{1};
            end

            session.Experiment = lastExperiment;
            session.Subject_Code = '000';
            session.Session_Code = 'Z';

            while(1)
                sessionDlg.Experiment = {possibleExperiments};
                sessionDlg.Experiment{1}{strcmp(possibleExperiments,lastExperiment)} = ['{'  lastExperiment '}'];

                sessionDlg.Subject_Code = session.Subject_Code;
                sessionDlg.Session_Code = session.Session_Code;

                session = StructDlg(sessionDlg, 'New Session');
                if ( isempty( session ) )
                    return
                end

                if ( ~ArumeCore.Session.IsValidSubjectCode(session.Subject_Code) )
                    uiwait(msgbox('The subject code is not valid', 'Error', 'Modal'));
                    continue;
                end
                if ( ~ArumeCore.Session.IsValidSessionCode(session.Session_Code) )
                    uiwait(msgbox('The session code is not valid', 'Error', 'Modal'));
                    continue;
                end

                % Check if session already exists
                if ( isempty(app.arumeController.currentProject.findSession( session.Subject_Code, session.Session_Code)))
                    break;
                else
                    uiwait(msgbox('There is already a session with app name/code', 'Error', 'Modal'));
                end
            end


            % Show the dialog for experiment options if necessary
            experiment = ArumeCore.ExperimentDesign.Create(session.Experiment);
            optionsDlg = experiment.GetExperimentOptionsDialog( 1 );
            if ( ~isempty( optionsDlg) )
                options = StructDlg(optionsDlg, 'Edit experiment options');
                if ( isempty( options ) )
                    return;
                end
            else
                options = [];
            end

            app.arumeController.importSession( session.Experiment, session.Subject_Code, session.Session_Code, options  );

            app.UpdateGui();
        end

        function [newSubjectCodes, newSessionCodes] = DlgNewSubjectAndSessionCodes(app)

            sessions = app.arumeController.selectedSessions;

            newSubjectCodes = cell(length(sessions),1);
            newSessionCodes = cell(length(sessions),1);
            for i=1:length(sessions)
                newSubjectCodes{i} = sessions(i).subjectCode;
                newSessionCodes{i} = sessions(i).sessionCode;
            end

            while(1)
                newNamesDlg = [];
                for i=1:length(sessions)
                    newNamesDlg.([sessions(i).name '_New_Subject_Code' ]) = newSubjectCodes{i};
                    newNamesDlg.([sessions(i).name '_New_Session_Code' ]) = newSessionCodes{i};
                end

                P = StructDlg(newNamesDlg);
                if ( isempty( P ) )
                    newSubjectCodes = {};
                    newSessionCodes = {};
                    return
                end

                newSubjectCodes = cell(length(sessions),1);
                newSessionCodes = cell(length(sessions),1);
                for i=1:length(sessions)
                    newSubjectCodes{i} = P.([sessions(i).name '_New_Subject_Code' ]);
                    newSessionCodes{i} = P.([sessions(i).name '_New_Session_Code' ]);
                end

                allgood = 1;
                for i=1:length(sessions)

                    if ( ~ArumeCore.Session.IsValidSubjectCode(newSubjectCodes{i}) || ~ArumeCore.Session.IsValidSessionCode(newSessionCodes{i}) )
                        allgood = 0;
                        break;
                    end

                    if ( ~isempty( app.arumeController.currentProject.findSession(newSubjectCodes{i}, newSessionCodes{i})) )
                        allgood = 0;
                        uiwait(msgbox(['One of the names is repeated ' newSubjectCodes{i} '-' newSessionCodes{i} '.'], 'Error', 'Modal'));
                        break;
                    end
                end

                if ( allgood)
                    break;
                end
            end
        end

        function CopySessions( app, ~, ~ )

            [newSubjectCodes, newSessionCodes] = app.DlgNewSubjectAndSessionCodes();

            if ( ~isempty( newSubjectCodes ) )
                app.arumeController.copySessions(app.arumeController.selectedSessions, newSubjectCodes, newSessionCodes);
                app.UpdateGui();
            end
        end

        function RenameSessions( app, ~, ~ )

            sessions = app.arumeController.selectedSessions;

            [newSubjectCodes, newSessionCodes] = app.DlgNewSubjectAndSessionCodes();

            if ( ~isempty( newSubjectCodes ) )
                for i=1:length(sessions)
                    app.arumeController.renameSession(sessions(i), newSubjectCodes{i}, newSessionCodes{i});
                end
                app.UpdateGui();
            end
        end

        function DeleteSessions( app, ~, ~ )
            choice = questdlg('Are you sure you want to delete the sessions?', ...
                'Closing', ...
                'Yes','No','No');
            switch choice
                case 'Yes'
                    app.arumeController.deleteSessions(app.arumeController.selectedSessions);
                    app.UpdateGui();
            end
        end

        function EditSessionSettings(app, ~, ~ )

            session = app.arumeController.currentSession;

            if ( session.isStarted )
                msgbox('This is session is already started, cannot change settings.');
                return;
            end

            % Show the dialog for experiment options if necessary
            experiment = ArumeCore.ExperimentDesign.Create(session.experimentDesign.Name);
            optionsDlg = experiment.GetExperimentOptionsDialog( );
            if ( ~isempty( optionsDlg) )
                options = StructDlg(optionsDlg,'Edit experiment options',session.experimentDesign.ExperimentOptions);
                if ( isempty( options ) )
                    return;
                end
            else
                options = [];
            end

            session.updateExperimentOptions( options );

            app.UpdateGui();
        end

        function SendDataToWorkspace( app, ~, ~ )
            if (~isempty(app.arumeController.currentSession))
                TrialDataTable = app.arumeController.currentSession.trialDataTable;
                SamplesDataTable = app.arumeController.currentSession.samplesDataTable;
                ProjectDataTable = app.arumeController.currentProject.GetDataTable();
                analysisResults = app.arumeController.currentSession.analysisResults;
                currentRun = app.arumeController.currentSession.currentRun;

                assignin('base','currentRun',currentRun);
                assignin('base','TrialDataTable',TrialDataTable);
                assignin('base','SamplesDataTable',SamplesDataTable);
                assignin('base','AnalysisResults',analysisResults);
                assignin('base','ProjectDataTable',ProjectDataTable);
            end
        end

        function NewSessionAggregate( app, ~, ~ )

            session.Session_Code = 'Z';

            while(1)
                sessionDlg.Session_Code = session.Session_Code;

                session = StructDlg(sessionDlg, 'New Session');
                if ( isempty( session ) )
                    return
                end

                if ( ~ArumeCore.Session.IsValidSessionCode(session.Session_Code) )
                    uiwait(msgbox('The session code is not valid', 'Error', 'Modal'));
                    continue;
                end

                % Check if session already exists
                if ( isempty(app.arumeController.currentProject.findSession( 'AGGREGATED', session.Session_Code)))
                    break;
                else
                    uiwait(msgbox('There is already a session with app name/code', 'Error', 'Modal'));
                end
            end


            % Show the dialog for experiment options if necessary

            app.arumeController.newAggregatedSession( session.Session_Code );

            app.UpdateGui();
        end

        % Run

        function startSession( app, ~, ~ )
            app.arumeController.runSession();
            app.UpdateGui();
        end

        function resumeSession( app, ~, ~ )
            app.arumeController.resumeSession();
            app.UpdateGui();
        end

        function resumeSessionFrom( app, ~, ~ )

            choice = questdlg('Are you sure you want to resume from a past point?', ...
                'Closing', ...
                'Yes','No','No');
            switch choice
                case 'Yes'
                    pastRunNumber = regexp(source.Text, '\[(?<runNumber>\d+)\]', 'names');
                    pastRunNumber = str2double(pastRunNumber.runNumber);

                    app.arumeController.resumeSession(pastRunNumber);
                    app.UpdateGui();
            end
        end

        function restartSession( app, ~, ~ )

            choice = questdlg('Are you sure you want to restart the sessions?', ...
                'Closing', ...
                'Yes','No','No');
            switch choice
                case 'Yes'
                    app.arumeController.restartSession();
                    app.UpdateGui();
            end
        end

        % Analysis

        function RunDataAnalyses( app, ~, ~ )

            optionsDlg = app.arumeController.getAnalysisOptions( );
            optionsLast = app.arumeController.getAnalysisOptionsLast();
            if ( ~isempty( optionsDlg) )
                if ( ~app.arumeController.SIMPLEMODE )
                    options = StructDlg(optionsDlg, 'Edit analysis options', optionsLast);
                    if ( isempty( options ) )
                        return;
                    end
                else
                    options = StructDlg(optionsDlg,'',[],[],'off');
                end
            else
                options = [];
            end

            app.arumeController.runDataAnalyses(options);
            app.UpdateGui();
        end

        function MarkData(app, ~, ~)
            options = app.arumeController.getAnalysisOptionsDefault( );
            app.arumeController.currentSession.experimentDesign.markData(options);
            app.UpdateGui();
        end

        function ClearAllPrepareAndAnalyses( app, ~, ~ )

            choice = questdlg(...
                ['Are you sure you want to delete the analsys data for the sessions?' ...
                'You will need to prepare everything again. You will not lose the raw data.'], ...
                'Clear data', ...
                'Yes','No','No');
            switch choice
                case 'Yes'
                    app.arumeController.clearAllData(app.arumeController.selectedSessions);
                    app.UpdateGui();
            end
        end

        % Plot

        function Plot( app, source, ~ )

            optionsDlg = app.arumeController.getPlotOptionsDlg(source.Label, app.arumeController.selectedSessions);
            optionsLast = app.arumeController.getPlotOptionsLast(source.Label);
            options = [];

            if ( ~isempty(optionsDlg) )
                options = StructDlg(optionsDlg, ['Plot ' source.Label ' options'], optionsLast);
                if ( isempty( options ) )
                    return
                end
            end

            app.arumeController.generatePlot(source.Label, options);

            app.UpdateGui();
        end

        function OpenProjectFolderInExplorer( app, ~, ~ )
            if ( ~isempty(app.arumeController.currentProject))
                winopen(app.arumeController.currentProject.path)
            end
        end

        % Comments
        function commentsTextBoxCallBack( app, ~, ~ )
            app.arumeController.currentSession.updateComment(app.commentsTextBox.Value);
        end

        % Tabs
        function tabRightPanelCallBack(app, ~, eventdata )
            switch(eventdata.NewValue.Title)
                case 'Session table'
                    if ( ~isempty( app.arumeController.currentProject ) )
                        app.sessionTable.Data = app.arumeController.currentProject.sessionsTable;
                    end
                case 'Trial table'
                    if ( ~isempty( app.arumeController.currentSession) && ~isempty(app.arumeController.currentSession.currentRun))
                        app.trialTable.Data = app.arumeController.currentSession.currentRun.pastTrialTable;
                    end
            end
        end

    end

    methods(Access=public)

        function updateSessionTree(app)
            % update session listbox
            if ( ~isempty( app.arumeController.currentProject ) && ~isempty(app.arumeController.currentProject.sessions) )

                % delete sessions that do not exist anymore and updte text
                % of existing ones
                for iSubj = length(app.sessionTree.Children):-1:1
                    subjNode = app.sessionTree.Children(iSubj);
                    for iSess = length(subjNode.Children):-1:1
                        sessNode = subjNode.Children(iSess);
                        session = app.arumeController.currentProject.findSessionByIDNumber( sessNode.NodeData );
                        if ( isempty( session ) )
                            delete(sessNode);
                        else
                            sessNode.Text = session.sessionCode;
                        end
                    end

                    % if the subject does not have children (sessions) delete it too
                    if ( isempty(subjNode.Children) )
                        delete(subjNode);
                    end
                end

                % add nodes for new sessions. Add subject node if necessary
                for i=1:length(app.arumeController.currentProject.sessions)
                    foundSession = 0;
                    foundSubject = 0;
                    session = app.arumeController.currentProject.sessions(i);
                    for iSubj = length(app.sessionTree.Children):-1:1
                        subjNode = app.sessionTree.Children(iSubj);
                        if ( strcmp(subjNode.Text, session.subjectCode ) )
                            foundSubject = iSubj;
                            for iSess = length(subjNode.Children):-1:1
                                sessNode = subjNode.Children(iSess);
                                if (sessNode.NodeData == session.sessionIDNumber)
                                    foundSession = 1;
                                    break;
                                end
                            end
                            break;
                        end
                    end

                    if ( ~foundSession )
                        if ( foundSubject > 0 )
                            newSubjNode = app.sessionTree.Children(foundSubject);
                        else
                            newSubjNode = uitreenode(app.sessionTree);
                            newSubjNode.Text = session.subjectCode;

                            % move to keep alphabetical sorting
                            for iSubj = 1:length(app.sessionTree.Children)
                                subjNode = app.sessionTree.Children(iSubj);
                                [~,j] = sort(upper({subjNode.Text, newSubjNode.Text}));
                                if ( j(1) > 1 )
                                    move(newSubjNode, subjNode, 'before');
                                    break;
                                end
                            end
                        end

                        newSessNode = uitreenode(newSubjNode);
                        newSessNode.Text = session.sessionCode;
                        newSessNode.NodeData = session.sessionIDNumber;

                        % move to keep alphabetical sorting
                        for iSess = 1:length(newSubjNode.Children)
                            sessNode = newSubjNode.Children(iSess);
                            [~,j] = sort(upper({sessNode.Text, newSessNode.Text}) );
                            if ( j(1) > 1 )
                                move(newSessNode, sessNode, 'before');
                                break;
                            end
                        end
                    end
                end

                % find the nodes corresponding with the selected sessions
                nodes = [];
                for i=1:length(app.arumeController.selectedSessions)
                    session = app.arumeController.selectedSessions(i);
                    for iSubj = length(app.sessionTree.Children):-1:1
                        subjNode = app.sessionTree.Children(iSubj);
                        if ( strcmp(subjNode.Text, session.subjectCode ) )
                            for iSess = length(subjNode.Children):-1:1
                                sessNode = subjNode.Children(iSess);
                                if (sessNode.NodeData == session.sessionIDNumber)
                                    nodes = cat(1,nodes, sessNode);
                                end
                                expand(subjNode);
                            end
                        end
                    end
                end
                app.sessionTree.SelectedNodes = nodes;
            else
                delete(app.sessionTree.Children);
                app.sessionTree.SelectedNodes = [];
            end
        end

        function row = GetInfoRow(app, varName, dataItm)
            optionClass = class(dataItm);
            switch(optionClass)
                case 'double'
                    if ( isscalar(dataItm))
                        row = sprintf('%-25s: %s\n', varName, num2str(dataItm));
                    else
                        row = sprintf('%-25s: %s\n', varName, evalc('disp(dataItm)'));
                    end
                case 'char'
                    if ( length(dataItm) <= 50 )
                        row = sprintf('%-25s: %s\n', varName, dataItm);
                    else
                        row = sprintf('%-25s: %s\n', varName, [dataItm(1:20) ' [...] ' dataItm(end-30:end)]);
                    end
                case 'string'
                    fieldText = char(dataItm);
                    if ( length(fieldText) > 50 )
                        fieldText = [fieldText(1:20) ' [...] ' fieldText(end-30:end)];
                    end
                    row = sprintf('%-25s: %s\n', varName, fieldText);
                case 'categorical'
                    row = sprintf('%-25s: %s\n', varName, string(dataItm));
                case 'cell'
                    row = '';
                    if ( length(size(dataItm))<=2 && min(size(dataItm))==1 && ischar(dataItm{1}) && ~isempty(dataItm) )
                        for j=1:length(dataItm)
                            row = [row app.GetInfoRow([varName num2str(j)], dataItm{j})];
                        end
                    else
                        row = sprintf('%-25s: %s\n', varName, 'CELL');
                    end
                case 'struct'
                    row = '';
                    fields = fieldnames(dataItm);
                    for i=1:length(fields)
                        row = [row app.GetInfoRow([ varName '.' fields{i}], dataItm.(fields{i}))];
                    end
                otherwise
                    row = sprintf('%-25s: %s\n', varName, '-');
            end
        end

        function updateInfoBox(app)
            if ( ~isempty( app.arumeController.currentSession ) && length(app.arumeController.selectedSessions)==1 )
                s = '';
                if ( ~isempty(app.arumeController.currentSession.sessionDataTable) )
                    dataTable = app.arumeController.currentSession.sessionDataTable;
                else
                    dataTable = app.arumeController.currentSession.GetBasicSessionDataTable();
                end

                for i=1:min(length(dataTable.Properties.VariableNames), 100)
                    s = [s app.GetInfoRow(dataTable.Properties.VariableNames{i}, dataTable{1,i})];
                end
                if ( length(dataTable.Properties.VariableNames) > 100 )
                    s = [s sprintf('...\n')];
                end


                app.infoBox.Value = s;
            elseif ( length(app.arumeController.selectedSessions) > 1 )
                sessions = app.arumeController.selectedSessions;
                sessionNames = cell(length(sessions),1);
                for i=1:length(sessions)
                    sessionNames{i} = [sessions(i).subjectCode ' __ ' sessions(i).sessionCode];
                end
                sessionNames = strcat(sessionNames,'\n');
                sessionNames = horzcat(sessionNames{:});
                app.infoBox.Value = sprintf(['\nSelected sessions: \n\n' sessionNames]);
            else
                app.infoBox.Value = '';
            end
        end
    end

end

