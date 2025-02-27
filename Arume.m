classdef Arume < handle
    % ARUME Is a GUI to control experiments and analyze their results.
    %
    %   Usage   : Arume, opens Arume GUI.
    %           : Arume( 'open', 'C:\path\to\project' ), opens a given project
    %
    % A project in Arume consists on multiple experimental SESSIONS and the
    % results ana analyses associted with them.
    %
    % A session is asociated with a given experimental paradigm selected
    % when the session is created.
    %
    % A session can be restarted, paused and resumed. Every time a new
    % experiment run is created containing the data related to each run.
    % That is, if you run the experiment, almost finish and the restart
    % over. All the data will be saved. For the first partial run and for
    % the second complete run.
    % 
    % A project can have sessions of different paradigms but a session will
    % have runs of one individual paradigm.
    %
    % The projects can be managed with the GUI but also with command line.
    
    properties (Constant)
        version_number = '1.1.20240226'
    end
    
    properties( SetAccess=private )
        gui                 % Current gui associated with the controller
        configuration       % Configuration options saved into a mat file in the Arume folder
        possiblePlots
        priv_currentProject      
        priv_selectedSessions   

        SIMPLEMODE = 0 ;
    end
    
    properties(Dependent=true)
        currentProject      % Current working project
        currentSession      % Current selected session (empty if none)
        selectedSessions    % Current selected sessions (if multiple selected enabled)
    end
    
    % methods for dependent properties
    
    methods
        function session = get.currentSession( this )
            if ( length(this.selectedSessions) >= 1 )
                session = this.selectedSessions(end);
            else
                session = [];
            end
        end
        
        function selectedSessions = get.selectedSessions( this )
            selectedSessions = this.priv_selectedSessions;
        end
        
        function set.selectedSessions(this, sessions)
            
            % Updates the current session selection
            this.priv_selectedSessions = sessions;
            
            % update plotlist every time the selection of sessions change
            this.UpdatePossiblePlots();
        end
        
        function currentProject = get.currentProject( this )
            currentProject = this.priv_currentProject;
        end
        
        function set.currentProject(this, project)
            
            % Updates the current session selection
            this.priv_currentProject = project;
            
            if ( ~isempty(this.currentProject) && ~isempty(this.currentProject.sessions) )
                this.selectedSessions = this.currentProject.sessions(1);
            else
                this.selectedSessions = [];
            end
            
            if (~isempty(project) )
                
                if ( ~isfield(this.configuration, 'recentProjects' ) )
                    this.configuration.recentProjects = {};
                end
                
                % remove the current file and keep it down to a max number
                maxNumberRecentProjects = 30;
                this.configuration.recentProjects = unique(this.configuration.recentProjects,'stable');
                if (~isempty( this.configuration.recentProjects ) )
                    this.configuration.recentProjects =  this.configuration.recentProjects(1:min(maxNumberRecentProjects,length(this.configuration.recentProjects)));
                end
                this.configuration.recentProjects(strcmp(this.configuration.recentProjects, project.path)) = [];
                % add it again at the top
                this.configuration.recentProjects = {project.path this.configuration.recentProjects{:}};
                
                this.saveConfiguration();
                
            end
        end
    end
    
    methods( Access=public )
        
        %
        % Main constructor
        %
        
        function arumeController = Arume(command, param)
            
            % Persistent variable to keep the singleton to make sure there is only one
            % arume controller loaded at any point in time. That way we can open the UI
            % and then also call arume in the command line to get a reference to the
            % controller and write scripts working with the current project.
            persistent arumeSingleton;
            
            if ( isempty( arumeSingleton ) )
                % The persistent variable gets deleted with clear all. However,
                % variables within the UI do not and stay until UI is
                % closed. So, we can search for the handle of the UI window
                % and get the controller from there. This way we avoid
                % problems if clear all is called with the UI open 
                % and then Arume is called again.
                h = findall(0,'tag','Arume');
                if ( ~isempty(h) )
                    arumeSingleton = h.UserData.arumeController;
                end
            end
            
            useGui = 1;
            loadlast = 0;
            
            % option to clear the singleton
            if ( exist('command','var') )
                switch (command )
                    case 'open'
                        if ( exist('param','var') )
                            projectPath = param;
                        end
                    case 'nogui'
                        if ( exist('param','var') )
                            projectPath = param;
                        end
                        useGui = 0;
                    case 'clear'
                        clear arumeSingleton;
                        clear arumeController;
                        return;
                    case 'last'
                        loadlast = 1;
                end
            end
            
            if isempty(arumeSingleton)
                % Initialization, object is created automatically
                % (this is the constructor) and then initialized
                
                arumeSingleton = arumeController;
                arumeSingleton.init();

                arumeSingleton.SIMPLEMODE = arumeSingleton.configuration.userMode=="Simple";
            end
            
            if ( loadlast )
                projectPath = arumeSingleton.configuration.recentProjects{1};
            end
            
            if ( exist('projectPath','var') )
                arumeSingleton.loadProject( projectPath );
            end
            
            if ( useGui )
                if ( isempty(arumeSingleton.gui) || ~arumeSingleton.gui.isvalid)
                    % Load the GUI
                    arumeSingleton.gui = ArumeCore.ArumeGui( arumeSingleton );
                end
                % make sure the Arume gui is on the front and update
                figure(arumeSingleton.gui.figureHandle)
                arumeSingleton.gui.UpdateGui();
            end
            
            arumeController = arumeSingleton;
            
            if nargout == 0
                clear arumeController
            end
        end
        
        function init( this )
            
            this.initConfiguration();
        
            % create a default folder for data
            if ( ~exist( this.configuration.defaultDataFolder, 'dir') )
                % find the folder of arume
                folder = fileparts(which('Arume'));
                if ( ~exist(fullfile(folder, 'ArumeData'), 'dir') )
                    mkdir(folder, 'ArumeData');
                end
            end
        end
        
        function initConfiguration( this )
            % find the folder of arume
            folder = fileparts(which('Arume'));
            configFile = fullfile(folder,'arumeconf.json');
            
            % find the configuration file
            if ( ~exist(configFile,'file'))
                if ( exist(fullfile(folder,'arumeconf.mat'),'file'))
                    confdata = load(fullfile(folder,'arumeconf.mat'));
                    conf = confdata.conf;
                else
                    conf = [];
                    conf.recentProjects = {};
                end
                this.configuration = conf;
                this.saveConfiguration();
            else
                % read the json file
                fid = fopen(configFile);
                raw = fread(fid,inf);
                str = char(raw');
                fclose(fid);
                conf = jsondecode(str);
            end
            
            % double check configuration fields
            if ( ~isfield( conf, 'recentProjects') )
                conf.recentProjects = {};
            end
            if ( length(conf.recentProjects)==1 && isempty(conf.recentProjects{1}) )
                conf.recentProjects = {};
            end
            if ( ~isfield( conf, 'defaultDataFolder') )
                conf.defaultDataFolder = fullfile(folder, 'ArumeData');
            end
            if ( ~iscell(conf.recentProjects) )
                conf.recentProjects = {conf.recentProjects};
            end


            if ( ~isfield( conf, 'userMode') )
                conf.userMode = 'Normal';
            end
            
            % save the updated configuration
            this.configuration = conf;
            this.saveConfiguration()
        end
        
        function saveConfiguration( this, conf )
            if ( ~exist('conf','var'))
                conf = this.configuration;
            end
            [folder, ~, ~] = fileparts(which('Arume'));
            configFile = fullfile(folder,'arumeconf.json');

            jconf = jsonencode(conf,'PrettyPrint',true);
            fid = fopen(configFile,'w');
            fprintf(fid,'%s',jconf);
            fclose(fid);
        end
        
        %
        % Managing projects
        %
        
        function newProject( this, parentPath, projectName )
            % Creates a new project
            this.currentProject = ArumeCore.Project.NewProject( parentPath, projectName);
        end
        
        function loadProject( this, folder )
            % Loads a project from a project folder
            
            if ( ~exist( folder, 'dir') )
                msgbox( 'The project folder does not exist.');
                return;
            end
            
            if ( ~isempty(this.currentProject) && strcmp(this.currentProject.path, folder))
                disp('Loading the same project folder that is currently loaded');
                return;
            end
            
            this.currentProject = ArumeCore.Project.LoadProject( folder );
        end
        
        function loadProjectBackup( this, file, parentPath )
            % Loads a project from a project file
            if ( ~exist( file, 'file') )
                msgbox( 'The project file does not exist.');
                return;
            end
            
            [~,projectName] = fileparts(file);
            
            if ( ~isempty(this.currentProject) && strcmp(this.currentProject.name, projectName))
                disp('Loading the same project file that is currently loaded');
                return;
            end
            
            this.currentProject = ArumeCore.Project.LoadProjectBackup( file, parentPath );
        end
        
        function saveProjectBackup(this, file)
            if ( exist( file, 'file') )
                msgbox( 'The file already exists.');
            end
            
            this.currentProject.backup(file);
        end
        
        function closeProject( this )
            % Closes the current project (always saves)
            this.currentProject.save();
            this.currentProject = [];
        end
        
        %
        % Managing sessions
        %
        
        function session = newSession( this, experiment, subjectCode, sessionCode, experimentOptions )
            % Crates a new session to start the experiment and collect data
            
            % check if session already exists with that subjectCode and
            % sessionCode
            if ( ~isempty(this.currentProject.findSession(subjectCode, sessionCode) ) )
                error( 'Arume: session already exists, use a diferent name' );
            end
            
            session = ArumeCore.Session.NewSession( this.currentProject.path, experiment, subjectCode, sessionCode, experimentOptions );
            this.currentProject.addSession(session);
            this.selectedSessions = session;
            this.currentProject.save();
        end
    
        function session = newAggregatedSession( this, sessionCode )
            options.sessions  = {};
            for i=1:length(this.selectedSessions)
                options.sessions{end+1} = this.selectedSessions(i).name;
            end
            session = this.newSession( 'AGGREGATED', 'AGGREGATED', sessionCode, options );
        end
        
        function session = importSession( this, experiment, subjectCode, sessionCode, options )
            % Imports a session from external files containing the data. It
            % will not be possible to run this session
            
            % check if session already exists with that subjectCode and
            % sessionCode
            if ( ~isempty(this.currentProject.findSession(subjectCode, sessionCode) ) )
                error( 'Arume: session already exists, use a diferent name' );
            end
            
            session = ArumeCore.Session.NewSession( this.currentProject.path, experiment, subjectCode, sessionCode, options, 1 );
            this.currentProject.addSession(session);
            this.selectedSessions = session;
            
            session.importSession();
            
            this.currentProject.save();
        end
        
        function renameSession( this, session, subjectCode, sessionCode)
            % Renames the current session
            
            for session1 = this.currentProject.sessions
                if ( isequal(subjectCode, session1.subjectCode) && isequal( sessionCode, session1.sessionCode) )
                    error( 'Arume: session already exists use a diferent name' );
                end
            end
            disp(['Renaming session' session.subjectCode ' - ' session.sessionCode ' to '  subjectCode ' - ' sessionCode]);
            
            [~, i] = this.currentProject.findSession(session.subjectCode, session.sessionCode);
            this.currentProject.sessions(i).rename(subjectCode, sessionCode);
            this.currentProject.save();
        end
        
        function newSessions = copySessions( this, sessions, newSubjectCodes, newSessionCodes)
            
            newSessions = [];
            for i =1:length(sessions)
                newSession = sessions(i).copy(newSubjectCodes{i}, newSessionCodes{i});
                this.currentProject.addSession(newSession);
                newSessions = cat(1,newSessions, newSession);
            end
            
            this.selectedSessions = newSessions;
            this.currentProject.save();
        end
        
        function deleteSessions( this, sessions )
            
            for i =1:length(sessions)
                this.currentProject.deleteSession(sessions(i));
            end
            
            this.selectedSessions = [];
            this.currentProject.save();
        end
        
        function options = getDefaultExperimentOptions(this, experiment, importing)
            if ( ~exist('importing','var') )
                importing = 0;
            end
            options = ArumeCore.ExperimentDesign.GetDefaultExperimentOptions(experiment, importing);
        end
        
        %
        % Running sessions
        %
        
        function runSession( this )
            % Start running the experimental session
            
            this.currentSession.start();
            this.currentProject.save();
        end
        
        function resumeSession( this )
            % Resumes running the experimental session
            
            this.currentSession.resume();
            this.currentProject.save();
        end
        
        function resumeSessionFrom( this, runNumber )
            % Resumes running the experimental session
            
            this.currentSession.resumeFrom(runNumber);
            this.currentProject.save();
        end
        
        function restartSession( this )
            % Restarts a session from the begining. Past data will be saved.
            
            this.currentSession.restart();
            this.currentProject.save();
        end
        
        %
        % Analyzing
        %
        function runDataAnalyses(this, options, sessions)
            useWaitBar = 0;
            
            if ( ~exist('sessions','var') )
                sessions = this.selectedSessions;
                useWaitBar = 1;
            end
            
            if ( exist( 'options','var') )
                this.configuration.LastAnalysisOptions = options;
                this.saveConfiguration();
            end
            
            n = length(sessions);
            
            if (useWaitBar)
                h = waitbar(0,'Please wait...');
            end
            
            for i =1:n
                try
                    cprintf('blue', '++ ARUME::running analyses for session %s\n', sessions(i).name);
                    session = sessions(i);
                    session.runAnalysis(options);
                    if ( useWaitBar )
                        waitbar(i/n,h)
                    end
                catch ex
                    
                    beep
                    cprintf('red', '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n')
                    cprintf('red', '!!!!!!!!!!!!! ARUME ERROR RUNNING ANALYSES: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n')
                    cprintf('red', '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n')
                    disp(ex.getReport);
                    cprintf('red', '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n')
                    cprintf('red', '!!!!!!!!!!!!! END ARUME ERROR: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n')
                    cprintf('red', '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n')
                end
            end
            
            this.currentProject.save();
            
            if (useWaitBar)
                close(h);
            end
        end
        
        function clearAllData(this, sessions)
            for i =1:length(sessions)
                sessions(i).clearData();
            end
        end
        
        function options = getAnalysisOptions(this, sessions)
            
            if ( ~exist('sessions','var') )
                sessions = this.selectedSessions;
            end
                        
            options = struct();
            options.Prepare_For_Analysis_And_Plots = { {'0','{1}'} };
            
            for session = sessions
                dlg1 = session.experimentDesign.GetAnalysisOptionsDialog();
                if ( ~isempty(dlg1) )
                    f1 = fields(options);
                    f2 = fields(dlg1);
                    for i=1:length(f2)
                        if ( ~any(contains(f1,f2{i})) )
                            options.(f2{i}) = dlg1.(f2{i});
                        end
                    end
                end
            end
        end
        
        function options = getAnalysisOptionsLast(this)
            options = [];
            if (isfield( this.configuration, 'LastAnalysisOptions') )
                options = this.configuration.LastAnalysisOptions;  
            end
        end
                
        function options = getAnalysisOptionsDefault(this, sessions)
            
            if ( ~exist('sessions','var') )
                sessions = this.selectedSessions;
            end
            options = StructDlg(this.getAnalysisOptions(sessions),'',[],[],'off');
        end
        
        %
        % Plotting
        %
        function generatePlot( this, plotName, options)
            
            try
                % save last options
                if ( exist( 'options','var') )
                    this.configuration.LastPlotOptions.(plotName) = options;
                    this.saveConfiguration();
                end
                
                plots = this.possiblePlots(strcmp(this.possiblePlots.PlotName, plotName),:);
                if ( height(plots) > 1 )
                    experimentClassName = class(this.currentSession.experimentDesign);
                    plots = plots(strcmp(this.possiblePlots.ExperimentClassName, experimentClassName),:);
                end
                plot = table2struct(plots(1,:));

                
                if ( ~plot.IsAggregate )
                    if (~isfield(options,'Combine_Sessions_Into_One_Figure') || ~options.Combine_Sessions_Into_One_Figure)
                        % Single sessions plot
                        for session = this.selectedSessions
                            experimentClassName = class(session.experimentDesign);
                            plot = table2struct(this.possiblePlots(strcmp(this.possiblePlots.PlotName, plotName) & strcmp(this.possiblePlots.ExperimentClassName, experimentClassName),:));
                            if ( ~plot.HasOptions)
                                session.experimentDesign.(plot.PlotMethodName)();
                            else
                                session.experimentDesign.(plot.PlotMethodName)(options);
                            end
                        end
                    else
                        nplot1 = [1 1 1 2 2 2 2 2 3 2 3 3 3 3 4 4 4 4 4 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5];
                        nplot2 = [1 2 3 2 3 3 4 4 3 5 4 4 5 5 4 5 5 5 5 5 5 5 5 5 6 6 6 6 6 6 7 7 7 7 7 7];
                        nSessions = length(this.selectedSessions);
                        p1 = nplot1(nSessions);
                        p2 = nplot2(nSessions);
                        iSession = 0;
                        
                        combinedFigure = figure;
                        
                        % Single sessions plot
                        for session = this.selectedSessions
                            experimentClassName = class(session.experimentDesign);
                            plot = table2struct(this.possiblePlots(strcmp(this.possiblePlots.PlotName, plotName) & strcmp(this.possiblePlots.ExperimentClassName, experimentClassName),:));
                            
                            iSession = iSession+1;
                            handles = get(0,'children');
                            
                            if ( ~exist( 'options','var') || numel(fieldnames((rmfield(options, 'Combine_Sessions_Into_One_Figure'))))==0)
                                session.experimentDesign.(plot.PlotMethodName)();
                            else
                                session.experimentDesign.(plot.PlotMethodName)(options);
                            end
                            
                            set(gcf,'name',strrep(session.shortName,'_',' '))
                            newhandles = get(0,'children');
                            for iplot = 1:(length(newhandles)-length(handles))
                                
                                idx = length(handles)+1;
                                axorig = get(newhandles(1),'children');
                                theTitle1 = get(gca,'title');
                                theName = strrep(get(newhandles(1),'name'),'_',' ');
                                theTitle = [theName ' - ' theTitle1.String];
                                if ( iSession > 1 )
                                    axcopy = copyobj(axorig(end), combinedFigure);
                                else
                                    % copy all including legend
                                    axcopy = copyobj(axorig(:), combinedFigure);
                                end
                                ax = subplot(p1,p2,iSession,axcopy(end));
                                title(ax,theTitle);
                            end
                            
                            close(setdiff( newhandles,handles))
                        end
                    end
                    
                else
                    experimentClassName = class(this.currentSession.experimentDesign);
                    plot = table2struct(this.possiblePlots(strcmp(this.possiblePlots.PlotName, plotName) & strcmp(this.possiblePlots.ExperimentClassName, experimentClassName),:));
                    
                    % Aggregate session plots
                    if ( ~plot.HasOptions)
                        this.currentSession.experimentDesign.(plot.PlotMethodName)(this.selectedSessions);
                    else
                        this.currentSession.experimentDesign.(plot.PlotMethodName)( this.selectedSessions, options );
                    end
                end
            catch err
                beep
                cprintf('red', '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n')
                cprintf('red', '!!!!!!!!!!!!! ARUME PLOT ERROR: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n')
                cprintf('red', '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n')
                cprintf('red', '\n')
                cprintf('red', 'Error ploting, try preparing the session first!\n')
                cprintf('red', '\n')
                disp(err.getReport);
                cprintf('red', '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n')
                cprintf('red', '!!!!!!!!!!!!! END PLOT ARUME ERROR: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n')
                cprintf('red', '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n')
            end
        end
        
        function [plotListNames, plotListNamesAgg] = GetPlotList( this )
            if ( ~isempty(this.possiblePlots) )
                plotListNames = unique(this.possiblePlots.PlotName(~this.possiblePlots.IsAggregate));
                plotListNamesAgg = unique(this.possiblePlots.PlotName(this.possiblePlots.IsAggregate==1));
            else
                plotListNames = {};
                plotListNamesAgg = {};
            end
        end
                
        function options = getPlotOptionsLast(this, plotName)
            options = [];
            if (isfield( this.configuration, 'LastPlotOptions') && isfield( this.configuration.LastPlotOptions, plotName) )
                options = this.configuration.LastPlotOptions.(plotName);  
            end
        end
        
        function options = getPlotOptionsDlg(this, plotName, sessions)
            plot = table2struct(this.possiblePlots(strcmp(this.possiblePlots.PlotName, plotName),:));
            
            options = [];
            
            if ( plot.HasOptions)
                [~,options] = this.currentSession.experimentDesign.(plot.PlotMethodName)('get_options');
            end
            if ( ~plot.IsAggregate && ~strcmp(plotName,'VOG_Data_Explorer') && numel(sessions) > 1)
                options.Combine_Sessions_Into_One_Figure = { {'0','{1}'} };
            end
        end
        
    end
    
    methods( Access = private )
        function UpdatePossiblePlots(this)
            this.possiblePlots = table();
            
            for j=1:length(this.selectedSessions )
                experimentDesign = this.selectedSessions (j).experimentDesign;
                experimentClassName = class(experimentDesign);
                if (~isempty(this.possiblePlots) && any(contains(this.possiblePlots.ExperimentClassName, experimentClassName)))
                    continue;
                end
                methodList = meta.class.fromName(experimentClassName).MethodList;
                for i=1:length(methodList)
                    methodName = methodList(i).Name;
                    HasOptions = false;
                    if ( contains( methodName, 'Plot_') )
                        IsAggregate = false;
                        if ( length(methodList(i).InputNames) == 2 && strcmp(methodList(i).InputNames{2},'options') )
                            HasOptions = true;
                        end
                        name = strrep(methodName, 'Plot_' ,'');
                        
                    elseif ( contains( methodName, 'PlotAggregate_') )
                        IsAggregate = true;
                        if ( length(methodList(i).InputNames) == 3 && strcmp(methodList(i).InputNames{3},'options') )
                            HasOptions = true;
                        end
                        name = strrep(methodName, 'PlotAggregate_' ,'');
                    else
                        continue;
                    end
                    newPlot = table(string(name),string(experimentClassName), string(methodName),HasOptions, IsAggregate, ...
                        'VariableNames', {'PlotName', 'ExperimentClassName', 'PlotMethodName', 'HasOptions', 'IsAggregate'});
                    this.possiblePlots = vertcat(this.possiblePlots, newPlot);
                end
            end
                
            if (~isempty(this.possiblePlots ) )
                if ( height(unique(this.possiblePlots(:,{'PlotName','ExperimentClassName','PlotMethodName'}))) ~= height(unique(this.possiblePlots(:,{'PlotName','ExperimentClassName'}))) )
                    error('There is a conflict in plotting methods. Same plot is aggreagete in one experiment design and not aggregate in another experiment design');
                end
            end
        end
        
    end
end

