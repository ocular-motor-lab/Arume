classdef Project < handle
    %PROJECT Class handingling Arume projects
    %
    %
    %
    
    properties( SetAccess = private )
        name            % Name of the project
        path            % Working path of the project
        
        sessions        % Sessions that belong to this project
        
        sessionsTable   % Table with information about the sessions. It
                        % corresponds with the concatenation of all the 
                        % "sessionDataTable" for each session and it will 
                        % be updated every time the project is saved.
    end
        
    methods(Access=public)
        %
        % Save project
        %
        function save( this )
            
            for i =1:length(this.sessions)
                session = this.sessions(i);
                if ( ~isempty(session) )
                    session.save();
                end
            end
            
            try
                % sort the sessions and update the sessions table.
                this.sortSessions();
                this.sessionsTable = this.GetDataTable();
                
                if (~isempty(this.sessionsTable) )
                    writetable(...
                        this.sessionsTable, ...
                        fullfile(this.path, ...
                        [this.name '_ArumeSessionTable.csv']));
                end
                cprintf('blue', '++ ARUME::Done saving project session table to disk.\n');
                disp('======= ARUME EXCEL DATA SAVED TO DISK ==============================')
            catch err
                disp('ERROR saving excel data');
                disp(err.getReport);
            end
            
             
            disp('======= ARUME PROJECT SAVED TO DISK REMEMBER TO BACKUP ==============================')
        end
        
        function backup(this, file)
            if (~isempty(file) )
                % create a backup of the last project file before
                % overriding it
                if ( exist(file,'file') )
                    copyfile(file, [file '.aruback']);
                end
                
                % compress project files
                zip(file , this.path);
            end
        end
        
        %
        % Other methods
        %
        function addSession( this, session)
            if ( ~isempty(this.findSession(session.subjectCode, session.sessionCode) ) )
                error( ['Arume: session already exists use a diferent name ' session.subjectCode '__' session.sessionCode ] );
            end
            
            this.sessions = horzcat(this.sessions, session);
        end
        
        function deleteSession( this, session )
            session.deleteFolders();
            this.sessions( this.sessions == session ) = [];
        end
        
        function [session, i] = findSessionByIDNumber( this, sessionIDNumber)
            
            for i=1:length(this.sessions)
                if ( this.sessions(i).sessionIDNumber ==sessionIDNumber )
                    session = this.sessions(i);
                    return;
                end
            end
            
            % if not found
            session = [];
            i = 0;
        end
        
        function [session, i] = findSession( this, subjectCode, sessionCode)
            
            for i=1:length(this.sessions)
                if ( exist('sessionCode','var') )
                    if ( strcmpi(this.sessions(i).subjectCode, subjectCode) &&  ...
                            strcmpi(this.sessions(i).sessionCode, sessionCode))
                        session = this.sessions(i);
                        return;
                    end
                else
                    if ( strcmpi(this.sessions(i).experimentDesign.Name, experimentName) &&  ...
                            strcmpi([upper(this.sessions(i).subjectCode) upper(this.sessions(i).sessionCode),], subjectCode))
                        session = this.sessions(i);
                        return;
                    end
                end
            end
            
            % if not found
            session = [];
            i = 0;
        end

        function [sessions, idx] = findSessionBySubjectAndExperiment( this, subjectCode, experiment)
            sessions = [];
            idx = [];
            for i=1:length(this.sessions)
                if ( strcmpi(this.sessions(i).subjectCode, subjectCode) &&  ...
                        strcmpi(this.sessions(i).experimentDesign.Name, experiment))
                    if ( isempty(sessions))
                        sessions = this.sessions(i);
                    else
                        sessions(end+1) = this.sessions(i);
                    end
                    idx(end+1) = i;
                end
            end
        end
        
        
        function sortSessions(this)
            
            sessionNames = cell(length(this.sessions),1);
            for i=1:length(this.sessions)
                sessionNames{i} = [this.sessions(i).subjectCode this.sessions(i).sessionCode];
            end
            [~, i] = sort(upper(sessionNames));
            this.sessions = this.sessions(i);
        end
        
        %
        % Analysis methods
        %
        function dataTable = GetDataTable(this, sessions)
            if( ~exist( 'sessions', 'var' ) )
                sessions =  this.sessions;
            end
         
            try
                dataTable = table();
                
                for isess=1:length(sessions)
                    session = sessions(isess);
                    
                    if ( ~isempty( session.sessionDataTable ) )
                        sessionRow = session.sessionDataTable;
                    else
                        sessionRow = session.GetBasicSessionDataTable();
                    end
                    
                    dataTable = VertCatTablesMissing(dataTable, sessionRow);
                end
                
                %disp(dataTable);
                assignin('base','ProjectTable',dataTable);
            catch
                disp('ERROR getting data table');
            end
        end
    end
    
    methods ( Static = true )
        
        %
        % Factory methods
        %
        function project = NewProject( parentPath, projectName )
            
            % check if parentFolder exists
            if ( ~exist( parentPath, 'dir' ) )
                error('Arume: parent folder does not exist');
            end
            
            if ( exist( fullfile(parentPath, projectName), 'dir' ) )
                error('Arume: project folder already exist');
            end
            
            % check if name is a valid name
            if ( ~ArumeCore.Project.IsValidProjectName( projectName ) )
                error('Arume: project name is not valid');
            end
            
            % create project object
            project = ArumeCore.Project();
            
            % Initializes a new project
            project.name               = projectName;
            project.path               = fullfile(parentPath, projectName);
            
            % prepare folder structure
            mkdir( parentPath, projectName );
            
            % save the project
            project.save();
        end
        
        function project = LoadProject( projectPath )
            
            % check if parentFolder exists
            if ( ~exist( projectPath, 'dir' ) )
                error('Arume: parent folder does not exist');
            end
            
            project = ArumeCore.Project();
            ArumeCore.Project.UpdateFileStructure(projectPath); % update old projects (.aruprj) to new structure
            
            
            % Initializes a project loading from a folder
            
            [~, projectName]    = fileparts(projectPath);
            project.name        = projectName;
            project.path     	= projectPath;
            
            % find the session folders
            sessionDirs = sortrows(struct2table(dir(projectPath)),'date');
            sessionDirs = sessionDirs(sessionDirs.isdir & ~strcmp(sessionDirs.name,'.') & ~strcmp(sessionDirs.name,'..'),:);
            
            % load sessions
            for i=1:length(sessionDirs.name)
                sessionName = sessionDirs.name{i};
                session     = ArumeCore.Session.LoadSession( fullfile(projectPath, sessionName) );
                if ( ~isempty(session) )
                    project.addSession(session);
                else
                    disp(['WARNING: session ' sessionName ' could not be loaded. May be an old result of corruption.']);
                    load(fullfile(fullfile(project.path, sessionName), 'sessionData.mat'));
                    sessionData = sessionDataBack;
                    
                        disp(sprintf('... updating %s ...',sessionName));
                        try
                            sessionData.currentRun = ArumeCore.Project.UpdateRun(sessionData.currentRun, sessionData.experimentName );
                            
                            newPastRuns = [];
                            for j=1:length(sessionData.pastRuns)
                                if (isempty( newPastRuns ) )
                                    newPastRuns = ArumeCore.Project.UpdateRun(sessionData.pastRuns(j), sessionData.experimentName );
                                else
                                    newPastRuns = cat(1,newPastRuns,  ArumeCore.Project.UpdateRun(sessionData.pastRuns(j), sessionData.experimentName ));
                                end
                            end
                            sessionData.pastRuns = newPastRuns;
                            
                            filename = fullfile( fullfile(project.path, sessionName), 'ArumeSession.mat');
                            save( filename, 'sessionData' );
                        catch err
                            disp('Could not update');
                            disp(err.getReport);
                        end
                end
            end
                     
            project.sortSessions();
            
            try
                project.sessionsTable = project.GetDataTable();
            catch
                disp('ERROR getting data table');
            end
        end
        
        function project = LoadProjectBackup(file, parentPath)
            
            % check if parentFolder exists
            if ( ~exist( file, 'file' ) )
                error('Arume: file does not exist');
            end
            
            % check if parentFolder exists
            if ( ~exist( parentPath, 'dir' ) )
                error('Arume: parent folder does not exist');
            end
            
            [~, projectName, ext] = fileparts(file);
            
            projectPath = fullfile(parentPath, projectName);
            
            mkdir(projectPath);
            
            % uncompress project file into temp folder
            if ( strcmp(ext, '.aruprj' ) )
                untar(file, parentPath);
            else
                unzip(file, parentPath);
            end
            
            ArumeCore.Project.UpdateFileStructure(projectPath);
            
            project = ArumeCore.Project.LoadProject(projectPath);
        end
        
        function UpdateFileStructure(path)
            [~, projectName] = fileparts(path);
            
            if ( exist(fullfile(path, 'project.mat'),'file') )
                disp('Updated file structure to new version of Arume ...');
                movefile(fullfile(path, 'project.mat'), fullfile(path, [projectName '_ArumeProject.mat']),'f');
                movefile(fullfile(fullfile(path,'dataAnalysis'),'*'), path,'f');
                movefile(fullfile(fullfile(path,'dataRaw'),'*'), path,'f');
                
                if ( exist( fullfile(path, 'analysis'),'dir') )
                    rmdir(fullfile(path, 'analysis'),'s');
                end
                if ( exist( fullfile(path, 'dataAnalysis'),'dir') )
                    rmdir(fullfile(path, 'dataAnalysis'),'s');
                end
                if ( exist( fullfile(path, 'dataRaw'),'dir') )
                    rmdir(fullfile(path, 'dataRaw'),'s');
                end
                if ( exist( fullfile(path, 'figures'),'dir') )
                    rmdir(fullfile(path, 'figures'),'s');
                end
                if ( exist( fullfile(path, 'stimuli'),'dir') )
                    rmdir(fullfile(path, 'stimuli'),'s');
                end
                
                
                projectMatFile = fullfile(path, [projectName '_ArumeProject.mat']);
                
                % load project data
                data = load( projectMatFile, 'data' );
                data = data.data;
                
                badSessions = {};
                for sessionData = data.sessions
                    sessionDataBack = sessionData;
                    try
                        
                        % TEMPORARY
                        oldSessionName = [sessionData.experimentName '_' sessionData.subjectCode sessionData.sessionCode];
                        if ( strcmp(sessionData.experimentName, 'MVSTorsion') )
                            sessionData.experimentName = 'EyeTracking';
                        end
                        newSessionName = [sessionData.experimentName '__' sessionData.subjectCode '__' sessionData.sessionCode];
                        if ( ~strcmp(oldSessionName, newSessionName) )
                            movefile(fullfile(path,oldSessionName) ,fullfile(path,newSessionName))
                        end
                        
                        
                        sessionName = [sessionData.experimentName '__' sessionData.subjectCode '__' sessionData.sessionCode];
                        
                        disp(sprintf('... updating %s ...',sessionName));
                        sessionData.currentRun = ArumeCore.Project.UpdateRun(sessionData.currentRun, sessionData.experimentName );
                        
                        newPastRuns = [];
                        for i=1:length(sessionData.pastRuns)
                            if (isempty( newPastRuns ) )
                                newPastRuns = ArumeCore.Project.UpdateRun(sessionData.pastRuns(i), sessionData.experimentName );
                            else
                                newPastRuns = cat(1,newPastRuns,  ArumeCore.Project.UpdateRun(sessionData.pastRuns(i), sessionData.experimentName ));
                            end
                        end
                        sessionData.pastRuns = newPastRuns;
                        
                        
                        
                        filename = fullfile( fullfile(path, sessionName), 'ArumeSession.mat');
                        save( filename, 'sessionData' );
                    catch
                        disp(['Error loading session ' sessionName]);
                        badSessions{end+1} = sessionName;
                        newSessionName = [sessionData.experimentName '__' sessionData.subjectCode '__' sessionData.sessionCode];
                        save(fullfile(fullfile(path,newSessionName),'sessionData.mat'),'sessionDataBack');
                    end
                end
                data = rmfield(data,'sessions');
                % TODO: maybe save the updated data without sessions.
                
                disp('... Done updating file structure.');
                
                for i=1:length(badSessions)
                    disp(['Error loading session ' badSessions{i} ' this session will not be in the new project.']);
                    if (exist(fullfile(path, badSessions{i}), 'dir'))
                        rmdir(fullfile(path, badSessions{i}));
                    end
                end
            end
        end
        
        function newRun = UpdateRun(runData,experimentName)
            
            experimentDesign = ArumeCore.ExperimentDesign.Create( experimentName );
            experimentDesign.init();
            
            if ( isempty( runData) )
                newRun = ArumeCore.ExperimentRun.SetUpNewRun( experimentDesign );
            else
                
                newRun = runData;
                
                futureConditions = runData.futureConditions;
                f2 = table();
                f2.Condition = futureConditions(:,1);
                f2.BlockNumber = futureConditions(:,2);
                f2.BlockSequenceNumber = futureConditions(:,3);
                f2.Session = ones(size(f2.Condition));
                
                t2 = table();
                for i=1:height(f2)
                    vars = experimentDesign.getVariablesCurrentCondition( f2.Condition(i) );
                    t2 = cat(1,t2,struct2table(vars,'AsArray',true));
                end
                
                newRun.futureTrialTable = [f2 t2];
                
                
                futureConditions = runData.originalFutureConditions;
                f2 = table();
                f2.Condition = futureConditions(:,1);
                f2.BlockNumber = futureConditions(:,2);
                f2.BlockSequenceNumber = futureConditions(:,3);
                f2.Session = ones(size(f2.Condition));
                
                
                t2 = table();
                for i=1:height(f2)
                    vars = experimentDesign.getVariablesCurrentCondition( f2.Condition(i) );
                    t2 = cat(1,t2,struct2table(vars,'AsArray',true));
                end
                
                newRun.originalFutureTrialTable = [f2 t2];
                
                
                pastConditions = runData.pastConditions;
                
                f2 = table();
                f2.TrialNumber = (1:length(pastConditions(:,1)))';
                f2.Session = pastConditions(:,5);
                f2.Condition = pastConditions(:,1);
                f2.BlockNumber = pastConditions(:,3);
                f2.BlockSequenceNumber = pastConditions(:,4);
                f2.Session = ones(size(f2.TrialNumber));
                
                t2 = table();
                for i=1:height(f2)
                    vars = experimentDesign.getVariablesCurrentCondition( f2.Condition(i) );
                    t2 = cat(1,t2,struct2table(vars,'AsArray',true));
                end
                f2 = [f2 t2];
                
                i=1;
                Enum = ArumeCore.ExperimentDesign.getEnum();
                Enum.Events.EYELINK_START_RECORDING     = i;i=i+1;
                Enum.Events.EYELINK_STOP_RECORDING      = i;i=i+1;
                Enum.Events.PRE_TRIAL_START             = i;i=i+1;
                Enum.Events.PRE_TRIAL_STOP              = i;i=i+1;
                Enum.Events.TRIAL_START                 = i;i=i+1;
                Enum.Events.TRIAL_STOP                  = i;i=i+1;
                Enum.Events.POST_TRIAL_START            = i;i=i+1;
                Enum.Events.POST_TRIAL_STOP             = i;i=i+1;
                Enum.Events.TRIAL_EVENT                 = i;i=i+1;
                ev = runData.Events;
                ev(ev(:,4)>height(f2),:) = []; % remove events for trials that are not in pastConditions
                
                f2.TimePreTrialStart = nan(size(f2.TrialNumber));
                f2.TimePreTrialStop = nan(size(f2.TrialNumber));
                f2.TimeTrialStart = nan(size(f2.TrialNumber));
                f2.TimeTrialStop = nan(size(f2.TrialNumber));
                f2.TimePostTrialStart = nan(size(f2.TrialNumber));
                f2.TimePostTrialStop = nan(size(f2.TrialNumber));
                
                f2.TimePreTrialStart(ev(ev(:,3)==Enum.Events.PRE_TRIAL_START ,4)) = ev(ev(:,3)==Enum.Events.PRE_TRIAL_START ,1);
                f2.TimePreTrialStop(ev(ev(:,3)==Enum.Events.PRE_TRIAL_STOP ,4)) = ev(ev(:,3)==Enum.Events.PRE_TRIAL_STOP ,1);
                f2.TimeTrialStart(ev(ev(:,3)==Enum.Events.TRIAL_START ,4)) = ev(ev(:,3)==Enum.Events.TRIAL_START ,1);
                f2.DateTimeTrialStart(ev(ev(:,3)==Enum.Events.TRIAL_START ,4),:) = datestr(ev(ev(:,3)==Enum.Events.TRIAL_START ,2));
                f2.TrialResult = Enum.trialResult.PossibleResults(pastConditions(:,2)+1);
                % from here on only if trialresult is correct or abort
                
                f2.TimeTrialStop(ev(ev(:,3)==Enum.Events.TRIAL_STOP ,4)) = ev(ev(:,3)==Enum.Events.TRIAL_STOP ,1);
                f2.TimePostTrialStart( ev(ev(:,3)==Enum.Events.POST_TRIAL_START ,4)) = ev(ev(:,3)==Enum.Events.POST_TRIAL_START ,1);
                f2.TimePostTrialStop(ev(ev(:,3)==Enum.Events.POST_TRIAL_STOP ,4)) = ev(ev(:,3)==Enum.Events.POST_TRIAL_STOP ,1);
                
                tout = table();
                for i=1:height(f2)
                    if ( isfield(runData.Data{i}, 'trialOutput' ) && ~isempty(runData.Data{i}.trialOutput) )
                        trialOutput = runData.Data{i}.trialOutput;
                        if ( isfield(trialOutput,'Response') && (trialOutput.Response == 'L' || trialOutput.Response == 'R') )
                            trialOutput.Response = categorical(cellstr(trialOutput.Response));
                        elseif ( isfield(trialOutput,'Response') )
                            trialOutput = rmfield(trialOutput,'Response');
                        end
                        if ( isfield(trialOutput,'ReactionTime') && (trialOutput.ReactionTime == -1 || isempty(trialOutput.ReactionTime)) )
                            trialOutput = rmfield(trialOutput,'ReactionTime');
                        end
                    else
                        trialOutput = struct();
                    end
                    
                    if ( ~isempty( tout ) )
                        trialOutput.HasDataOutput = true;
                        trialOutputTable = struct2table(trialOutput,'AsArray',true);
                        tout = VertCatTablesMissing(tout,trialOutputTable);
                    else
                        trialOutput.HasDataOutput = false;
                        tout = struct2table(trialOutput,'AsArray',true);
                    end
                    
                end
                
                if ( ~isempty(tout) )
                    f2 = [f2 tout];
                end
                
                newRun.pastTrialTable = f2;
                
                newRun;
            end
        end
        
        %
        % Other methods
        %
        function result = IsValidProjectName( name )
            result = ~isempty(regexp(name,'^[_a-zA-Z0-9]+$','ONCE') );
        end
    end
end

