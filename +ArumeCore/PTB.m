classdef PTB < handle
    %PTB this class is a wrapper of psychtoolbox for Arume
    properties
        screens = [];
        selectedScreen = 1;
        
        window = [];
        wRect = [];
        
        black = [];
        white = [];
        
        dlgTextColor = [255 255 255];
        dlgBackgroundScreenColor = [0 0 0];
        
        frameRate = [];
        nominalFrameRate = [];
        
        
        reportedmmWidth = [];
        reportedmmHeight = [];
        
        pxWidth = [];
        pxHeight = [];
        
        windiwInfo = [];
        
        mmWidth = [];
        mmHeight = [];
        
        distanceToMonitor = [];
        
        NumFlips = 0;
        NumSlowFlips = 0;
        NumSuperSlowFlips = 0;
    end
    
    properties(Access=private)
        lastfliptime = 0;
    end

    methods (Static=true)
        function Test( debugMode)

            if ( ~exist('debugMode','var'))
                debugMode = 1;
            end

            Graph = ArumeCore.PTB(debugMode );

            Graph.DlgHitKey( 'Hit a key to continue',[],[]);

            result = Graph.DlgSelect( ...
                'Choose an option:', ...
                { 'n' 'c' 'q'}, ...
                { 'Choose n' 'Choose c' 'quit'} , [],[])

            Graph.Clear();
        end

        % Gets the options for 
        function displayOptions = GetDisplayOptions()
            displayOptions.ForegroundColor      = 0;
            displayOptions.BackgroundColor      = 127;
            displayOptions.ScreenWidth          = { 142.8 '* (cm)' [1 3000] };
            displayOptions.ScreenHeight         = { 80 '* (cm)' [1 3000] };
            displayOptions.ScreenDistance       = { 85 '* (cm)' [1 3000] };
            displayOptions.ShowTrialTable       = { {'0','{1}'} };
            displayOptions.PlaySound            = { {'0','{1}'} };
            displayOptions.StereoMode           = { 0 '* (mode)' [0 9] }; % SR added, 0 should be the default
            displayOptions.SelectedScreen       = { 2 '* (screen)' [0 5] }; % SR added, screen 2 should perhaps be the default
        end

        function displayOptions = GetDisplayOptionsDefault()
            displayOptions = ArumeCore.PTB.GetOptions();
            displayOptions = StructDlg(displayOptions,'',[],[],'off');
        end
    end
    
    methods
        
        %% Display
        function graph = PTB( debugMode, displayOptions)
            
            if ( ~exist('displayOptions','var') )
                displayOptions = ArumeCore.PTB.GetDisplayOptionsDefault();
            end

            selectedScreenFromOptions = displayOptions.SelectedScreen;
            stereoMode = displayOptions.StereoMode;

            foregroundColor = displayOptions.ForegroundColor;
            backgroundColor = displayOptions.BackgroundColor;
            mmWidthOptions    = displayOptions.ScreenWidth;
            mmHeightOptions          = displayOptions.ScreenHeight;
            distanceToMonitorOptions = displayOptions.ScreenDistance;



            % -- GRAPHICS KEYBOARD and MOUSE
            %-- hide the mouse cursor during the experiment
            if ( ~debugMode)
                HideCursor;
                ListenChar(2);
            else
                ListenChar(1);
            end
            
            %           Screen('Preference', 'VisualDebugLevel', 3);
            Screen('Preference', 'SkipSyncTests', 1);
            Screen('Preference', 'VisualDebugLevel', 0);
            
            %-- screens
            
            graph.screens = Screen('Screens');
            if selectedScreenFromOptions > max(graph.screens) % if the defaul screen is 2 but you don't have screen number 2
                graph.selectedScreen=max(graph.screens);
            else
                graph.selectedScreen= selectedScreenFromOptions; %max(graph.screens);
            end
            % graph.selectedScreen=1;

            % fix issue with mac multisampling, especially with custom shaders?
            if ismac 
                multisample = [];
            else
                multisample = 10;
                Screen('Preference', 'ConserveVRAM', 64);
            end

            %-- window
            % Screen('Preference', 'ConserveVRAM', 64);
            if (~debugMode)
                [graph.window, graph.wRect] = Screen('OpenWindow', graph.selectedScreen, backgroundColor, [], [], [], stereoMode, multisample);%SR changed stereomode to 4 from 0
            else
                [graph.window, graph.wRect] = Screen('OpenWindow', graph.selectedScreen, backgroundColor, [10 10 900 600], [], [], stereoMode, multisample); %SR changed stereomode to 4 from 0
            end
            
            %-- color
            
            graph.black = BlackIndex( graph.window );
            graph.white = WhiteIndex( graph.window );
            graph.dlgTextColor = foregroundColor;
            graph.dlgBackgroundScreenColor = backgroundColor;
            
            % FOR OKN
            AssertOpenGL;
            Screen('BlendFunction', graph.window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            Priority(MaxPriority(graph.window));
            
            
            %-- font
            Screen('TextSize', graph.window, 18);
            
            
            %-- frame rate
            graph.frameRate         = Screen('FrameRate', graph.selectedScreen);
            ifi                                             = Screen('GetFlipInterval', graph.window);
            if graph.frameRate==0
                graph.frameRate=1/ifi;
            end
            graph.nominalFrameRate  = Screen('NominalFrameRate', graph.selectedScreen);
            
            %-- size
            [graph.reportedmmWidth, graph.reportedmmHeight] = Screen('DisplaySize', graph.selectedScreen);
            [graph.pxWidth, graph.pxHeight]                 = Screen('WindowSize', graph.window);
            graph.windiwInfo                                = Screen('GetWindowInfo',graph.window);           
            %TODO: force resolution and refresh rate
            

            %-- physical dimensions
            graph.mmWidth           = mmWidthOptions;
            graph.mmHeight          = mmHeightOptions;
            graph.distanceToMonitor = distanceToMonitorOptions;
        end
        
        function Clear(graph)
            ShowCursor;
            ListenChar(0);
            Priority(0);
            
            Screen('CloseAll');
        end
        
        function ResetBackground( this )
            Screen('FillRect', this.window, this.dlgBackgroundScreenColor);
            Screen('Flip', graph.window);
        end
        
        function ResetFlipTimes(this)
            this.NumFlips = 0;
            this.NumSlowFlips = 0;
            this.NumSuperSlowFlips = 0;
        end
        
        %% Flip
        %--------------------------------------------------------------------------
        function fliptime = Flip( this, exper, thisTrialData, secondsRemaining )
            
            Enum = ArumeCore.ExperimentDesign.getEnum();
            
            if ( nargin == 4)
                if ( exper.ExperimentOptions.Debug.DebugMode )
                    Screen('DrawText', this.window, sprintf('%i seconds remaining...', round(secondsRemaining)), 20, 50, this.white);
                    currentline = 50 + 25;
                    vNames = thisTrialData.Properties.VariableNames;
                    for iVar = 1:length(vNames)
                        if ( ischar(thisTrialData.(vNames{iVar})) )
                            s = sprintf( '%s = %s',vNames{iVar},thisTrialData.(vNames{iVar}) );
                        elseif ( isnumeric(thisTrialData.(vNames{iVar})) )
                            s = sprintf( '%s = %s',vNames{iVar},num2str(thisTrialData.(vNames{iVar})) );
                        else
                            s = sprintf( '%s = -',vNames{iVar});
                        end
                        Screen('DrawText', this.window, s, 20, currentline, this.white);
                        
                        currentline = currentline + 25;
                    end
                    %
                    %                             if ( ~isempty( this.EyeTracker ) )
                    %                                 draweye( this.EyeTracker.eyelink, graph)
                    %                             end
                end
            end
            
            fliptime = Screen('Flip', this.window);
            this.NumFlips = this.NumFlips +1;
            if ( this.lastfliptime>0 && fliptime-this.lastfliptime > 1.5/this.frameRate)
                this.NumSlowFlips = this.NumSlowFlips + 1;
            end
            if ( this.lastfliptime>0 && fliptime-this.lastfliptime > 10/this.frameRate)
                this.NumSuperSlowFlips = this.NumSuperSlowFlips + 1;
            end
            this.lastfliptime = fliptime;
            
            %-- Check for keyboard press
            [keyIsDown,secs,keyCode] = KbCheck;
            if keyCode(Enum.keys.ESCAPE)

                % in the special case where a movie is playing, close the
                % movie to stop audio playback (and subsequent interference
                % with other movies if the user continues the experiment)
                Screen('CloseMovie');

                if nargin >1
                    exper.abortExperiment();
                else
                    throw(MException('PSYCORTEX:USERQUIT', ''));
                end
            end
        end
        
        %% dva2pix
        %--------------------------------------------------------------------------
        function pix = dva2pix( this, dva )
            
            horPixPerDva = this.pxWidth/2 / (atan(this.mmWidth/2/this.distanceToMonitor)*180/pi);
            %             verPixPerDva = this.pxHeight/2 / (atan(this.mmHeight/2/this.distanceToMonitor)*180/pi);
            
            pix = round( horPixPerDva * dva );
            
            % TODO: improve
            
            % function pix = psyCortex_dva2pixExact( graph, poit1, point2 )
            %
            % distanceToCenter = min
            %
            % horPixPerDva = graph.pxWidth / atan(graph.mmWidth/2/graph.distanceToMonitor)*180/pi;
            % verPixPerDva = graph.pxHeight / atan(graph.mmHeight/2/graph.distanceToMonitor)*180/pi;
        end
        
        %% pix2dva
        function dva = pix2dva( this, pix )
            
            horPixPerDva = this.pxWidth/2 / (atan(this.mmWidth/2/this.distanceToMonitor)*180/pi);
            %             verPixPerDva = this.pxHeight/2 / (atan(this.mmHeight/2/this.distanceToMonitor)*180/pi);
            
            %dont need to round dva
            dva =    pix/ horPixPerDva ;
            
            % TODO: improve
            
            % function pix = psyCortex_dva2pixExact( graph, poit1, point2 )
            %
            % distanceToCenter = min
            %
            % horPixPerDva = graph.pxWidth / atan(graph.mmWidth/2/graph.distanceToMonitor)*180/pi;
            % verPixPerDva = graph.pxHeight / atan(graph.mmHeight/2/graph.distanceToMonitor)*180/pi;
        end
        
        %% rotatePointCenter
        %--------------------------------------------------------------------------
        function [x y] = rotatePointCenter( graph, point, angle )
            
            [mx, my] = RectCenter(graph.wRect);
            
            
            p = rotatePoint( point, angle/180*pi, [mx my]);
            
            x = p(1);
            y = p(2);
        end
        
        %------------------------------------------------------------------
        %% Dialog Functions  ----------------------------------------------
        %------------------------------------------------------------------
        
        %% DlgHitKey
        function result = DlgHitKey( this, message, varargin )
            % DlgHitKey(window, message, [, sx][, sy][, color][, wrapat][, flipHorizontal][, flipVertical])
            
            
            if nargin < 1
                error('DlgHitKey: Must provide at least the first argument.');
            end
            
            if ( nargin < 2 || isempty(message) )
                message = 'Hit a key to continue';
            end
            
            if ( nargin < 3 || isempty(varargin{1}) )
                varargin{1} = 'center';
            end
            if ( nargin < 4 || isempty(varargin{2}) )
                varargin{2} = 'center';
            end
            
            oldDefaultColor = Screen( 'TextColor', this.window); % recover previous default color
            Screen( 'TextColor', this.window, this.dlgTextColor);
            % relevant keycodes
            ESCAPE = 27;
            
            % remove previous key presses
            FlushEvents('keyDown');
            
            DrawFormattedText( this.window, message, varargin{:} );
            Screen('Flip', this.window);
            cprintf('blue','\n---------------------------------------------------------\n');
            cprintf('blue','---------------------------------------------------------\n');
            cprintf('blue',[ message '\n']);
            cprintf('blue','---------------------------------------------------------\n');
            cprintf('blue','---------------------------------------------------------\n');
            
            KeyNotDown = 0;
            while(1)
                
                try
                    % g = ArumeHardware.GamePad();
                    % [ ~, ~, ~, a, b, x, y] = g.Query;
                    % 
                    % if ( a || b || x || y)
                    %     result = char('a');
                    %     break;
                    % end
                    
                    [~,~,buttons] = GetMouse();
                    
                    if buttons(2) % wait for release
                        
                        result = char('a');
                        break
                    end
                    
                catch
                end
                
                [keyIsDown, ~, keyCode, ~] = KbCheck();
                if ( keyIsDown )
                    if ( KeyNotDown )
                        keys = find(keyCode);
                        result = keyCode(keys(1));
                        break;
                    end
                else
                    % at some point it has to be not clicked before it is
                    % clicked
                    KeyNotDown = 1;
                end
                
                
            end
            
            %             char = GetChar;
            %             switch(char)
            %
            %                 case ESCAPE
            %                     result = 0;
            %
            %                 otherwise
            %                     result = char;
            %             end
            %
            
            Screen( 'TextColor', this.window, oldDefaultColor); % recover previous default color
        end
        
        
        %% DlgHitMouse
        function result = DlgHitMouse( this, message, varargin )
            % DlgHitMouse(window, message, [, sx][, sy][, color][, wrapat][, flipHorizontal][, flipVertical])
            
            
            if nargin < 1
                error('DlgHitMouse: Must provide at least the first argument.');
            end
            
            if ( nargin < 2 || isempty(message) )
                message = 'Click to continue';
            end
            
            if ( nargin < 3 || isempty(varargin{1}) )
                varargin{1} = 'center';
            end
            if ( nargin < 4 || isempty(varargin{2}) )
                varargin{2} = 'center';
            end
            
            oldDefaultColor = Screen( 'TextColor', this.window); % recover previous default color
            Screen( 'TextColor', this.window, this.dlgTextColor);
            % relevant keycodes
            ESCAPE = 27;
            
            % remove previous key presses
            FlushEvents('keyDown');
            
            DrawFormattedText( this.window, message, varargin{:} );
            Screen('Flip', this.window);
            buttons(1) = 0;
            
            while(~buttons(1))
                [x,y,buttons] = GetMouse;
                result = buttons(1);
            end
            
            Screen( 'TextColor', this.window, oldDefaultColor); % recover previous default color
        end
        
        %% DlgYesNo
        function result = DlgYesNo( this, message, yesText, noText, varargin )
            % DlgYesNo(window, message, yesText, noText, [, sx][, sy][, color][, wrapat][, flipHorizontal][, flipVertical])
            
            
            if nargin < 2
                error('DlgYesNo: Must provide at least the first two arguments.');
            end
            
            if ( nargin < 3 || isempty(yesText) )
                yesText = 'Yes';
            end
            
            if ( nargin < 4 || isempty(noText) )
                noText = 'No';
            end
            
            oldDefaultColor = Screen( 'TextColor', this.window); % recover previous default color
            Screen( 'TextColor', this.window, this.dlgTextColor);
            % possible results
            YES = 1;
            NO  = 0;
            
            % relevant keycodes
            ESCAPE  = 27;
            ENTER   = {13,3,10};
            
            % remove previous key presses
            FlushEvents('keyDown');
            
            DrawFormattedText(this.window, [message ' ' yesText ' (enter), ' noText ' (escape)'], varargin{:});
            Screen('Flip', this.window);
            
            cprintf('blue','\n---------------------------------------------------------\n');
            cprintf('blue','---------------------------------------------------------\n');
            cprintf('blue', [ message ' ' yesText ' (enter), ' noText ' (escape)\n']);
            cprintf('blue','---------------------------------------------------------\n');
            cprintf('blue','---------------------------------------------------------\n');
            
            
            while(1)
                char = GetChar;
                switch(char)
                    
                    case ENTER
                        result = YES;
                        break;
                        
                    case ESCAPE
                        result = NO;
                        break;
                end
            end
            
            Screen( 'TextColor', this.window, oldDefaultColor); % recover previous default color
        end
        
        
        %% DlgTimer
        function result = DlgTimer( this, message, maxTime, varargin )
            % DlgTimer(window, message [, maxTime][, sx][, sy][, color][, wrapat][, flipHorizontal][, flipVertical])
            
            result = 0;
            
            if nargin < 1
                error('DlgHitKey: Must provide at least the first argument.');
            end
            
            if ( nargin < 2 || isempty(message) )
                message = 'Hit a key to continue';
            end
            
            if ( nargin < 3 || isempty(maxTime) )
                maxTime = 90;
            end
            
            oldDefaultColor = Screen( 'TextColor', this.window); % recover previous default color
            Screen( 'TextColor', this.window, 255);
            
            % relevant keycodes
            ESCAPE = 27;
            
            % remove previous key presses
            FlushEvents('keyDown');
            
            tini = GetSecs;
            while(1)
                t = GetSecs-tini;
                DrawFormattedText( this.window, sprintf('%s - %4.1f seconds',message,maxTime-t), varargin{:} );
                
                % draw a fixation spot in the center;
                [mx, my] = RectCenter(this.wRect);
                fixRect = [0 0 10 10];
                fixRect = CenterRectOnPointd( fixRect, mx, my );
                Screen('FillOval', this.window,  255, fixRect);
                fliptime = Screen('Flip', this.window);
                
                Screen('Flip', this.window);
                
                if ( CharAvail )
                    char = GetChar;
                    switch(char)
                        
                        case ESCAPE
                            result = -1;
                            break;
                    end
                end
                if ( maxTime > 0 && (GetSecs-tini> maxTime ) )
                    break
                end
            end
            
            Screen( 'TextColor', this.window, oldDefaultColor); % recover previous default color
        end
        
        %% DlgSelect
        function result = DlgSelect( this, message, optionLetters, optionDescriptions, varargin )
            
            %DlgInput(window, message, optionLetters, optionDescriptions, tstring [, sx][, sy][, color][, wrapat][, flipHorizontal][, flipVertical])
            if nargin < 3
                error('DlgSelect: Must provide at least the first three arguments.');
            end
            
            if ( nargin < 4 || isempty(optionDescriptions) )
                optionDescriptions = optionLetters;
            end
            
            if ( length(optionLetters) ~= length(optionDescriptions) )
                error('DlgSelect: the number of options does not match the number of letters.');
            end
            
            oldDefaultColor = Screen( 'TextColor', this.window); % recover previous default color
            Screen( 'TextColor', this.window, this.dlgTextColor);
            
            ESCAPE = 27;
            ENTER   = {13,3,10};
            DOWN = 40;
            UP = 38;
            
            % remove previous key presses
            FlushEvents('keyDown');
            
            % draw options
            text = message;
            for i=1:length(optionLetters)
                text = [text '\n\n( ' optionLetters{i} ' ) ' optionDescriptions{i}];
            end
            
            selection = 0;
            DrawFormattedText( this.window, text, varargin{:} );
            Screen('Flip', this.window);
            
            cprintf('blue','\n---------------------------------------------------------\n');
            cprintf('blue','---------------------------------------------------------\n');
            cprintf('blue', [  text '\n']);
            cprintf('blue','---------------------------------------------------------\n');
            cprintf('blue','---------------------------------------------------------\n');
            
            while(1) % while no valid key is pressed
                
                c = GetChar;
                
                switch(c)
                    
                    case ESCAPE
                        result = 0;
                        break;
                        
                    case ENTER
                        if ( selection > 0 )
                            result = optionLetters{1};
                            break;
                        else
                            continue;
                        end
                    case {'a' 'z'}
                        if ( c=='a' )
                            selection = mod(selection-1-1,length(optionLetters))+1;
                        else
                            selection = mod(selection+1-1,length(optionLetters))+1;
                        end
                        text = message;
                        for i=1:length(optionLetters)
                            if ( i==selection )
                                text = [text '\n\n ->( ' optionLetters{i} ' ) ' optionDescriptions{i}];
                            else
                                text = [text '\n\n ( ' optionLetters{i} ' ) ' optionDescriptions{i}];
                            end
                        end
                        
                        DrawFormattedText( this.window, text, varargin{:} );
                        Screen('Flip', this.window);
                        
                    otherwise
                        if ( ~isempty( intersect( upper(optionLetters), upper( char(c) ) ) ) )
                            
                            result = optionLetters( streq( upper(optionLetters), upper( char(c) ) ) );
                            if ( iscell(result) )
                                result = result{1};
                            end
                            break;
                        end
                end
            end
            Screen( 'TextColor', this.window, oldDefaultColor); % recover previous default color
        end
        
        
        %% DlgInput
        function answer = DlgInput( this, message, varargin )
            
            %DlgInput(win, tstring [, sx][, sy][, color][, wrapat][, flipHorizontal][, flipVertical])
            
            if nargin < 1
                error('DlgHitKey: Must provide at least the first argument.');
            end
            
            if ( nargin < 2 || isempty(message) )
                message = '';
            end
            
            oldDefaultColor = Screen( 'TextColor', this.window); % recover previous default color
            Screen( 'TextColor', this.window, this.dlgTextColor);
            
            ESCAPE = 27;
            ENTER = {13,3,10};
            DELETE = 8;
            
            answer = '';
            
            FlushEvents('keyDown');
            
            while(1)
                text = [message ' ' answer ];
                
                DrawFormattedText( this.window, text, varargin{:} );
                Screen('Flip', this.window);
                
                
                char=GetChar;
                switch(abs(char))
                    
                    case ENTER,	% <return> or <enter>
                        break;
                        
                    case ESCAPE, % <scape>
                        answer  = '';
                        break;
                        
                    case DELETE,			% <delete>
                        if ~isempty(answer)
                            answer(end) = [];
                        end
                        
                    otherwise,
                        answer = [answer char];
                end
            end
            
            Screen( 'TextColor', this.window, oldDefaultColor); % recover previous default color
            
        end
        
    end
end


