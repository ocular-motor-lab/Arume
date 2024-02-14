classdef EyeTrackerEyelink  < handle
    %VOG Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        graph
        el
        
    end
    
    methods
        function result = Connect(this, graph, ip, port, openirispath)

            result = 0;

            % do some basic eyelink setup
            this.el = EyelinkInitDefaults(graph.window);

            % change the calibration point properties
            this.el.calibrationtargetsize = 1;
            this.el.calibrationtargetwidth = .25;
            this.el.backgroundcolour = 0;
            this.el.imgtitlecolour = 1;
            this.el.foregroundcolour = 1;
            this.el.msgfontcolour = 255;
            this.el.calibrationtargetcolour = [255,255,255];
            EyelinkUpdateDefaults(this.el);

            % Initialization of the connection with the Eyelink Gazetracker.
            % exit program if this fails.
            if ~EyelinkInit()
                fprintf('Eyelink Init aborted.\n');
                Eyelink('Shutdown');  % cleanup function
                return;
            end

            % set up some basic calibration stuff after we initialize the eyelink
            % (we cannot do this before EyelinkInit)
            Eyelink('Command','calibration_area_proportion = 0.4 0.4'); % should be about 32 deg across!!!
            Eyelink('Command','validation_area_proportion = 0.4 0.4');

            % make sure that we get gaze data from the Eyelink
            Eyelink('command', 'link_sample_data = LEFT,RIGHT,GAZE,AREA,INPUT');
            Eyelink('command', 'file_sample_data = LEFT,RIGHT,GAZE,AREA,INPUT');

            % we can also extract event data if we like
            Eyelink('command', 'link_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,FIXUPDATE,INPUT');
            Eyelink('command', 'file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,FIXUPDATE,INPUT');

            % open file to record data to
            this.el.edfFile = sprintf('ArumeTmp.edf'); % TODO: maybe change this
            Eyelink('openfile', this.el.edfFile);


            result = 1;
        end
        
        function result = IsRecording(this)
            result = 0;
            if ( ~isempty( this.el) )
                error=Eyelink('CheckRecording');
                if(error==0)
                    result = 1;
                end
            end
        end
        
        function SetSessionName(this, sessionName)
            if ( ~isempty( this.el) )
            end
        end
        
        function StartRecording(this)
            if ( ~isempty( this.el) )
                Eyelink('StartRecording');
            end
        end
        
        function StopRecording(this)
            if ( ~isempty( this.el) )
                Eyelink('StopRecording');
            end
        end
        
        function [frameNumber, timestamp] = RecordEvent(this, message)
            frameNumber = nan;
            timestamp = nan;
            if ( ~isempty( this.el) )
                timestamp=EyelinkGetTime(this.el); % [, maxwait]) % TODO: this will be a timestamp not a frame number
                Eyelink('Message',sprintf('ELtime=%d PTBtime=%d    %s',timestamp, GetSecs, message))
                timestamp = timestamp/1000; % convert to seconds for Arume
            end
        end
        
        function evt = GetCurrentData(this, message)
            % data =[];
            if ( ~isempty( this.el) )

                eye_used = Eyelink('EyeAvailable'); % get eye that's tracked
                if eye_used == this.el.BINOCULAR % if both eyes are tracked
                    eye_used = this.el.LEFT_EYE; % use left eye?
                end
                
                % get all gaze pos and pupil data 
                if Eyelink('NewFloatSampleAvailable') > 0

                    % get the sample in the form of an event structure
                    evt = Eyelink('NewestFloatSample');

                    if eye_used ~= -1 % do we know which eye to use yet?

                        % if we do, get current gaze position from sample
                        x = evt.gx(eye_used+1); % +1 as we're accessing MATLAB array
                        y = evt.gy(eye_used+1);
                        % do we have valid data and is the pupil visible?
                        if x~=this.el.MISSING_DATA && y~=this.el.MISSING_DATA && evt.pa(eye_used+1) > 0
                            evt.mx=x;
                            evt.my=y;
                        end
                    end
                end

                if exist('message','var')
                    Eyelink('Message',message);
                end
                % data = this.el.GetCurrentData();
            end
        end

        function calibrationSuccessful = Calibrate(this)
            result = EyelinkDoTrackerSetup(this.el);

            if ( result == 0 )
                calibrationSuccessful = 1;
            end
        end

        function Disconnect(this)

            Eyelink('Shutdown');
        end
        
        function [files]= DownloadFile(this, path, newFileName)
            files = {};
            if ( isempty( this.el) )
                return;
            end
            if (~exist('path','var'))
                path = '';
            end

            if (~exist('newFileName','var'))
                newFileName = this.el.edfFile;
            end

            % download data file
            try
                % close EL file and copy file over to local directory
                Eyelink('CloseFile');

                fprintf('Receiving data file ''%s''\n', newFileName);

                % do not overwrite existing file
                [~,name,ext] = fileparts(newFileName);
                ctr = 1;
                while exist(fullfile(path, newFileName),'file')
                    ctr = ctr+1;
                    newFileName = [name, sprintf('%02d',ctr), ext];
                end

                % send that data over boi!
                status=Eyelink('ReceiveFile',this.el.edfFile, ...
                    convertStringsToChars(fullfile(path, newFileName)));
                if status > 0
                    fprintf('ReceiveFile status %d\n', status);
                end
                
                files = {newFileName};

            catch ex
                getReport(ex)
                cprintf('red', sprintf('++ EYELINK :: Problem receiving data file ''%s''\n', this.el.edfFile));
                files = {};
            end
        end
    end
    
    methods(Static = true)
        
    end
    
end



