function this = eyelinkSetupCustomCalib(this)

    % do some basic eyelink setup
    this.el=EyelinkInitDefaults(this.Graph.window);
    
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
    this.el.edfFile = sprintf('%s.edf',this.ExperimentOptions.ObserverID);
    Eyelink('openfile', this.el.edfFile);

end