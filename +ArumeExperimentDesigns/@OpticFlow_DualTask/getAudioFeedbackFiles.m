function this = getAudioFeedbackFiles(this)

    InitializePsychSound

    % For some reason AudioRead doesn't like the relative path, so we find
    % the absolute path
    currpath = mfilename('fullpath');
    [filepath,~,~] = fileparts(currpath);
    rootdir = [filepath,'/responseSounds/'];
    correctsound = [rootdir,'success-1-6297.mp3'];
    incorrectsound = [rootdir,'failure-drum-sound-effect-2-7184.mp3'];
    % repetitions = 1;
    correctvol = .25;

    [ycorrect,freqcorrect] = audioread(correctsound);
    wavedatacorrect = (ycorrect').*correctvol;
    nrchannelscorrect = size(wavedatacorrect,1); % Number of rows == number of channels.

    [yincorrect,freqincorrect] = audioread(incorrectsound);
    wavedataincorrect = yincorrect';
    nrchannelsincorrect = size(wavedataincorrect,1); % Number of rows == number of channels.

    % Make sure we have always 2 channels stereo output.
    % Why? Because some low-end and embedded soundcards
    % only support 2 channels, not 1 channel, and we want
    % to be robust in our demos.
    if nrchannelscorrect < 2
        wavedatacorrect = [wavedatacorrect; wavedatacorrect];
        nrchannelscorrect = 2;
    end

    if nrchannelsincorrect < 2
        wavedataincorrect = [wavedataincorrect; wavedataincorrect];
        nrchannelsincorrect = 2;
    end

    this.audio.wavedatacorrect = wavedatacorrect;
    this.audio.wavedataincorrect = wavedataincorrect;

    device = [];

    % Open the  audio device, with default mode [] (==Only playback),
    % and a required latencyclass of zero 0 == no low-latency mode, as well as
    % a frequency of freq and nrchannels sound channels.
    % This returns a handle to the audio device:
    try
        % Try with the 'freq'uency we wanted:
        this.audio.pahandlecorrect = PsychPortAudio('Open', device, [], 0, freqcorrect, nrchannelscorrect);
    catch
        % Failed. Retry with default frequency as suggested by device:
        fprintf('\nCould not open device at wanted playback frequency of %i Hz. Will retry with device default frequency.\n', freq);
        fprintf('Sound may sound a bit out of tune, ...\n\n');
    
        psychlasterror('reset');
        this.audio.pahandlecorrect = PsychPortAudio('Open', device, [], 0, [], nrchannelscorrect);
    end

    try
        this.audio.pahandleincorrect = PsychPortAudio('Open', device, [], 0, freqincorrect, nrchannelsincorrect);
    catch
        fprintf('\nCould not open device at wanted playback frequency of %i Hz. Will retry with device default frequency.\n', freq);
        fprintf('Sound may sound a bit out of tune, ...\n\n');
    
        psychlasterror('reset');
        this.audio.pahandleincorrect = PsychPortAudio('Open', device, [], 0, [], nrchannelsincorrect);
    end

    % Fill the audio playback buffer with the audio data 'wavedata':
    PsychPortAudio('FillBuffer', this.audio.pahandlecorrect, this.audio.wavedatacorrect);
    PsychPortAudio('FillBuffer', this.audio.pahandleincorrect, this.audio.wavedataincorrect);

end