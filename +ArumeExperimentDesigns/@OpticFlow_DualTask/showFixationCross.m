function showFixationCross(this, thisTrialData)

    % add a half-second trial buffer
    startt = GetSecs;
    while (GetSecs-startt) < .5
        this.Graph.Flip(this, thisTrialData);
    end

    % add a response prompt
    switch thisTrialData.Task
        case 'Visual Search'
            msg = sprintf( 'Search for the %s',thisTrialData.SearchTarget);
        case 'Heading'
            msg = 'Estimate the heading';
        case 'Both'
            msg = sprintf( 'Search for the %s, and estimate the heading',thisTrialData.SearchTarget);
    end

    startt = GetSecs;
    while (GetSecs-startt) < 2
        DrawFormattedText(this.Graph.window, msg, 'center', 'center', [255,255,255]);
        this.Graph.Flip(this, thisTrialData);
    end

    % check for response
    keyIsDown = false;
    while ~keyIsDown
        % show fixation until any button press
        Screen('FillRect', this.Graph.window, this.camera.fixcol, this.uicomponents.fixrects)
        this.Graph.Flip(this, thisTrialData);
        [keyIsDown, ~, ~, ~] = KbCheck();
    end

end