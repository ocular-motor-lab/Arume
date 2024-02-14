function [this, thisTrialData]  = getVisualSearchResponse(this, thisTrialData)

    noresp = true;
    startt = GetSecs;
    resp = nan;

    while noresp
    
        % we flip the heading so that the labels are more intuitive
        msg = sprintf( 'Were there any %s?',thisTrialData.SearchTarget);
        DrawFormattedText(this.Graph.window, msg, 'center', 'center', [255,255,255]);
        this.Graph.Flip(this, thisTrialData);
    
        % enter accepts the user's response
        [keyIsDown, ~, keyCode, ~] = KbCheck();
        if keyIsDown
            keys = find(keyCode);
            for i=1:length(keys)
                KbName(keys(i));
    
                switch(KbName(keys(i)))
                    case 'y'
                        resp = 1;
                        noresp = false;
    
                    case 'n'
                        resp = 0;
                        noresp = false;

                end
            end
        end
    
    end

    % only end trial once key is up
    while keyIsDown
        [keyIsDown, ~, ~, ~] = KbCheck();
    end


    thisTrialData.SearchResp = resp;
    thisTrialData.ResponseTime = GetSecs - startt;

    if this.ExperimentOptions.AuditoryFeedback

        if thisTrialData.SearchResp == thisTrialData.TargetPresent
            PsychPortAudio('Start', this.audio.pahandlecorrect, 1, 0, 1);
        else
            PsychPortAudio('Start', this.audio.pahandleincorrect, 1, 0, 1);
        end

    end

end