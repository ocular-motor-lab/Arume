function [this, thisTrialData]  = getHeadingResponse(this, thisTrialData)

    % random angle starting position
    respang = (rand(1))*pi;
    startt = GetSecs;
    
    noresp = true;
    
    while noresp
    
        % convert angle into x-y position
        resplinex2 = this.uicomponents.resplinex1+cos(respang)*this.uicomponents.linelen;
        respliney2 = this.uicomponents.respliney1-sin(respang)*this.uicomponents.linelen;
    
        Screen('DrawLines', this.Graph.window, this.uicomponents.outsidelinexy, 3, [127,127,127], [], 2);
        Screen('DrawLines', this.Graph.window, this.uicomponents.rulerlinexy, 3, [127,127,127], [], 2);
        Screen('DrawLines', this.Graph.window, [[this.uicomponents.resplinex1;this.uicomponents.respliney1]...
            ,[resplinex2;respliney2]], 3, [255,255,255], [], 2);
    
        % we flip the heading so that the labels are more intuitive
        DrawFormattedText(this.Graph.window, num2str(-round(rad2deg(respang)-90)), 'center', this.uicomponents.texty, [255,255,255]);
        this.Graph.Flip(this, thisTrialData);
    
        % enter accepts the user's response
        [keyIsDown, ~, keyCode, ~] = KbCheck();
        if keyIsDown
            keys = find(keyCode);
            for i=1:length(keys)
                KbName(keys(i));
    
                switch(KbName(keys(i)))
                    case 'y'
                        noresp = false;
    
                    case 'LeftArrow'
                        respang = min(respang+this.uicomponents.anginc,pi);

                    case 'RightArrow'
                        respang = max(respang-this.uicomponents.anginc,0);

                end
            end
        end
    
    end

    % only end trial once key is up
    while keyIsDown
        [keyIsDown, ~, ~, ~] = KbCheck();
    end

    % note that we do not negate the headings when saving responses, because
    % the raw angles and heading directions are already aligned (we negate them
    % only for showing the angle in text to observers). Negative = right,
    % positive = left
    thisTrialData.HeadingResponse = round(rad2deg(respang)-90);
    thisTrialData.ResponseTime = GetSecs - startt;

    if this.ExperimentOptions.AuditoryFeedback
        
        % set some sort of threshold error in degrees?
        if abs(thisTrialData.HeadingResponse - thisTrialData.HeadingChange) < 10
            PsychPortAudio('Start', this.audio.pahandlecorrect, 1, 0, 1);
        else
            PsychPortAudio('Start', this.audio.pahandleincorrect, 1, 0, 1);
        end

    end

end