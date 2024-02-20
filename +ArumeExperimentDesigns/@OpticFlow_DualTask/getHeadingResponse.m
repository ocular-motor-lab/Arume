function [this, thisTrialData]  = getHeadingResponse(this, thisTrialData)

    % random angle starting position
    respang = (rand(1))*pi;
    startt = GetSecs;
    
    noresp = true;
    
    while noresp

        % show final frame of stimulus, but convert to grey (so they can't
        % continue searching). 
        Screen('FillOval', this.Graph.window, [100,100,100], this.shapes.lastframecirclerects)
        Screen('FillRect', this.Graph.window, [100,100,100], this.shapes.lastframesquarerects)
        if thisTrialData.TargetPresent
            Screen(this.shapes.targetdrawtype, this.Graph.window, [100,100,100], this.shapes.lastframetargetrects)
        end
    
        % convert angle into x-y position
        resplinex2 = this.uicomponents.resplinex1+cos(respang)*this.uicomponents.linelen;
        respliney2 = this.uicomponents.respliney1-sin(respang)*this.uicomponents.linelen;
    
        Screen('DrawLines', this.Graph.window, this.uicomponents.outsidelinexy, 3, [127,127,127], [], 2);
        Screen('DrawLines', this.Graph.window, this.uicomponents.rulerlinexy, 3, [127,127,127], [], 2);
        Screen('DrawLines', this.Graph.window, [[this.uicomponents.resplinex1;this.uicomponents.respliney1]...
            ,[resplinex2;respliney2]], 3, [255,255,255], [], 2);
    
        % we flip the heading so that the labels are more intuitive
        DrawFormattedText(this.Graph.window, sprintf('Your Answer: %i',-round(rad2deg(respang)-90)), 'center', this.uicomponents.texty, [255,255,255]);
        this.Graph.Flip(this, thisTrialData);
    
        % enter accepts the user's response
        [keyIsDown, ~, keyCode, ~] = KbCheck();
        if keyIsDown
            keys = find(keyCode);
            for i=1:length(keys)
                KbName(keys(i));
    
                switch(KbName(keys(i)))
                    case 'Return'
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
        % if abs(thisTrialData.HeadingResponse - thisTrialData.HeadingChange) < 10
        %     PsychPortAudio('Start', this.audio.pahandlecorrect, 1, 0, 1);
        % else
        %     PsychPortAudio('Start', this.audio.pahandleincorrect, 1, 0, 1);
        % end

        % equivalent of auditory feedback is showing what the ground-truth
        % angle was
        gtang = deg2rad(thisTrialData.HeadingChange+90);

        resplinex2 = this.uicomponents.resplinex1+cos(respang)*this.uicomponents.linelen;
        respliney2 = this.uicomponents.respliney1-sin(respang)*this.uicomponents.linelen;

        resplinex3 = this.uicomponents.resplinex1+cos(gtang)*this.uicomponents.linelen;
        respliney3 = this.uicomponents.respliney1-sin(gtang)*this.uicomponents.linelen;
        
        currt = GetSecs;
        while (GetSecs-currt)<2

            Screen('DrawLines', this.Graph.window, this.uicomponents.outsidelinexy, 3, [127,127,127], [], 2);
            Screen('DrawLines', this.Graph.window, this.uicomponents.rulerlinexy, 3, [127,127,127], [], 2);
            Screen('DrawLines', this.Graph.window, [[this.uicomponents.resplinex1;this.uicomponents.respliney1]...
                ,[resplinex2;respliney2]], 3, [255,255,255], [], 2);
    
            Screen('DrawLines', this.Graph.window, [[this.uicomponents.resplinex1;this.uicomponents.respliney1]...
                ,[resplinex3;respliney3]], 3, [50,255,50], [], 2);
        
            % we flip the heading so that the labels are more intuitive
            DrawFormattedText(this.Graph.window, sprintf('Your Answer: %i',-round(rad2deg(respang)-90)), 'center', this.uicomponents.texty, [255,255,255]);
            DrawFormattedText(this.Graph.window, sprintf('Correct Answer: %i',-round(rad2deg(gtang)-90)), 'center', this.uicomponents.texty+30, [50,255,50]);
            this.Graph.Flip(this, thisTrialData);
        end


    end

end