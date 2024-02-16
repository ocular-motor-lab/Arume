function endOfTrialSequence(this,thisTrialData)

    % just exit if end of experiment
    if thisTrialData.Trial ~= this.ExperimentOptions.nTrialsTotal
        
        % end of block = break
        if thisTrialData.BlockChange

            msg = 'Switching Task. then press any key to continue...';
            noresp = true;

            while noresp

                DrawFormattedText(this.Graph.window,msg,'center','center',[255,255,255])
                this.Graph.Flip(this, thisTrialData);

                % enter accepts the user's response
                [keyIsDown, ~, ~, ~] = KbCheck();
                if keyIsDown
                    % just start next block
                    noresp = false;
                end
            end

        end
    end
end