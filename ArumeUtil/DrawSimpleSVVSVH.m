
clear screen
%           Screen('Preference', 'VisualDebugLevel', 3);
Screen('Preference', 'SkipSyncTests', 1);
Screen('Preference', 'VisualDebugLevel', 0);

%-- screens

graph.screens = Screen('Screens');
graph.selectedScreen=max(graph.screens);
graph.selectedScreen=1;

%-- window
Screen('Preference', 'ConserveVRAM', 64);
[graph.window, graph.wRect] = Screen('OpenWindow', graph.selectedScreen, 0, [], [], [], 0, 10);


whichTask = 'SVV';
angle = randn(1)*5;
targetPix = 25;
linelength = 200;



% draw a fixation spot in the center;
[mx, my] = RectCenter(graph.wRect);
fixRect = [0 0 10 10];
fixRect = CenterRectOnPointd( fixRect, mx, my );
Screen('FillOval', graph.window,  100, fixRect);

KbName('UnifyKeyNames')
while(1)
    if ( strcmp(whichTask,'SVV'))
        Screen('DrawLine', graph.window, 100, mx+linelength*sind(angle), my-linelength*cosd(angle), mx-linelength*sind(angle), my+linelength*cosd(angle), 2);
    else
        angleh = angle + 90;
        Screen('DrawLine', graph.window, 100, mx+linelength*sind(angleh), my-linelength*cosd(angleh), mx-linelength*sind(angleh), my+linelength*cosd(angleh), 2);
    end

    fixRect = [0 0 targetPix/2 targetPix/2];
    fixRect = CenterRectOnPointd( fixRect, mx, my);
    Screen('FillOval', graph.window, [255 255 255], fixRect);

    %-- Check for keyboard press
    [keyIsDown,secs,keyCode] = KbCheck;
    if keyCode(KbName('ESCAPE'))
        break;
    end
    if keyCode(KbName('RightArrow'))
        angle = angle + 0.1;
    end
    if keyCode(KbName('LeftArrow'))
        angle = angle - 0.1;
    end
    if keyCode(KbName('UpArrow'))
            whichTask = 'SVV';
%             angle = randn(1)*5;
    end
    if keyCode(KbName('DownArrow'))
        whichTask = 'SVH';
%         angle = randn(1)*5;
    end


    fl0iptime = Screen('Flip', graph.window);

end

angle

clear Screen