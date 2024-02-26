     



  




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




     % draw a fixation spot in the center;
            [mx, my] = RectCenter(graph.wRect);
            fixRect = [0 0 10 10];
            fixRect = CenterRectOnPointd( fixRect, mx, my );
            Screen('FillOval', graph.window,  100, fixRect);
            
            Screen('DrawLine', graph.window, 100, mx, my-2500, mx, my+2500, 1);
            Screen('DrawLine', graph.window, 100, mx-2500, my, mx+2500, my, 1);

            fl0iptime = Screen('Flip', graph.window);
            
            pause

            clear Screen