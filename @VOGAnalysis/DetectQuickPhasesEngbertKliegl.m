function [data, info] = DetectQuickPhasesEngbertKliegl(data, params)

info = [];

if ( ischar(data) )
    command = data;
    switch( command)
        case 'get_options'
            optionsDlg = [];

            optionsDlg.Components = {{'H','{HV}','HVT'}};
            optionsDlg.MINDUR = 3;
            optionsDlg.VFAC = 6;
            optionsDlg.Remove_Monoculars = { {'0','{1}'} };
            optionsDlg.Binocular_Minimum_Overlap = { {'0','{1}'} };
            optionsDlg.Remove_Overshoots = { {'0','{1}'} };
            optionsDlg.Overshoot_Interval = 30;
            optionsDlg.Refine_Begining_And_End  = { {'{0}','1'} };
            optionsDlg.Recover_Monoculars       = { {'0','{1}'} };
            optionsDlg.Recover_Monoculars_Threshold = 10;

            data = optionsDlg;
            return;
        case 'get_defaults'
            optionsDlg = VOGAnalysis.DetectQuickPhasesEngbertKliegl('get_options');
            defaultOptions = StructDlg(optionsDlg,'',[],[],'off');
            data = defaultOptions;
            return;
    end
end

if ( ~exist('params','var') )
    params.Detection.Engbert = VOGAnalysis.DetectQuickPhasesEngbertKliegl('get_defaults');
end


samplerate = data.Properties.UserData.sampleRate;

eyes = {'Left' 'Right'};
eyeSignals = {'X','Y', 'T'};
LEFT = true;
RIGHT = true;

if ( isfield( data.Properties.UserData, 'Eyes') )
    eyes = data.Properties.UserData.Eyes;
end

if ( isfield( data.Properties.UserData, 'EyeSignals') )
    eyeSignals = data.Properties.UserData.EyeSignals;
end
if ( isfield( data.Properties.UserData, 'LEFT') )
    LEFT = data.Properties.UserData.LEFT;
end
if ( isfield( data.Properties.UserData, 'RIGHT') )
    RIGHT = data.Properties.UserData.RIGHT;
end

%% FIND SACCADES in each component
textprogressbar('++ DetectQuickPhaseEngbertKliegl :: Detecting quick phases using E&K method modified: ');
Nprogsteps = 5/100;
tic

try

    switch(params.Detection.Engbert.Components)
        case 'H'
            componetsToUse = {'X'};
        case 'HV'
            componetsToUse = {'X','Y'};
        case 'HVT'
            componetsToUse = {'X','Y', 'T'};
    end

    % get the velocity for each of the componets
    for k=1:length(eyes)
        for j=1:length(eyeSignals)
            data.([eyes{k} 'Vel' eyeSignals{j}]) = engbert_vecvel(data.([eyes{k} eyeSignals{j}]), samplerate, 2);
            data.([eyes{k} 'Accel' eyeSignals{j}]) = engbert_vecvel(data.([eyes{k} 'Vel' eyeSignals{j}]), samplerate, 2);
            data.([eyes{k} 'Jerk' eyeSignals{j}]) = engbert_vecvel(data.([eyes{k} 'Accel' eyeSignals{j}]), samplerate, 2);
        end
    end

    textprogressbar(1/Nprogsteps);

    lrsac = [];
    xy = cell(size(eyes));
    radEyes = {};
    for k=1:length(eyes)

        xy{k} = data{:,{[eyes{k} 'X'], [eyes{k} 'Y'], [eyes{k} 'T']}};
        v = data{:,{[eyes{k} 'VelX'], [eyes{k} 'VelY'], [eyes{k} 'VelT']}};

        isInTrial = ones(size(data.Time)); % TODO!!!
        blinkYesNo = isnan(v(:,1)) | isnan(v(:,2));

        [sac, rad] = engbert_microsacc(xy{k}, v ,params.Detection.Engbert.VFAC,params.Detection.Engbert.MINDUR, blinkYesNo, componetsToUse);
        %--------------------------------------------------------------------
        % OUTPUT
        %  sac(1:num,1)   onset of saccade
        %  sac(1:num,2)   end of saccade
        %  sac(1:num,3)   peak velocity of saccade (vpeak)
        %  sac(1:num,4)   horizontal component     (dx)
        %  sac(1:num,5)   vertical component       (dy)
        %  sac(1:num,6)   horizontal amplitude     (dX)
        %  sac(1:num,7)   vertical amplitude       (dY)
        %  sac(1:num,9)   peak test variable       (dY)


        %% -- Filter bad microsaccades --------------------------------------------
        %-- remove usacc in intertrials and in blinks
        good_samples = isInTrial & ~ blinkYesNo;

        if ( ~isempty( sac ) )
            is_good = zeros(size(sac,1),1);
            for j=1:length(is_good)
                is_good(j) = sum( ~good_samples( sac(j,1):sac(j,2)) ) == 0;
            end
            sac = sac( is_good==1, :);
        end

        lrsac.(eyes{k}) = sac;
        radEyes{k} = rad;
    end

    textprogressbar(2/Nprogsteps);

    if ( (LEFT && ~isempty(lrsac.Left)) || (RIGHT && ~isempty(lrsac.Right)) ) % some saccades detected


        %% - Remove monoculars (overlap threshold)
        if ( LEFT && RIGHT && ~isempty(lrsac.Left) && ~isempty(lrsac.Right) && params.Detection.Engbert.Remove_Monoculars )
            [left_monoculars_idx, right_monoculars_idx] = FindMonoculars( lrsac.Left, lrsac.Right, params.Detection.Engbert.Binocular_Minimum_Overlap, params.Detection.Engbert.Recover_Monoculars, params.Detection.Engbert.Recover_Monoculars_Threshold );
            lrsac.Left(left_monoculars_idx,:) = [];
            lrsac.Right(right_monoculars_idx,:)  = [];
        end

        textprogressbar(3/Nprogsteps);

        %% - Remove overshoots
        if ( params.Detection.Engbert.Remove_Overshoots )
            if ( LEFT && RIGHT )
                if ( params.Detection.Engbert.Remove_Monoculars && (length(lrsac.Left) == length(lrsac.Right) && ~params.Detection.Engbert.Recover_Monoculars) )
                    % if monoculars have been removed we can find overshoots
                    % binocularly
                    [overshoots_idx] = FindOvershoots( xy{1}, lrsac.Left, xy{2}, lrsac.Right, samplerate,  params.Detection.Engbert.Overshoot_Interval );
                    left_overshoots_idx = overshoots_idx;
                    right_overshoots_idx = overshoots_idx;
                else
                    % if monoculars are still present we find the overshoots
                    % independently in each eye
                    left_overshoots_idx = FindOvershoots( xy{1}, lrsac.Left, [], [], samplerate, params.Detection.Engbert.Overshoot_Interval );
                    right_overshoots_idx = FindOvershoots( [], [], xy{2}, lrsac.Right, samplerate, params.Detection.Engbert.Overshoot_Interval );
                end
                lrsac.Left(left_overshoots_idx,:) = [];
                lrsac.Right(right_overshoots_idx,:) = [];
            elseif ( LEFT )
                [left_overshoots_idx] = FindOvershoots( xy{1}, lrsac.Left, [], [], samplerate, params.Detection.Engbert.Overshoot_Interval );
                lrsac.Left(left_overshoots_idx,:) = [];
            elseif ( RIGHT )
                [right_overshoots_idx] = FindOvershoots( [], [], xy{1}, lrsac.Right, samplerate, params.Detection.Engbert.Overshoot_Interval );
                lrsac.Right(right_overshoots_idx,:) = [];
            end
        end

        textprogressbar(4/Nprogsteps);

        % up to this point because we may have tried to recover monoculars, it is
        % possible that the left and right eyes don't have the same number of
        % saccades
        %
        % here we will condense them into a single binocular sac variable that has
        % the begining and ends according to the single eye begining and end
        if ( LEFT && RIGHT  && (length(lrsac.Left) == length(lrsac.Right) && ~params.Detection.Engbert.Recover_Monoculars) )
            % if we truly have matching binocular saccades in the left and right
            % eye we consider the begining the first begining between the two eyes
            % and the end the last ending between the two eyes
            sac = [min([lrsac.Left(:,1) lrsac.Right(:,1)],[],2) max([lrsac.Left(:,2) lrsac.Right(:,2)],[],2)];
        elseif ( LEFT && RIGHT)
            % if we are not the same, we can build a yes/no vector for all the
            % samples for each eye to see if they belong to a saccade or not and
            % then do the logical OR between the two eyes and then get the begining
            % and end of whatever periods with "yes" we get
            lu = zeros(size(data.Time,1)+1,1);
            ru = zeros(size(data.Time,1)+1,1);
            lu(lrsac.Left(:,1)) = 1;
            lu(lrsac.Left(:,2)+1) = -1;
            lu = cumsum(lu); % yes/no saccade in the left eye
            ru(lrsac.Right(:,1)) = 1;
            ru(lrsac.Right(:,2)+1) = -1;
            ru = cumsum(ru); % yes/no saccade in the right eye
            sac = [];
            u = double(lu(1:end-1) | ru(1:end-1)); % yes/no left OR right eye
            sac(:,1) = find(diff([0;u])>0);
            sac(:,2) = find(diff([0;u])<0);
        elseif (LEFT)
            sac = lrsac.Left;
        elseif (RIGHT)
            sac = lrsac.Right;
        end


        textprogressbar(5/Nprogsteps);


    end

catch ME
    timeElapsed = toc;

    textprogressbar(sprintf('NOT DONE DUE TO ERROR  after %0.2f seconds.', timeElapsed));
    rethrow(ME)
end
timeElapsed = toc;

textprogressbar(sprintf('Done in %0.2f seconds.', timeElapsed));

if ( (LEFT && ~isempty(lrsac.Left)) || (RIGHT && ~isempty(lrsac.Right)) ) % some saccades detected


    %% refine beginings and ends

    if ( params.Detection.Engbert.Refine_Begining_And_End )
        textprogressbar('++ DetectQuickPhaseEngbertKliegl :: Adjusting begining and end of quick phases: ');
        Nprogsteps = size(sac,1)/100;
        tic

        newqp = sac;
        qp =sac;
        for i=1:size(qp,1)
            if ( mod(i,100) == 0 )
                textprogressbar(i/Nprogsteps);
            end
            %%
            idx = qp(i,1):qp(i,2);
            idx2 = max(1,(qp(i,1)-20)):min(height(data),(qp(i,2)+20)); % indices of the window of data in the table
            idx = idx-idx2(1)+1; % indices of the qp inside the window of data

            idxmaxqp = nan(length(eyes), length(eyeSignals));
            accelPeakIdx = nan(length(eyes), length(eyeSignals));
            breakPeakIdx = nan(length(eyes), length(eyeSignals));
            start = nan(length(eyes), length(eyeSignals));
            finish = nan(length(eyes), length(eyeSignals));

            for k=1:length(eyes)
                for j=1:length(eyeSignals)

                    vel = data{idx2,[eyes{k} 'Vel' eyeSignals{j}]};
                    accel = data{idx2,[eyes{k} 'Accel' eyeSignals{j}]};
                    jerk = data{idx2,[eyes{k} 'Jerk' eyeSignals{j}]};

                    % find velocity peak
                    [~,idxmaxqp(k,j)] = max(abs(vel(idx)));
                    idxmaxqp(k,j) = idxmaxqp(k,j) + idx(1) - 1; % index of the peak within the window

                    [t1, t2, tacc, tbreak] = VOGAnalysis.FindBeginEnd(vel, accel, jerk, idxmaxqp(k,j));
                    if ( isnan(t1) && isnan(t2) )
                        continue;
                    end

                    start(k,j) = t1;
                    finish(k,j) = t2;
                    accelPeakIdx(k,j) = tacc;
                    breakPeakIdx(k,j) = tbreak;
                end
            end

            newqp(i,1) = nanmin(start(:)) + idx2(1) - 1;
            newqp(i,2) = nanmax(finish(:)) + idx2(1) - 1;
            %%
            if (0)
                %%
                f =figure
                comps = {'', 'Vel','Accel', 'Jerk'};
                for j=1:numel(eyeSignals)
                    for icomp=1:numel(comps)
                        nplot = (icomp-1)*numel(eyeSignals) + j;
                        subplot(numel(comps),numel(eyeSignals),nplot,'nextplot','add')
                        for k=1:length(eyes)
                            x = data{idx2,[eyes{k} comps{icomp} eyeSignals{j}]};
                            peakVel = nan(size(x));
                            peakVel(idxmaxqp(k,j)) = x(idxmaxqp(k,j));
                            peak1Accel = nan(size(x));
                            peak1Accel(accelPeakIdx(k,j)) = x(accelPeakIdx(k,j));
                            peak1Accel(breakPeakIdx(k,j)) = x(breakPeakIdx(k,j));
                            qpIdx = start(k,j):finish(k,j);
                            qpBin = nan(size(x));
                            qpBin(qpIdx) = x(qpIdx);

                            qpIdx2 = [newqp(i,1) newqp(i,1) nan newqp(i,2) newqp(i,2)]-(+ idx2(1) - 1);
                            qpBin2 = [min(x) max(x) nan min(x) max(x)];

                            plot(qpBin, 'g','linewidth',3);
                            plot(x)
                            plot(qpIdx2, qpBin2, 'r','linewidth',2);
                            plot(peakVel,'rv');
                            plot(peak1Accel, 'ks');
                            if ( icomp > 1)
                                line(get(gca,'xlim'),[0 0])
                            end
                            set(gca,'xlim',[0 80])
                        end
                    end
                end
                %%
                %             pause;
                close(f);
            end
        end
        sac = newqp;

        timeElapsed = toc;
        textprogressbar(sprintf('Done in %0.2f seconds.', timeElapsed));
    end



    % build the zero and one vector

    l = height(data);
    starts = sac(:,1);
    stops = sac(:,2);
    yesNo = zeros(l,1);
    [us ius] = unique(starts);
    yesNo(us) = yesNo(us)+ diff([ius;length(starts)+1]);
    [us ius] = unique(stops);
    yesNo(us) = yesNo(us)- diff([ius;length(stops)+1]);
    yesNo = cumsum(yesNo)>0;

    data.QuickPhase = yesNo;
else

    data.QuickPhase = zeros(height(data),1);
end

info.Engbert_Vth = rad;

for k=1:length(eyes)
    cprintf('blue', sprintf('++ VOGAnalysis :: Engbert thresholds %s (H V T deg/s) %0.1f %0.1f %0.1f\n', ...
        eyes{k}, radEyes{k}(1), radEyes{k}(2), radEyes{k}(3) ));
end



end

function v = engbert_vecvel(xx,SAMPLING,TYPE)
%------------------------------------------------------------
%
%  VELOCITY MEASURES
%  - EyeLink documentation, p. 345-361
%  - Engbert, R. & Kliegl, R. (2003) Binocular coordination in
%    microsaccades. In:  J. Hyönä, R. Radach & H. Deubel (eds.)
%    The Mind's Eyes: Cognitive and Applied Aspects of Eye Movements.
%    (Elsevier, Oxford, pp. 103-117)
%
%  (Version 1.2, 01 JUL 05)
%-------------------------------------------------------------
%
%  INPUT:
%
%  xy(1:N,1:2)     raw data, x- and y-components of the time series
%  SAMPLING        sampling rate
%
%  OUTPUT:
%
%  v(1:N,1:2)     velocity, x- and y-components
%
%-------------------------------------------------------------
N = length(xx(:,1));            % length of the time series
M = length(xx(1,:));
v = zeros(N,M);

if ( SAMPLING < 1000 )
    switch TYPE
        case 1
            v(2:N-1,:) = SAMPLING/2*[xx(3:end,:) - xx(1:end-2,:)];
        case 2
            v(3:N-2,:)	= SAMPLING/6 * [xx(5:end,:) + xx(4:end-1,:) - xx(2:end-3,:) - xx(1:end-4,:)];
            v(2,:)		= SAMPLING/2 * [xx(3,:) - xx(1,:)];
            v(N-1,:)	= SAMPLING/2 * [xx(end,:) - xx(end-2,:)];
    end
else

    v(9:end-8,:)	= SAMPLING/24 * [xx(17:end,:) + xx(13:end-4,:) - xx(5:end-12,:) - xx(1:end-16,:)];

end
end

function [sac, radius] = engbert_microsacc(x,vel,VFAC,MINDUR, blink, comps)
%-------------------------------------------------------------------
%
%  FUNCTION microsacc.m
%
%  (Version 1.0, 22 FEB 01)
%  (Version 2.0, 18 JUL 05)
%  (Version 2.1, 03 OCT 05)
%
%-------------------------------------------------------------------
%
%  INPUT:
%
%  x(:,1:2)         position vector
%  vel(:,1:2)       velocity vector
%  VFAC             relative velocity threshold
%  MINDUR           minimal saccade duration
%
%  OUTPUT:
%
%  sac(1:num,1)   onset of saccade
%  sac(1:num,2)   end of saccade
%  sac(1:num,3)   peak velocity of saccade (vpeak)
%  sac(1:num,4)   horizontal component     (dx)
%  sac(1:num,5)   vertical component       (dy)
%  sac(1:num,6)   horizontal amplitude     (dX)
%  sac(1:num,7)   vertical amplitude       (dY)
%  ADDED BY JORGE
%  sac(1:num,8)   mean velocity of the saccade
%
%---------------------------------------------------------------------
vel2 = vel(~blink,:);warning off MATLAB:divideByZero
%vel2 =vel;

% compute threshold
medx = median(vel2(:,1),'omitnan');
msdx = sqrt( median((vel2(:,1)-medx).^2,'omitnan') );
medy = median(vel2(:,2),'omitnan');
msdy = sqrt( median((vel2(:,2)-medy).^2,'omitnan') );
medt = median(vel2(:,3),'omitnan');
msdt = sqrt( median((vel2(:,3)-medt).^2,'omitnan') );

if msdx<realmin
    msdx = sqrt( mean(vel2(:,1).^2) - (mean(vel2(:,1)))^2 );
end
if msdy<realmin
    msdy = sqrt( mean(vel2(:,2).^2) - (mean(vel2(:,2)))^2 );
end
if msdt<realmin
    msdt = sqrt( mean(vel2(:,3).^2) - (mean(vel2(:,3)))^2 );
end
radiusx = VFAC*msdx;
radiusy = VFAC*msdy;
radiust = VFAC*msdt;
radius	= [radiusx radiusy radiust];

% compute test criterion: ellipse equation
if ( radiusy == 0 || length(comps) == 1)
    test = (vel(:,1)/radiusx).^2 ;
elseif ( radiusx == 0 )
    test =  (vel(:,2)/radiusy).^2;
elseif (length(comps) == 2)
    test = (vel(:,1)/radiusx).^2 + (vel(:,2)/radiusy).^2;
elseif (length(comps) == 3)
    test = (vel(:,1)/radiusx).^2 + (vel(:,2)/radiusy).^2 + (vel(:,3)/radiust).^2;
end

indx = find(test>1 & ~isnan(test));

% determine saccades
N = length(indx);
sac = zeros(N, 8);
nsac = 0;
dur = 1;
a = 1;
k = 1;
while k<N
    if ( indx(k+1)-indx(k) )==1
        dur = dur + 1;
    else
        if dur>=MINDUR
            nsac = nsac + 1;
            b = k;
            sac(nsac,1:2) = [indx(a) indx(b)];
        end
        a = k+1;
        dur = 1;
    end
    k = k + 1;
end

% check for minimum duration
if dur>=MINDUR
    nsac = nsac + 1;
    b = k;
    sac(nsac,1:2) = [indx(a) indx(b)];
end
sac = sac(1:nsac, :);

% compute peak velocity, horizonal and vertical components
for s=1:nsac
    % onset and offset
    a = sac(s,1);
    b = sac(s,2);
    % saccade peak velocity (vpeak)
    vpeak = max( sqrt( vel(a:b,1).^2 + vel(a:b,2).^2 ) );
    sac(s,3) = vpeak;
    % saccade vector (dx,dy)
    dx = x(b,1)-x(a,1);
    dy = x(b,2)-x(a,2);
    sac(s,4) = dx;
    sac(s,5) = dy;
    % saccade amplitude (dX,dY)
    i = sac(s,1):sac(s,2);
    [minx, ix1] = min(x(i,1));
    [maxx, ix2] = max(x(i,1));
    [miny, iy1] = min(x(i,2));
    [maxy, iy2] = max(x(i,2));
    dX = sign(ix2-ix1)*(maxx-minx);
    dY = sign(iy2-iy1)*(maxy-miny);
    sac(s,6:7) = [dX dY];

    sac(s, 8) = mean( sqrt( vel(a:b,1).^2 + vel(a:b,2).^2 ) );
    testsac = sort(test(a:b),'descend');
    sac(s, 9) = testsac(3);
end
end

function [left_monoculars_idx, right_monoculars_idx] = FindMonoculars( lsac, rsac, Minimum_Overlap, Recover_Monoculars, Recover_Monoculars_Threshold)
ls = lsac(:,1);
le = lsac(:,2);
rs = rsac(:,1);
re = rsac(:,2);

ltest = lsac(:,9);
rtest = rsac(:,9);

binocs2 = [];
% for each microsaccade in the left eye
for i=1:length(ls)
    e = ones(size(re))*le(i);
    s = ones(size(rs))*ls(i);
    row = max(min( re, e) - max(rs, s),0);
    % row has the overlap of this left microsaccade with all the right
    % microsaccades
    if( sum(row) > 0 )
        % we get the right microsaccade that overlaps more with the
        % left one and we check that the left one is also the one that
        % overlaps more with the right one.
        [~, max_index_r] = max( row );
        e = ones(size(le))*re(max_index_r);
        s = ones(size(ls))*rs(max_index_r);
        col = max( min( le, e) - max(ls, s), 0 );
        [~, max_index_l] = max( col );
        if (max_index_l == i )
            binocs2(end+1, : ) = [max_index_l;max_index_r];
        end
    end
end
left_monoculars_idx = setdiff(1:size(lsac,1),binocs2(:,1));
right_monoculars_idx = setdiff(1:size(rsac,1),binocs2(:,2));

if ( Recover_Monoculars )
    left_monoculars_idx(ltest(left_monoculars_idx)>Recover_Monoculars_Threshold) = [];
    right_monoculars_idx(rtest(right_monoculars_idx)>Recover_Monoculars_Threshold) = [];
end
end

function overshoots_idx = FindOvershoots( lxy, lsac, rxy, rsac, samplerate, Overshoot_Time )
MIN_ISI = Overshoot_Time*samplerate/1000;
if ( ~isempty(lsac) && ~isempty( rsac) )
    lmags = getmagnitude( lxy, lsac );
    rmags = getmagnitude( rxy, rsac );
elseif( ~isempty(lsac) )
    lmags = getmagnitude( lxy, lsac );
    rmags = lmags;
    rsac = lsac;
elseif (~isempty(rsac) )
    rmags = getmagnitude( rxy, rsac );
    lmags = rmags;
    lsac = rsac;
else
    % no saccades in either eye
    overshoots_idx = [];
    return;
end

mean_mag = mean([lmags rmags],2);
bad_index =[];
too_close_group = [];
for j = 1:length(mean_mag(:,1))
    % finds the microsaccades that are too close to the current
    % one
    too_close_idx = find(...
        (lsac(:,1) - lsac(j,2) <= MIN_ISI & ...
        lsac(:,1) - lsac(j,2) > 0) | ...
        (rsac(:,1) - rsac(j,2) <= MIN_ISI & ...
        rsac(:,1) - rsac(j,2) > 0) );

    if ( ~isempty(too_close_idx) )
        % if there is any microsaccade too close
        if( isempty(too_close_group) )
            % if there is no current group, start a new one
            too_close_group =  union(j,too_close_idx);
        else
            % if there is already a group, add these
            % microsaccades to the group
            too_close_group =  union(too_close_group,too_close_idx);
        end
    else
        % if the current microsaccade is not too close and
        % there is a current group, deal with it
        if ( ~isempty(too_close_group) )
            % keep only the largest microsaccade
            [max_mag,max_idx] = max(mean_mag(too_close_group));
            bad_index = union(bad_index, too_close_group(setdiff(1:length(too_close_group),max_idx)));

            % initialize the current group
            too_close_group = [];
        end
    end
end
overshoots_idx = bad_index;
end

function magnitude = getmagnitude( xy, sacc  )
magnitude = zeros(size(sacc,1),1);
for i = 1:size(sacc,1)
    idx = sacc(i,1):sacc(i,2);
    % saccade magnitude (amplitude)
    dX = max(xy(idx,1)) - min(xy(idx,1));
    dY = max(xy(idx,2)) - min(xy(idx,2));
    magnitude(i)  = sqrt( dX.^2 + dY.^2);
end
end


function textprogressbar(c)
% This function creates a text progress bar. It should be called with a
% STRING argument to initialize and terminate. Otherwise the number correspoding
% to progress in % should be supplied.
% INPUTS:   C   Either: Text string to initialize or terminate
%                       Percentage number to show progress
% OUTPUTS:  N/A
% Example:  Please refer to demo_textprogressbar.m

% Author: Paul Proteus (e-mail: proteus.paul (at) yahoo (dot) com)
% Version: 1.0
% Changes tracker:  29.06.2010  - First version

% Inspired by: http://blogs.mathworks.com/loren/2007/08/01/monitoring-progress-of-a-calculation/

%% Initialization
persistent strCR;           %   Carriage return pesistent variable

% Vizualization parameters
strPercentageLength = 10;   %   Length of percentage string (must be >5)
strDotsMaximum      = 10;   %   The total number of dots in a progress bar

%% Main

if isempty(strCR) && ~ischar(c),
    % Progress bar must be initialized with a string
    error('The text progress must be initialized with a string');
elseif isempty(strCR) && ischar(c),
    % Progress bar - initialization
    fprintf('%s',c);
    strCR = -1;
elseif ~isempty(strCR) && ischar(c),
    % Progress bar  - termination
    strCR = [];
    fprintf([c '\n']);
elseif isnumeric(c)
    % Progress bar - normal progress
    c = floor(c);
    percentageOut = [num2str(c) '%%'];
    percentageOut = [percentageOut repmat(' ',1,strPercentageLength-length(percentageOut)-1)];
    nDots = floor(c/100*strDotsMaximum);
    dotOut = ['[' repmat('.',1,nDots) repmat(' ',1,strDotsMaximum-nDots) ']'];
    strOut = [percentageOut dotOut];

    % Print it on the screen
    if strCR == -1,
        % Don't do carriage return during first run
        fprintf(strOut);
    else
        % Do it during all the other runs
        fprintf([strCR strOut]);
    end

    % Update carriage return
    strCR = repmat('\b',1,length(strOut)-1);

else
    % Any other unexpected input
    error('Unsupported argument type');
end
end
