function this = createLocomotionTrajectory(this,thisTrialData)

    % walking = 1.31, jogging = 3.25, running = 5.76, driving = 13.41
    zinc = thisTrialData.WalkingSpeed/this.camera.fps; % gives the delta-z per frame
    this.camera.deltaz = repmat(zinc,this.camera.fps*this.ExperimentOptions.TrialDuration,1);

    % the target angle is specified, per trial, in degrees
    ang = deg2rad(thisTrialData.HeadingChange);

    % tangent is opposite/adjacent = ang, which we rearrange to opposite = adjacent*theta
    xinc = tan(ang)*zinc;

    % first specify the interval where we are just travelling straight
    xstraight = zeros(round(this.camera.fps*thisTrialData.HeadingChangeOnsetTime),1);

    if strcmp(this.ExperimentOptions.Smoothing,'Gaussian')
        % gaussian ramp
        n = (this.camera.fps*this.ExperimentOptions.HeadingChangeDuration);
        mu = n/2;
        sigma = n/7;
        xcurve = (normcdf(1:n,mu,sigma).*xinc)';
    elseif strcmp(this.ExperimentOptions.Smoothing,'Linear')
        % linear x-velocity ramp
        xcurve = linspace(0,xinc,this.camera.fps*this.ExperimentOptions.HeadingChangeDuration)';
    elseif  strcmp(this.ExperimentOptions.Smoothing,'None')
        % sudden jump in x velocity
        xcurve = [zeros(this.camera.fps*this.ExperimentOptions.HeadingChangeDuration/2,1); ...
            repmat(xinc,this.camera.fps*this.ExperimentOptions.HeadingChangeDuration/2,1)];
    else
        error('Unknown Smoothing Type')
    end

    % then the final straight-ahead direction
    xstraight2 = repmat(xinc,round(this.camera.fps*(this.ExperimentOptions.TrialDuration-thisTrialData.HeadingChangeOnsetTime-this.ExperimentOptions.HeadingChangeDuration)),1);
    this.camera.deltax = [xstraight;xcurve;xstraight2];

end