function this = setUpShapeAppearanceAndCameraProjectionMatrix(this)

    this.camera.itemcol1 = uint8([1.0000    0.4118    0.1608] * 255);
    this.camera.itemcol2 = uint8([0.3922    0.8314    0.0745] * 255);
    this.camera.fixcol = [255, 255, 255];
    this.camera.hdotsz = round((this.ExperimentOptions.ShapeSz)/2);

    % area of circle versus square
    this.camera.dotarea = pi*(this.camera.hdotsz)^2;
    this.camera.hsquaresz = round(sqrt(this.camera.dotarea)/2);

    % how much additional space to have on the x-axis to account for
    % diagonal motion trajectories
    this.camera.widthfactor = 1.5;
    
    % Define camera parameters
    this.camera.ar = this.Graph.pxHeight/this.Graph.pxWidth;
    this.camera.hfov = this.ExperimentOptions.HFOV;
    this.camera.vfov = this.camera.hfov*this.camera.ar;

    % Set the max distance of the targets (so they aren't predominantly
    % hidden), and the max possibly frame-delay for the onset of the targets
    % (we do not want them to all appear at once). 
    this.camera.fcp = this.ExperimentOptions.FCP;
    this.camera.maxtargetzthreshold = min(75,this.camera.fcp); 
    this.camera.ntargets = this.ExperimentOptions.NumberTargets;

    % dot coordinate grid resolution
    this.camera.nptsX = 500;
    this.camera.nptsZ = 500; 
    
    % set near and far clipping plane in meters. The near-clipping plane will
    % actually be constrained by the view frustum later on
    this.camera.ncp = 0.001;
    
    % standard perspective projection matrix
    this.camera.rhfov = deg2rad(this.camera.hfov / 2);
    this.camera.rvfov = deg2rad(this.camera.vfov / 2);
    this.camera.projection = [
        [1/tan(this.camera.rhfov), 0, 0, 0];
        [0, 1/tan(this.camera.rvfov), 0, 0];
        [0, 0, (-this.camera.ncp-this.camera.fcp)/(this.camera.ncp-this.camera.fcp), -1];
        [0, 0, (2*this.camera.fcp*this.camera.ncp)/(this.camera.ncp-this.camera.fcp), 0]
    ];
    
    % we will later convert our NDUs to pixels on the display
    this.camera.z_offset = [this.Graph.pxWidth/2, this.Graph.pxHeight/2];
    
    % shape lifetime in frames.
    this.camera.fps = this.Graph.frameRate;
    this.camera.nshapeframes = round(this.camera.fps*this.ExperimentOptions.ShapeLifetime);

    % trial duration in frames
    this.camera.ntrialframes = round(this.camera.fps*this.ExperimentOptions.TrialDuration);
    
    % initialize our dots. treat z=depth distribution differently based on desired geometry
    this.camera.maxx = this.camera.fcp*tan(this.camera.rhfov);

    % if we are using the size cue, set up the distance-based scaling
    % factor
    if this.ExperimentOptions.ShapeSizeCue
        this.camera.minshapeszpx = 1.5;
        this.camera.shapescalfac = 1/this.camera.fcp;
    end

end