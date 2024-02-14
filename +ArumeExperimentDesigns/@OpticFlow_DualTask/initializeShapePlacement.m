function this = initializeShapePlacement(this, thisTrialData)

    % render dot locations based on 2d or 3d uniformity condition
    switch thisTrialData.Density
        case 'Uniform'
            error('Only non-uniform condition available currently (stay tuned)')
        case 'NonUniform'
    
            % define a linearly spaced grid of points in x-y space and identify all points that fall outside clipping plane
            this.camera.minzWorld = 1/tan(this.camera.rvfov)*this.ExperimentOptions.ObserverHeight;
            zWorld = repmat(linspace(this.camera.minzWorld, this.camera.fcp, this.camera.nptsZ)', 1, this.camera.nptsX);
    
            % note, if maxx is multiplied by a scalar, the number of dots WILL
            % differ from the specified number above. It's just that it would
            % be pretty complicated to figure out just how many we would need, and how
            % widely distributed in X, when we are co-varying locomotion
            xWorld = repmat(linspace(-this.camera.maxx*this.camera.widthfactor,this.camera.maxx*this.camera.widthfactor,this.camera.nptsX),this.camera.nptsZ,1); 
            xWorldFlatten = reshape(xWorld,1,[]);
            zWorldFlatten = reshape(zWorld,1,[]);
    
    end
    
    yWorldFlatten = zeros(size(xWorldFlatten));
    this.shapes.allValidWorldCoords = [xWorldFlatten',yWorldFlatten',zWorldFlatten',ones(size(yWorldFlatten))'];

    sampleidxs = randi(size(this.shapes.allValidWorldCoords,1), 1, this.ExperimentOptions.NumShapes);
    this.shapes.currentWorldCoords = this.shapes.allValidWorldCoords(sampleidxs,:);
    this.shapes.nValidLocations = size(this.shapes.allValidWorldCoords,1);

end