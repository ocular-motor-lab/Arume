%
%PAL_ndhisto  Bin multivariable values into multidimensional grid
%
%   syntax: [B edges bin] = PAL_ndhisto(samples)
%
%Internal function
%
%Input:
%
%   samples: Numeric MxN array containing M samples drawn from
%       N-dimensional distribution. 
%
%Optional arguments:
%
%   PAL_ndhisto(samples, gridsize) where 'gridsize' is an integer scalar 
%       uses 'gridsize' bins for each variable in 'samples'. Beware of
%       combinatorial explosion as output 'B' will contain gridsize^N
%       values. Default: 64. Different gridsizes for different dimensions
%       is not supported.
%
%   PAL_ndhisto(samples, edges1, edges2,...), where each 'edges[x]' is a
%       vector, uses for each variable in 'samples' the bin edges specified
%       in 'edges[x]'. If used, edges must be supplied for all N variables 
%       'samples' (i.e., edges1, edges2, ..., edgesN). All must be equal
%       length.
%
%Output: 
%   
%   B: N-dimensional array with side-lengths equal to gridsize (default: 
%       64) containing counts of values in samples that lie between values 
%       given in output edges.  
%
%   edges: (gridsize + 1) x N array containing edges of bins for each
%       dimension of 'B'
%
%   bin: M x N array listing the combination of bin indeces that each value 
%       in 'samples' is contained in.
%
%Example:
%
%   %Draw samples from bivariate normal distribution, bin values in 2-D
%   %grid and create surface plot.
%   samples = randn(200000,2);
%   B = PAL_ndhisto(samples);
%   surface(B)
%
%Introduced: Palamedes version 1.11.13 (NP)
 
function [B edges bin] = PAL_ndhisto(samples,varargin)

gridsize = 64;

edgesSupplied = false;

if ~isempty(varargin)
    if isscalar(varargin{1})
        gridsize = varargin{1};
    elseif length(varargin) == size(samples,2)
        if any(diff(cellfun(@numel,varargin)))
            warning('PALAMEDES:invalidOption','Edge vectors must be equal length. Ignored.');
        else
            edgesSupplied = true;
            for I = 1:size(samples,2)
                edges(:,I) = varargin{I};
            end
        end
    end
end

bin = zeros(size(samples));

for column = 1:size(samples,2)
    if ~edgesSupplied
        edges(:,column) = linspace(min(samples(:,column)),max(samples(:,column)),gridsize+1);
    end
    [trash trash bin(:,column)] = histcounts(samples(:,column),edges(:,column));
end

bin = bin(min(bin,[],2)>0,:);

B = accumarray(bin,1);