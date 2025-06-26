function im_out = alphaShape2im_fast(shp,center_coord,im_size,scaling)
%ALPHASHAPE2IM Convert alphaShape object to binary image
%   im_out = alphaShape2im(shp,center_coord,im_size,scaling) creates a
%   binary image with size im_size that represents a given alphaShape
%   "shp". "alphaShape2im" maps the center voxel(s) to the reference center
%   coordinate "center_coord". Each voxel represents an isotropic volume
%   which can be scaled with a factor "scaling"; higher scaling factors
%   correspond to higher resolutions.
%
%   "scaling" is an optional argument (default 1 - native resolution
%   relative to coordinate system).
%
%   Nade Sritanyartatana
%   October 31, 2014
%   v1.0

if ~exist('scaling','var')
    scaling = 1;
end

% Map image voxel locations to coordinates
coords = mapVox2coords(im_size,center_coord,scaling);

% Initialize output
im_out = false(im_size);

% Create spheres of a given radius along every voxel that corresponds with
% a point in shp.
im_query = false(im_size);
p = shp.Points;
subs = mapPoints2subs(im_size,center_coord,p,scaling);
subs = subs(~isnan(subs));
subs = reshape(subs,[],3);
idx  = sub2ind(im_size,subs(:,1),subs(:,2),subs(:,3));
im_query(idx) = true;
% r = shp.Alpha; % Radius of spheres
r = scaling*shp.Alpha*2;
im_query = binaryspheres(im_query,r);

% Render the binary image along the queried coordinates
query_coords = cellfun(@(x)x(im_query),coords,'UniformOutput',false);
im_out(im_query) = shp.inShape(query_coords{:});

% Fill the gaps/holes using connected components and binary image
% operations. Check what value to fill in each connected component by
% querying a single point using "inShape"
CC = bwconncomp(im_out);
for i=1:length(CC.PixelIdxList)
    query_pt = CC.PixelIdxList{i}(1);
    if ~im_out(query_pt)
        query_coord = cellfun(@(x)x(query_pt),coords,'UniformOutput',false);
        val = shp.inShape(query_coord{:});
        im_out(CC.PixelIdxList{i}) = val;
    end
end

end

function coords = mapVox2coords(im_size,center_coord,scaling)
% Input: image voxel locations, center_coord, scaling factor
% Output: Cartesian coordinates (same size as image voxel locations)

% Create 1D arrays that describe grid in each direction
mesh_inputs = [ones([1, length(im_size)]); im_size];
mesh_inputs = mat2cell(mesh_inputs, 2, ones([1, length(im_size)]) );
mesh_inputs = cellfun(@(x)x(1):x(end), ...
                      mesh_inputs, ...
                      'UniformOutput', false);

voxlocs = cell(1,length(mesh_inputs)); % Initialize voxel locations
% Calculate image voxel locations based on im_size
if length(im_size)==2
    [voxlocs{1} voxlocs{2}] = meshgrid(mesh_inputs{:});
elseif length(im_size)==3
    [voxlocs{1} voxlocs{2} voxlocs{3}] = meshgrid(mesh_inputs{:});
else
    error('MATLAB:alphaShape2Im:im_size','im_size can have only 2 or 3 elements');
end

coords = cell(1,length(im_size)); % Initialize Cartesian coordinates

% Find a transformation function to map voxlocs to coords
% For each dimension of im_size:
%  If even, then map left center voxel 0.5 to the left of center_coord
%  If odd, then map center voxel onto center_coord
for i=1:length(coords)
    coords{i} = voxlocs{i} - im_size(i)/2 - 0.5; % "Center" the locations to 0
    coords{i} = coords{i}/scaling; % Scale the coordinates
    coords{i} = coords{i} + center_coord(i); % Offset the locations by center_coord
end

end

function subs = mapPoints2subs(im_size,center_coord,points,scaling)
% MAPPOINTS2SUBS Map points to subscripts
% subs = mapPoints2subs(im_size,center_coord,points,scaling) linearly
% transforms "points" to subscript indices.
    points = bsxfun(@minus,points,center_coord);
    points = points*scaling;
    points = round(bsxfun(@plus,points,im_size/2) + 0.5);
    SZ = repmat(im_size,[size(points,1) 1]);
    points( any(points>SZ | points<1, 2), : ) = nan;
    subs = points;
end

function im = binaryspheres(bwim,r)
% BINARYSPHERES Generate spheres of N-D binary image
% im = binaryspheres(bwim,r) is a simple function that utilizes the
% distance transform "bwdist" to generate an N-D binary image where there
% are spheres of a radius r centered at every non-zero voxel.
    im = bwdist(bwim);
    im = im <= r;
end