function [point, pos, isInOn] = intersectLineTriFAST(line, triangle, varargin)
%INTERSECTLINETRIANGLE3D Intersection point of a 3D line and a 3D triangle
% lines is an 6xN with each coumn containing l0 and nhat

%% Default values
[r,N] = size(line);

if r ~= 6
    error('line should have 6 rows')
end
[r,num_tri] = size(triangle);
if r ~= 9
    error('each triangle requires 9 elements')
end

point = NaN(3,N*num_tri);
pos = NaN(N,num_tri);
isInOn = false(N,num_tri);

tol = 1e-13;
if ~isempty(varargin)
    tol = varargin{1};
end


%% Process inputs

% triangle is given as a 3-by-3 array
t0  = triangle(1:3, :);
u   = triangle(4:6, :) - t0;
v   = triangle(7:9, :) - t0;

% triangle normal
n   = cross(u, v);
% test for degenerate case of flat triangle %% FIX THIS
if vectorNorm3d(n) < tol
    return;
end

%% Compute intersection with plane
% % line direction vector
% dir = line(4:6,:);

% vector between triangle origin and line origin
% w0 = line(1:3,:) - t0;
w0x = bsxfun(@minus,line(1,:)',t0(1,:));
w0y = bsxfun(@minus,line(2,:)',t0(2,:));
w0z = bsxfun(@minus,line(3,:)',t0(3,:));

bx = bsxfun(@times,line(4,:)',n(1,:));
by = bsxfun(@times,line(5,:)',n(2,:));
bz = bsxfun(@times,line(6,:)',n(3,:));

% compute projection of each vector on the plane normal
% a = -dot(repmat(n,1,N), w0); 
a = -(repmat(n(1,:),N,1).*w0x + repmat(n(2,:),N,1).*w0y...
    +repmat(n(3,:),N,1).*w0z);
% b = dot(repmat(n,1,N), line(4:6,:));

b = bx+by+bz;

% test case of line parallel to the triangle and ignore these lines
% inds = abs(b) > tol;        % maybe don't need this check as divide by 0 will give NaN??
% NN = sum(inds(:));
% compute intersection point of line with supporting plane
% If r < 0: point before ray
% If r > 1: point after edge
pos = a./ b;


% coordinates of intersection point
% point(:,inds) = line(1:3,inds) + repmat(pos(inds),3,1) .* line(4:6,inds);
pointx = repmat(line(1,:)',1,num_tri) + pos .* repmat(line(4,:)',1,num_tri);
pointy = repmat(line(2,:)',1,num_tri) + pos .* repmat(line(5,:)',1,num_tri);
pointz = repmat(line(3,:)',1,num_tri) + pos .* repmat(line(6,:)',1,num_tri);

point(1,:) = pointx(:)';
point(2,:) = pointy(:)';
point(3,:) = pointz(:)';

%% test if intersection point is inside triangle

% normalize direction vectors of triangle edges
uu  = dot(u, u);
uv  = dot(u, v);
vv  = dot(v, v);

T0 = reshape(repmat(t0,N,1),3,[]);

w = point - T0;
% coordinates of vector v in triangle basis
U = reshape(repmat(u,N,1),3,[]);
V = reshape(repmat(v,N,1),3,[]);
wu = dot(U,w);
wv = dot(V,w);

% normalization constant
D = uv.^2 - uu .* vv;

DD = reshape(repmat(D,N,1),1,N*num_tri);
UU = reshape(repmat(uu,N,1),1,N*num_tri);
UV = reshape(repmat(uv,N,1),1,N*num_tri);
VV = reshape(repmat(vv,N,1),1,N*num_tri);

% test first coordinate
s = (UV .* wv - VV .* wu) ./ DD;
in1 = (s >= -eps) .* (s <=1.0+eps);

% test second coordinate, and third triangle edge
t = (UV .* wu - UU .* wv) ./ DD;
in2 = (t >=-eps).*(s+t <= 1.0+eps);


isInOn = logical(in1.*in2)';
point(:,~isInOn) = NaN;
pos(~isInOn) = NaN;pos = pos(:);

end

