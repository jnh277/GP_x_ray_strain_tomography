function [entry,exit,nhat,L,yInds] = find_intersects_3D(shape,lines)
%[entry,exit,nhat,L] = find_intersects_3D(shape,lines)
%   shape is a structure with
%       N: normals to each face, N(i,:) is the normal to the ith face
%       C: centers of each face, C(i,:) is the center of the ith face
%       F: Indices of vertices corresponding to each face, F(i,:) are the
%       indices of the vertices corresponding to face i
%       V: vertices of the shape, V(k,:) = [x,y,z] of the kth vertices
%   lines is a 2D array containing the starting and direction of each line,
%   line(:,j) = [x,y,z,nhat_x,nhat_y,nhat_z]' is the start point and
%   direction of the jth line

[nn,N] = size(lines);
assert(nn==6,'line(:,k) must constain [x,y,znhat_x,nhat_y,nhat_z] transposed');

[num_faces,~] = size(shape.F);

Point = NaN(3,num_faces*N);
Pos = NaN(1,N*num_faces);
IsInOn = NaN(1,N*num_faces);
rayInd = repmat(1:N,1,num_faces);
for i = 1:num_faces
    tri = shape.V(shape.F(i,:),:)';
    [Point(:,N*i-(N-1):N*i), Pos(N*i-(N-1):N*i), IsInOn(N*i-(N-1):N*i)] = intersectLineTri(lines, tri);
end
% remove NaNs
Point(:,~IsInOn) = [];
Pos(~IsInOn) =[];
rayInd(~IsInOn) =[];

% remove duplicates
[C,IA,IC] = unique(Point','rows','stable');
%
Point = C';
Pos = Pos(IA);
rayInd = rayInd(IA);

% sort by distance along ray
[s, sInd] = sort(Pos);
Point = Point(:,sInd);
rayInd = rayInd(sInd);

% sort by ray ind, this must be done after sorting by distance
[rayInd,rInd] = sort(rayInd);
Point = Point(:,rInd);
s = s(rInd);

% % check number of intersections per ray
intsPerRay = hist(rayInd,1:N);
% uneven number of intersects are caught later

% % this stuff only works if all rays have 0 or 2 intersections
magicInd = rayInd*max(intsPerRay)-(max(intsPerRay)-1) + [0 ~diff(rayInd)];

tmpx = NaN(max(intsPerRay),N);
tmpy = NaN(max(intsPerRay),N);
tmpz = NaN(max(intsPerRay),N);

tmpx(magicInd) = Point(1,:)';
tmpy(magicInd) = Point(2,:)';
tmpz(magicInd) = Point(3,:)';

entry = NaN(ceil(max(intsPerRay)/2)*3,N);
exit = NaN(ceil(max(intsPerRay)/2)*3,N);

entry(1:3:end,:) = tmpx(1:2:end,:);
entry(2:3:end,:) = tmpy(1:2:end,:);
entry(3:3:end,:) = tmpz(1:2:end,:);

exit(1:3:end,:) = tmpx(2:2:end,:);
exit(2:3:end,:) = tmpy(2:2:end,:);
exit(3:3:end,:) = tmpz(2:2:end,:);

del = logical((intsPerRay > 4) + ~intsPerRay + mod(intsPerRay,2));
yInds = find(~del);
entry(:,del) = [];
exit(:,del) = [];

L = sqrt(sum((exit-entry).^2));       % this is dodgy but there isnt a 3d hypot function
nhat = (exit-entry)./L;


end

