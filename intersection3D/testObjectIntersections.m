% developing code for intersects with an object defined by face and vertex
% data
clear all
clc

cube = load('cube.mat');


%% define some lines
N = 5;
lines = NaN(6,N);
lines(:,1) = [0 1e-3 0 1 0 0]';
lines(:,2) = [1e-3,0 0 0 0 1]';
lines(:,3) = [0 0 0 1 0 0]';
lines(:,4) = [1 0 0 0 0 1]';
lines(:,5) = [cube.V(1,:) -1 0 0]';


[entry,exit,nhat,L,yInds] = find_intersects_3D(cube,lines);
% [num_faces,~] = size(cube.F);
% 
% Point = NaN(3,num_faces*N);
% Pos = NaN(1,N*num_faces);
% IsInOn = NaN(1,N*num_faces);
% rayInd = repmat(1:N,1,num_faces);
% for i = 1:num_faces
%     tri = cube.V(cube.F(i,:),:)';
%     [Point(:,N*i-(N-1):N*i), Pos(N*i-(N-1):N*i), IsInOn(N*i-(N-1):N*i)] = intersectLineTri(lines, tri);
% end
% % remove NaNs
% Point(:,~IsInOn) = [];
% Pos(~IsInOn) =[];
% rayInd(~IsInOn) =[];
% 
% % remove duplicates
% [C,IA,IC] = unique(Point','rows','stable');
% %
% Point = C';
% Pos = Pos(IA);
% rayInd = rayInd(IA);
% 
% % sort by distance along ray
% [s, sInd] = sort(Pos);
% Point = Point(:,sInd);
% rayInd = rayInd(sInd);
% 
% % sort by ray ind, this must be done after sorting by distance
% [rayInd,rInd] = sort(rayInd);
% Point = Point(:,rInd);
% s = s(rInd);
% 
% % % check number of intersections per ray
% intsPerRay = hist(rayInd,1:N);
% % % this stuff only works if all rays have 0 or 2 intersections
% magicInd = rayInd*max(intsPerRay)-(max(intsPerRay)-1) + [0 ~diff(rayInd)];
% 
% tmpx = NaN(max(intsPerRay),N);
% tmpy = NaN(max(intsPerRay),N);
% tmpz = NaN(max(intsPerRay),N);
% 
% tmpx(magicInd) = Point(1,:)';
% tmpy(magicInd) = Point(2,:)';
% tmpz(magicInd) = Point(3,:)';
% 
% entry = NaN(max(intsPerRay)/2*3,N);
% exit = NaN(max(intsPerRay)/2*3,N);
% 
% entry(1:3:end,:) = tmpx(1:2:end,:);
% entry(2:3:end,:) = tmpy(1:2:end,:);
% entry(3:3:end,:) = tmpz(1:2:end,:);
% 
% exit(1:3:end,:) = tmpx(2:2:end,:);
% exit(2:3:end,:) = tmpy(2:2:end,:);
% exit(3:3:end,:) = tmpz(2:2:end,:);
% 
% del = logical((intsPerRay > 4) + ~intsPerRay + mod(intsPerRay,2));
% yInds = find(~del);
% entry(:,del) = [];
% exit(:,del) = [];
%%
sc = 1e-2;
figure(1)
clf
p = patch('Faces',cube.F,'Vertices',cube.V);
p.FaceAlpha =0.1;
axis equal
hold on

for i = 1:N
    plot3([lines(1,i)-lines(4,i)*sc lines(1,i)+lines(4,i)*sc],...
        [lines(2,i)-lines(5,i)*sc lines(2,i)+lines(5,i)*sc],...
        [lines(3,i)-lines(6,i)*sc lines(3,i)+lines(6,i)*sc])   
end
% plot3(Point(1,:),Point(2,:),Point(3,:),'go')
plot3(entry(1,:),entry(2,:),entry(3,:),'go','MarkerFaceColor','g')
plot3(exit(1,:),exit(2,:),exit(3,:),'ro','MarkerFaceColor','r')
hold off
xlim([min(cube.V(:,1))*1.1 max(cube.V(:,1))*1.1])
ylim([min(cube.V(:,2))*1.1 max(cube.V(:,2))*1.1])
zlim([min(cube.V(:,3))*1.1 max(cube.V(:,3))*1.1])
view(3)

