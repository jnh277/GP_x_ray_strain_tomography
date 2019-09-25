addpath(genpath('./geom3d'))

%%

p1 = [0 0 0]';
p2 = [1 0 0]';
p3 = [0.5 0.5 0]';

tri = [p1 p2 p3];
TRI = tri';

p1 = [0 0 0.05]';
p2 = [1 0 0.3]';
p3 = [0.5 0.5 0.1]';
tri2 = [p1 p2 p3];

N = 5;
lines = NaN(6,N);
lines(:,1) = [0.5,0.25,1,0,0,-1];
lines(:,2) = [0.25,0.25,0.1,0.1,0,-0.1];
lines(:,3) = [0,0,1,1,0,0];
lines(:,4) = [2,2,1,0,0,-1];
lines(:,5) = [0 0 0.5,0,0,-1];




% tic
% for i = 1:length(LINES)
% [POINT POS ISINSIDE] = intersectLineTriangle3d(LINES(:,i)', TRI);
% end
% t1 = toc
% 

TRI = repmat([tri(:) tri2(:)],1,1000);
LINES = repmat(lines,1,5000);
tic
point1 = NaN(3,25000,2000);
for i = 1:length(TRI)
[point1(:,:,i), pos, isInOn] = intersectLineTri(LINES, TRI(:,i));
end
t2 =toc             % this is still the fastest           
tic

tic
[point2, pos, isInOn] = intersectLineTriFAST(LINES, TRI);
t3 = toc

[point, pos, isInOn] = intersectLineTri(lines, tri);
[point2, pos, isInOn2] = intersectLineTri(lines, tri2);
point = [point,point2];
isInOn = [isInOn,isInOn2];

% TriMulti = [tri(:) tri2(:)];
% [point, pos, isInOn] = intersectLineTriFAST(lines, TriMulti);

%%

figure(1)
clf
plot3([tri(1,:) tri(1,1)],[tri(2,:) tri(2,1)],[tri(3,:) tri(3,1)])
hold on
plot3([tri2(1,:) tri2(1,1)],[tri2(2,:) tri2(2,1)],[tri2(3,:) tri2(3,1)])
sc =10;
for i = 1:N
    plot3([lines(1,i) lines(1,i)+lines(4,i)*sc],...
        [lines(2,i) lines(2,i)+lines(5,i)*sc],[lines(3,i) lines(3,i)+lines(6,i)*sc])   
end
plot3(point(1,isInOn),point(2,isInOn),point(3,isInOn),'go')
hold off
axis equal
xlim([-0.1 1.1]);
ylim([-0.1 1.1])
zlim([-0.1 0.5])
