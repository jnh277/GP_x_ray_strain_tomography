addpath('STLRead')

[~,V,~] = stlread('10mm_test_cube.stl');

V = V*1e-3;

DT = delaunayTriangulation(V);
[K,v] = convexHull(DT);

F = K;
V = DT.Points;
V = V-mean(V);

cx = mean(reshape(V(F',1),3,[]));
cy = mean(reshape(V(F',2),3,[]));
cz = mean(reshape(V(F',3),3,[]));
C = [cx',cy',cz'];

TR = triangulation(F,V);
N =faceNormal(TR);

figure(4)
clf
p = patch('Faces',F,'Vertices',V);
p.FaceAlpha =0.1;
axis equal
hold on
plot3(C(:,1),C(:,2),C(:,3),'ro')
quiver3(C(:,1),C(:,2),C(:,3),N(:,1),N(:,2),N(:,3))
hold off

save('cube.mat','F','V','C','N');


