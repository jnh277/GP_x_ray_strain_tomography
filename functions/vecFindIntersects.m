function [mStructs, numM] = vecFindIntersects(lines, FV, dVertices, numFaces, raylength)
% function vecFindIntersects(lines, FV, dVertices, numFaces)
% Johannes Hendriks, 2016
% The University of Newcastle, Australia
% This function takes in ray and sample data and determines intersect
% locations. Stores data into an mStructs structure. 
% This function utilises the inbuilt function polyxpoly in order to
% vectorise the finding of intersects
% Version history
%   parFindIntersects (Johannes Henriks, 2015).
%   as_parFindIntersects (Alex Gregg, 2016).


% Identify how many unique rays have been projected.
[~, ~, numRays] = size(lines);

%%Initialise structures for storing measurement data.
% Initialise the measurement structure template
mStruct.xj = NaN(2,1);                  % exit intersection in assembly coords
mStruct.xjDef = NaN(2,1);               % exit intersection with deformed particle in assembly coords    
mStruct.jFace = NaN;                    % facet with which exit intersection occured
mStruct.xi = NaN(2,1);                  % entry intersectoin in assembly coords
mStruct.xiDef = NaN(2,1);               % entry deformed intersection in assembly coords
mStruct.iFace = NaN;                    % facet with which entry intersection occured
mStruct.nHat = NaN(2,1);                % unit normal of measurement ray
mStruct.pixel = NaN(2,1);               % position of measurement pixel
mStruct.jN = NaN;                       % exit shape function
mStruct.iN = NaN;                       % entry shape function
mStruct.L = NaN;                        % measurement lengths
mStruct.y = NaN;                        % measurement value
mStruct.nPerp = NaN(2,1);
mStruct.A=NaN;
mStruct.Yw=NaN;
mStruct.ray=NaN;

% get ordered indices of all vertices around the object border
i_verts = [FV.faces(:,1); 1];

x_verts = dVertices(i_verts,1);     % x coordinates of ordered deformed boundary vertices
y_verts = dVertices(i_verts,2);     % y coordinates of ordered deformed boundary vertices

% either use provided or calculate max needed ray length
if nargin < 5      % required ray length not provided
    d = max(sqrt(x_verts.^2 + y_verts.^2)); % max distance from origin ('defined center of pixel grid') to a vertex
    raylength = 2*d*1.5;                   % twice radius like distance and then factor of safety
end

nhats = squeeze(lines(:,2,:) - lines(:,1,:));

% build X and Y vectors for the rays
X = [squeeze(lines(1,1,:))' - nhats(1,:)*raylength/2;                       % x coord of start of segment
    squeeze(lines(1,1,:))' + nhats(1,:)*raylength/2;                        % x coord of end of segment
    NaN(1,numRays)];                                                        % breaks segments up
Y = [squeeze(lines(2,1,:))'- nhats(2,:)*raylength/2;                        % y coord of start of segment
    squeeze(lines(2,1,:))'+ nhats(2,:)*raylength/2;                         % y coord of end of segment
    NaN(1,numRays)];                                                        % breaks segments up

X = X(:);                               % reorders into single line where the NaNs represent breaks
Y = Y(:);                       


[xi,yi,ii] = polyxpoly(X(:),Y(:),x_verts,y_verts);
ii(:,1) = (ii(:,1)+2)/3;        % remove the NaN segments from the index numbering

nhati = nhats(:,ii(:,1))';       % direction of ray that created each intersect
pixeli = squeeze(lines(:,1,ii(:,1)));   % pixel location that created each intersect
% check that an even number of intersections was found
if mod(length(ii(:,1)),2)
    warning('An uneven number of intersections was found')
end

% sort ii according to first column, this colum contains the index of which
% ray the intersection came from
[~,Is] = sort(ii(:,1));
pairs = NaN(floor(length(ii(:,1))/2),2);          % store the exit entry pair indeces for each ray
i = 1;
numM = 0;
while(i < length(ii))
    % check that there is a pair of intersections for the ray
    if ii(Is(i),1) == ii(Is(i+1),1)
        numM = numM+1;
        pairs(numM,:) = [Is(i), Is(i+1)];
        i = i+2;
    else
        warning('ray did not have entry and exit intersect, intersect discarded')
        i = i+1;
    end
end

% discard non used pair storage
pairs = pairs(1:numM,:);         % these are pairs telling you what rows of ii belong to the same measurement ray

% Initialise the measurement structure (matrix form)
mStructs = repmat(mStruct,numM,1);          % using the number of measurement rays with two intersects


Ints1 = [xi(pairs(:,1)), yi(pairs(:,1))];       % first intersect in each pair
Ints2 = [xi(pairs(:,2)), yi(pairs(:,2))];       % second intersect in each pair

ray = ii(Is);
rayCells = num2cell(ray(1:2:end));


d1 = sum(Ints1.*nhati(pairs(:,1),:),2);           % distance first intersects are along nhat direction
d2 = sum(Ints2.*nhati(pairs(:,1),:),2);           % distance second intersects are along nhat direction

% now if d1(i) > d2(i) then Ints1(i) is the exit and ints2(i) is the entry
sel = d1 > d2;

jFace = num2cell(sel.*ii(pairs(:,1),2) + (~sel) .* ii(pairs(:,2),2));       % use sel to determine jFace and store as cell
iFace = num2cell((~sel).*ii(pairs(:,1),2) + sel .* ii(pairs(:,2),2));       % use sel to determine iFace and store as cell

xjDef = repmat(sel,1,2).*Ints1 + repmat(~sel,1,2).*Ints2;        % use sel to select correct intersect for xjDef
xiDef = repmat(~sel,1,2).*Ints1 + repmat(sel,1,2).*Ints2;        % use sel to select correct intersect for xiDef

xjDefCells = mat2cell(xjDef',2,ones(numM,1));                   % store as cells for easy putting into mStructs
xiDefCells = mat2cell(xiDef',2,ones(numM,1));

nHatCells = mat2cell(nhati(pairs(:,1),:)',2,ones(numM,1));      % get the nHats and pixel points for each ray
pixelCells = mat2cell(pixeli(:,pairs(:,1)),2,ones(numM,1));     % and store as cells
 



% mStructs
% each row of the structure has fields: 
% s - vector of s co-ords of intersections [0, s1, s2 ... sn]
% nsegs
% iface -> vector, col 1 is first i face, col 2 is second i face
% jface -> vector, col 1 is first j face, col 2 is second j face
% xidef -> matrix, each column is an intersect, row 1 is x, row 2 is y
% xjdef -> matrix, each column is an intersect, row 1 is x, row 2 is y




[mStructs.jFace] = jFace{:};
[mStructs.iFace] = iFace{:};
[mStructs.nHat] = nHatCells{:};
[mStructs.pixel] = pixelCells{:};
[mStructs.xjDef] = xjDefCells{:};
[mStructs.xiDef] = xiDefCells{:};
[mStructs.ray] = rayCells{:};

end