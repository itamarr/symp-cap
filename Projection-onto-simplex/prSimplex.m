function [faceArray,alpha,projectionData,distMatrix] = prSimplex(n,simplex,path,tolerance)
%PRSIMPLEX Projects path on simplex boundary onto the separate faces
%   INPUT:
%
%   n = half of dimension
%
%   simplex = collection of simplex vertiecs written in the colums of a
%   matrix.
%
%   path = path on the surface of the simplex. Should be either in a 
%
%   tolerance = tolerance of the computation, should have a small value,
%   sensible value would for example be 1e-8.
%
%   OUTPUT:
%   
%   faceArray = row vector assigning every connecting vector between two
%   path vectors 
%
%   alpha = angle between the characteristic foliation vector and the
%   connecting vector between two path vectors
%
%   projectionData = cell array containing the projections: 2n+1 rows, each row i
%   encoding the projection of the path onto the face opposite the vertex
%   i.
%   First column contains the translation vector v_tr of the coordinate
%   transformation. Corresponds to one of the 
%   Second column contains the (orthonormal) base transformation matrix B
%   (base vectors in the columns). Last base vector is orthogonal to the
%   face.
%   So the new vector v_new and v_old have the relation
%       v_new=(B^-1)*(v_old-v_tr)
%   where the last coordinate of v_new is zero and is omitted, so v_new is
%   (2n-1)-dimensional.
%
%   distMatrix = distances of the path vectors to the faces, distance of
%   vector j from face i is stored in distMatrix(i,j)

%%  Some preparations: Cheching that the input is sensible
if rem(n,1)~=0 || n<=0
    error('Invalid value of "n"');
end
dim=2*n;

%determine path length m
if size(path,1)==1
    m=size(path,2)/dim;
elseif size(path,2)==1
    m=size(path,1)/dim;
else
    error('Wrong dimensions of "path"!!!');
end

if rem(m,1)~=0
    error('Wrong dimensions of "path"!!!');
end

% if size(simplex,1)~=dim || size(simplex,2)~=(dim+1)
%     if size(simplex,2)==dim && size(simplex,2)==(dim+1)
%         simplex=simplex';
%     else
%         error('Wrong dimensions of "simplex"!!!');
%     end
% end

if tolerance<=0
    error('Unacceptable value for the tolerance!');
elseif tolerance>0.1
    disp('Tolerance seeps pretty big...');
end

barycenter = sum(simplex,2)/size(simplex,2);


%%  Looking at the convex polytope: determine the base change info

%check that the simplex is nondegenerate
S=simplexSetCornerToZero(simplex,1);
if abs(det(S))<tolerance
    error('Simplex is (almost) degenerate!');
end

numOfFaces=size(simplex,2);

%create empty data cell array
projectionData = cell(numOfFaces,2);

%Determine the base change information for every face

%store normal vectors for every face separately

normalVecs=zeros(dim,numOfFaces);

for j=1:numOfFaces
    %Determine linear subspace corresponding to face
    %k: determines which corner is set to zero
    if j~=1
        k=j-1;
    else
        k=1;
    end
    S=simplex;
    S(:,j)=[];
    barycenterS=sum(S,2)/size(S,2);
    projectionData{j,1}=S(:,k); %store v_tr
    S=simplexSetCornerToZero(S,k);
    %determine vector orthogonal to this
    N=null(S');
    if mean(size(N)==[dim,1])~=1
        error('unexpected dimensions while using "null"');
    end
    %normalize N
    N=N/norm(N);
    %correct the sign of N by looking at vector from barycenter to
    %barycenter of a face
    ddd=barycenterS-barycenter;
    if dot(ddd,N)<0
        N=-N;
    end
    normalVecs(:,j)=N;
    %From the above computed matrix S and vector n, determine the base
    %change matrix 'pMat'
    pMat=[gramSchmidt(S),N];
    shouldBeZero=pMat'*pMat-eye(dim);
    if norm(shouldBeZero(:))>tolerance
        error('Base transform matrix is not orthogonal... Check your code, man!');
    else
        projectionData{j,2}=pMat;
    end
end
%%  Assign the points to a face using the foliation
%Calculate foliation vectors
J = [zeros(n),eye(n);-eye(n),zeros(n)];
foliationVecs=J*normalVecs;

%store path in columns of matrix
pathVecs=zeros(dim,m);
pathVecs(:)=path(:);

%Calculate 'connecting vectors' between path vectors
y=zeros(dim,m);
for i=1:(m-1)
    y(:,i)=pathVecs(:,i+1)-pathVecs(:,i);
end
y(:,m)=pathVecs(:,1)-pathVecs(:,m);

% dotMatrix stores the inner products of a connecting vector with each
% foliation vector as a column.

% here size(simplex,2) means the number of vertices in P. In the case of a
% simplex this is equal to the number of faces. For a general polytope this
% is not true. We will need additional input to the functionto account for 
% this in the future.

dotMatrix=zeros(size(simplex,2),m);
for i=1:m
    for j=1:(dim+1)
        dotMatrix(j,i)=abs(dot(y(:,i),foliationVecs(:,j)));
    end
end

%
[M,faceArray]=max(dotMatrix);

alpha=zeros(1,m);
for i=1:m
    yy=y(:,i);
    ff=foliationVecs(:,faceArray(i));
    alpha(i)=dot(yy,ff)/norm(yy);
end
alpha=acos(alpha);


%%  Calculate distances of path vectors from faces, and if n<=2, plot

distMatrix=zeros(numOfFaces,m);

if n==2
    %first column will contain path vectors, second column will contain
    %indices from the set {1,...,m} corresponding to path vectors
    projectedVectors=cell(numOfFaces,2);
end

for k=1:numOfFaces
    B=projectionData{k,2};
    vTr=projectionData{k,1};
    for i=1:m
        ProjectedVec=B\(pathVecs(:,i)-vTr); %B\v calls a more efficient
        %Matlab algorithm to compute (B^-1)*v
        dist=abs(ProjectedVec(dim));
        %if n==2 && dist<tolerance
        %Above line commented out because this doesn't work!
        if n==2 && faceArray(i)==k
            %Information for drawing and labeling
            projectedVectors{k,1}=[projectedVectors{k,1},ProjectedVec(1:(dim-1))];
            projectedVectors{k,2}=[projectedVectors{k,2},i];
        end
        distMatrix(k,i)=dist;
    end
end

if n==1
    drawSimplex=[simplex,simplex(:,1)];
    figure
    title('trajectory on simplex')
    plot(drawSimplex(1,:),drawSimplex(2,:))
    hold on
    plot(pathVecs(1,:),pathVecs(2,:))
    hold on
    labels=cellstr(num2str([1:m]'))';
    text(pathVecs(1,:),pathVecs(2,:),labels)
elseif n==2
    for j=1:numOfFaces
        S=simplex;
        S(:,j)=[];
        B=projectionData{j,2};
        vTr=projectionData{j,1};
        S=B\(S-repmat(vTr,1,size(S,2)));
        if norm(S(dim,:))>tolerance
            error('Projection of simplex corners not working as expected..');
        end
        S=S(1:(dim-1),:);
        drawSimplexFace=[S,S(:,[1,3,2,4])];
        ttl=strcat('face ',num2str(j));
        figure('name',ttl)
        title(ttl)
        plot3(drawSimplexFace(1,:),drawSimplexFace(2,:),drawSimplexFace(3,:))
        hold on
        PP=projectedVectors{j,1};
        plot3(PP(1,:),PP(2,:),PP(3,:))
        hold on
        labels=cellstr(num2str(projectedVectors{j,2}'))';
        text(PP(1,:),PP(2,:),PP(3,:),labels)
    end
end

end
