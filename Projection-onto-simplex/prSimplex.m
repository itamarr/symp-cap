function [faceArray,alpha,projectionData] = prSimplex(n,simplex,path,tolerance)
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

%%  Some preparations: Cheching that the input is sensible
if rem(n,1)~=0
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

%create empty data cell array
projectionData = cell(size(simplex,2),2);

%Determine the base change information for every face

%store normal vectors for every face separately
numOfFaces=size(simplex,2);
normalVecs=zeros(dim,numOfFaces);

for j=1:numOfFaces
    %Determine linear subspace corresponding to face
    %k: determines which corner is set to zero
    if j~=1
        k=j-1;
    else
        k=1;
    end
    projectionData{j,1}=simplex(:,k); %store v_tr
    S=simplex;
    S(:,j)=[];
    barycenterS=sum(S,2)/size(S,2);
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
        j
        N=-N;
    end
    normalVecs(:,j)=N;
    %From the above computed matrix S and vector n, determine the base
    %change matrix 'pMat'
    pMat=[gramSchmidt(S),N];
    if norm(pMat'*pMat-eye(dim))>tolerance
        error('Base transform matrix is not orthogonal... Check your code, man!');
    else
        projectionData{j,2}=pMat;
    end
end
%%  Assign the points to a face
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

% dotArray stores the inner products of a connecting vector with each
% foliation vector as a column.

% here size(simplex,2) means the number of vertices in P. In the case of a
% simplex this is equal to the number of faces. For a general polytope this
% is not true. We will need additional input to the functionto account for 
% this in the future.

dotArray=zeros(size(simplex,2),m);
for i=1:m
    for j=1:(dim+1)
        dotArray(j,i)=abs(dot(y(:,i),foliationVecs(:,j)));
    end
end

%
[M,faceArray]=max(dotArray);

alpha=zeros(1,m);
for i=1:m
    yy=y(:,i);
    ff=foliationVecs(:,faceArray(i));
    alpha(i)=dot(yy,ff)/norm(yy);
end
alpha=acos(alpha);

end
