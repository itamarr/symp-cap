function [faceArray,projectionData] = prSimplex(n,simplex,path,tolerance)
%PRSIMPLEX Projects path on simplex boundary onto the separate faces
%   INPUT:
%
%   n = half of dimension
%
%   simplex = collection of simplex corners written either in a single
%   vector or in columns of a matrix.
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

%determine m
if size(path,1)==1
    m=size(path,2)/dim;
elseif size(path,2)==1
    m=size(path,1)/dim;
else
    error('Worng dimensions of "path"!!!');
end

if rem(m,1)~=0
    error('Worng dimensions of "path"!!!');
end

if size(simplex,1)~=dim || size(simplex,2)~=(dim+1)
    if size(simplex,2)==dim && size(simplex,2)==(dim+1)
        simplex=simplex';
    else
        error('Wrong dimensions of "simplex"!!!');
    end
end

if tolerance<=0
    error('Unacceptable value for the tolerance!');
elseif tolerance>0.1
    disp('Tolerance seeps pretty big...');
end

%%  Looking at the simplex: determine the base change info
%store simplex in columns of matrix
simplexVecs=zeros(dim,dim+1);
simplexVecs(:)=simplex(:);

%check that the simplex is nondegenerate
S=simplexSetCornerToZero(simplexVecs,1);
if abs(det(S))<tolerance
    error('Simplex is (almost) degenerate!');
end

%create empty data cell array
projectionData = cell(dim+1,2);

%Determine the base change information for every face

%store normal vectors for every face separately
normalVecs=zeros(dim,dim+1);

for j=1:(dim+1)
    %Determine linear subspace corresponding to face
    %k: determines which corner is set to zero
    if j~=1
        k=j-1;
    else
        k=dim; %should be 'k=dim+1' but because a column disappears index shifts
    end
    projectionData{j,1}=simplexVecs(:,k); %store v_tr
    S=simplexVecs;
    S(:,j)=[];
    S=simplexSetCornerToZero(S,k);
    %determine vector orthogonal to this
    n=null(S');
    if mean(size(n)==[dim,1])~=1
        error('unexpected dimensions while using "null"');
    end
    %normalize n
    n=n/norm(n);
    normalVecs(:,j)=n;
    %From the above computed matrix S and vector n, determine the base
    %change matrix 'pMat'
    pMat=[gramSchmidt(S),n];
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

dotArray=zeros(dim+1,m);
for i=1:m
    for j=1:(dim+1)
        dotArray(j,i)=abs(dot(y(:,i),foliationVecs(:,j)));
    end
end
[M,faceArray]=max(dotArray);

end

