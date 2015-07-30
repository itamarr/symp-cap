function [data] = prSimplex(n,simplex,path,tolerance)
%PRSIMPLEX Projects path on simplex boundary onto the separate faces
%   INPUT:
%
%   n = half of dimension
%
%   simplex = collection of simplex corners written either in a single
%   vector or in columns of a matrix.
%   
%   path = 2n*m-dimensional vector representing the path of m points
%
%   tolerance = tolerance of the computation, should have a small value,
%   sensible value would for example be 1e-8.
%
%   OUTPUT:
%   
%   data = cell array containing the projections: 2n+1 rows, each row i
%   encoding the projection of the path onto the face opposite the vertex
%   i.
%   First column contains the projected points as a matrix: vectors written
%   in the columns, have the dimension 2n-1. For details, see comments
%   about the third and fourth column.
%   Second column contains the indices of the points as a row vector,
%   showing which position they have in the path. (Values between 1 and m)
%   Third column contains the translation vector v_tr of the coordinate
%   transformation. Corresponds to one of the 
%   Fourth column contains the (orthonormal) base transformation matrix B
%   (base vectors in the columns). Last base vector is orthogonal to the
%   face.
%   So the new vector v_new and v_old have the relation
%       v_new=(B^-1)*(v_old-v_tr)
%   where the last coordinate of v_new is zero and is omitted, so v_new is
%   (2n-1)-dimensional.

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

%store path in columns of matrix
pathVecs=zeros(dim,m);
pathVecs(:)=path(:);

%store simplex in columns of matrix
simplexVecs=zeros(dim,dim+1);
simplexVecs(:)=simplex(:);

%check that the simplex is nondegenerate
S=simplexSetCornerToZero(simplexVecs,1);
if abs(det(S))<tolerance
    error('Simplex is (almost) degenerate!');
end

%create empty data cell array
data = cell(dim+1,4);

%for every face i opposite vertex i, v_tr is set to the corner i-1.

%determine normal vectors for every face
normalVecs=zeros(dim,dim+1);
for j=1:(dim+1)
    %Determine linear subspace corresponding to face
    %k: determines which corner is set to zero
    if j~=1
        k=j-1;
    else
        k=dim; %should be 'k=dim+1' but because a column disappears index shifts
    end
    data{j,3}=simplexVecs(:,k); %store v_tr
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
    %From the above computed matrix S and vector n, determine the base
    %change matrix B
    data{j,4}=[gramSchmidt(S),n];
    if norm(data{j,4}'*data{j,4}-eye(dim))>tolerance
        error('Base transform matrix is not orthogonal... Check your code, man!');
    end
    normalVecs(:,j)=n;
end

%determine which parts of the path are on which face and store in 'data'

xxx=zeros(dim+1,m);
for i=1:(dim+1)
    for j=1:m
        xxx(i,j)=abs(dot(pathVecs(:,j),normalVecs(:,i)));
        if abs(dot(pathVecs(:,j),normalVecs(:,i)))<tolerance
            data{i,2}=[data{i,2},j];
        end
    end
end
xxx %displaying distances of points from sides of simplex

for i=1:(dim+1)
    vNew=data{i,4}\(pathVecs(:,data{i,2})-repmat(data{i,3},1,size(data{i,2},2)));
    if norm(vNew(dim,:))>tolerance
        disp(strcat('Something is wrong with the projection onto plane',num2str(i),' .. Last coordinate is not zero!!'));
    end
    vNew(dim,:)=[];
    data{i,1}=vNew;
end

end

