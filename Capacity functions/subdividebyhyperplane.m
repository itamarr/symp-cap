function [data]=subdividebyhyperplane(A,B,C,D,E,F)
P = [1,0,0,0;0,1,A,0;0,0,B,D;0,F/D,C,E;0,0,0,0];
%barycenter = sum(P)/size(P,1)
%P = P - repmat(barycenter,size(P,1),1)
P = P - repmat(min(P),size(P,1),1);

n = size(P,2)/2;


M = max(P);
centerpt = M.*rand(1,2*n);
while not(inhull(centerpt,P));
    centerpt = M.*rand(1,2*n);
end   


%figure
%hold on
%scatter(P(:,1),P(:,2))
%scatter(newv(1),newv(2))

hpnormal = rand(1,2*n);
abovevertices = [];
belowvertices = [];

for i = 1:size(P,1)
vert = P(i,:);
if (vert - centerpt)*hpnormal' > 0
abovevertices(end+1,:) = vert;
else
belowvertices(end+1,:) = vert;   
end
end
   
newvertices = [];

for i = 1:size(abovevertices,1)
for j = 1:size(belowvertices,1)
    
abovevert = abovevertices(i,:);
belowvert = belowvertices(j,:);

t = (centerpt*hpnormal' - belowvert*hpnormal')/(abovevert*hpnormal' - belowvert*hpnormal');
newvertices(end+1,:) = t*abovevert + (1-t)*belowvert;

end
end

aboveP = abovevertices;
for i = 1:size(newvertices,1)
aboveP(end+1,:) = newvertices(i,:);
end

belowP = belowvertices;
for i = 1:size(newvertices,1)
belowP(end+1,:) = newvertices(i,:);   
end

%figure
%scatter(aboveP(:,1),aboveP(:,2))
%figure
%scatter(belowP(:,1),belowP(:,2))

[c1,char1,udot1] = Capacity(aboveP,n);
[c2,char2,udot2] = Capacity(belowP,n);
[c,char,udot] = Capacity(P,n);

c
additiveCap = c1 + c2

comment='No comment.';

if(additiveCap==c)
    comment='equal. Interesting...';
elseif(additiveCap<c)
    comment='counterexample!!!';
end

disp(comment);

data=cell(1,12);
data{1}=[A,B,C,D,E,F];
data{2}=[c1,c2,c];
data{3}=char1;
data{4}=char2;
data{5}=char;
data{6}=udot1;
data{7}=udot2;
data{8}=udot;
data{9}=comment;
data{10}=aboveP;
data{11}=belowP;
data{12}=P;
end
