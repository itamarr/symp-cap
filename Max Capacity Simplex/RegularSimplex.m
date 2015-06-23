%for i=1:10
%c = rand(1)/4;
reg = [zeros(1,4);eye(4)];% - c;
[R,v] = convhulln(reg);
reg = v ^(-1/4) * reg;
avg = reg'*ones(5,1) / 5;
reg = reg - repmat(avg',5,1);
Capacity(reg,2)
%[T, cap] = evalc('Capacity(reg,2);');
%disp('shift, capcity')
%end