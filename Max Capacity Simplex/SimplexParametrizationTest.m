A = GenerateSimplexParamertization(2);

% plug in the numbers for regular polygon in R^4
reg = subs(A, 'a1_2',1);
reg = subs(reg, 'a1_3',0);
reg = subs(reg, 'a1_4',0);
reg = subs(reg, 'a2_3',0);
reg = subs(reg, 'a2_4',0);
reg = subs(reg, 'a3_4',1);
