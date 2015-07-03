A = GenerateSimplexParamertization(2);

% plug in the numbers for regular polygon in R^4
reg = subs(A, 't1_2',1);
reg = subs(reg, 't1_3',0);
reg = subs(reg, 't1_4',0);
reg = subs(reg, 't2_3',0);
reg = subs(reg, 't2_4',0);
reg = subs(reg, 't3_4',1);