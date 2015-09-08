This results folder contains data files of a script that calculated the ratio c(K-K)/c(K) for K the standard simplex of dimension 2*n, for n=2,...,13.

Besides the workspace and the log file, there's a a csv file with the following structure for each row:
n,cK,length of udotK,udotK,cDiff,length of udotDiff,udotDiff

where udotK and udotDiff are the vectors of the derivatives of the minimizers of the dual action (returned from the Capacity function), for K and K-K respectively.