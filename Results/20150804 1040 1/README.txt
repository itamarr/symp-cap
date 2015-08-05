The file in this folder contains the capacities of simplices with these vertices:
[ 1,    0, -a2_3/a1_2,                                 -a2_4/a1_2]
[ 0, a1_2,       a1_3,                                       a1_4]
[ 0,    0,          1,                                          0]
[ 0,    0,          0, a3_4 - (a1_3*a2_4)/a1_2 + (a1_4*a2_3)/a1_2]
(The vertices are the columns of this matrix)
Where the coordinates are: p1,p2,q1,q2.
We fix ai_j to be 1, except for a1_2, a1_4 which run on -5 to 5 with steps 0.5

To analyze the results, run either bars_generation.m or AnalyzeResultsScript.m.
