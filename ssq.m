function t = ssq(A)
%SSQ	SSQ(A) is the sum of squares of the elements of matrix A.
t = sum(sum(A.^2));
