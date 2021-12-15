A=[1 2 3;4 5 6;7 8 9];
B=triu(A);
sum_value = sum(sum(B));
B=B/sum_value;
A=tril(B',-1)+B;