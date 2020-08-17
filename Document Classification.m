clear all
clc
load('news_posts.mat');
Data=double(full(documents'));
%PART A
%Applying PCA 
Data=Data-mean(Data);
[u s v]=svd(Data);
d=diag(s).^2;
disp('Variance explained by the first 3 PCs are:')
Var_exp=[d(1) d(2) d(3)]
%PART B
%Finding the first sparse PCA
X=Data;
num_comp1=1;
%cardinality=18 gives the adj_var value which is almost equal to 75per of
%first PC
card1=10;
[F1,adj_var1,cum_var1] = sparsePCA(X, card1, num_comp1);
disp('Variance explained by the first Sparse PC')
adj_var1
disp('Number of non-zero components in first sparse PC are:')
nnz(F1)
F1=F1';
words=string(wordlist);
r1=find(F1);
first_words1=[];
for i=1:size(r1,2)
    first_words1=[first_words1; words(r1(i))];
end
disp('Here is the list of words first Sparse PC contains:')
first_words1
%PART C
%Finding the second sparse PC
X=Data;
card2=11;
num_comp2=2;
[F2,adj_var2,cum_var2] = sparsePCA(X, card2, num_comp2);
disp('Cummulative frequency explained by the first 2 Sparse PCs:')
cum_var2(2)
Sparse2=F2(:,2)';
disp('Number of non-zero components in 2nd Sparse PC are:')
nnz(Sparse2)
r2=find(Sparse2);
first_words2=[];
for i=1:size(r2,2)
    first_words2=[first_words2; words(r2(i))];
end
disp('The list of words present in second Sparse PC')
first_words2