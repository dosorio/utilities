# function idx = knnsearch_same(X,k)
# % same as idx = knnsearch(X,[],k), but faster
# aa = sum(X.*X);
# n = size(aa,2);
# D = repmat(aa',[1 n]) + repmat(aa,[n 1]) - 2*X'*X;
# [~,idx] = sort(D);
# idx = idx(2:k+1,:)';

knnSearch <- function(X, k){
  aa = colSums(X * X)
  n = length(aa)
  D = matrix(rep(aa, each=n), ncol = n, byrow = TRUE) + matrix(rep(aa, each=n), ncol = n, byrow = FALSE) - (2 * t(X) %*% X)
  D <- apply(D,2,order)
  return(D[2:(k+1),])
}
