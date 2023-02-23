M = 3;
N = 5;
oldMat = rand(M, N);

clear newMat1
for ii=1:M
  newMat1(ii,:,:) = diag(oldMat(ii,:));
end
newMat1
reshape(newMat1,N,M,N)

clear newMat2
newMat2 = zeros(M,M,N);
%newMat2(bsxfun(@plus,(1:M)', (0:N-1)*(M*N+M))) = oldMat;
newMat2(reshape(find(repmat(logical(eye(M)),[1,1,N])),M,N)) = oldMat;
newMat2

isequal(newMat1,newMat2)

