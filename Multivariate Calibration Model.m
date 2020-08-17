%loading the data
load('Inorfull.mat');
%---------------------------------------------------------------------------
%A) part//////////////////////////////////////////////////////////////////
%---------------------------------------------------------------------------
%Finding the max wavelength and their indices so as to select only those
%indices in absorbance matrix.
[PureCoMax, i] = max(PureCo);
[PureNiMax, j] = max(PureNi);
[PureCrMax, k] = max(PureCr);
disp('part A answers:')
disp('lambda max values are:')
LamdaMax =[PureCoMax PureNiMax PureCrMax]
C=[];
A=[];
%Selecting the first measurement out of the replicates. So each has 26 rows
for t=1:5:130
    C=[C;CONC(t,:)];
    A=[A;DATA(t,:)];
end
%Building an OLS model
Amax=[A(:,i) A(:,j) A(:,k)];
disp('The OLS MODEL matrix is:') 
Model_mat=inv(Amax'*Amax)*Amax'*C
%Selecting the i,j,k columns of data matrix and building a model at the
%same time cross validating using LOOCV.
C_estimated=[];
for m=1:26
    LOOA=[A(1:(m-1),:);A((m+1):26,:)];
    LOOC=[C(1:(m-1),:);C((m+1):26,:)];
    LOOAmax=[LOOA(:,i) LOOA(:,j) LOOA(:,k)];
    B=inv(LOOAmax'*LOOAmax)*LOOAmax'*LOOC;
    c_estm=[A(m,i) A(m,j) A(m,k)]*B;
    C_estimated=[C_estimated; c_estm];
end
Error_mat=C-C_estimated;
disp('RMSE value is:')
RMSE=sqrt(sum(sum(Error_mat.^2))/(26*3))
%--------------------------------------------------------------------------
%B) part///////////////////////////////////////////////////////////////////
%---------------------------------------------------------------------------
disp('part B answers')
%mean shifting
 for a=1:176
     me_an=sum(A(:,a))/26;
     As(:,a)=A(:,a)-me_an;
 end
 As=A;
[u s v]=svd(As);
T=As*v;
%Multivariate calibration model built using the scores obtained by applying PCA
M=(inv(T'*T))*T'*C;
RMSE_PCR=[];
%Selecting the number of PC's
for n=1:5
    C_estimatedpca=[];
    for l=1:26
    LOOAs_pca=[As(1:(l-1),:);As((l+1):26,:)];
    LOOC_pca=[C(1:(l-1),:);C((l+1):26,:)];
    [u1 s1 v1]=svd(LOOAs_pca);
    Tpca=LOOAs_pca*v1(:,1:n);
    Mpca=(inv(Tpca'*Tpca))*Tpca'*LOOC_pca;
    c_estmpca=[As(l,:)*v1(:,1:n)]*Mpca;
    C_estimatedpca=[C_estimatedpca; c_estmpca];
    end
    Error_pca=C-C_estimatedpca;
    RMSE_PCR=[RMSE_PCR;sqrt(sum(sum(Error_pca.^2))/(26*3))];
end
disp('RMSE values for different choices of PCs')
fprintf('no.of PCs -- RMSE\n');
for n=1:5
    fprintf(' %d -- %f\n',n,RMSE_PCR(n));
end
figure(1)
plot(RMSE_PCR,'*-');
xlabel('No. of Principal Components');
ylabel('RMSE');
%--------------------------------------------------------------------------
%C) part///////////////////////////////////////////////////////////////////
%---------------------------------------------------------------------------
disp('part C answers')
Std=[];
for n=1:5:130
    Std=[Std;stdDATA(n,:)];
end
L=diag(mean(Std));
%Scaling the data
Ascaled=As*inv(L);
figure(2)
plot(mean(Std,2),'*-y');
xlabel('Mixtures');
ylabel('Standard deviation values');
figure(3)
plot(mean(Std,1),'*-g');
xlabel('Wavelengths');
ylabel('Standard deviation values');
%From graphs, Error std deviation values vary significantly w.r.t
%wavelength
RMSE_MLPCR=[];
for n=1:5
    C_estimatedmlpca=[];
    for l=1:26
    LOOAs_mlpca=[Ascaled(1:(l-1),:);Ascaled((l+1):26,:)];
    LOOC_mlpca=[C(1:(l-1),:);C((l+1):26,:)];
    [u2 s2 v2]=svd(LOOAs_mlpca);
    Tmlpca=LOOAs_mlpca*v2(:,1:n);
    Mmlpca=(inv(Tmlpca'*Tmlpca))*Tmlpca'*LOOC_mlpca;
    c_estmmlpca=[Ascaled(l,:)*v2(:,1:n)]*Mmlpca;
    C_estimatedmlpca=[C_estimatedmlpca; c_estmmlpca];
    end
    Error_mlpca=C-C_estimatedmlpca;
    RMSE_MLPCR=[RMSE_MLPCR;sqrt(sum(sum(Error_mlpca.^2))/(26*3))];
end
disp('RMSE values for different choices of PCs')
fprintf('no.of PCs -- RMSE\n');
for n=1:5
    fprintf(' %d -- %f\n',n,RMSE_MLPCR(n));
end
figure(4)
plot(RMSE_MLPCR,'*-');
xlabel('No. of Principal Components');
ylabel('RMSE(MLPCR)');
%--------------------------------------------------------------------------
%D) part///////////////////////////////////////////////////////////////////
%--------------------------------------------------------------------------
RMSE_IPCR=[];
Eigenanalysis=[];
%n---no. of PCs selection
for n=1:10
    flag=1;
    sumtest=0;
    %Applying PCA to use the model as starting point in IPCA
    [U S V]=svd(As);
    Model=[];
    for k = n+1:176
      Model(k-n,:) = V(:,k)';
    end
    %Here starts the IPCA method, iterations go on till the target
    %tolerance is met.
    while(flag)
       [L2] = stdest(Model,As');
       L2=diag(L2);
       As2=As*inv(L2);
       [u3 s3 v3] = svd(As2);
       for k = n+1:176
          Model(k-n,:) = v3(:,k)';
       end
       Model=Model*inv(L2);
       s3_diag=diag(s3);
       sumtestnew=sum(s3_diag);
       if(sumtestnew-sumtest<=0.0001)
           flag=0;
       else
           sumtest=sumtestnew;
       end
    end
    Ascaled2=As*inv(L2);
    [u5 s5 v5]=svd(Ascaled2);
    s5_diag=diag(s5);
    Eigenanalysis=[Eigenanalysis s5_diag];
    C_estimatedipca=[];
    for w=1:26    
    LOOAs_ipca=[Ascaled2(1:(w-1),:);Ascaled2((w+1):26,:)];
    LOOC_ipca=[C(1:(w-1),:);C((w+1):26,:)];
    [u4 s4 v4]=svd(LOOAs_ipca);
    Tipca=LOOAs_ipca*v4(:,1:n);
    Mipca=(inv(Tipca'*Tipca))*Tipca'*LOOC_ipca;
    c_estmipca=[Ascaled2(w,:)*v4(:,1:n)]*Mipca;
    C_estimatedipca=[C_estimatedipca; c_estmipca];
    end
    Error_ipca=C-C_estimatedipca;
    RMSE_IPCR=[RMSE_IPCR;sqrt(sum(sum(Error_ipca.^2))/(26*3))];
end
disp('part D answers')
disp('RMSE values for different choices of PCs')
fprintf('no.of PCs -- RMSE\n');
for n=1:10
    fprintf(' %d -- %f\n',n,RMSE_IPCR(n));
end
Eigenanalysis;
figure(5)
plot(log(RMSE_IPCR),'*-r');
xlabel('No. of Principal Components');
ylabel('RMSE(IPCR)');