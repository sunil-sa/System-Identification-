clear all
clc
%loading the data
load('flowdata3.mat');
Fs=zeros(1000,5);
%mean shifting
for i=1:5
    m=sum(Fmeas(:,i))/1000;
    Fs(:,i)=Fmeas(:,i)-m;
end
%Fs=Fmeas;
%If asked to use MLPCA
 %Std=std(1:5);
 %Fs = Fmeas*inv(diag(Std))/sqrt(1000);
%--------------------------------------------------------------------------
% a) part Applying PCA
%--------------------------------------------------------------------------
%singular value deomposition
[U S V]=svd(Fs,0);
%Amat is the constraint matrix
for k = 3:5
    Amat(k-2,:) = V(:,k)';
end
 %Amat=Amat*inv(diag(Std));
sval = diag(S);
eigval=sval.^2;
% subspace angle
theta_pca = 180*subspace(Atrue', Amat')/pi;
Amatd = [Amat(:,1:2) Amat(:,4)];
Amati = [Amat(:,3) Amat(:,5)];
EstimReg = - inv(Amatd)*Amati;
Atrued = [Atrue(:,1:2) Atrue(:,4)];
Atruei = [Atrue(:,3) Atrue(:,5)];
TrueReg = - inv(Atrued)*Atruei;
disp('a part answers');
eigval
maxdiff=max(abs(TrueReg-EstimReg),[],'all')
theta_pca
%--------------------------------------------------------------------------
%b) Applying IPCA to determine error covariance matrix
%--------------------------------------------------------------------------
Amat1=Amat;
flag=1;
iter=0;
Fs=Fmeas;
 sumtest=0;
 while(flag)
     iter=iter+1;
     [ErcovEst1] = stdest(Amat1,Fs');
     L1=diag(ErcovEst1);
     Fs1=Fs*inv(L1)/sqrt(1000);
     [u1 s1 v1]=svd(Fs1,0);
     for z = 3:5
        Amat1(z-2,:) = v1(:,z)';
     end
     Amat1=Amat1*inv(L1);
     theta_pca1 = 180*subspace(Atrue', Amat1')/pi;
     eigval1=diag(s1).^2;
     s1_diag=diag(s1);
     sumtestnew=sum(s1_diag);
     if(abs(sumtestnew-sumtest)<=0.000001)
           flag=0;
     else
           sumtest=sumtestnew;
     end
end
AmatdIPCA = [Amat1(:,1:2) Amat1(:,4)];
AmatiIPCA = [Amat1(:,3) Amat1(:,5)];
EstimRegIPCA = - inv(AmatdIPCA)*AmatiIPCA;
disp('part b answers')
disp('Final Estimated Error Variances are:')
diag(L1*L1')
eigval1
maxdiffIPCA=max(abs(TrueReg-EstimRegIPCA),[],'all')
theta_pca1
%--------------------------------------------------------------------------
% %c part Applying IPCA without prior knowledge of no. of independent
%--------------------------------------------------------------------------
%variables.
for k = 2:5
    Amatc(k-1,:) = V(:,k)';
end
Amat2=Amatc;
flag=1;
sumtest=0;
iter1=0;
while(flag)
    iter1=iter1+1;
   [ErcovEst2] = stdest(Amat2,Fs');
    L2=diag(ErcovEst2);
    Fs2=Fs*inv(L2)/sqrt(1000);
   [u2 s2 v2]=svd(Fs2,0);
    for z = 2:5
      Amat2(z-1,:) = v2(:,z)';
    end
    Amat2=Amat2*inv(L2);
    eigval2=diag(s2).^2;
    s2_diag=diag(s2);
    sumtestnew1=sum(s2_diag);
    if(iter1==6)
           flag=0;
    else
           sumtest=sumtestnew;
    end
end
diag(s2);
disp('part c answers')
disp('Estimated Error Variances are:')
diag(L2*L2')
eigval2
diag(s2)
disp('As the last 4 singular values are not close to 1, we can say that our assumption is wrong.')
%--------------------------------------------------------------------------
% %d part finding the best and worse set of independent variables using IPCA
%--------------------------------------------------------------------------
%constarint model
%Using condition number as a measure
%Total 10 combinations
condition=zeros(10,1);
for p=1:4
    for x=p+1:5
        if p==1
            o=p+x-2;
        elseif p==2
             o=p+x;
            else
             o=p+x+1;
        end
a=sort(setdiff([1 2 3 4 5],[p x]));
AmatDd = [Amat1(:,a(1)) Amat1(:,a(2)) Amat1(:,a(3))];
AmatId = [Amat1(:,p) Amat1(:,x)];
EstimRegd = - inv(AmatDd)*AmatId;
condition(o)=cond(EstimRegd,2);
    end
end
%Checking the consistency
condition1=zeros(10,1);
for p=1:4
    for x=p+1:5
        if p==1
            o1=p+x-2;
        elseif p==2
             o1=p+x;
            else
             o1=p+x+1;
        end
a1=sort(setdiff([1 2 3 4 5],[p x]));
AmatDd1 = [Atrue(:,a(1)) Atrue(:,a(2)) Atrue(:,a(3))];
AmatId1 = [Atrue(:,p) Atrue(:,x)];
EstimRegd1 = - inv(AmatDd1)*AmatId1;
condition1(o1)=cond(EstimRegd1,2);
    end
end