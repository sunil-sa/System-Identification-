clear all
clc
%Storing all the images in the imds
imds=imageDatastore('yalefaces','IncludeSubfolders',1,'LabelSource','foldernames');
D=[];
for m=0:14
for j=1:6
i=readimage(imds,m*6+j);
i=im2double(i);
D=[D i(:)];
end
end
S=D;
Sr=zeros(15,77760); % Represenatative image matrix
StdS = zeros(15,77760);
%Compressing 15 Image to 1 image
z=[];
for i=0:14
 S1 = S(:,(i*6+1):(i*6+6));
 M2=mean(S1,2); %Mean of each pixel
 M = repmat(M2,[1,6]);
 S1= 255*im2double(S1);
 S1s=S1-M; %Mean shifted data
 CovS1 = cov(S1s); %Covariance matrix
 [V,D] = eig(CovS1);
 Vlarge=V(:,6); %eigen vector corresponding to largest eigen value
 X=reshape(S1s*Vlarge,[243,320]);
 z=[z X];
 Sr(i+1,:)=(S1s*Vlarge+ M2)';
end
%Finding Standard deviation of error
for i=0:14
 for j=1:6
 StdS(i+1,:)=StdS(i+1,:)+(Sr(i+1,:)-255*im2double(S(:,i*6+j))').^2;
 end
 StdS(i+1,:)= StdS(i+1,:)/5;
 u=unique(StdS(i+1,:));
 m=u(2);
 for k=1:77760
 if StdS(i+1,k)==0
 %weighting won't be possible for if stdS == 0
 %hence 0 is replaced with min of StdS exept 0
 StdS(i+1,k)=m;
 end
 end
end
Edis=zeros(15,1); % Euclidian distance vector
correct_match = 0; % No.of correct matches
%testing model
for i=0:14
 for j=1:6
 for k=0:14
 Edis(k+1)=sum(((Sr(k+1,:)-255*im2double(S(:,i*6+j))').^2)./StdS(i+1,:));
 end
 [H,J]=min(Edis);
 if J == (i+1)
 correct_match = correct_match + 1;
 end
 end
end
fprintf('correct matches = %d \n', correct_match);
fprintf('Accuracy = %f %%\n', correct_match*100/90);
% for j=0:14
% figure
% imshow(uint8(z(:,320*j+1:320*(j+1))))
% hold on
% end