function bands = jily_NMFW( Y,p,n )

%normalized matched filter weight (NMFW), band selection method for
%hyperspectral target detection

%Ji, L., L. Wang, and X. Geng. 2019. "An Automatic Bad Band Pre-Removal Method for Hyperspectral Imagery."  
%IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing:1-10. 
%doi: 10.1109/JSTARS.2019.2944930.

%Input: Y: the data matrix, L*N, L is the number of bands, N is the number
%of pixels
%       p: the percentage of the selected target pixels in the total pixels, 
%        1-100,generally set to 2~10;
%       n:the number of bad bands, n should be n<=L

%Output: bands: the band indices for the n bad bands.

% Luyan Ji, jily@mail.ustc.edu.cn
%2023.7.18


p=p/100;
[nb,nN]=size(Y);

  for i=1:nb
     Yi=Y(i,:);
     ind=find(Yi ~=0);
     nind=length(ind);
     if nind ==0 
         Yi=rand(nN,1);
         i
     end
     a=mean(Yi);
     Yi=Yi-a;
     Y(i,:)=Yi/norm(Yi); 
  end
 
%首先计算自相关矩阵

R=Y*Y'/nN; 

num=round(p*nN);
ind=randi(nN,1,num);
% ind=indgt(ind);
wcem=zeros(nb,1);
for i=1:num
    %i
    di=Y(:,ind(i));
    wcemi=jily_CEM_R(R,di);
    wcem =wcem+abs(wcemi); 
end
wcem=wcem/num;
[~,ind]=sort(wcem,'descend');
% bands=sort(ind(1:nb-n));
bands=sort(ind(nb-n+1:nb));

end

