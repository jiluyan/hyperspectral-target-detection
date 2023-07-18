function BandIndex=jily_FVGBS(X,n)
%Fast Volume-Gradient-Based Band Selection (FVGBS), band selection method for
%hyperspectral target detection

%Ji, L., L. Zhu, L. Wang, Y. Xi, K. Yu, and X. Geng. 2021. "FastVGBS: A Fast Version of the Volume-Gradient-Based 
%Band Selection Method for Hyperspectral Imagery."  IEEE Geoscience and Remote Sensing Letters 18 (3):514-7. 
%doi: 10.1109/LGRS.2020.2980108.

%Input: X: the data matrix, L*N, L is the number of bands, N is the number
%of pixels
%       n:the number of selected bands, n should be n<=L

%Output: BandIndex: the band indices for the n selected bands.


% Luyan Ji, jily@mail.ustc.edu.cn
%2023.7.18

L=size(X,2);
X=X-mean(X,1);
invK=inv(X'*X);
BandIndex=1:L;

for k=n+1:L 
    Nd=diag(invK);
    ind=find(Nd ~=0);
    [~,delete]=max(Nd(ind));
    delete=ind(delete);
    BandIndex(delete)=0;
    c=invK(:,delete);
    c1=c(delete);
    c(delete)=0;
    invK(delete,:)=0;
    invK(:,delete)=0;
    invK=invK-c*c'/c1; 
end
BandIndex=(find(BandIndex));

end