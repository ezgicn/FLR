N=50;k=1; ss=100;
x1=normrnd(0,1,ss,1); % x'in birinci sütunu
x2=normrnd(0,1,ss,1); % x'in ikinci sütunu

X=[ones(ss,1) x1 x2]; x=[x1 x2]; Xint=ones(ss,1); Xyeniint=ones(ss,1);
y=binornd(1, 0.10, ss, 1); 
beta0=0; beta1=0; beta2=0;

beta=[beta0;beta1;beta2];

TPB=1./(1+exp(-X*beta));

tempbiasLR=zeros(N,3);
tempbiasfirth=zeros(N,3);

topprobbiasLR=zeros(1,N);
topprobbiasfirth=zeros(1,N);

topprobbiasLRkare=zeros(1,N);
topprobbiasfirthkare=zeros(1,N);

while k<N+1
betaLR=beta;
betaLR_old=betaLR+1;

betafirth=beta;
betafirth_old=betafirth+1;

%%%PARAMETRE TAHMÝNLERÝ, Newton-Raphson yöntemi
%%%Lojistik regresyon yöntemi (LR)
while (max(abs(betaLR-betaLR_old))>0.001)
u=X*betaLR;
p = 1./(1+exp(-u));
W=p.*(1-p);
z=u+(y-p)./W;
betaLR_old=betaLR;
betaLR=inv(X'*diag(W)*X)*X'*diag(W)*z;
end
betaLR1(k)=betaLR(1); betaLR2(k)=betaLR(2); betaLR3(k)=betaLR(3);


%%%Firth's Lojistik regresyon yöntemi (FLR)
while(max(abs(betafirth-betafirth_old))>0.001) 
uf=X*betafirth;
pfirth = 1./(1+exp(-uf));
Wfirth=pfirth.*(1-pfirth);
Hfirth=sqrt(diag(Wfirth))*X*inv(X'*diag(Wfirth)*X)*X'*sqrt(diag(Wfirth));
ufirth=diag(Hfirth).*(pfirth-1/2);  
Uscore=(X'*(y-pfirth))-(X'*ufirth);
betafirth_old=betafirth;
betafirth=betafirth_old+(inv(X'*diag(Wfirth)*X))*Uscore;
end
betafirth1(k)=betafirth(1);
betafirth2(k)=betafirth(2);
betafirth3(k)=betafirth(3);

%%%%LR için parametre yanlýlýðý 
tempbiasLR(k,1)=betaLR(1)-beta0;
tempbiasLR(k,2)=betaLR(2)-beta1;
tempbiasLR(k,3)=betaLR(3)-beta2;
%%%%FLR için parametre yanlýlýðý 
tempbiasfirth(k,1)=betafirth(1)-beta0;
tempbiasfirth(k,2)=betafirth(2)-beta1;
tempbiasfirth(k,3)=betafirth(3)-beta2;

%%%LR için kestirilen olasýlýk (predicted probability)
ppbLR=1./(1+exp(-X*betaLR));

%%%FLR için kestirilen olasýlýk 
ppbfirth=1./(1+exp(-X*betafirth));

%%%LR için kestirilen olasýlýk yan (Predicted probability bias)
probbiasLR=ppbLR-TPB; probbiasLRkare=(ppbLR-TPB).^2;
topprobbiasLR(k)=sum(probbiasLR); topprobbiasLRkare(k)=sum(probbiasLRkare);

%%%FLR için kestirilen olasýlýk yan (Predicted probability bias)
probbiasfirth=ppbfirth-TPB; probbiasfirthkare=(ppbfirth-TPB).^2;
topprobbiasfirth(k)=sum(probbiasfirth); topprobbiasfirthkare(k)=sum(probbiasfirthkare);


check1 = isnan(tempbiasLR(k,1)); check2 = isnan(tempbiasfirth(k,1));
 if check1 == 1| check2 == 1;
            k=k;
        else
            k=k + 1
        end
        if k==0; k=1;
         check1 = 0; check2 = 0;
        end
end

%%%%LR ortalama parametre deðerleri
ortparLR=[mean(betaLR1); mean(betaLR2);mean(betaLR3)]
ortparfirth=[mean(betafirth1); mean(betafirth2);mean(betafirth3)]

%%%LR için kestirilen parametre yan (predicted parameter bias)
meanbiasLR= mean(tempbiasLR);

%%%FLR için kestirilen parametre yan (predicted parameter bias)
meanbiasfirth= mean(tempbiasfirth);

%%%LR için parametre yanlar (estimated parameter bias)
parametreyanbeta0=[meanbiasLR(1); meanbiasfirth(1)]                                    
parametreyanbeta1=[meanbiasLR(2); meanbiasfirth(2)]                                     
parametreyanbeta2=[meanbiasLR(3); meanbiasfirth(3)]

%%%%LR için kestirilen olasýlýk yan (predicted probability bias) ve RMSE
sonprobbiasLR=sum(topprobbiasLR)/ss*N; 
stdprobbiasLR=std(topprobbiasLR); 
RMSELR=sqrt((sum(topprobbiasLRkare))/ss*N);

%%%%FLR için kestirilen olasýlýk yan (predicted probability bias) ve RMSE
sonprobbiasfirth=sum(topprobbiasfirth)/ss*N;
stdprobbiasfirth=std(topprobbiasfirth);
RMSEfirth=sqrt((sum(topprobbiasfirthkare))/ss*N);

