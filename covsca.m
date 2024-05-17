function  [Z, C,fp,dys, func] = covsca(S,L,Q,Cstart,conC,nanal)
                        
% ALS algorithm for Misca2, possibly with maximally L rows of C fixed
% to be binary (1/0), with maximally one fixed 1 per column of C
%
% To partially identify solution: Rescale C columnwise, such that max(C)=1,
% and compensate in Zar(:,:,l) 
%   
%
% INPUT:
% S         =   (J,JxK) K Observed covariance matrices S_k, k=1,...,K, positioned next to each other
% L         =   (1x1) # clusters
% Q         =   (1x1) or (Lx1) # components per cluster
%
% Cstart    =   EITHER:  Cstart = (1x1) arbitrary scalar: initial C will be
%                 rowwise random permutation of [eye;zeros];
%               OR:     (KxL) matrix with binary indicators, per row maximally 1 weight of 1, 0 otherwise:  
%                       Rows with Weight of 1 will be kept fixed during analysis, 
%                       Rows with only zeros will be estimated freely; (E.g., Cstart = [1 0 0;0 0 0; 0 0 0; 0 0 1], 
%                       C = [1 0 0;x x x;x x x;0 0 1], with x freely estimated)
%                       
%              
% conC      =   (1x1) 1 = no row sum constraint on C; 
%                     2 = 1'C'=1 (i.e., sigma_l c_kl=1, k=1,...,K)
% 
% nanal     =   (1x1) # requested randomly started analyses
% [Z,C,fp,dys] = covsca(S,L,Q,1,1,1000)
%
%
% OUTPUT:
% Z   =   (Jxsum(Q)) loading matrices Z_l, l=1,...,L, positioned next to
% each other
% C     =   (KxL) matrix with weights for L clusters for K covariance
% matrices; constraint on C: c_ij>=0, all i,j; 
% bestFp=   (1x1) fit%
% out   =   (runsx1) ssq(Zsup-ZEsup) for all requested runs

% uses: ed.m ssq.m lsi.m fastnnls.m
% Marieke Timmerman, 01-09-2011

dys=[];Cdiag=[];Zsol=[];Csol=[];
conv=1e-6;
[J,n]=size(S);
K=n/J;
U=orth(null(ones(L,1)'));   % (L x L-1 basis matrix)


Cbasis=[eye(L);zeros(K-L,L)];


if length(Q)==1, Q=ones(L,1)*Q; end;
Qcum=zeros(L+1,1);
for l=1:L, Qcum(l+1,1)=sum(Q(1:l,1));end
Z=zeros(J,sum(Q));

for ana=1:nanal,

    disp(sprintf('Anlysis %i of %i',ana,nanal))
    %Random initialisation of C (random L rows of C fixed at row of eye(L)
    %without replacement, and zeros otherwise) 
    Fixk=[];
    if  numel(Cstart)==1,
        Ri=randperm(K);
        C=Cbasis(Ri,:);
        Fixk=zeros(K,1);
    else
        Fixk=sum(Cstart,2);
        [iF,jF]=find(Cstart);
        C=Cstart;
        if sum(Fixk)<L,
            Ri=randperm(K)';
            [iR,jR]=find(Ri==jF);
            Ri(iR,1)=Ri(iF,1);
            Ri(iF,1)=jF;
            C=Cbasis(Ri,:);
        end
     end
    
    %Initialise Zar
    for l=1:L,
         k=find(C(:,l));
         Stemp=S(:,(k-1)*J+1:k*J);
        [KK,LL]=ed(Stemp);
    
        Z(:,(Qcum(l,1)+1):(Qcum(l+1,1)))=KK(:,1:Q(l))*(LL(1:Q(l),1:Q(l)).^0.5);
    end
    
   f=0;  %compute function value
   for l=1:L, Cdiag(:,(Qcum(l,1)+1):(Qcum(l+1,1)))=C(:,l)*ones(1,Q(l)); end %Full diagonal C
   for k=1:K,
      f=f+ ssq(S(:,(k-1)*J+1:k*J) - Z*diag(Cdiag(k,:)) *Z');
   end;
   
   
fold=f+2*conv*f;
iter=1;   
f1=0;f2=0;stat=[];

 while fold-f>conv*f %| iter<1000
  fold=f;
  stat(iter,:) = [iter f1 f2 f];

  iter=iter+1;
    
  %Update %May yield dysmonotony, Sigma_k||S_k - Zsup2 Cdiag Zsup1'||, only
  %Zsup1 is being updated. Zsup2 is being replaced by Zsup1 (which needs
  %not be optimal)
  B=zeros(sum(Q),sum(Q));
  BB=zeros(J,sum(Q));
  for k=1:K, B=B+diag(Cdiag(k,:))*Z'*Z*diag(Cdiag(k,:)); end;
  for k=1:K,
    BB=BB+S(:,(k-1)*J+1:k*J)*Z*diag(Cdiag(k,:));
  end;
  Z=BB/B;
  
  f=0;  %compute function value
   for l=1:L, Cdiag(:,(Qcum(l,1)+1):(Qcum(l+1,1)))=C(:,l)*ones(1,Q(l)); end %Full diagonal C
   for k=1:K,
      f=f+ ssq(S(:,(k-1)*J+1:k*J) - Z*diag(Cdiag(k,:)) *Z');
   end; f1=f;
  
  %Update C, for all zero elements of Cstart
  
  
  VecZZ=[];
  for l=1:L,
  Zl=Z(:,(Qcum(l,1)+1):(Qcum(l+1,1)))*(Z(:,(Qcum(l,1)+1):(Qcum(l+1,1))))';
  VecZZ=[VecZZ Zl(:)];  
  end
  VecZZhat=inv(VecZZ'*VecZZ)*VecZZ';
  VecZVecZ=(VecZZ'*VecZZ);
  
  
   for k=1:K,
       if Fixk(k)==0, %update of non-fixed rows of C only
        Sk=S(:,(k-1)*J+1:k*J);
              
            if conC==1, 
                %C(k,:)=(VecZZhat*Sk(:))';%unconstrained update of row k of C 
                C(k,:)=fastnnls(VecZVecZ,VecZZ'*Sk(:)); %constrained update of row k of C, subject to C(k,:)>=0                
            end
        
            if conC==2, %constrained update of row k of C, 1'C'=1, C>=0
              f=Sk(:)- VecZZ*ones(L,1)*(1/L);
              E=VecZZ*U;
              G=U;
              h=-ones(L,1)/L;
              x=lsi(E,f,G,h);  %  minimizes || Ex - f || subject to Gx >= h
              p=U*x+ones(L,1)/L;
              C(k,:)=p';
        end
         
       end
        
   end
   
  
  for l=1:L, Cdiag(:,(Qcum(l,1)+1):(Qcum(l+1,1)))=C(:,l)*ones(1,Q(l)); end %Full diagonal C
  
    
     
   f=0;  %compute function value
     for l=1:L, Cdiag(:,(Qcum(l,1)+1):(Qcum(l+1,1)))=C(:,l)*ones(1,Q(l)); end %Full diagonal C
   for k=1:K,
      f=f+ ssq(S(:,(k-1)*J+1:k*J) - Z*diag(Cdiag(k,:)) *Z');
   end;
   f2=f;
 
 
    if f>fold
   %disp(' dysmonotony ')
    dys=[dys; ana];
   end

 end
 
 func(ana)=f;
 Zsol(:,ana)=Z(:);
 Csol(:,ana)=C(:);
 
end
 
%disp(' Minimal function value is ');
[f,mi]=min(func);
C(:)=Csol(:,mi);
Z(:)=Zsol(:,mi);

 %rescale C such that mean(C)=1
  Cr=C./(ones(K,1)*mean(C));
  
  Zr=[];
  for l=1:L,
       Zt=Z(:,(Qcum(l,1)+1):(Qcum(l+1,1)));
       Zr=[Zr Zt*(sqrt(mean(C(:,l))))];
  end
  Z=Zr;C=Cr;
  
   f=0;  %compute function value
     for l=1:L, Cdiag(:,(Qcum(l,1)+1):(Qcum(l+1,1)))=C(:,l)*ones(1,Q(l)); end %Full diagonal C
   for k=1:K,
      f=f+ ssq(S(:,(k-1)*J+1:k*J) - Z*diag(Cdiag(k,:)) *Z');
   end;
   f2=f;
   
  fp=100*(ssq(S)-f)/ssq(S);
  
%   
% %Imponi constraint ortogonalità   su Zl
% 
% 
% sumat = [];
% for c = 1 : length(Q)
%     
%     R = Q(c);
%     sumat = [sumat R*ones(1,R);];
%     
% end
% 
% ZO = [];
% for c = 1 : length(Q)
% 
%     
%     z = Z(:, find(sumat==c));
%     
%     %Fai svd
%     [U,S,V] = svd(z);
%     
%     %Ridefinisco Z
%     zo = U*S;
%     clear U S V z
%     
%     ZO = [ZO zo];
%     clear zo
%     
% end

  
  
 
