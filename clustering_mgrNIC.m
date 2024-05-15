function  [P,S,B,F,H,err1]=clustering_NIC(X,M,k,k2,alpha,lamta,beta)
%%%%The model----------------
% min{P[l],S[l],S} sumTr(P[l]'X[l]L[l]X[l]'P[l])+alpha
% sum||S[l]||^2+sum(\|S[l]-B[l]F\|^2+\|M[l]-H[l]F\|^2)
% You only need to provide the above three inputs.
% Notation:
% X[l] ... (dl x n) gene-cell data 
% n  ... number of cells            
% k ... number of features (the number of cluster)
% k2 ...the number of cluster)
% alpha ... Regularization parameter
% iter ... The maximum number of iterations
%%%
% Coder Wenming Wu Email: wenmingwu55 at 163.com
%
% version --September/2021
    iter = 100;
    err1 = zeros(iter,1);err2 = zeros(iter,1);
    m = length(X);
    [d{1},n] = size(X{1});
    [d{2},n] = size(X{2});
%     P{1}=rand(d{1},k);
%     P{2}=rand(d{2},k);
    for l=1:m
%         options = [];
%         option.Metric = 'Cosine';%娆у紡璺濈
%         options.NeighborMode = 'KNN';%KNN
%         options.k = 5;%5涓渶杩戦偦
%         options.WeightMode = 'Cosine';%鏉冨?间负0鎴?1 鍙洿鎹负 'HeatKernel', 'Cosine' 
% 
%         S{l} = constructW(X{l}',options);%%%%%%X涓鸿緭鍏ョ煩闃?
        S{l} = M{l};
        Ds{l} = diag(sum(S{l},2));
	    Ls{l} = Ds{l}-S{l};
        
        J{l}=S{l};
        [U,V,D] = svds(S{l},k2);
        B{l} = abs(U*sqrt(V));
        
        fff{l} = abs(sqrt(V)*D');
        [U1,V1,D1] = svds(M{l},k2);
        H{l} = abs(U1*sqrt(V1));
        
        Y = X{l}*Ls{l}*X{l}';
        [Vy,Dy] = eig(Y);
        dy = diag(Dy);
        [t,v] = sort(dy,'ascend');
        P{l} = Vy(:,v(1:k));
    end
    F = (fff{1}+fff{2})/2;
for o = 1:iter
%%%%%--------------Update variables P{l} by iteration------------
    for l = 1:m
        %%========P{l}=========
        Y = X{l}*Ls{l}*X{l}';
        [Vy,Dy] = eig(Y); %特征值分解
        dy = diag(Dy); %Y的特征值
        [t,v] = sort(dy,'ascend'); %升序排列
        P{l} = Vy(:,v(1:k)); %选取特征值最低的k个特征向量
        
    %%%%%====================S{l}======================
    wW = P{l};xX = X{l};
        for i = 1:n
            for j = 1:n
                f(i,j) = norm(wW'*xX(:,i)-wW'*xX(:,j),'fro')^2;
            end
        end
    S{l} = S{l}.*((B{l}*F)./((alpha+1)*S{l}+f'+eps));
    S{l} = mapminmax(S{l}, 0, 1);
    end
     %%%%%--------------Update variables B{l},H{l},F by iteration------------

    for l = 1:m
        ssS{l} = (S{l}+S{l}')/2;
       B{l} = B{l}.*((ssS{l}*F'+lamta*J{l}*B{l})./(B{l}*F*F'+lamta*Ds{l}*B{l}+eps)); 
    end
   % for l = 1:m
    %    ssS{2} = (S{2}+S{2}')/2;
    %   B{2} = B{2}.*((ssS{2}*F'+2*J{2}'*B{1})./(B{2}*F*F'+2*D{2}'*B{1})); 
    %end
  %for l = 1:m
      %  ssS{1} = (S{1}+S{1}')/2;
     %  B{1} = B{1}.*((ssS{1}*F'+2*J{1}*B{2})./(B{1}*F*F'+2*D{1}*B{2})); 
  %end


   for l = 1:m
       H{l} = H{l}.*((M{l}*F'+beta*J{l}*B{l})./(H{l}*F*F'+beta*Ds{l}*H{l}+eps)); 
    end
    
    ff1 = zeros(k2,n); ff2=zeros(k2,n);
     [t1, t2] = size(F*F');
    E = eye(t1, t2);
    for l = 1:m
        ssS{l} = (S{l}+S{l}')/2;
        ff1 = ff1+B{l}'*ssS{l}+H{l}'*M{l};
        ff2 = ff2+B{l}'*B{l}*F+H{l}'*H{l}*F+eps;
    end
    F=F.*((ff1)./(ff2));

%%%%%%%%%%%%%%%-------------Error-----------------------
    ee =(norm(S{1}-B{1}*F,'fro')+norm(S{2}-B{2}*F,'fro')+norm(M{1}-H{1}*F,'fro'));
    err1(o,1)=ee;
     disp([' 迭代次数 ' num2str(o) ' temp1 ' num2str(ee)]);

    if ee < 1.000000000000000e-15
        break;
    else 
        P = P;
        S = S;
    end
end

end