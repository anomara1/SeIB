% Script 2 for block update
% Complete solution for sparse code update
% ReCheck with Adaptive Lambda (1 Sep 2020) & Alpha consistent with
% Manuscript
function [Dk,Xk,Err,Diff,AvgB] = Robust_BlockDL_DR2_ANO_128(Y,Dini,Xk,noIt,lambda,alpha,delta,verb)
    if ~exist('verb','var');    verb = 0;   end
    % function to calculate weights
    wei = @(aa,e) (exp(-aa./2.*(e.^2)));
    BlStruct = ones(1,size(Dini,2));
    nB = max(BlStruct);
    [n,N] = size(Y);    len = n*N;
    K = size(Dini,2);   fact = max(1,ceil(delta*n));
    [Diff,Err,AvgB] = deal(zeros(1,noIt));
    [Dk,Dold] = deal(Dini);
    %Dk = Dini;
    %Dold = Dini;
    for it = 1:noIt 
        % Using delta info, finding the optimal alpha
        % only once every iteration
        if (it == 1 || ~alpha)
            W = ones(n,N);
        elseif alpha 
            EE = Y - Dk*Xk;
            EE_ = sort(EE.^2,'descend');
            tt = EE_(fact,:);
            alpha_ = alpha*(-2*log(0.5)./tt);
            W = wei(alpha_,EE);
        end
        % block update begins
        for bb = 1:nB     % only a single block to update for digit recognition  
            Indd = BlStruct == bb;   BSize = nnz(Indd); 
            D = Dk(:,Indd);    % Current Updating Block
            X = Xk(Indd,:);            
            
            Ek = Y - Dk*Xk + D*X; 
            if it == 1 
                LAMBDAS = lambda.*ones(1,N).*sqrt(BSize);
            else
                LAMBDAS = lambda./max(diag(X'*X),0.01).*sqrt(BSize);
            end
            dj = 0;  jj = 0;    
            while norm(dj - D,'fro') > 10^-2 && jj < 2
%                 fprintf('IIter: %d, Block: %d, Conv: %0.4f\n',jj+1,bb,norm(dj - D,'fro'));
                dj = D;     jj = jj + 1;    %[k,jj]
                % Approximation solution
%                 gg = dj'*(W.*Ek);
%                 gg2 = sqrt(diag(gg'*gg)');
%                 X = max(0,1-lambda./gg2).*gg;
                % Complete/Exact solution
                for NN = 1:N
                    [U,S] = eig(D'*(W(:,NN).*D));
                    U = real(U);
                    S = real(S);
                    y_i = U'*D'*(W(:,NN).*Ek(:,NN));
                     LAMBDAS = lambda./sqrt(y_i'*y_i)*sqrt(BSize);
                    if norm(y_i) <= LAMBDAS%(NN)
                        X(:,NN) = 0;
                    else
                        z_norm = LAMBDAS^2*((y_i'*y_i) / (y_i'*S*y_i));%LAMBDAS(NN)^2*((y_i'*y_i) / (y_i'*S*y_i));
                        z = (S + LAMBDAS/z_norm*eye(BSize)) \ y_i;%(S + LAMBDAS(NN)/z_norm*eye(BSize)) \ y_i;
                        X(:,NN) = U*real(z);
                    end
                end
                
                temp = (Ek.*W)*X';
                for nn = 1:n
%                     D(nn,:) = Ek(nn,:)*(diag(W(nn,:)))*X' / (X*(diag(W(nn,:)))*X');
                    D(nn,:) = temp(nn,:) / ((X.*W(nn,:))*X');
                end
                D = normc(D);
                stop_ = 0;
            end
            
            Dk(:,Indd) = D;        Xk(Indd,:) = X;
%             if verb
%                 tt = norm(Dk(:,Indd) - Dold(:,Indd),'fro');
%                 fprintf('Iteration: %d, Block: %d, Conv: %0.3f\n',it,bb,tt);
%             end            
        end
        Diff(it) = norm(Dk - Dold,'fro')/sqrt(K);
        AvgB(it) = nnz(Xk)/(BSize*N);
        if verb   
            Err(it) = DispError(YOrig,Dk,Xk,0);
            fprintf('Iteration: %2d, Error: %0.3f, DConv: %0.3f, BSpar: %0.3f\n',...
                it,Err(it),Diff(it),AvgB(it));            
        end
        if Diff(it) < 0.05 %|| AvgB(it) < 1.9
            break;
        end
        Dold = Dk;
    end     
    stop = 0;
end