% we will use 50 images (321001 ->321050) to get 256*256 dictionaries
function Start_PAP7_SISM_v_4_EXP4_get_Dics()
imageSize = 128;
img_double = [];
for jj = 310001:1:310050%321050
    imExt ='.jpg';
    fName = strcat('(',int2str(jj),')',imExt);
    img_pixels = imread(fName);
    img_double_ = im2double(img_pixels);
    No_ATOMS = [];
    img_double = [img_double,img_double_(:,:,1)];
end
a = ones(imageSize,1);
dc_atom = normc(a);
for z=1:1:size(img_double,2)
    img_double(:,z)= img_double(:,z)-dot(dc_atom,img_double(:,z))*dc_atom;
    cccc = 0;
end
No_ATOMS = 128;

%RBDL 2023
alg = 'rksvd';%'rbdl';%'odl';'ksvd';'mod';'bssdl';'bksvd'
switch alg
          case 'rksvd'
            for j = 1:1:1
                
                RBDL_Lambda = 0.15;
                mydelta = [4,8,16,32,64];
                xinit = [32, 16, 8,4,2];
                %Dini = Dict_Learn_DR(training_feats,H_train,DSize,noIt4ini,K,Algo,RBDL_Lambda)
                %Dic{j} = Dict_Learn_DR_ANO_128(img_double,1,No_ATOMS(j),5,5,'bksvd',RBDL_Lambda,mydelta(j));
                p = 128;         % problem dimension
                n = No_ATOMS;         % number of atoms in the dictionary
                m = 50*128;        % number of training signals
                s = 4;          % sparsity constraint
                reg = 0.1;      % regularization
                vanish = 1;     % regularization vanishing factor
                regstop = 31;   % cancel regularization term starting with this iteration
                iters = 50;     % DL iterations
                %%-------------------------------------------------------------------------
                % Path to OMP implementation by Ron Rubinstein
                % (http://www.cs.technion.ac.il/~ronrubin/Software/ompbox10.zip)
                addpath('ompbox10');
                %%-------------------------------------------------------------------------
                %Y = randn(p,m);
                Y = img_double;
                D0 = normc(randn(p,n));
                
                params = {'reg', reg, 'vanish', vanish, 'regstop', regstop};
                [Df, Xf, errs] = DL(Y, D0, s, iters, params{:});
                %plot(1:iters, errs);
                Dic = Df;
                
                
            end
          save(strcat('EXP4_RKSVD_DIC_',int2str(No_ATOMS)),'Dic');            
          stop = 0;


    case 'rdl'
        for j =1:1:1
            K = 128;
            Ys = img_double;%rand(256,256);
            Y = Ys;
            [n,N] = size(Ys);
            imgsize = imageSize;
            %% Learning process
            Dini = normc(abs(randn(n,K)));
            Xini = zeros(K,N);
            % Xini = pinv(Dini)*Ys;    %zeros(K,N);
            alpha = 1;  lambda = 0.02;   delta = 0.1;    verb = 1;
            noIt = 40;  tic;
            [Dic,X] = OWN_RDL(Ys,Dini,Xini,noIt,lambda,alpha,delta,verb);
        end
        save(strcat('EXP4_RBDL_DIC_',int2str(No_ATOMS)),'Dic');
        stop = 0;
    case 'rdl1'
        for j = 1:1:1
            
            RBDL_Lambda = 0.15;
            delta = 50;
            xinit = [26,13,9,7,5];
            %Dini = Dict_Learn_DR(training_feats,H_train,DSize,noIt4ini,K,Algo,RBDL_Lambda)
            Dic = Dict_Learn_DR_ANO(img_double,1,No_ATOMS(j),5,5,'rdl1',RBDL_Lambda,delta(j));
            x = 0;
        end
        save(strcat('EXP4_RDL1_DIC_',int2str(No_ATOMS)),'Dic');

    
    case 'bksvd'
        for j = 1:1:1
            
            RBDL_Lambda = 0.15;
            delta = 50;
            xinit = [26,13,9,7,5];
            %Dini = Dict_Learn_DR(training_feats,H_train,DSize,noIt4ini,K,Algo,RBDL_Lambda)
            Dic = Dict_Learn_DR_ANO(img_double,1,No_ATOMS(j),5,5,'bksvd',RBDL_Lambda,delta(j));
            x = 0;
        end
        save(strcat('EXP4_BKSVD_DIC_',int2str(No_ATOMS)),'Dic');
        
        
    case 'bssdl'
        for j = 1:1:1
            
            RBDL_Lambda = 0.15;
            delta = 50;
            xinit = [26,13,9,7,5];
            %Dini = Dict_Learn_DR(training_feats,H_train,DSize,noIt4ini,K,Algo,RBDL_Lambda)
            Dic = Dict_Learn_DR_ANO(img_double,1,No_ATOMS(j),5,5,'bssdl',RBDL_Lambda,delta(j));
            x = 0;
        end
        save(strcat('EXP4_BSSDL_DIC_',int2str(No_ATOMS)),'Dic');
        
    case 'rbdl'
        for j = 1:1:1
            
            RBDL_Lambda = 0.001;
            delta = 4;
            xinit = [26,13,9,7,5];
            %Dini = Dict_Learn_DR(training_feats,H_train,DSize,noIt4ini,K,Algo,RBDL_Lambda)
            Dic = Dict_Learn_DR_ANO(img_double,1,No_ATOMS(j),5,5,'rbdl',RBDL_Lambda,delta(j));
            x = 0;
        end
        save(strcat('EXP4_RBDL_DIC_',int2str(No_ATOMS)),'Dic');
        xx = 0;
    case 'odl'
        % ODL
        for j = 1:1:1
            lambda =0.1;
            [Dic, X] = ODL(img_double, No_ATOMS(j), lambda);%, opts, method)
            
            
        end
        save(strcat('EXP4_ODL_DIC_',int2str(No_ATOMS)),'Dic');
        
    case 'ksvd'
        %1    KSVD
        for j = 1:1:1
            
            Dic = getDic_new(size(img_double,1),1,img_double,No_ATOMS(j));
        end
        save(strcat('EXP4_KSVD_DIC_',int2str(No_ATOMS)),'Dic');
    case 'mod'
        for j = 1:1:1
            
            %2    MOD
            Dic = getDic_new(size(img_double,1),2,img_double,No_ATOMS(j));
        end
        save(strcat('EXP4_MOD_DIC_',int2str(No_ATOMS)),'Dic');
end

end




function [Dic] = getDic_new(vecLength,LDIC_TYPE,img_double,No_ATOMS)
if LDIC_TYPE == 1 | LDIC_TYPE == 2
    Dic = getLearnedAtoms_new(LDIC_TYPE,img_double,No_ATOMS);
    switch LDIC_TYPE
        case 'ksvd'
            %dd = load('PAP7_KSVD_LEARNED_256_512.mat');
        case 'mod'
            %dd = load('PAP7_MOD_LEARNED_256_512.mat');
        case 'randomized'
            %dd = load('PAP7_RAND_LEARNED_256_512.mat');
        case 'modortho'
            %dd = load('PAP7_MODO_LEARNED_256_512.mat');
    end
    %dd = load('MOD_LEARNED_256_512.mat');
    %Dic = dd.Dic;
else
    
    [Dic_r,nbvect] = wmpdictionary(vecLength,'lstcpt',{{'wphaar',2},'dct',});%{'wphaar',2},,{'wpsym4',2}}'dct'
    aa = ones(vecLength,1);
    Dic = full(Dic_r); %[aa_n full(Dic_r)];
    
end
end
function Dic = getLearnedAtoms_new(DIC_LEAR_ALG,img_double_,No_ATOMS)
%file = 220001:1:220015;
k = 1;
for i = 1:1:1
    Vector = img_double_;
    
end
if DIC_LEAR_ALG == 1
    options.learning_method = 'ksvd';
    options.niter_learning = 40;
    options.K = No_ATOMS;
    options.nbr_max_atoms = 30;
    [Dic1,X,E] = perform_dictionary_learning(Vector, options);
    %[Dic2,nbvect] = wmpdictionary(size(Vector,1),'lstcpt',{{'wphaar',2},'dct',});
    %Dic = [full(Dic2), Dic1];
    Dic = Dic1;
    
    
elseif DIC_LEAR_ALG == 2
    options.learning_method = 'mod';
    options.niter_learning = 40;
    options.K = No_ATOMS;
    options.nbr_max_atoms = 30;
    [Dic1,X,E] = perform_dictionary_learning(Vector, options);
    %[Dic2,nbvect] = wmpdictionary(size(Vector,1),'lstcpt',{{'wphaar',2},'dct',});
    %Dic = [full(Dic2), Dic1];
    Dic = Dic1;
end


end
function [batchMatrix, batchPerRow,batchPerCol] = getBatches_(singleColorLevelMatrix,batchSize)

batchPerRow = size(singleColorLevelMatrix,1)/batchSize;
batchPerCol = size(singleColorLevelMatrix,2)/batchSize;

kr = 0:1:batchPerRow;
kc = 0:1:batchPerCol;
for i=1:1:floor(batchPerRow)
    for j=1:1:floor(batchPerCol)
        batchMatrix(i,j,:,:)= singleColorLevelMatrix(kr(i)*batchSize+1:kr(i)*batchSize+batchSize,kc(j)*batchSize+1:kc(j)*batchSize+batchSize);
    end
    
end

end

