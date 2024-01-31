% we need to know 26, 13,9,7,5 atoms
function Start_PAP7_SISM_v_4_EXP2_get_Dics(file_)
for jj = file_:file_
    imExt ='.jpg';
    fName = strcat('(',int2str(jj),')',imExt);
    img_pixels = imread(fName);
    img_double_ = im2double(img_pixels);
    No_ATOMS = [];
    img_double = img_double_(:,:,1);
    for i = 0.05:0.05:0.3
        No_ATOMS = [No_ATOMS ceil(size(img_double,2)*i)];
    end
    No_ATOMS = [26, 13, 9,7,5];
    %RBDL 2023
    alg = 'bksvd';%'rbdl';%'odl';'ksvd';'mod';'bssdl';'bksvd'
    switch alg

        case 'bksvd'
            for j = 1:1:5

                RBDL_Lambda = 0.15;
                delta = [10,20,30,40,50];
                xinit = [26,13,9,7,5];
                %Dini = Dict_Learn_DR(training_feats,H_train,DSize,noIt4ini,K,Algo,RBDL_Lambda)
                Dic{j} = Dict_Learn_DR_ANO(img_double,1,No_ATOMS(j),5,5,'bksvd',RBDL_Lambda,delta(j));
                x = 0;
            end
            save(strcat('(',int2str(jj),')_BKSVD'),'Dic');


        case 'bssdl'
            for j = 1:1:5

                RBDL_Lambda = 0.15;
                delta = [10,20,30,40,50];
                xinit = [26,13,9,7,5];
                %Dini = Dict_Learn_DR(training_feats,H_train,DSize,noIt4ini,K,Algo,RBDL_Lambda)
                Dic{j} = Dict_Learn_DR_ANO(img_double,1,No_ATOMS(j),5,5,'bssdl',RBDL_Lambda,delta(j));
                x = 0;
            end
            save(strcat('(',int2str(jj),')_BSSDL'),'Dic');
        
        case 'rbdl'
            for j = 1:1:5

                RBDL_Lambda = 0.15;
                delta = [10,20,30,40,50];
                xinit = [26,13,9,7,5];
                %Dini = Dict_Learn_DR(training_feats,H_train,DSize,noIt4ini,K,Algo,RBDL_Lambda)
                Dic{j} = Dict_Learn_DR_ANO(img_double,1,No_ATOMS(j),5,5,'rbdl',RBDL_Lambda,delta(j));
                x = 0;
            end
            save(strcat('(',int2str(jj),')_RBDL'),'Dic');
        case 'odl'
            % ODL
            for j = 1:1:5
                lambda =0.1;
                [Dic{j}, X] = ODL(img_double, No_ATOMS(j), lambda);%, opts, method)


            end
            save(strcat('(',int2str(jj),')_ODL'),'Dic');

        case 'ksvd'
            %1    KSVD
            for j = 1:1:5

                Dic{j} = getDic_new(size(img_double,1),1,img_double,No_ATOMS(j));
            end
            save(strcat('(',int2str(jj),')_KSVD'),'Dic')
        case 'mod'
            for j = 1:1:5

                %2    MOD
                Dic{j} = getDic_new(size(img_double,1),2,img_double,No_ATOMS(j));
            end
            save(strcat('(',int2str(jj),')_MOD'),'Dic')
    end

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

