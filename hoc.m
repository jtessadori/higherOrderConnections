classdef hoc < handle
    % January 2020, Jacopo Tessadori
    % The idea is to verify whether higher order connections are useful in
    % discriminating between healthy and RR subjects: in theory, if AB, BC
    % and CA are often active togther in a particular subject but one of
    % the connections is systematically missing in the average of healthy
    % controls, it is likely that that subject is re-routing a connection
    % to compensate for a missing link
    
    properties
        dataPath;
        FC;
        lbls;
        nROIs;
        feats;
        tripletIdx;
        tripletCoords;
        featIdx;
        relFeatIdx;
        subjAge;
    end
    
    methods
        function this=hoc(path)
            % Constructor for hoc class.
            % Only argument (required) is absolut path of data file.
            this.dataPath=path;
        end
        
        function loadDataHC_RR(this)
            % Recovers dynamical functional connectivity matrices as well
            % as lbls from the file specified as dataPath. Task is
            % classification between HC and RR subjects
            inData=load(this.dataPath);
            
            % Old datasets may have data split in different ways
            if isfield(inData,'CM')
                if length(inData.CM)==28
                    % This is a very specific test dataset, define lbls
                    % here
                    this.lbls=zeros(1,28);
                    this.lbls(15:end)=1;
                    this.FC=inData.CM;
                end
            else
                % Define labels and remove unlabeled data
                this.lbls=cat(1,zeros(length(inData.idHC),1),ones(length(inData.idRR),1));
                inData.cm_all_subj_corr=inData.cm_all_subj_corr(union(inData.idHC,inData.idRR));
                this.FC=inData.cm_all_subj_corr;
            end
        end
        
        function loadDataHC_MS(this)
            % Recovers dynamical functional connectivity matrices as well
            % as lbls from the file specified as dataPath. Task is
            % classification between HC and union of SPMS and PPMS subjects
            inData=load(this.dataPath);
            
            % Old datasets may have data split in different ways
            if isfield(inData,'CM')
                if length(inData.CM)==28
                    % This is a very specific test dataset, define lbls
                    % here
                    this.lbls=zeros(1,28);
                    this.lbls(15:end)=1;
                    this.FC=inData.CM;
                end
            else
                % Define labels and remove unlabeled data
                try
%                     HCidx=find(inData.HC_All_lab==0);
%                     MSidx=find(ismember(inData.HC_All_lab,[2,3]));
%                     this.lbls=cat(1,zeros(length(HCidx),1),ones(length(MSidx),1));
%                     inData.cm_all_subj_corr=inData.cm_all_subj_corr(union(HCidx,MSidx));
%                     this.FC=inData.cm_all_subj_corr;
                    featMat=cat(3,inData.HC,inData.MS);
                    featMat=permute(featMat,[3,1,2]);
                    this.FC=num2cell(featMat,[2,3]);
                    this.FC=cellfun(@(x)squeeze(x),this.FC,'UniformOutput',false);
                    this.lbls=cat(1,zeros(size(inData.HC,3),1),ones(size(inData.MS,3),1));
                catch
                    warning('Data path file does not contain lbl information. Trying to recover it from xls file');
                    num=xlsread(inData.xlsFile);
                    relLbls=num(:,8);
                    relLbls(isnan(relLbls))=0;
                    this.FC=inData.FCs(relLbls~=1); % Excluding RR and bening
                    relLbls(relLbls==1)=[];
                    this.lbls=ismember(relLbls,[2,3]);
                end
            end
        end
        
        function loadData(this)
            % Recovers dynamical functional connectivity matrices as well
            % as lbls from the file specified as dataPath
            inData=load(this.dataPath);
            
            % Old datasets may have data split in different ways
            if isfield(inData,'CM')
                if length(inData.CM)==28
                    % This is a very specific test dataset, define lbls
                    % here
                    this.lbls=zeros(1,28);
                    this.lbls(15:end)=1;
                    this.FC=inData.CM;
                end
            else
                % Define labels and remove unlabeled data
                this.lbls=cat(1,zeros(length(inData.idHC),1),ones(length(inData.idRR),1));
                inData.cm_all_subj_corr=inData.cm_all_subj_corr(union(inData.idHC,inData.idRR));
                this.FC=inData.cm_all_subj_corr;
            end
        end
        
        function loadDataImprovement(this)
            % Recovers static functional connectivity matrices as well
            % as lbls from the file specified as dataPath. Task is
            % classification between stable/improving subjects and
            % worsening ones
            inData=load(this.dataPath);
            
            this.lbls=inData.lbls;
            this.FC=inData.FCs;
        end
        
        function computeOutlyingConnections(this)
            % Initialize large variables
            this.tripletIdx=zeros(this.nROIs^3,length(this.lbls),3,'single');
            this.feats=zeros(this.nROIs^3,length(this.lbls),'single');
            
            for currSubj=1:length(this.lbls)
                if ~isnan(this.lbls(currSubj))
                    % Compute three-way coactivations for current subject
                    sData=this.FC{currSubj};
                    sa=repmat(permute(sData,[4 1 2 3]),this.nROIs,1,1,1);
                    sb=repmat(permute(sData,[1 4 2 3]),1,this.nROIs,1,1);
                    sc=repmat(permute(sData,[1 2 4 3]),1,1,this.nROIs,1);
                    sHOCavg=squeeze(mean(sa.*sb.*sc.*(sa>0).*(sb>0).*(sc>0),4));
                    
                    % Fix threshold for "often active together" triplets
                    th=.125;
                    
                    % Recover data for healthy subjects, excluding current one, if
                    % relevant
                    hSubjsIdx=find(this.lbls==0);
                    hData=cat(3,this.FC{setdiff(hSubjsIdx,currSubj)});
                    
                    % Define average 3-way coactivations above a
                    % certain threshold as "often active together"
                    actIdx=find(sHOCavg>th);
                    [idx1,idx2,idx3]=ind2sub(size(sHOCavg),actIdx);
                    
                    % Remove all entries with a duplicate coord
                    actSubs=[idx1,idx2,idx3];
                    actSubs(cellfun(@(x)length(unique(x))~=3,num2cell(actSubs,2)),:)=[]; % Removing entries with duplicate coords
                    actSubs=sort(unique(sort(actSubs,2),'rows'),2); % Can do this, as order is irrelevant (matrix is VERY symmetric)
                    
                    % Determine average FC for current subject and healthy
                    % controls
                    hMat=squeeze(median(hData,3));
                    sMat=squeeze(median(sData,3));
                    
                    % For each strong triplet, save difference between
                    % strongest and second strongest connections on healthy
                    % subjects
                    for currSub=1:size(actSubs,1)
                        % Recover involved connection strengths for current
                        % subject and corresponding strengths for healthy ones
                        hVals=[hMat(actSubs(currSub,1),actSubs(currSub,2)),hMat(actSubs(currSub,1),actSubs(currSub,3)),hMat(actSubs(currSub,2),actSubs(currSub,3))];
                        sVals=[sMat(actSubs(currSub,1),actSubs(currSub,2)),sMat(actSubs(currSub,1),actSubs(currSub,3)),sMat(actSubs(currSub,2),actSubs(currSub,3))];
                        
                        %                     % Reorder indexes so that strongest connection is first
                        %                     switch (sVals==max(sVals))*(1:3)' % Case 1 requires no changes
                        %                         case 2 % i.e. 1-3
                        %                             actSubs(currSub,:)=actSubs(currSub,[1 3 2]);
                        %                         case 3 % i.e. 2-3
                        %                             actSubs(currSub,:)=actSubs(currSub,[2 3 1]);
                        %                     end
                        this.feats(actSubs(currSub,:)*90.^(2:-1:0)',currSubj)=median(sVals)-median(hVals);
                        this.tripletIdx(actSubs(currSub,:)*90.^(2:-1:0)',currSubj,:)=actSubs(currSub,:);
                    end
                    
                    % Print progress
                    fprintf('%d/%d subj lbl: %d, # 3-way connections: %d\n',currSubj,length(this.lbls),this.lbls(currSubj),sum(this.feats(:,currSubj)~=0));
                end
            end
            
            % Remove feature that are zeros for all subjects
            this.tripletIdx=this.tripletIdx(sum(this.feats~=0,2)>0,:,:);
            this.feats=this.feats(sum(this.feats~=0,2)>0,:);
        end
        
        function plotDistr(this)
            featsMask=this.feats~=0;
            hMaskDistr=sum(featsMask(:,this.lbls==0),2);
            hMaskDistr=sort(hMaskDistr,'descend')/sum(hMaskDistr);
            semilogx(hMaskDistr)
            sMaskDistr=sum(featsMask(:,this.lbls==1),2);
            sMaskDistr=sort(sMaskDistr,'descend')/sum(sMaskDistr);hold on
            plot(sMaskDistr)
            set(findall(gca,'type','line'),'linewidth',1.5);
            xlabel('Sorted feature #');
            ylabel('Fraction of total triplets');
            title('Distribution of "strong" triplets in HC and RR population');
            legend({'HC','RR'});
        end
        
        function showFeatureMask(this)
            figure;
            imagesc([this.feats(:,this.lbls==0)>.2,(this.feats(:,this.lbls==1)>.2)*2]')
            xlabel('Features');
            ylabel('Subjects');
            title('Thresholded representation of "strong" triplets');
        end
        
        function BAcc=crossvalClassify(this,nReps,nFeats)
            if ~exist('nReps','var')
                nReps=100;
            end
            if ~exist('nFeats','var')
                nFeats=80;
            end
            featSelFun=@(x,y)rankfeatures(x,y,'CRITERION','ttest','CCWEIGHTING',1,'NUMBEROFINDICES',nFeats);
%             featSelFun=@(x,y)1:length(x);
            pdf=histcounts(this.lbls(~isnan(this.lbls)),length(unique(this.lbls(~isnan(this.lbls)))),'Normalization','probability');
            costMat=[0,1/pdf(1);1/pdf(2),0];
%             classifierFun=@(x,y)fitcsvm(x,y,'Cost',costMat,'KernelFunction','polynomial','PolynomialOrder',2,'KernelScale','auto','BoxConstraint',1,'Standardize',true);
            classifierFun=@(x,y)fitcknn(x,y,'Cost',costMat,'Distance','euclidean','DistanceWeight','inverse','NumNeighbors',3,'Standardize',true);
            BAcc=zeros(nReps,1);
            for currRep=1:nReps
                lblsEst=crossVal(featSelFun,classifierFun,this.feats,this.lbls);
                BAcc(currRep)=computeBAcc(this.lbls,lblsEst);
                fprintf('%d/%d, last BAcc: %0.2f, mean BAcc: %0.2f\n',currRep,nReps,BAcc(currRep),mean(BAcc(1:currRep)));
            end
        end
        
        function [relFeatsCoordsLog,relFeatsIdxLog]=recoverMostRelevantTriplets(this,nFeats,nReps)
            nFolds=5;
            relFeatsCoordsLog=cell(0);
            relFeatsIdxLog=cell(0);
            for currRep=1:nReps
                % Determine folds
                foldID=ceil(linspace(.1,nFolds-.1,length(this.lbls)));
                foldID=foldID(randperm(length(foldID)));
                
                % Cycle through folds and perform feature selection
                for currFold=1:nFolds
                    testIdx=find(foldID==currFold);
                    trainIdx=setdiff(1:length(this.lbls),testIdx);
                    relFeatsIdx=rankfeatures(this.feats(:,trainIdx),this.lbls(trainIdx),'CRITERION','ttest','CCWEIGHTING',1,'NUMBEROFINDICES',nFeats);
                    relFeatsCoords=this.tripletCoords(relFeatsIdx,:);
                    relFeatsCoordsLog=cat(1,relFeatsCoordsLog,relFeatsCoords);
                    relFeatsIdxLog=cat(1,relFeatsIdxLog,relFeatsIdx);
                end
                fprintf('%d/%d\n',currRep,nReps);
            end
        end
                
        function BAcc=testClassifier(this)
            % Recover labels proportions
            pdf=histcounts(this.lbls,length(unique(this.lbls)),'Normalization','probability');
            costMat=[0,1/pdf(1);1/pdf(2),0];
            
%             % Lasso (or elastic net, depending on Alpha parameter)
%             if isempty(this.B) % Skip lasso if it has already been performed here
%                 [this.B,this.FitInfo] = lassoglm(this.feats',this.lbls,'binomial','CV',5,'Alpha',.5,'link','logit');
%             end
            
            % Determine best classifier structure, using half data for
            % training and half for validation (this causes some degree of
            % double dipping with the next step)
%             BAccFun=@(Yreal,Yest)((sum((Yreal==0).*(Yest==0))/sum(Yreal==0))+(sum((Yreal==1).*(Yest==1))/sum(Yreal==1)))/2;
            
%             %% Random forest (or is it?)
%             errFun=@(x)1-BAccFun(this.lbls(2:2:end),...
%                 predict(fitcensemble(this.feats(featsIdx,1:2:end)',this.lbls(1:2:end)',...
%                 'Cost',costMat,'NumLearningCycles',x.nCycles,'Method',char(x.Method),'Cost',costMat,...
%                 'Learners',templateTree('MaxNumSplits',x.MaxNumSplits,'MinLeafSize',x.nLeaves)),this.feats(featsIdx,2:2:end)'));
%             optCycles=optimizableVariable('nCycles',[20 500],'Type','integer');
%             optLeaves=optimizableVariable('nLeaves',[1 100],'Type','integer');
%             optMethod=optimizableVariable('Method',{'Bag', 'GentleBoost', 'LogitBoost', 'RUSBoost'});
%             optMaxSplits=optimizableVariable('MaxNumSplits',[2 100],'Type','integer');
%             results=bayesopt(errFun,[optCycles,optLeaves,optMethod,optMaxSplits],...
%                 'AcquisitionFunctionName','expected-improvement-plus','MaxObjectiveEvaluations',50);%,'Verbose',0,'PlotFcn',{});
%                 
%             % Train cross-validated classifier
%             ens=fitcensemble(this.feats(featsIdx,:)',this.lbls','KFold',2,'Cost',costMat,'Method',char(results.XAtMinEstimatedObjective.Method),...
%                 'NumLearningCycles',results.XAtMinEstimatedObjective.nCycles,'Learners',templateTree('MaxNumSplits',...
%                 results.XAtMinEstimatedObjective.MaxNumSplits,'MinLeafSize',results.XAtMinEstimatedObjective.nLeaves));
%             
%             % Recover predictions
%             lblsEst=ens.kfoldPredict;
            
%             %%  Gauss SVM
%             errFun=@(x)1-BAccFun(this.lbls(2:2:end),...
%                 predict(fitcsvm(this.feats(featsIdx,1:2:end)',this.lbls(1:2:end)',...
%                 'Cost',costMat,'BoxConstraint',x.BC,'KernelScale',x.KS,'Cost',costMat,...
%                 'KernelFunction','gaussian'),this.feats(featsIdx,2:2:end)'));
%             optBC=optimizableVariable('BC',[1e-7,1e7],'Transform','log');
%             optKS=optimizableVariable('KS',[1e-3,1e3],'Transform','log');
%             results=bayesopt(errFun,[optBC,optKS],...
%                 'AcquisitionFunctionName','expected-improvement-plus','MaxObjectiveEvaluations',50);%,'Verbose',0,'PlotFcn',{});
%             
%             % Train cross-validated classifier
%             svmMdl=fitcsvm(this.feats(featsIdx,:)',this.lbls','KFold',5,'Cost',costMat,'KernelFunction','gaussian',...
%                 'KernelScale',results.XAtMinEstimatedObjective.KS,...
%                 'BoxConstraint',results.XAtMinEstimatedObjective.BC);
%             
%             % Recover predictions
%             lblsEst=svmMdl.kfoldPredict;

%             %%  Linear SVM (best so far)
%             errFun=@(x)1-BAccFun(this.lbls(2:2:end),...
%                 predict(fitcsvm(this.feats(featsIdx,1:2:end)',this.lbls(1:2:end)',...
%                 'Cost',costMat,'BoxConstraint',x.BC,'Cost',costMat,...
%                 'KernelFunction','linear'),this.feats(featsIdx,2:2:end)'));
%             optBC=optimizableVariable('BC',[1e-7,1e7],'Transform','log');
%             results=bayesopt(errFun,optBC,...
%                 'AcquisitionFunctionName','expected-improvement-plus','MaxObjectiveEvaluations',20);%,'Verbose',0,'PlotFcn',{});
%             
%             % Train cross-validated classifier
%             svmMdl=fitcsvm(this.feats(featsIdx,:)',this.lbls','KFold',5,'Cost',costMat,'KernelFunction','linear',...
%                 'BoxConstraint',results.XAtMinEstimatedObjective.BC);
%             
%             % Recover predictions
%             lblsEst=svmMdl.kfoldPredict;
            
%             %%  Naive Bayes
%             errFun=@(x)1-BAccFun(this.lbls(2:2:end),...
%                 predict(fitcnb(this.feats(featsIdx,1:2:end)',this.lbls(1:2:end)',...
%                 'Cost',costMat,'Width',x.K,'Cost',costMat,...
%                 'DistributionNames','kernel'),this.feats(featsIdx,2:2:end)'));
%             optK=optimizableVariable('K',[1e-4,1e4],'Transform','log');
%             results=bayesopt(errFun,optK,...
%                 'AcquisitionFunctionName','expected-improvement-plus','MaxObjectiveEvaluations',30);%,'Verbose',0,'PlotFcn',{});
%             
%             % Train cross-validated classifier
%             nbMdl=fitcnb(this.feats(featsIdx,:)',this.lbls','KFold',5,'Cost',costMat,'DistributionNames','kernel',...
%                 'Width',results.XAtMinEstimatedObjective.K);
%             
%             % Recover predictions
%             lblsEst=nbMdl.kfoldPredict;
            
%             %% Poly-SVM
%             errFun=@(x)1-BAccFun(this.lbls(2:2:end),...
%                 predict(fitcsvm(this.feats(featsIdx,1:2:end)',this.lbls(1:2:end)',...
%                 'Cost',costMat,'BoxConstraint',x.BC,'Cost',costMat,...
%                 'PolynomialOrder',x.PO,'KernelFunction','polynomial'),this.feats(featsIdx,2:2:end)'));
%             optBC=optimizableVariable('BC',[1e-4,1e4],'Transform','log');
%             optPO=optimizableVariable('PO',[2 4],'Type','integer');
%             results=bayesopt(errFun,[optBC,optPO],...
%                 'AcquisitionFunctionName','expected-improvement-plus','MaxObjectiveEvaluations',50);%,'Verbose',0,'PlotFcn',{});
%             
%             % Train cross-validated classifier
%             svmMdl=fitcsvm(this.feats(featsIdx,:)',this.lbls','KFold',5,'Cost',costMat,'KernelFunction','polynomial',...
%                 'BoxConstraint',results.XAtMinEstimatedObjective.BC,...
%                 'PolynomialOrder',results.XAtMinEstimatedObjective.PO);
%             
%             % Recover predictions
%             lblsEst=svmMdl.kfoldPredict;
% 
%             % Train cross-validated classifier
%             svmMdl=fitcsvm(this.feats(featsIdx,:)',this.lbls','KFold',5,'Cost',costMat,'KernelFunction','polynomial',...
%                 'BoxConstraint',50,...
%                 'PolynomialOrder',2);
%             
%             % Recover predictions
%             lblsEst=svmMdl.kfoldPredict;
            
            % Train cross-validated classifier
            C=cvpartition(length(this.lbls),'kfold',5);
            lblsEst=zeros(size(this.lbls));
            for currPart=1:C.NumTestSets
                % Recover train and test sets
                trainData=this.feats(:,C.training(currPart));
                trainLbls=this.lbls(C.training(currPart));
                testData=this.feats(:,C.test(currPart));
                
%                 % Ignore "bad" features
                featsIdx=fisherScore(trainData',trainLbls)>.25;
%                 featsIdx=1:size(this.feats,1);
                
%                 %%  Gauss SVM
%                 errFun=@(x)1-BAccFun(trainLbls(2:2:end),...
%                     predict(fitcsvm(trainData(featsIdx,1:2:end)',trainLbls(1:2:end)',...
%                     'Cost',costMat,'BoxConstraint',x.BC,'KernelScale',x.KS,'Cost',costMat,...
%                     'KernelFunction','gaussian'),trainData(featsIdx,2:2:end)'));
%                 optBC=optimizableVariable('BC',[1e-7,1e7],'Transform','log');
%                 optKS=optimizableVariable('KS',[1e-3,1e3],'Transform','log');
%                 results=bayesopt(errFun,[optBC,optKS],...
%                     'AcquisitionFunctionName','expected-improvement-plus','MaxObjectiveEvaluations',25);%,'Verbose',0,'PlotFcn',{});
%                 
%                 % Train cross-validated classifier
%                 svmMdl=fitcsvm(trainData(featsIdx,:)',trainLbls,'Cost',costMat,'KernelFunction','gaussian',...
%                     'KernelScale',results.XAtMinEstimatedObjective.KS,...
%                     'BoxConstraint',results.XAtMinEstimatedObjective.BC);
                
                %%  Linear SVM (best so far)
                errFun=@(x)1-computeBAcc(trainLbls(2:2:end),...
                    predict(fitcsvm(trainData(featsIdx,1:2:end)',trainLbls(1:2:end)',...
                    'Cost',costMat,'BoxConstraint',x.BC,'Cost',costMat,...
                    'KernelFunction','linear'),trainData(featsIdx,2:2:end)'));
                optBC=optimizableVariable('BC',[1e-6,1e6],'Transform','log');
                results=bayesopt(errFun,optBC,...
                    'AcquisitionFunctionName','expected-improvement-plus','MaxObjectiveEvaluations',20);%,'Verbose',0,'PlotFcn',{});
                
                % Train cross-validated classifier
                svmMdl=fitcsvm(trainData(featsIdx,:)',trainLbls','Cost',costMat,'KernelFunction','linear',...
                    'BoxConstraint',results.XAtMinEstimatedObjective.BC);

                % Recover predictions
                lblsEst(C.test(currPart))=svmMdl.predict(testData(featsIdx,:)');
            end
            
            % Compute balanced accuracy
            BAcc=computeBAcc(this.lbls,lblsEst);
        end
        
        function BAcc=repTestClassifier(this,nReps)
            BAcc=nan(nReps,1);
            for currRep=1:nReps
                % Recover labels proportions
                pdf=histcounts(this.lbls,length(unique(this.lbls)),'Normalization','probability');
                costMat=[0,1/pdf(1);1/pdf(2),0];
                BAccFun=@(Yreal,Yest)((sum((Yreal==0).*(Yest==0))/sum(Yreal==0))+(sum((Yreal==1).*(Yest==1))/sum(Yreal==1)))/2;
                % Train cross-validated classifier
                C=cvpartition(length(this.lbls),'kfold',5);
                lblsEst=zeros(size(this.lbls));
                for currPart=1:C.NumTestSets
                    % Recover train and test sets
                    trainIdx=find(C.training(currPart));
%                     valIdx=trainIdx(randperm(length(trainIdx),round(length(trainIdx)/5)));
%                     trainIdx=setdiff(trainIdx,valIdx);
                    trainData=this.feats(:,trainIdx);
                    trainLbls=this.lbls(trainIdx);
%                     valData=this.feats(:,valIdx);
%                     valLbls=this.lbls(valIdx);
                    testData=this.feats(:,C.test(currPart));
%                     
%                     % Train cross-validated classifier
%                     errFun=@(x)1-BAccFun(valLbls,...
%                         predict(fitcsvm(trainData',trainLbls,...
%                         'Cost',costMat,'BoxConstraint',x.BC,'Cost',costMat,...
%                         'KernelFunction','linear'),valData'));
%                     optBC=optimizableVariable('BC',[1e-4,1e4],'Transform','log');
%                     results=bayesopt(errFun,optBC,...
%                         'AcquisitionFunctionName','expected-improvement-plus','MaxObjectiveEvaluations',20,'Verbose',0,'PlotFcn',{});
%                     
%                     % Train cross-validated classifier
%                     svmMdl=fitcsvm(trainData',trainLbls','Cost',costMat,'KernelFunction','linear',...
%                         'BoxConstraint',results.XAtMinEstimatedObjective.BC);
                                        
                    % Train cross-validated classifier
                    svmMdl=fitcsvm(trainData',trainLbls','Cost',costMat,'KernelFunction','linear',...
                        'BoxConstraint',25);
                    
                    % Recover predictions
                    lblsEst(C.test(currPart))=svmMdl.predict(testData');
                end
                BAcc(currRep)=BAccFun(this.lbls,lblsEst);
                fprintf('%d/%d; BAcc: %0.2f\n',currRep,nReps,BAcc(currRep));
            end
        end
        
        function n=get.nROIs(this)
            n=size(this.FC{1},1);
        end
        
        function tCoords=get.tripletCoords(this)
            tCoords=this.tripletIdx;
            tCoords(tCoords==0)=NaN;
            tCoords=squeeze(nanmean(tCoords,2));
        end
        
        function fIdx=get.featIdx(this)
            fIdx=this.tripletIdx;
            fIdx(fIdx==0)=NaN;
            fIdx=squeeze(nanmean(fIdx,2));
            fIdx=fIdx*[this.nROIs^2;this.nROIs;1];
        end
        
        function fIdx=get.relFeatIdx(this)
            fIdx=zeros(this.nROIs^3,1);
            tempFeatIdx=this.featIdx;
            for currRelFeatIdx=1:length(tempFeatIdx)
                fIdx(tempFeatIdx(currRelFeatIdx))=currRelFeatIdx;
            end
        end
    end
    
    methods (Static)
        function prepareDataFile(fileName,xlsFile,dataFolder)
            % Recover patient data from xls file
            [~,~,RAW]=xlsread(xlsFile);
            f=RAW(1,:);
            f(cellfun(@(x)~ischar(x),f))=[];
            s=cell2struct(RAW(2:end,1:length(f)),f,2);
            
            % Use recovered data to load patients data
            lbls=zeros(size(s));
            FCs=cell(length(s),1);
            D=dir(dataFolder);
            for currSubj=1:length(s)
                for currD=3:length(D)
                    D2=dir([dataFolder,'\',D(currD).name]);
                    T=struct2table(D2);
                    subjIdx=find(cellfun(@(x)strcmp(x,sprintf('%04d',s(currSubj).IDsoggetto)),T.name));
                    if ~isempty(subjIdx)
                        subjData=dlmread(sprintf('%s\\%s\\%s\\resting_conn\\AAL116_gm_ROISignals.txt',dataFolder,D(currD).name,D2(subjIdx).name));
                        
                        % Only use first 90 ROIs
                        subjData=subjData(:,1:90);
                        
                        % Occasionally, one region may be missing from a
                        % subject. Put white noise there
                        if sum(sum(isnan(subjData)))
                            subjData(isnan(subjData))=randn(size(subjData(isnan(subjData))));
                            fprintf('ROI missing data, substituting with white noise\n');
                        end
                        
                        % Compute covariance matrices
                        subjCorr=cov(subjData);
                        break;
                    end
                end
                FCs{currSubj}=subjCorr;
                
                % Recover lables
                if isempty(s(currSubj).Dis_Prog_2)
                    lbls(currSubj)=nan;
                else
                    lbls(currSubj)=s(currSubj).Dis_Prog_2;
                end
                fprintf('%d/%d\n',currSubj,length(s));
            end
            fprintf('Saving results...\n');
            save(fileName,'FCs','-v7.3');
            save(fileName,'lbls','-append');
            save(fileName,'xlsFile','-append');
            save(fileName,'dataFolder','-append');
        end
        
        function plotTripletsHist(relFeatsLog)
            % Put together all iterations and order triplets
            commonFeats=cat(1,relFeatsLog{:});
            ordrdFeats=sort(commonFeats,2);
            
            % Compute triplet indexes and count istances
            tripletID=ordrdFeats*[90^2;90;1];
            cnts=histcounts(tripletID,.5:max(tripletID)+.5);
            
            % Plot triplets "histogram"
            figure;
            for currBin=1:length(cnts)
                if cnts(currBin)>0
                    featIdxs=find(tripletID==currBin);
                    currTriplet=ordrdFeats(featIdxs(1),:);
                    plot3(currTriplet(1),currTriplet(2),currTriplet(3),'.','MarkerSize',cnts(currBin)*3);
                    hold on
                end
            end
            title('Triplets most commonly selected as significant');
            xlabel('ROI 1');
            ylabel('ROI 2');
            zlabel('ROI 3');
            h=findall(gcf,'Type','Line');
            set(h,'LineWidth',1.5)
            set(gca,'LineWidth',1.5,'TickDir','out')
            
            % Plot ROIs histogram
            figure;
            [ROIcnts,binLims]=histcounts(reshape(commonFeats,[],1),'BinMethod','integers','Normalization','probability');
            bar(binLims(1:end-1)+(binLims(2)-binLims(1))*.5,ROIcnts);
            title('ROIs most commonly selected as significant')
            xlabel('ROI');
            ylabel('Relative frequence');
            h=findall(gcf,'Type','Line');
            set(h,'LineWidth',1.5)
            set(gca,'LineWidth',1.5,'TickDir','out')
        end
    end
end