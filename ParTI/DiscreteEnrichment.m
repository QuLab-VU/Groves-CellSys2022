function [table, pval, PoverQ, Pmaxim, isSignificantAfterFDR] = DiscreteEnrichment(DataPointsInd,EnMatDis,binSize,DiscFeatName,OutputFileName,evalPmax)

%Inputs
% 1. data point indexes sorted according to distances from archetypes (the 
%    first output of sortDataByDistane(DataPoints,Archetypes).
% 2. EnMatDis is a boolean matrix with each Discrete feature as a column, each row is a
%  datapoint (order of datapoints and their features should match, order of columns and feature labels should match).
% 3. binSize is the value of the percent each bin should have in the
%    enrichment calculation
% 4. DiscFeatName contains the (human-readable) names of the discrete features
% 5. OutputFileName is a tag that will be used to name the file in which
% the figure produced by this function will be saved.
% 6. evalPmax determines whether to compute the probability that enrichment is maximal in the first bin. Setting this parameter to 0 skips this computation to save time.
% Initialazing
% significance threshold after FDR (Benjamini Hochberg). 
ThreshHoldBH = 1;

[Numarchs, numDataPoints] = size(DataPointsInd);
[numDataPoints2, numFeatures] = size(EnMatDis);
if(numDataPoints2 ~= numDataPoints)
    table = NaN;
    return; 
end
numOfBins = round(1 / binSize);

fprintf('Your data was divided into %d bins.\n',numOfBins);
%%
% numDataPoints = 105; numOfBins = 5; %just for testing
breakPoints = floor(linspace(0.5, numDataPoints + 0.5, numOfBins+1));
% diff(breakPoints) %just for testing 
numPointInBin = diff(breakPoints);
breakPoints = breakPoints(2:end);
%%
pval    = zeros(Numarchs,numFeatures);
Pmaxim  = ones(Numarchs,numFeatures);
PoverQ  = cell(Numarchs,1);


%On Each Archetype
for arch = 1:Numarchs
   %Divide to bins
   tempEnrich =  cumsum(EnMatDis(DataPointsInd(arch,:),:)); %returns list of cumulative sums for boolean for each feature
   %e.g. 1 2 3 4 4 4 4 5 6 
   %(would correspond to boolean EnMatDis of 1 1 1 1 0 0 0 1 1 
   binnedEnrichment = diff([zeros(1,numFeatures) ; tempEnrich(breakPoints,:) ]); %breaks the cumulative sums at breakpoints and finds the difference between each (I think)
   % for above example with  breakpoints at 4 and 9, would return
   % diff([0,4,6]) = [4,2] in the bins
   
   %EnGen - calculate p-val by hypergeometric test for each feature
   % hyperGeometricTest
   pval(arch,:) = 1 - hygecdf(binnedEnrichment(1,:)-1,numDataPoints,tempEnrich(end,:),numPointInBin(1));
   % The following is a more precise solution, but works only in Matlab 2013
   % pval(arch,:) = hygecdf(binnedEnrichment(1,:)-1,numDataPoints,tempEnrich(end,:),diff(breakPoints),'upper');

   % Pval of first bin Vs. all data 
           
   %EnPQ - Prepare enrichment graphs
   P = (bsxfun(@rdivide,(binnedEnrichment'),numPointInBin))'; %%round(mean(numPointInBin)));
   % getting percentage (+) in each bin
   % in example above, this is [4/4,2/5] = [1, .2]
   Q = tempEnrich(end,:)./numDataPoints;
   %getting total percentage (+) by dividing last index of tempEnrich
   %(cumsum, i.e. # of +) by total number of data points
   %example, Q = 6/9
   PoverQ{arch} = bsxfun(@rdivide, P,Q);
   
   
   %PMax
   for j = 2:numOfBins
		if evalPmax > 0
			Pmaxim(arch,:) = Pmaxim(arch,:).*HGfuncReg(binnedEnrichment(1,:),...
				numPointInBin(1) - binnedEnrichment(1,:),binnedEnrichment(j,:),...
				numPointInBin(j) - binnedEnrichment(j,:));
		else
			[~, ind] =  max(PoverQ{arch},[],1);
			Pmaxim(arch,:) = (ind == 1); %if index (ind) of the maximum bin is 1 (i.e. closest to arch is highest), return True/1.
		end
   end
   
end



%FDR
FDRs = mafdr(pval(:),'BHFDR',true);
isSignificantAfterFDR = (FDRs<=ThreshHoldBH); 
%WriteTable
%Cols: ArchetypeNum, Features, pvals,Significant After FDR ,Pmaxim
table =  [repmat((1:Numarchs)',numFeatures,1),... %arch number
          reshape(ones(Numarchs,1)*(1:numFeatures),Numarchs*numFeatures,1),... %feat ID
          pval(:),... %p-value
          isSignificantAfterFDR,... %is significant
          Pmaxim(:)]; %Pmaximal
          
 style={'-r','-g','-b','-m','-y','-c','-k','--r','--b','--g','--k','--m','--c','--y'};
   bins = linspace(0, 1, numOfBins);
   subRows = 4; 
   subCols = 4;
   numFeatPerPlot = subRows*subCols;
   %at most 10 figures   
   maxNumOfPlot = 10; 
   
   featuresToDraw = unique( table( (table(:,4)> 0) ,2));
   numOfPlot = ceil(length(featuresToDraw)/numFeatPerPlot);
   if numOfPlot < maxNumOfPlot
    for fig = 1:numOfPlot
        %    for feat = 1:numFeatures
        h = figure;
        for subp = ((fig-1)*numFeatPerPlot+1):min(fig*numFeatPerPlot,length(featuresToDraw))
%            subplot(subRows,subCols,mod(subp,numFeatPerPlot)+1);
            subplot(subRows,subCols,mod(subp-1,numFeatPerPlot)+1);

            hold on
            leg = cell(1,Numarchs);
            for ar = 1:Numarchs
                plot(bins,PoverQ{ar}(:,featuresToDraw(subp)),style{mod(ar,14)},'LineWidth',2);
                title(sprintf('feature: %s',DiscFeatName{featuresToDraw(subp)}),'FontWeight','bold');
                leg{ar} = ['archetype:', num2str(ar)];
            end
            xlabel('bin # / number of bins');ylabel('enrichment');
            hold off
        end
        legend(leg);
        set(gcf,'units','normalized','outerposition',[0, 0 , 1, 1]);
        %figname = [OutputFileName,num2str(fig + 4),'_discreteEnrichment',num2str(fig)]; 
        figname = sprintf('%s_discreteEnrichment%d', OutputFileName, fig);
        if exist('savefig')
            savefig(h, figname);
        end
    end
   else
       fprintf('Skipping enrichment plots because showing them would open %.0f figures.\n', numOfPlot);
   end


