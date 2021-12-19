% Author: Sarah Maddox Groves
% Date: December 10, 2021
% This is a script based on the mlx file with the same name. 
% It is to be run on a cluster that supports MATLAB.
% File structure in Accre:

% - this script
% - data  
%     - downsampled_X_magic_for_ParTI
%     - X_magic_for_ParTI_GeneNames
%     - X_magic_for_ParTI_Attr
% - ParTI (package files)
% - out
%     - 3
%     - 4
%     - 5

clear all
% Shown to be best number of archetypes by elbow plot from py_pcha
nArchetypes = 3;

global ForceNArchetypes; ForceNArchetypes = nArchetypes;

origPath = pwd;
inputPath = strcat(origPath, '/data/');

geneExpression = dlmread(strcat(inputPath,'downsampled_X_magic_for_ParTI.csv'), ',');
% geneExpression = geneExpression'
% The file is formated as samples (i.e. patients) x genes.
% We load gene names.
geneNames = importdata(strcat(inputPath,'X_magic_for_ParTI_GeneNames.csv'), ',');

cd ./ParTI
[discrAttrNames, discrAttr] = ...
    read_enriched_csv(strcat(inputPath,'X_magic_for_ParTI_Attr.csv'), ',');
contIdcs = [1 3 4 5 6 7 9 10];
contAttrNames = discrAttrNames(:,contIdcs);
contAttr = discrAttr(:,contIdcs);

contAttr = str2double(contAttr);
discIdcs = [2 8];
discrAttrNames = discrAttrNames(:,discIdcs);
discrAttr = discrAttr(:,discIdcs);

[GOExpression,GONames,~,GOcat2Genes] = MakeGOMatrix(geneExpression, geneNames, ...
                {'MSigDB/c5.go.bp.v7.2.symbols.gmt'}, ...
                10);
contAttrNames = [contAttrNames, GONames];
contAttr = [contAttr, GOExpression];
discrAttrNames = regexprep(discrAttrNames, '_', ' ');
contAttrNames = regexprep(contAttrNames, '_', ' ');

dim = 12;
binSize = 0.1 ;
for nArchetypes = [3 4 5]
    algNum = 5
    ForceNArchetypes = nArchetypes;
    random_int = randi(1000000);
    outputFileName = sprintf('%d',random_int);
    outPath = strcat(origPath, sprintf('/out/%d',nArchetypes));
    [arc, arcFull, pc, errs, pval, coefs, P, tRatioRand, tRatioReal] = ParTI(geneExpression, algNum, dim, discrAttrNames, ...
    discrAttr, 0, contAttrNames, contAttr,GOcat2Genes, binSize, strcat(outPath,'/', outputFileName));
    fid = fopen(strcat(origPath,'/params.txt'), 'a+');
    fprintf(fid, '%d, %d, %d, %d, %d,%d, %d\n', random_int,algNum,dim,binSize, nArchetypes, pval, tRatioReal);
    fclose(fid);
    csvwrite(strcat(outPath,sprintf('/%s_arc_full.csv',outputFileName)), arcFull);
    csvwrite(strcat(outPath,sprintf('/%s_errs.csv',outputFileName)), errs);
    csvwrite(strcat(outPath,sprintf('/%s_coefs.csv',outputFileName)), coefs);
    csvwrite(strcat(outPath,sprintf('/%s_tRatioRand.csv',outputFileName)), tRatioRand);
    csvwrite(strcat(outPath,sprintf('/%s_arc.csv',outputFileName)), arc);
    csvwrite(strcat(outPath,sprintf('/%s_pc.csv',outputFileName)), pc);
    clear arc arcFull pc errs pval coefs P tRatioRand tRatioReal
end