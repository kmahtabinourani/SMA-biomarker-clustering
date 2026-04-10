%% 01_compute_coexpression_by_biomarker
% This script computes Pearson correlation values between each gene and
% three reference SMA biomarkers in both Control and SMA groups.
%
% Workflow:
% 1. Load the filtered gene expression table
% 2. Compute Pearson correlation for each gene relative to each biomarker
% 3. Save full correlation tables for Control and SMA groups
% 4. Extract co-expressed genes using a correlation cutoff
% 5. Save overlap gene sets between Control and SMA for each biomarker
%
% Expected input:
%   data/Filtered_genes.xlsx
%
% Output files:
%   outputs/Control_CrossCorr.csv
%   outputs/SMA_CrossCorr.csv
%   outputs/Control_Coexpressed_<Biomarker>.csv
%   outputs/SMA_Coexpressed_<Biomarker>.csv
%   outputs/Overlap_<Biomarker>.csv

clc;
clear;
close all;

%% ===================== USER SETTINGS =====================
inputFile = fullfile('data', 'Filtered_genes.xlsx');
outputDir = 'outputs';

% Correlation cutoff for defining co-expressed genes
cutoff = 0.7;

% Row indices of biomarkers in the filtered table
% Order: NfL, HSPA7, SMN2
biomarkerIndices = [8361, 7976, 7746];
biomarkerNames = {'NfL', 'HSPA7', 'SMN2'};

%% ===================== SETUP =====================
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

%% ===================== LOAD INPUT TABLE =====================
filteredTable = readtable(inputFile);

% Convert gene IDs to string format
geneIDs = string(filteredTable.Gene_ID);

% Extract biomarker expression values
% Columns 2:6 = Control samples
% Columns 7:11 = SMA samples
controlBiomarkers = table2array(filteredTable(biomarkerIndices, 2:6));
smaBiomarkers     = table2array(filteredTable(biomarkerIndices, 7:11));

numGenes = height(filteredTable);
numBiomarkers = numel(biomarkerIndices);

% Preallocate correlation matrices
controlCrossCorr = zeros(numGenes, numBiomarkers);
smaCrossCorr     = zeros(numGenes, numBiomarkers);

%% ===================== COMPUTE PEARSON CORRELATIONS =====================
for g = 1:numGenes
    controlVals = table2array(filteredTable(g, 2:6));
    smaVals     = table2array(filteredTable(g, 7:11));

    for b = 1:numBiomarkers
        controlBio = controlBiomarkers(b, :);
        smaBio     = smaBiomarkers(b, :);

        controlCrossCorr(g, b) = corr(controlVals(:), controlBio(:), 'Type', 'Pearson');
        smaCrossCorr(g, b)     = corr(smaVals(:), smaBio(:), 'Type', 'Pearson');
    end
end

%% ===================== SAVE FULL CORRELATION TABLES =====================
controlTable = array2table(controlCrossCorr, 'VariableNames', biomarkerNames);
controlTable.Gene_ID = geneIDs;
controlTable = movevars(controlTable, 'Gene_ID', 'Before', 1);

smaTable = array2table(smaCrossCorr, 'VariableNames', biomarkerNames);
smaTable.Gene_ID = geneIDs;
smaTable = movevars(smaTable, 'Gene_ID', 'Before', 1);

writetable(controlTable, fullfile(outputDir, 'Control_CrossCorr.csv'));
writetable(smaTable,     fullfile(outputDir, 'SMA_CrossCorr.csv'));

%% ===================== EXTRACT CO-EXPRESSED GENES =====================
for b = 1:numBiomarkers
    biomarkerName = biomarkerNames{b};

    % ----- Control -----
    controlR = controlCrossCorr(:, b);
    controlIdx = find(controlR > cutoff | controlR < -cutoff);

    controlCoexpressed = table( ...
        geneIDs(controlIdx), controlR(controlIdx), ...
        'VariableNames', {'Gene_ID', 'Pearsons_r'});

    writetable(controlCoexpressed, ...
        fullfile(outputDir, ['Control_Coexpressed_' biomarkerName '.csv']));

    % ----- SMA -----
    smaR = smaCrossCorr(:, b);
    smaIdx = find(smaR > cutoff | smaR < -cutoff);

    smaCoexpressed = table( ...
        geneIDs(smaIdx), smaR(smaIdx), ...
        'VariableNames', {'Gene_ID', 'Pearsons_r'});

    writetable(smaCoexpressed, ...
        fullfile(outputDir, ['SMA_Coexpressed_' biomarkerName '.csv']));

    fprintf('Control %s co-expressed genes: %d\n', biomarkerName, height(controlCoexpressed));
    fprintf('SMA %s co-expressed genes: %d\n', biomarkerName, height(smaCoexpressed));
end

%% ===================== DETECT OVERLAP BETWEEN CONTROL AND SMA =====================
for b = 1:numBiomarkers
    biomarkerName = biomarkerNames{b};

    controlGenes = controlTable{ ...
        controlCrossCorr(:, b) > cutoff | controlCrossCorr(:, b) < -cutoff, ...
        'Gene_ID'};

    smaGenes = smaTable{ ...
        smaCrossCorr(:, b) > cutoff | smaCrossCorr(:, b) < -cutoff, ...
        'Gene_ID'};

    overlapGenes = intersect(controlGenes, smaGenes);

    fprintf('%s - Control vs SMA overlap: %d genes\n', biomarkerName, numel(overlapGenes));

    overlapTable = table(overlapGenes, 'VariableNames', {'Gene_ID'});
    writetable(overlapTable, ...
        fullfile(outputDir, ['Overlap_' biomarkerName '.csv']));
end

disp('Co-expression analysis completed successfully.');
