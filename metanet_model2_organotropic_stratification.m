% Predict survival-associated metastatic risk
% ordinal model

clear
clc
close all

addpath /Users/bjiang/Documents/MATLAB/orca/src/Utils
addpath /Users/bjiang/Documents/MATLAB/orca/src/Measures
addpath /Users/bjiang/Documents/MATLAB/orca/src/Algorithms

if (exist ('OCTAVE_VERSION', 'builtin') > 0)
  pkg load statistics
  try
    graphics_toolkit ('gnuplot')
  catch
    error('This code uses gnuplot for plotting. Please install gnuplot and restart Octave to run this code.')
  end
end

% Disable MATLAB warnings
warning('off','MATLAB:nearlySingularMatrix')
warning('off','stats:mnrfit:IterOrEvalLimit')

%% load data

load input.mat

%% Apply POM model
% Create the POM object
algorithmObj = POM();

% Train POM
POMinfo = algorithmObj.fitpredict(train,test);

%% TCGA Validation
model = POMinfo.model;
[vprob,vdlow,vdhi] = mnrval([-model.thresholds'; model.projection],...
    scaledXval,model.Stats,'model','ordinal','interactions','off',...
    'Link','logit');
[~,vpredlabel] = max(vprob,[],2);
ystar = scaledXval*model.projection;

[p, fh, stats] = svval(vpredlabel,clinic,false);
pairwisepvalue = table(stats.ParwiseName, [stats.ParwiseStats(1).p_MC;stats.ParwiseStats(2).p_MC;stats.ParwiseStats(3).p_MC]);
disp(pairwisepvalue)
