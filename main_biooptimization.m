addpath(genpath('/home/paula/Documents/matRad/DKFZCollab/RemoFig1OptimizationWithBio/matRad_files4bioOpt'));

% Water phantom
ctDim = [160,160,160]; % x,y,z dimensions
ctResolution = [3,3,3]; % x,y,z the same here!

builder = matRad_PhantomBuilder(ctDim,ctResolution,1);

builder.addBoxTarget('Target',[20 20 20],'HU', 0,'objectives',struct(DoseObjectives.matRad_SquaredDeviation(200,60)));
builder.addBoxOAR('OAR',[20 11 20],'offset',[0 20 0],'HU',0 ,'objectives',struct(DoseObjectives.matRad_SquaredDeviation(5,0)));
builder.addBoxOAR('BODY',[80 80 80],'HU',0);

[ct,cst] = builder.getctcst();
cst{1,5}.Priority = 1;
cst{2,5}.Priority = 1;
cst{3,5}.Priority = 2;

%% Plan configuration
pln.numOfFractions  = 30;
pln.radiationMode   = 'protons';
pln.machine         = 'Generic';
pln.propOpt.quantityOpt = 'effect';

% beam geometry settings
pln.propStf.bixelWidth      = 5;
pln.propStf.gantryAngles    = [0,90]; % [°] ;
pln.propStf.couchAngles     = zeros(numel(pln.propStf.gantryAngles),1);
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);

% optimization settings
pln.propDoseCalc.calcLET = true;

pln.propOpt.spatioTemp      = 0;

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 5; % [mm]

% retrieve bio model parameters
pln.bioModel = matRad_Wedenberg();

% robustness
robSetting = 'STOCH';

%% Scenario settings
% Uncertainty model 
shiftSD = [2.25 2.25 2.25]; %Standard deviation for x/y/z shifts in mm
rangeRelSD = 3.5; % Standard deviation in % of relative range error
rangeAbsSD = 1;   % Standard deviation in mm of absolute range error
alphaRelSD = 20;  % percent
betaRelSD = 20; % percent
numOptScen = 36;

% Random Scenarios for Optimization with physical shifts
optScenarios = matRad_RandomScenarios(ct);
optScenarios.shiftSD = shiftSD;
optScenarios.rangeRelSD = rangeRelSD;
optScenarios.rangeAbsSD = rangeAbsSD;
optScenarios.nSamples = numOptScen;
optScenarios.listAllScenarios();

% Random Scenarios for Optimization with physical and biological shifts
optbioScenarios = matRad_RandomScenariosWithBio(ct);
optbioScenarios.shiftSD = shiftSD;
optbioScenarios.rangeRelSD = rangeRelSD;
optbioScenarios.rangeAbsSD = rangeAbsSD;
optbioScenarios.bioRelSD = [alphaRelSD, betaRelSD];
optbioScenarios.nSamples = numOptScen;
optbioScenarios.listAllScenarios();

%% generate steering file 
stf = matRad_generateStf(ct,cst,pln);

%% dose calculation for robust optimization
% only physical shifts
pln_phys = pln;
pln_phys.multScen = optScenarios;
dij_phys = matRad_calcDoseInfluence(ct, cst, stf, pln_phys);

% physical + biological shifts
pln_bio = pln;
pln_bio.multScen = optbioScenarios;
dij_bio = matRad_calcDoseInfluence(ct, cst, stf, pln_bio);

linIxDIJ = find(~cellfun(@isempty,dij_bio.physicalDose(1,:,:)))';

for i=1:numOptScen
    ax(i) = (1+optbioScenarios.alphaShift(i)).*cst{1,5}.alphaX;
    bx(i) = (1+optbioScenarios.betaShift(i)).*cst{1,5}.betaX;
end

for i=1:numel(linIxDIJ)
    scen = linIxDIJ(i);
    physicalDose = dij_bio.physicalDose(scen);
    LET = dij_bio.mLETDose(scen);

    bixel.vABratio = ax(i)./bx(i);
    bixel.LET = LET{1};
    [RBEmin, RBEmax] = pln.bioModel.getRBEminMax(bixel);

    dij_bio.mAlphaDose(scen) = {RBEmax .* ax(i) .* physicalDose{1}};
    dij_bio.mSqrtBetaDose(scen) = {RBEmin .* sqrt(bx(i)) .* physicalDose{1}};
end 

%% inverse planning for imrt
%Enable robustness in objectives
for i = 1:size(cst,1)
    for j = 1:numel(cst{i,6})
        cst{i,6}{j}.robustness = robSetting;
    end
end
resultGUI_phys  = matRad_fluenceOptimization(dij_phys,cst,pln_phys);
resultGUI_bio  = matRad_fluenceOptimization(dij_bio,cst,pln_bio);

%% analysis
% Random Scenarios 4 analysis
anardnScenarios = matRad_RandomScenarios(ct);
anardnScenarios.shiftSD = shiftSD;
anardnScenarios.rangeRelSD = rangeRelSD;
anardnScenarios.rangeAbsSD = rangeAbsSD;
anardnScenarios.nSamples = 100;
anardnScenarios.listAllScenarios();

pln_phys_ana = pln_phys;
pln_phys_ana.multScen = anardnScenarios;
resultGUI_phys_rdnanalysis = matRad_calcDoseForward(ct,cst,stf,pln_phys_ana,resultGUI_phys.w);

pln_bio_ana = pln_bio;
pln_bio_ana.multScen = anardnScenarios;
resultGUI_bio_rdnanalysis = matRad_calcDoseForward(ct,cst,stf,pln_bio_ana,resultGUI_bio.w);

% bio shifts
truncSigma = 2;
Sigma = diag([alphaRelSD/100,betaRelSD/100].^2);
[cholcs,cholp] = chol(Sigma);
d = size(Sigma,1);
numScen = anardnScenarios.nSamples;
pd = makedist('Normal','mu',0,'sigma',1);
t = truncate(pd,-truncSigma,truncSigma);
bioshift = random(t,numScen,d) * cholcs;

alphaX= (1+bioshift(:,1)).*cst{1,5}.alphaX;
betaX = (1+bioshift(:,2)).*cst{1,5}.betaX;

% physical shifts
for i=1:numScen
    fnameRBExDose = sprintf('RBExDose_scen%d', i);
    fnamephysDose = sprintf('physicalDose_scen%d', i);
    fnameeffect = sprintf('effect_scen%d', i);
    fnameLET = sprintf('LET_scen%d',i);

    bixel.LET = resultGUI_phys_rdnanalysis.(fnameLET);

    bixel.vABratio = alphaX(i)/betaX(i);

    [RBEmin, RBEmax] = pln_phys.bioModel.getRBEminMax(bixel);

    alphaDose = (RBEmax.*alphaX(i)).*resultGUI_phys_rdnanalysis.(fnamephysDose);
    sqrtBetaDose = (RBEmin.*sqrt(betaX(i))).*resultGUI_phys_rdnanalysis.(fnamephysDose);

    resultGUI_phys_rdnanalysis_bio.(fnameeffect) = alphaDose + sqrtBetaDose.*sqrtBetaDose;
    resultGUI_phys_rdnanalysis_bio.(fnameRBExDose) = (sqrt(alphaX(i).^2 + 4 .* betaX(i) .* resultGUI_phys_rdnanalysis_bio.(fnameeffect)) - alphaX(i))./(2.*betaX(i));
end

% physical + biological shifts
for i=1:numScen
    fnameRBExDose = sprintf('RBExDose_scen%d', i);
    fnamephysDose = sprintf('physicalDose_scen%d', i);
    fnameeffect = sprintf('effect_scen%d', i);
    fnameLET = sprintf('LET_scen%d',i);

    bixel.LET = resultGUI_bio_rdnanalysis.(fnameLET);

    bixel.vABratio = alphaX(i)/betaX(i);

    [RBEmin, RBEmax] = pln_bio.bioModel.getRBEminMax(bixel);

    alphaDose = (RBEmax.*alphaX(i)).*resultGUI_bio_rdnanalysis.(fnamephysDose);
    sqrtBetaDose = (RBEmin.*sqrt(betaX(i))).*resultGUI_bio_rdnanalysis.(fnamephysDose);

    resultGUI_bio_rdnanalysis_bio.(fnameeffect) = alphaDose + sqrtBetaDose.*sqrtBetaDose;
    resultGUI_bio_rdnanalysis_bio.(fnameRBExDose) = (sqrt(alphaX(i).^2 + 4 .* betaX(i) .* resultGUI_bio_rdnanalysis_bio.(fnameeffect)) - alphaX(i))./(2.*betaX(i));
end

DVHPlot_ComprobustOptnominal(resultGUI_bio_rdnanalysis_bio,cst,'RBExDose',anardnScenarios.nSamples,5,95,25,75,1);
