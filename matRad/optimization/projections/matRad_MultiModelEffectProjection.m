classdef matRad_MultiModelEffectProjection < matRad_BackProjection
% matRad_EffectProjection class for effect-based optimization
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2019-2026 the matRad development team.
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    

    properties
        bioModel;
        alphaX;
        betaX;
    end

    methods
        function obj = matRad_MultiModelEffectProjection()
        end
    end
    
    methods
        function RBExD = computeSingleScenario(obj,pD, LET, alphaX, betaX)
            % Need:
                % physical dose distribution (dij.physicalDose{scen}*w);
                % let distribution (dij.mLETDose{scen}*w); % -> ! This is LET * dose
                % Need to divide by dose distribution (check for zeros in
                % physicaldose):
                    % nonZeroIdx = physicalDose ~= 0
                    % LET = zeros(size(physicalDose));
                    % LET(nonZeroIdx) = LETdose/dose;
                % Now that we have dose and LET we can compute RBEmin
                % RBEmax, then alpha*dose and sqrtBeta*dose (distributions)
                
                % Then effect and RBExDose

                bixel.vABratio = alphaX/betaX;
                bixel.LET = LET;

                [RBEmin, RBEmax] = obj.bioModel.getRBEminMax(bixel);

                alphaDose    = RBEmax .* alphaX .* pD;
                sqrtBetaDose = RBEmin .* sqrt(betaX) .* pD;

                effect  = alphaDose + sqrtBetaDose.^2;

                gamma = alphaX/(2*betaX);
                RBExD = sqrt(effect/betaX + gamma^2) - gamma;

        end
        
        function wGrad = projectSingleScenarioGradient(~,w)
            % Need to figure out the derivative of RBExD wrt w (can check
            % the VariableRBEProjection) and compute it through the model
            wGrad = ones(size(w));
            
        end
        
        function [eExp,dOmegaV] = computeSingleScenarioProb(~,dij,scen,w)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispError('Probabilistic Backprojection not availabel yet...\n');
           
        end
        
        function wGrad = projectSingleScenarioGradientProb(~,dij,dExpGrad,dOmegaVgrad,scen,~)
                    matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispError('Probabilistic Backprojection not availabel yet...\n');

        end

         function d = computeResult(obj,dij,w)

             physicalScenarios = find(~cellfun(@isempty,dij.physicalDose(:,:,:)))';
             nBioScenarios      = numel(obj.alphaX);


             % loop over physical scenarios
                % compute dose and LEt distribution
                % Loop over alphaX/betaX
                    % computeSingleScenario
             scenCounter = 0;
             for physScenIdx=physicalScenarios
                physicalDose = dij.physicalDose{physScenIdx}*w;

                nonZeroIdx = physicalDose ~= 0;
                LET = zeros(size(physicalDose));
                LET(nonZeroIdx) = dij.mLETDose{physScenIdx}(nonZeroIdx,:) * w ./ physicalDose(nonZeroIdx);

                for bioIdx=1:nBioScenarios
                    scenCounter = scenCounter +1;
                    d{scenCounter} = obj.computeSingleScenario(physicalDose, LET, obj.alphaX(bioIdx), obj.betaX(bioIdx));
                end
             end

         end

         function wGrad = projectGradient(obj,dij,doseGrad,w)
             physicalScenarios = find(~cellfun(@isempty,dij.physicalDose(:,:,:)))';
             nBioScenarios      = numel(obj.alphaX);

             scenCounter = 0;
             for physScenIdx=physicalScenarios
                
                for bioIdx=1:nBioScenarios
                    scenCounter = scenCounter +1;
                    wGrad{scenCounter} = obj.projectSingleScenarioGradient(w);
                end
             end
             % wGrad = arrayfun(@(scen) projectSingleScenarioGradient(obj,dij,doseGrad,scen,w),obj.scenarios,'UniformOutput',false);

              % Extract the bio scenarios from wGrad

        end
    end
    
    methods (Static)
        function optiFunc = setBiologicalDosePrescriptions(optiFunc,alphaX,betaX)
            doses = optiFunc.getDoseParameters();    
            effect = alphaX*doses + betaX*doses.^2;        
            optiFunc = optiFunc.setDoseParameters(effect);
        end
    end
end
