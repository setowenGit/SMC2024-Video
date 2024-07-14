function modelPD = spd_EM_TPGMM(x, modelPD, nbData, nbSamples, in, out)

nbIterEM = 10; % Number of iteration for the EM algorithm
nbIter = 10; % Number of iteration for the Gauss Newton algorithm (Riemannian manifold)

% EM for SPD matrices manifold
for f = 1:modelPD.nbFrames
    L = zeros(modelPD.nbStates, nbData*nbSamples);% 6*400
    xts = zeros(modelPD.nbVarVec, nbData*nbSamples, modelPD.nbStates);% 5*500*6
    for nb=1:nbIterEM
        fprintf('.');
        % E-step
        for i=1:modelPD.nbStates
            xts(in,:,i) = x(in,:,f)-repmat(modelPD.MuMan(in,i,f),1,nbData*nbSamples);
            xts(out,:,i) = logmap_vec(x(out,:,f), modelPD.MuMan(out,i,f));
            L(i,:) = modelPD.Priors(i,f) * gaussPDF(xts(:,:,i), modelPD.Mu(:,i,f), modelPD.Sigma(:,:,i,f));
        end
        % Responsibilities
        GAMMA = L ./ repmat(sum(L,1)+realmin, modelPD.nbStates, 1);
        H = GAMMA ./ repmat(sum(GAMMA,2)+realmin, 1, nbData*nbSamples);
        
        % M-step
        for i=1:modelPD.nbStates
            % Update Priors
            modelPD.Priors(i,f) = sum(GAMMA(i,:)) / (nbData*nbSamples);
            % Update MuMan
            for n=1:nbIter
                % Update on the tangent space
                uTmp = zeros(modelPD.nbVarVec,nbData*nbSamples);
                uTmp(in,:) = x(in,:,f) - repmat(modelPD.MuMan(in,i,f),1,nbData*nbSamples);
                uTmp(out,:) = logmap_vec(x(out,:,f), modelPD.MuMan(out,i,f));
                uTmpTot = sum(uTmp.*repmat(H(i,:),modelPD.nbVarVec,1),2);
                % Update on the manifold
                modelPD.MuMan(in,i,f) = uTmpTot(in,:) + modelPD.MuMan(in,i,f);
                modelPD.MuMan(out,i,f) = expmap_vec(uTmpTot(out,:), modelPD.MuMan(out,i,f));
            end
            % Update Sigma
            modelPD.Sigma(:,:,i,f) = uTmp * diag(H(i,:)) * uTmp' + eye(modelPD.nbVarVec) .* modelPD.params_diagRegFact;
        end
    end
    fprintf('\n');
end

% Eigendecomposition of Sigma
for f = 1:modelPD.nbFrames
    for i=1:modelPD.nbStates
        [V, D] = eig(modelPD.Sigma(:,:,i,f));
        modelPD.V(:,:,i,f) = real(V);
        modelPD.D(:,:,i,f) = real(D);
        modelPD.D_vec(:,i,f) = diag(modelPD.D(:,:,i,f));
    end
end

end
