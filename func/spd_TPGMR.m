function newmodelPD = spd_TPGMR(modelPD, xIn, in, out)

nbData = size(xIn,2);
nbVarOut = length(out);
outMat = out(1):modelPD.nbVar;
nbIter = 2;

for f = 1:modelPD.nbFrames
    uhat = zeros(nbVarOut,nbData); 
    xhat = zeros(nbVarOut,nbData); 
    uOut = zeros(nbVarOut,modelPD.nbStates,nbData); 
    expSigma = zeros(nbVarOut,nbVarOut,nbData); 
    H = [];
    H_raw = [];
    for t = 1:nbData
        for i=1:modelPD.nbStates
            H(i,t) = modelPD.Priors(i,f) * gaussPDF(xIn(:,t,f)-modelPD.MuMan(in,i,f),modelPD.Mu(in,i,f), modelPD.Sigma(in,in,i,f));
        end
        H_raw(:,t) = H(:,t);
        H(:,t) = H(:,t) / sum(H(:,t)+realmin);
        if t==1
            [~,id] = max(H(:,t));
            xhat(:,t) = modelPD.MuMan(out,id,f);
        else
            xhat(:,t) = xhat(:,t-1);
        end
        for n=1:nbIter
            uhat(:,t) = zeros(nbVarOut,1);
            for i=1:modelPD.nbStates
                S1 = vec2symmat(modelPD.MuMan(out,i,f));
                S2 = vec2symmat(xhat(:,t));
                Ac = blkdiag(eye(length(in)), transp_operator(S1,S2));
                for j = 1:size(modelPD.V,2)
                    vMat(:,:,j,i) = blkdiag(blkdiag(diag(modelPD.V(in,j,i,f))),vec2symmat(modelPD.V(out,j,i,f)));
                    pvMat(:,:,j,i) = Ac * modelPD.D(j,j,i,f)^.5 * vMat(:,:,j,i) * Ac';
                    pV(:,j,i) = [diag(pvMat(in,in,j,i)); symmat2vec(pvMat(outMat,outMat,j,i))];
                end
                pSigma(:,:,i) = pV(:,:,i)*pV(:,:,i)';
                uOut(:,i,t) = logmap_vec(modelPD.MuMan(out,i,f), xhat(:,t)) + pSigma(out,in,i)/pSigma(in,in,i)*(xIn(:,t,f)-modelPD.MuMan(in,i,f));%公式（26）
                uhat(:,t) = uhat(:,t) + uOut(:,i,t) * H(i,t);
            end
            xhat(:,t) = expmap_vec(uhat(:,t), xhat(:,t));
        end 
        for i=1:modelPD.nbStates
            SigmaOutTmp = pSigma(out,out,i) - pSigma(out,in,i)/pSigma(in,in,i)*pSigma(in,out,i);
            expSigma(:,:,t) = expSigma(:,:,t) + H(i,t) * (SigmaOutTmp + uOut(:,i,t)*uOut(:,i,t)');
        end
        expSigma(:,:,t) = expSigma(:,:,t) - uhat(:,t)*uhat(:,t)';
    end
    newmodelPD.Mu(:,:,f) = xhat;
    newmodelPD.Sigma(:,:,:,f) = expSigma;
    newmodelPD.H(:,:,f) = H;
    newmodelPD.H_raw(:,:,f) = H_raw;
end
end