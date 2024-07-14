%  Teleoperation-oriented Stiffness-adaptive GMM/GMR
%  demo
%  Author: Liwen Situ
%  Date: 24.07.24
%% env init
clear;
clc;
addpath('.\func\');
load_path = 'demonstration';
%% demonstration
nbSamples = 5;
demos = [];
for i = 1:nbSamples
    load(['.\' load_path '\demo_pos_stiff_' num2str(i) '.mat']);    
    demos = [demos, demo_pos_stiff];
end
nbData = size(demos,2)/nbSamples;
%% model parameters
modelPD.nbStates = 6;
modelPD.nbVar = 6;
modelPD.nbVarOut = 3;
modelPD.nbVarOutVec = modelPD.nbVarOut + modelPD.nbVarOut*(modelPD.nbVarOut-1)/2;
modelPD.nbVarVec = modelPD.nbVar - modelPD.nbVarOut + modelPD.nbVarOutVec;
modelPD.params_diagRegFact = 1E-4;
modelPD.nbFrames = 2;
%% process data
disp('Loading demonstration data');
X = zeros(modelPD.nbVar,modelPD.nbVar,nbData*nbSamples,modelPD.nbFrames); 
for n=1:nbSamples
    s(n).Data=[];
    dTmp = demos(1:3, 1+(n-1)*nbData:n*nbData);
    s(n).Data = [s(n).Data; dTmp];
    s(n).p(1).A = rotate_matrix_3d(0,0,0);
    s(n).p(1).b = s(n).Data(:, 1);
    s(n).p(2).A = rotate_matrix_3d(0,0,0);
    s(n).p(2).b = s(n).Data(:, end);
    for f = 1:modelPD.nbFrames
        fTmp = s(n).p(f).A \ (s(n).Data - repmat(s(n).p(f).b, 1, nbData));
        s(n).x(f).xx = fTmp;
        X(1,1,(n-1)*nbData+1:n*nbData,f) = reshape(fTmp(1, :), 1,1,nbData);
        X(2,2,(n-1)*nbData+1:n*nbData,f) = reshape(fTmp(2, :), 1,1,nbData);
        X(3,3,(n-1)*nbData+1:n*nbData,f) = reshape(fTmp(3, :), 1,1,nbData);
    end
    for f = 1:modelPD.nbFrames
        for num = 1:nbData
            X(4:6,4:6,(n-1)*nbData+num,f) = (s(n).p(f).A)' * vec2symmat(demos(4:9,(n-1)*nbData+num)) * s(n).p(f).A;
        end
    end
end
for f = 1:modelPD.nbFrames
    x(:,:,f) = [reshape(X(1,1,:,f),1,nbData*nbSamples); reshape(X(2,2,:,f),1,nbData*nbSamples);reshape(X(3,3,:,f),1,nbData*nbSamples);symmat2vec(X(4:6,4:6,:,f))]; %symmat2vec为用户自定义的，基于Mandel标记用于将matrix转为vector
end
%% Plot demonstration
is3Dplot = false;
clrmap = lines(nbSamples);
if is3Dplot == true
    figure('position',[10 10 1500 300],'color',[1 1 1]);
    for n=1:nbSamples
        subplot(1,nbSamples,n);hold on;
        for t=round(linspace(1,nbData,25))
            plotGMM3D([s(n).Data(1,t);s(n).Data(2,t);s(n).Data(3,t)], 4E-2*vec2symmat(demos(4:9,t+(n-1)*nbData)), clrmap(n,:), .4);
        end
        axis square;
        view(-5,85);
        xlabel('x');
        ylabel('y');
        zlabel('z');
    end
else
    figure('position',[10 10 1500 300],'color',[1 1 1]);
    for n=1:nbSamples
        subplot(1,nbSamples,n);
        hold on;
        for t=round(linspace(1,nbData,40))
            stiff_x_y = vec2symmat(demos(4:9,t+(n-1)*nbData));
            plotGMM([s(n).Data(1,t);s(n).Data(2,t)], 4E-2*stiff_x_y(1:2,1:2), clrmap(n,:), .4);
        end
        plot(s(n).Data(1,:),s(n).Data(2,:),'color', [0.5 0.5 0.5], 'Linewidth', 2)
        axis square;
        xlabel('x');
        ylabel('y');
    end

    figure('position',[10 10 1000 800],'color',[1 1 1]);
    for n=1:nbSamples
        subplot(nbSamples,1,n);
        hold on;
        for t = 1:30:nbData
            stiff_x_y = vec2symmat(demos(4:9,t+(n-1)*nbData));
            plotGMM([t*2E-3;0], 5E-2*stiff_x_y(1:2,1:2), clrmap(n,:), .4);
        end
        xlabel('t');
        ylabel('stiffness');
    end
    
end
figure;
hold on;
for n=1:nbSamples
    for t=round(linspace(1,nbData,20))
        scatter3(s(n).Data(1,:),s(n).Data(2,:),s(n).Data(3,:),'color',clrmap(n,:));
    end
end
grid on;
view(-37.5,30);
xlabel('x');
ylabel('y');
zlabel('z');
%% TOSA-GMM initial
disp('TOSA-GMM init');
in = 1:3;
outMat = 4:modelPD.nbVar;
out = 4:modelPD.nbVarVec;
modelPD = spd_init_TPGMM_kbins(x, modelPD, nbSamples, out);
%% TOSA-GMM EM
disp('TOSA-GMM EM');
modelPD.Mu = zeros(size(modelPD.MuMan));
modelPD = spd_EM_TPGMM(x, modelPD, nbData, nbSamples, in, out);
%% TOSA-GMR process
disp('TOSA-GMR process');
repro_n = '6';
load(['.\' load_path '\demo_pos_stiff_' repro_n '.mat'])
roate_x = 0;
roate_y = 0;
roate_z = 0;
dTmp = demo_pos_stiff(1:3,:);
xd = dTmp;
xInp(1).A = rotate_matrix_3d(roate_x,roate_y,roate_z);
xInp(1).b = dTmp(:,1);
xInp(2).A = rotate_matrix_3d(roate_x,roate_y,roate_z);
xInp(2).b = dTmp(:,end);
for f = 1:modelPD.nbFrames
    xIn(:, :, f) = xInp(f).A \ (dTmp - repmat(xInp(f).b, 1, nbData));
end
in = 1:3;
out = 4:modelPD.nbVarVec;
newmodelPD = spd_TPGMR(modelPD, xIn, in, out);
for f = 1:modelPD.nbFrames
    for num = 1:size(xd,2)
        newmodelPD.Mu_trans(:,:,num,f) = xInp(f).A * vec2symmat(newmodelPD.Mu(:, num, f)) * (xInp(f).A)';
        newmodelPD.Mu_trans_vec(:,num,f) = symmat2vec(newmodelPD.Mu_trans(:,:,num,f));
    end
end
nbIt = 10;
datanum = size(newmodelPD.Mu, 2);
for i = 1:datanum
    cov1 = trace(inv(newmodelPD.Sigma(:,:,i,1)));
    cov2 = trace(inv(newmodelPD.Sigma(:,:,i,2)));
    sum_weight = cov1 + cov2;
    newmodelPD.weight(1,i) = cov1 / sum_weight;
    newmodelPD.weight(2,i) = cov2 / sum_weight;
end

for i = 1:datanum
    M = newmodelPD.Mu_trans(:, :, i, 1);
    for n = 1:nbIt
        L = zeros(size(M,1), size(M,2));
        for f = 1:modelPD.nbFrames
            L = L + newmodelPD.weight(f,i) * logm(M^-.5 * newmodelPD.Mu_trans(:,:,i,f) * M^-.5);
        end
        M = M^.5 * expm(L) * M^.5;
    end
    newmodelPD.Result_mat(:,:,i) = M;
    newmodelPD.Result_vec(:,i) = symmat2vec(M);
end

%% Plot TOSA-GMM
disp('Finish and plot');
figure('position',[10 10 1500 700],'color',[1 1 1]);
subplot(1,2,1); hold on;
title('\fontsize{12}Stiffness TOSA-GMM means in Frame1');
clrmap = lines(modelPD.nbStates);
for n=1:nbSamples
	for t=round(linspace(1,nbData,40))
		plotGMM([X(1,1,t+(n-1)*nbData,1);X(2,2,t+(n-1)*nbData,1)], 4E-2*X(4:5,4:5,t+(n-1)*nbData,1), [.6 .6 .6], .1); % Scaled matrix!
        plot(X(1,1,t+(n-1)*nbData,1),X(2,2,t+(n-1)*nbData,1), 'color', [0.5 0.5 0.5], 'Linewidth', 2);
    end
    
end

for i=1:modelPD.nbStates
    sigmm = vec2symmat(modelPD.MuMan(out,i,1));
	plotGMM([modelPD.MuMan(1:2,i,1)], 4E-2*sigmm(1:2,1:2), clrmap(i,:), .3);
end
axis square;
set(gca, 'FontSize', 20)
xlabel('$x$', 'Fontsize', 28, 'Interpreter', 'Latex');
ylabel('$y$', 'Fontsize', 28, 'Interpreter', 'Latex');

subplot(1,2,2); hold on;
title('\fontsize{12}Stiffness TOSA-GMM means in Frame2');
clrmap = lines(modelPD.nbStates);
for n=1:nbSamples
	for t=round(linspace(1,nbData,40))
		plotGMM([X(1,1,t+(n-1)*nbData,2);X(2,2,t+(n-1)*nbData,2)], 4E-2*X(4:5,4:5,t+(n-1)*nbData,2), [.6 .6 .6], .1); % Scaled matrix!
        plot(X(1,1,t+(n-1)*nbData,2),X(2,2,t+(n-1)*nbData,2), 'color', [0.5 0.5 0.5], 'Linewidth', 2);
	end
end

for i=1:modelPD.nbStates
    sigmm = vec2symmat(modelPD.MuMan(out,i,2));
	plotGMM([modelPD.MuMan(1:2,i,2)], 4E-2*sigmm(1:2,1:2), clrmap(i,:), .3);
end
axis square;
set(gca, 'FontSize', 20)
xlabel('$x$', 'Fontsize', 28, 'Interpreter', 'Latex');
ylabel('$y$', 'Fontsize', 28, 'Interpreter', 'Latex');

%% Plot TOSA-GMR in frame
figure('position',[10 10 1800 600],'color',[1 1 1]);
subplot(1,2,1); hold on;
title('\fontsize{12}Reproduced stiffness in Frame1');
for t=1:15:nbData
    sigmm = vec2symmat(newmodelPD.Mu(:,t,1));
    plotGMM(xIn(1:2,t,1), 4E-2*sigmm(1:2,1:2), [0.2 0.8 0.2], .5, '-.', 2, 1);
end
axis square;
set(gca, 'FontSize', 20)
xlabel('$x$', 'Fontsize', 28, 'Interpreter', 'Latex');
ylabel('$y$', 'Fontsize', 28, 'Interpreter', 'Latex');

subplot(1,2,2); hold on;
title('\fontsize{12}Reproduced stiffness in Frame2');
for t=1:15:nbData
    sigmm = vec2symmat(newmodelPD.Mu(:,t,2));
    plotGMM(xIn(1:2,t,2), 4E-2*sigmm(1:2,1:2), [0.2 0.8 0.2], .5, '-.', 2, 1); % Scaled matrix!
end
axis square;
set(gca, 'FontSize', 20)
xlabel('$y$', 'Fontsize', 28, 'Interpreter', 'Latex');
ylabel('$z$', 'Fontsize', 28, 'Interpreter', 'Latex');

%% Plot result
figure('position',[10 10 1500 600],'color',[1 1 1]);
subplot(1,2,1);hold on;
title('\fontsize{12}Desired reproduction in 3D');
for t=1:15:nbData
    sigmm = vec2symmat(newmodelPD.Result_vec(:,t));
    plotGMM3D(xd(:,t), 4E-2*sigmm, [0.2 0.8 0.2], .4);
end
axis square;
view(-5,85);
set(gca, 'FontSize', 20)
xlabel('$x$', 'Fontsize', 28, 'Interpreter', 'Latex');
ylabel('$y$', 'Fontsize', 28, 'Interpreter', 'Latex');
zlabel('$z$', 'Fontsize', 28, 'Interpreter', 'Latex');

subplot(1,2,2);hold on;
title('\fontsize{12}Desired reproduction in 2D');
for t=1:15:nbData
    sigmm = vec2symmat(newmodelPD.Result_vec(:,t));
    plotGMM(xd(1:2,t), 4E-2*sigmm(1:2,1:2), [0.2 0.8 0.2], .5, '-.', 2, 1);
end
axis square;
set(gca, 'FontSize', 20)
xlabel('$x$', 'Fontsize', 28, 'Interpreter', 'Latex');
ylabel('$y$', 'Fontsize', 28, 'Interpreter', 'Latex');
%%
function R = rotate_matrix(theta)
    rad = theta * pi / 180;
    R = [cos(rad), -1*sin(rad); sin(rad), cos(rad)];
end

function R = rotate_matrix_3d(alpha, beta, gamma)
    rad = alpha * pi / 180;
    R = [1 0 0; 0 cos(rad) -sin(rad); 0 sin(rad) cos(rad)];
    rad = beta * pi /180;
    R = [cos(rad) 0 sin(rad); 0 1 0; -sin(rad) 0 cos(rad)] * R;
    rad = gamma * pi /180;
    R = [cos(rad) -sin(rad) 0; sin(rad) cos(rad) 0; 0 0 1] * R;
end