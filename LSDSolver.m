function LSDSolver(modelFEM_init,modelLS_init,modelResponse_init,modelPara)
% An ODE-driven Level-Set Density Method for Topology Optimization of Shell Structures
% Zuyu Li, Email:lizuyu0091@163.com; Yang Liu, Email:285548029@qq.com
modelFEM = modelFEM_init;
modelLS = modelLS_init;
modelResponse = modelResponse_init;
%% Algorithm parameter
disp(['Initial Structure, Comp: ',num2str(modelResponse.compliance),', Vol: ',num2str(modelLS.vol)]);
objComp = zeros(1+modelPara.nLoop,1); objComp(1,1) = 2*modelResponse.compliance;
constrVol = zeros(1+modelPara.nLoop,1); constrVol(1,1) = modelLS.vol/modelFEM_init.vol_total;
gLambda = zeros(modelPara.nLoop,2);
lambda = 0;
%% Level-set evolution
for iLoop = 1:modelPara.nLoop
    structTD =  modelResponse.eleSen;
    c = sum(modelFEM.eleVol(modelLS.freeDesEle))/sum(abs(structTD(modelLS.freeDesEle)));
    structTDN = zeros(modelFEM.nNode,1);
    structTDNCoeff = zeros(modelFEM.nNode,1);
    for iEle = 1:modelFEM.nEle
        eleNode = modelFEM.elementNode(iEle,:);
        structTDN(eleNode,1) = structTDN(eleNode,1)+structTD(iEle,1)/modelFEM.eleVol(iEle,1);
        structTDNCoeff(eleNode,1) = structTDNCoeff(eleNode,1)+1;
    end
    structTDN = structTDN./structTDNCoeff;
    %% Begain optimization
    figure(1); clf;
    Plot3DSurf(modelLS.nodeCoor,modelLS.valLS,[0,0],modelLS.elementNode,modelLS.axisLim);
    %Volume constraint
    Gmax = modelLS.vol_max/modelFEM_init.vol_total+(1-modelLS.vol_max/modelFEM_init.vol_total)*max(0,1-iLoop/modelPara.nVolRelax);
    lambda0 = -2e2/c; lambda1 = 2e2/c; iBisec = 0;
    while iBisec < 100
        lambda = (lambda0+lambda1)/2;
        modelLS0.valLS = modelLS.valLS+modelPara.dt*(structTDN-lambda);%Level-set update
        vol0 = (sum(modelLS0.valLS>=0,1))/modelFEM.nNode;
        if abs((vol0-Gmax)/Gmax)<=1e-3
            break;
        end
        if (vol0-Gmax)/Gmax>0
            lambda0 = lambda;
        else
            lambda1 = lambda;
        end
        iBisec = iBisec+1;
    end
    gLambda(iLoop,1) = lambda*c;
    gLambda(iLoop,2) = lambda*median(modelFEM.eleVol(modelLS.freeDesEle)./abs(structTD(modelLS.freeDesEle)));
    %Update level-set
    if iLoop>modelPara.nVolRelax && mod(iLoop,modelPara.nReg)==0 
        modelLS.valLS = modelPara.reg*modelLS0.valLS;%Regularization
    else
        modelLS.valLS = modelLS0.valLS;
    end 
    %Update structural response
    modelResponse = FEMSolver(modelFEM,modelLS,modelPara);
    modelLS.vol = modelResponse.vol;
    %Update optimization result
    disp(['Evolution No: ',num2str(iLoop),', Comp: ',num2str(modelResponse.compliance),', Vol: ',num2str(modelLS.vol)]);
    objComp(1+iLoop,1) = 2*modelResponse.compliance;
    constrVol(1+iLoop,1) = modelLS.vol/modelFEM_init.vol_total;
end
end