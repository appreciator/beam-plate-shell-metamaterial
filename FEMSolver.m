function modelResponse = FEMSolver(modelFEM,modelLS,modelPara)
% An ODE-driven Level-Set Density Method for Topology Optimization of Shell Structures
% Zuyu Li, Email:lizuyu0091@163.com; Yang Liu, Email:285548029@qq.com
%% Global stiffness matrix
structDisp = zeros(modelFEM.nDof,1);
[structK,eleVol] = StructKAssembly(modelFEM,modelLS,modelPara);
%% FEA solving
modelFEM.nodeForce = modelFEM.fixExternalForce;
%% Structural response
structDisp(modelFEM.freeDof,1) = structK(modelFEM.freeDof,modelFEM.freeDof)\modelFEM.nodeForce(modelFEM.freeDof,1);
[structComp,structSen] = StructResponse(modelFEM,modelLS,structDisp,modelPara);
modelResponse = struct('disp',structDisp,...
    'compliance',sum(structComp),'eleComp',structComp,'eleSen',structSen,'eleVol',eleVol,'vol',sum(eleVol));
end