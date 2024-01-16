function [structComp,structSen] = StructResponse(modelFEM,modelLS,structDisp,modelPara)
% An ODE-driven Level-Set Density Method for Topology Optimization of Shell Structures
% Zuyu Li, Email:lizuyu0091@163.com; Yang Liu, Email:285548029@qq.com
if isempty(gcp('nocreate'))
    parpool('local')
end
penal = modelPara.penal;
valLS =  modelLS.valLS;
nEle = size(modelFEM.elementNode,1);
structComp = zeros(nEle,1);%Strain energy density
structSen = zeros(nEle,1);%Strain energy density
for iEle = 1:nEle
    eleNode = modelFEM.elementNode(iEle,:);
    eleNodeCoor = modelFEM.nodeCoor(eleNode,:);
    eleMat = modelFEM.material{modelFEM.section(iEle,1)};
    eleProfile = modelFEM.profile{modelFEM.section(iEle,1)};
    eleDof = bsxfun(@plus,kron(eleNode',6*ones(6,1)),repmat((-5:1:0)',size(eleNode,2),1));
    eleNodeDisp = structDisp(eleDof,1);
    eleNodeLS =  valLS(eleNode,1);
    [structComp(iEle,1),structSen(iEle,1)] = ElementResponse_Density(penal,eleNodeDisp,eleNodeCoor,eleNodeLS,eleMat,eleProfile);
end
end