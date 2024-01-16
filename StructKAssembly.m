function [structK,eleVol] = StructKAssembly(modelFEM,modelLS,modelPara)
% An ODE-driven Level-Set Density Method for Topology Optimization of Shell Structures
% Zuyu Li, Email:lizuyu0091@163.com; Yang Liu, Email:285548029@qq.com
%% Global stiffness matrix assembly
if isempty(gcp('nocreate'))
    parpool('local')
end
penal = modelPara.penal;
valLS =  modelLS.valLS;
nEle = size(modelFEM.elementNode,1);
eleVol = zeros(nEle,1);
nEDofK = (size(modelFEM.elementNode,2)*modelFEM.nNodeDof)^2;
elementNode = modelFEM.elementNode;
structDof = bsxfun(@plus,kron(elementNode',6*ones(6,1)),repmat((-5:1:0)',size(elementNode,2),1));
nNodeDof = 6;
nEleNode = size(elementNode,2);%Node number
nEleDof = nEleNode*nNodeDof;%Node freedom
ii = repmat((1:1:nEleDof)',nEleDof,1);
jj = repmat((1:1:nEleDof),nEleDof,1); jj = jj(:);
dofI = structDof(ii,:);
dofJ = structDof(jj,:);
dofK = zeros(nEDofK,nEle);
for iEle = 1:nEle
    eleNode = modelFEM.elementNode(iEle,:);
    eleNodeCoor = modelFEM.nodeCoor(eleNode,:);
    eleMat = modelFEM.material{modelFEM.section(iEle,1)};
    eleProfile = modelFEM.profile{modelFEM.section(iEle,1)};
    eleNodeLS =  valLS(eleNode,1);
    [eleK,eleVol(iEle,1)] = ElementK_Density(penal,eleNodeCoor,eleNodeLS,eleMat,eleProfile);
    dofK(:,iEle) = eleK(:);
end
structK = sparse(dofI,dofJ,dofK);
end