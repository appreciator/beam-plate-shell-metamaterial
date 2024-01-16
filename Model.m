function [modelFEM,modelLS] = Model(modelPara)
% An ODE-driven Level-Set Density Method for Topology Optimization of Shell Structures
% Zuyu Li, Email:lizuyu0091@163.com; Yang Liu, Email:285548029@qq.com
%% Model definition
switch modelPara.Model
    case 1
        timenum = 8/10;
        length = 4;
        width = 2;
        thick = 1;
        %% Mesh
        nEX = 100;
        nEY = 50;
        nNodeDof = 6;
        eLen = length/nEX;
        eWidth = width/nEY;
        [coorX,coorY] = meshgrid(-length/2:eLen:length/2,-width/2:eWidth:width/2);
        coorZ = -timenum*coorX.^2/4-timenum*coorY.^2/4;
        nNode = (nEX+1)*(nEY+1);
        nEle = nEX*nEY;
        nDof = nNode*nNodeDof;
        nodeCoor = [coorX(:),coorY(:),coorZ(:)];
        eleN1 = reshape(1:nNode,nEY+1,nEX+1);
        eleN1 = eleN1(1:end-1,1:end-1);
        eleN1 = eleN1(:);
        eleN2 = eleN1+(nEY+1);
        eleN3 = eleN2+1;
        eleN4 = eleN1+1;
        elementNode = [eleN1,eleN2,eleN3,eleN4];
        %% Element volume and barycenter of the element
        [eleVol,elementCenter] = ElementVol(elementNode,nodeCoor,'barycenter' );
        vol_total = sum(eleVol);%Total volume
        %% Material property
        section = ones(nEle,1);
        material{1}.youngE = 1;%Young's modulus for solid material 
        material{1}.youngEmin = material{1}.youngE*1e-4;%Young's modulus for weak material 
        material{1}.E =  material{1}.youngE;
        material{1}.possionMu = 0.3;%Poisson's ratio
        material{1}.density = 1;%Material density for solid material 
        material{1}.densitymin = material{1}.density*1e-8;%Material density for void material 
        profile{1}.thick = thick;%Thickness
        %% Displacement boundary
        fixNode = find(nodeCoor(:,1)==-length/2);
        fixDof = repmat(6*fixNode,1,6)+repmat(0:-1:-5,size(fixNode,1),1); fixDof = fixDof(:);
        freeDof = setdiff(1:nDof,fixDof);
        %% Force boundary
        fixExternalForce = zeros(nDof,1);
        forceNode = find((nodeCoor(:,2)==0) & nodeCoor(:,1)==length/2);
        fixExternalForce(6*forceNode-3,1) = -1;
        %% FEA model
        modelFEM = struct('nNodeDof',nNodeDof,'nNode',nNode,'nEle',nEle,'nDof',nDof,'thick',thick,...
            'nodeCoor',nodeCoor,'elementNode',elementNode,'eleVol',eleVol,'vol_total',vol_total,'section',section,...
            'fixDof',fixDof,'freeDof',freeDof,'fixExternalForce',fixExternalForce,'elementCenter',elementCenter);
        modelFEM.material = material;
        modelFEM.profile = profile;
        %% Design domain
        fixDesEle = [];%Non-designable region
        freeDesEle = setdiff(1:nEle,fixDesEle);
        freeDesNode = unique(elementNode(freeDesEle,:));
        freeDesNodeIndex = zeros(nNode,1);
        freeDesNodeIndex(freeDesNode,1) = (1:numel(freeDesNode))';
        freeDesNodeCoor = nodeCoor(freeDesNode,:);
        freeDesEleNode = freeDesNodeIndex(elementNode(freeDesEle,:));
        %% Level-set function
        axisLim = [-length/2,length/2,-width/2,width/2,-0.5,0];
        charL = sqrt(length^2+width^2);
        valLS = zeros(nNode,1);%Initial level-set function
        vol = sum(eleVol);
        vol_max = 0.5*vol_total;%Prescribed material volume
        modelLS = struct('axisLim',axisLim,'nodeCoor',nodeCoor,'elementNode',elementNode,...
            'fixDesEle',fixDesEle,'freeDesEle',freeDesEle,'freeDesNode',freeDesNode,'freeDesNodeCoor',freeDesNodeCoor,...
            'vol_total',vol_total,'vol_max',vol_max,'vol_init',vol,'vol',vol,...
            'charL',charL,'freeDesEleNode',freeDesEleNode);
        modelLS.valLS_init = valLS; modelLS.valLS = valLS;
end
end