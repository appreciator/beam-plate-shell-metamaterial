function [eleComp,eleSen] = ElementResponse(eleNodeDisp,eleNodeCoor,eleMat,sen,eleProfile)
% An ODE-driven Level-Set Density Method for Topology Optimization of Shell Structures
% Zuyu Li, Email:lizuyu0091@163.com; Yang Liu, Email:285548029@qq.com
%% DKQ+Q4 with drilling dof
%Transformation from global coordinate to local coordinate
a = sum(eleNodeCoor([2;3],:)-eleNodeCoor([1;4],:),1)/2;
b = sum(eleNodeCoor([3;4],:))/2-sum(eleNodeCoor([1;4],:))/2;
cosX = a/sqrt(sum(a.^2));%[cos_xx', cos_xy', cos_xz']
cosZ = [a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3),a(1)*b(2)-a(2)*b(1)];
cosZ = cosZ/sqrt(sum(cosZ.^2));%[cos_zx', cos_zy', cos_zz']
cosY = [cosZ(2)*cosX(3)-cosX(2)*cosZ(3),cosZ(3)*cosX(1)-cosX(3)*cosZ(1),cosZ(1)*cosX(2)-cosX(1)*cosZ(2)];%[cos_yx', cos_yy', cos_yz']
transCos = [cosX;cosY;cosZ];
transT1 = blkdiag(transCos,transCos,transCos,transCos);
transCos = [transCos,zeros(3,3);zeros(3,3),transCos];
transT2 = blkdiag(transCos,transCos,transCos,transCos);
eleNodeCoor = eleNodeCoor';
eleNodeCoor = reshape(transT1*eleNodeCoor(:),3,4)';
h1 = (eleNodeCoor(2,3)+eleNodeCoor(4,3)-eleNodeCoor(1,3)-eleNodeCoor(3,3))/4;
h = [h1;-h1;h1;-h1];
eleNodeCoor(:,3) = eleNodeCoor(:,3)+h;
transW1 = diag(ones(6,1)); transW1(1,5)=h(1); transW1(2,4)=-h(1);
transW2 = diag(ones(6,1)); transW2(1,5)=h(2); transW2(2,4)=-h(2);
transW3 = diag(ones(6,1)); transW3(1,5)=h(3); transW3(2,4)=-h(3);
transW4 = diag(ones(6,1)); transW4(1,5)=h(4); transW4(2,4)=-h(4);
transW = blkdiag(transW1,transW2,transW3,transW4);
eleNodeDisp = transW*transT2*eleNodeDisp;
eleK = zeros(24,24);
eleK_ps = zeros(12,12);
eleK_bp = zeros(12,12);
eleS = zeros(24,24);
eleS_ps = zeros(12,12);
eleS_bp = zeros(12,12);
E = eleMat.E;
S = sen;
mu = eleMat.possionMu;
t = eleProfile.thick;%Thickness
D = t*E/(1-mu^2)*[1,mu,0;mu,1,0;0,0,(1-mu)/2];
Db = t^3/12*E/(1-mu^2)*[1,mu,0;mu,1,0;0,0,(1-mu)/2];
DS = t*S/(1-mu^2)*[1,mu,0;mu,1,0;0,0,(1-mu)/2];
DSb = t^3/12*S/(1-mu^2)*[1,mu,0;mu,1,0;0,0,(1-mu)/2];
gaussP = [-1/sqrt(3),1/sqrt(3)];
gaussW = [1,1];
for iGauss = 1:2
    for jGauss = 1:2
        ksi = gaussP(iGauss);
        eta = gaussP(jGauss);
        [B,J] = EleBJacob([ksi,eta],eleNodeCoor(:,1:2),2);
        eleK_ps = eleK_ps + gaussW(iGauss)*gaussW(jGauss)*B'*D*B*det(J);%Plane stress
        eleS_ps = eleS_ps + gaussW(iGauss)*gaussW(jGauss)*B'*DS*B*det(J);
        [B,J] = EleBJacob([ksi,eta],eleNodeCoor(:,1:2),1);
        eleK_bp = eleK_bp + gaussW(iGauss)*gaussW(jGauss)*B'*Db*B*det(J);%Bending of thin plate
        eleS_bp = eleS_bp + gaussW(iGauss)*gaussW(jGauss)*B'*DSb*B*det(J);
    end
end
index_ps = repmat([1;2;6;7;8;12;13;14;18;19;20;24;[1;2;6;7;8;12;13;14;18;19;20;24]+24;[1;2;6;7;8;12;13;14;18;19;20;24]+120],1,4)+kron(144*ones(36,1),[0,1,2,3]);
eleK(index_ps) = eleK_ps(:);
eleS(index_ps) = eleS_ps(:);
index_bp = repmat([51;52;53;57;58;59;63;64;65;69;70;71],1,12)+kron(24*ones(12,1),[0,1,2,6,7,8,12,13,14,18,19,20]);
eleK(index_bp) = eleK_bp(:);%Element stiffness under local coordinate
eleS(index_bp) = eleS_bp(:);
eleComp = eleNodeDisp'*eleK*eleNodeDisp/2;%Strain energy density
eleSen = eleNodeDisp'*eleS*eleNodeDisp/2;
end