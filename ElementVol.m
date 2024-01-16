function varargout = ElementVol(elementNode,nodeCoor,varargin)
% An ODE-driven Level-Set Density Method for Topology Optimization of Shell Structures
% Zuyu Li, Email:lizuyu0091@163.com; Yang Liu, Email:285548029@qq.com
%% Element volume
nEle = size(elementNode,1);
vol = zeros(nEle,1);
elementCenter = zeros(nEle,3);%Element barycenter
for iEle = 1:nEle
    eleNode = elementNode(iEle,:);
    eleNodeCoor = nodeCoor(eleNode,:);
    %Transformation from global coordinate to local coordinate
    a = sum(eleNodeCoor([2;3],:)-eleNodeCoor([1;4],:),1)/2;
    b = sum(eleNodeCoor([3;4],:))/2-sum(eleNodeCoor([1;4],:))/2;
    cosX = a/sqrt(sum(a.^2));%[cos_xx', cos_xy', cos_xz']
    cosZ = [a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3),a(1)*b(2)-a(2)*b(1)];
    cosZ = cosZ/sqrt(sum(cosZ.^2));%[cos_zx', cos_zy', cos_zz']
    cosY = [cosZ(2)*cosX(3)-cosX(2)*cosZ(3),cosZ(3)*cosX(1)-cosX(3)*cosZ(1),cosZ(1)*cosX(2)-cosX(1)*cosZ(2)];%[cos_yx', cos_yy', cos_yz']
    transCos = [cosX;cosY;cosZ];
    transT1 = blkdiag(transCos,transCos,transCos,transCos);
    eleNodeCoor = eleNodeCoor';
    eleNodeCoor = reshape(transT1*eleNodeCoor(:),3,4)';
    eleCenter_local = [PolygonBaryCenter(eleNodeCoor(:,1:2)),mean(eleNodeCoor(:,3))];
    elementCenter(iEle,:) = transCos'*eleCenter_local';
    eleNodeCoor = eleNodeCoor(:,1:2);
    %Area calculation under plane coordinate
    a = eleNodeCoor(2,:)-eleNodeCoor(1,:);
    b = eleNodeCoor(4,:)-eleNodeCoor(1,:);
    c = eleNodeCoor(2,:)-eleNodeCoor(3,:);
    d = eleNodeCoor(4,:)-eleNodeCoor(3,:);
    vol(iEle,1) = abs(a(:,1).*b(:,2)-a(:,2).*b(:,1))/2+abs(c(:,1).*d(:,2)-c(:,2).*d(:,1))/2;
end
varargout{1} = vol;
varargout{2} = elementCenter;
end