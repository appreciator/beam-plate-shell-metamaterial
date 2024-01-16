function [centerPoint,polyArea] = PolygonBaryCenter(polyPoint)
% An ODE-driven Level-Set Density Method for Topology Optimization of Shell Structures
% Zuyu Li, Email:lizuyu0091@163.com; Yang Liu, Email:285548029@qq.com
polyPoint2 = [polyPoint(2:end,:);polyPoint(1,:)];
detXY = polyPoint(:,1).*polyPoint2(:,2)-polyPoint(:,2).*polyPoint2(:,1);
facX = sum((polyPoint(:,1)+polyPoint2(:,1)).*detXY);
facY = sum((polyPoint(:,2)+polyPoint2(:,2)).*detXY);
facDet = sum(detXY);
centerPoint = [facX/facDet/3,facY/facDet/3];
polyArea = abs(facDet)/2;
end