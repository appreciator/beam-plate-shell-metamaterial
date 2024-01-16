function [eleComp,eleSen] = ElementResponse_Density(penal,eleNodeDisp,eleNodeCoor,eleNodeLS,eleMat,eleProfile)
% An ODE-driven Level-Set Density Method for Topology Optimization of Shell Structures
% Zuyu Li, Email:lizuyu0091@163.com; Yang Liu, Email:285548029@qq.com
%% Transformation from global coordinate to local coordinate
eleNodeCoor0 = eleNodeCoor;
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
eleNodeCoor = eleNodeCoor(:,1:2);
if min(eleNodeLS)>=0
    eleMat.E = eleMat.youngE;%Solid
    sen = penal*(eleMat.youngE-eleMat.youngEmin);
elseif max(eleNodeLS)<=0 && any(eleNodeLS<0)
    eleMat.E = eleMat.youngEmin;%Void
    sen = 0;
else
    [~,indexLS] = sort(eleNodeLS);
    if sum(eleNodeLS>0)==1 && any(eleNodeLS<0)
        %Triangle
        iP1 = indexLS(end);
        if iP1==1
            iP2 = 2; iP3 = 4;
        elseif iP1==2
            iP2 = 3; iP3 = 1;
        elseif iP1==3
            iP2 = 4; iP3 = 2;
        elseif iP1==4
            iP2 = 1; iP3 = 3;
        end
        %Local coordinate of polygon
        if eleNodeLS(iP1)>eleNodeLS(iP2)
            point2 = eleNodeCoor(iP1,:)-(eleNodeCoor(iP1,:)-eleNodeCoor(iP2,:))*eleNodeLS(iP1)/(eleNodeLS(iP1)-eleNodeLS(iP2));
        else
            point2 = eleNodeCoor(iP2,:)-(eleNodeCoor(iP2,:)-eleNodeCoor(iP1,:))*eleNodeLS(iP2)/(eleNodeLS(iP2)-eleNodeLS(iP1));
        end
        if eleNodeLS(iP1)>eleNodeLS(iP3)
            point3 = eleNodeCoor(iP1,:)-(eleNodeCoor(iP1,:)-eleNodeCoor(iP3,:))*eleNodeLS(iP1)/(eleNodeLS(iP1)-eleNodeLS(iP3));
        else
            point3 = eleNodeCoor(iP3,:)-(eleNodeCoor(iP3,:)-eleNodeCoor(iP1,:))*eleNodeLS(iP3)/(eleNodeLS(iP3)-eleNodeLS(iP1));
        end
        polygonPoint_solid = [eleNodeCoor(iP1,:);point2;point3];
    elseif sum(eleNodeLS>0)==2 && abs(indexLS(end)-indexLS(end-1))~=2 && any(eleNodeLS<0)
        %Quadrangle
        iP1 = min(indexLS(end),indexLS(end-1)); iP2 = max(indexLS(end),indexLS(end-1));
        if iP1==1&&iP2==2
            iP3 = 3; iP4 = 4;
        elseif iP1==2&&iP2==3
            iP3 = 4; iP4 = 1;
        elseif iP1==3&&iP2==4
            iP3 = 1; iP4 = 2;
        elseif iP1==1&&iP2==4
            iP1 = 4; iP2 = 1; iP3 = 2; iP4 = 3;
        end
        %Global coordinate of polygon
        if eleNodeLS(iP2)>eleNodeLS(iP3)
            point3 = eleNodeCoor(iP2,:)-(eleNodeCoor(iP2,:)-eleNodeCoor(iP3,:))*eleNodeLS(iP2)/(eleNodeLS(iP2)-eleNodeLS(iP3));
        else
            point3 = eleNodeCoor(iP3,:)-(eleNodeCoor(iP3,:)-eleNodeCoor(iP2,:))*eleNodeLS(iP3)/(eleNodeLS(iP3)-eleNodeLS(iP2));
        end
        if eleNodeLS(iP1)>eleNodeLS(iP4)
            point4 = eleNodeCoor(iP1,:)-(eleNodeCoor(iP1,:)-eleNodeCoor(iP4,:))*eleNodeLS(iP1)/(eleNodeLS(iP1)-eleNodeLS(iP4));
        else
            point4 = eleNodeCoor(iP4,:)-(eleNodeCoor(iP4,:)-eleNodeCoor(iP1,:))*eleNodeLS(iP4)/(eleNodeLS(iP4)-eleNodeLS(iP1));
        end
        polygonPoint_solid = [eleNodeCoor(iP1,:);eleNodeCoor(iP2,:);point3;point4];
    elseif (sum(eleNodeLS>0)==3 || (sum(eleNodeLS>0)==2 && abs(indexLS(end)-indexLS(end-1))==2)) && sum(eleNodeLS<0)==1
        %Pentagon
        iP1 = indexLS(1);
        if iP1==1
            iP2 = 2; iP3 = 4;
        elseif iP1==2
            iP2 = 3; iP3 = 1;
        elseif iP1==3
            iP2 = 4; iP3 = 2;
        elseif iP1==4
            iP2 = 1; iP3 = 3;
        end
        iP4 = setdiff(1:4,[iP1,iP2,iP3]);
        %Local coordinate of polygon
        if eleNodeLS(iP1)>eleNodeLS(iP2)
            point2 = eleNodeCoor(iP1,:)-(eleNodeCoor(iP1,:)-eleNodeCoor(iP2,:))*eleNodeLS(iP1)/(eleNodeLS(iP1)-eleNodeLS(iP2));
        else
            point2 = eleNodeCoor(iP2,:)-(eleNodeCoor(iP2,:)-eleNodeCoor(iP1,:))*eleNodeLS(iP2)/(eleNodeLS(iP2)-eleNodeLS(iP1));
        end
        if eleNodeLS(iP1)>eleNodeLS(iP3)
            point3 = eleNodeCoor(iP1,:)-(eleNodeCoor(iP1,:)-eleNodeCoor(iP3,:))*eleNodeLS(iP1)/(eleNodeLS(iP1)-eleNodeLS(iP3));
        else
            point3 = eleNodeCoor(iP3,:)-(eleNodeCoor(iP3,:)-eleNodeCoor(iP1,:))*eleNodeLS(iP3)/(eleNodeLS(iP3)-eleNodeLS(iP1));
        end
        polygonPoint_solid = [point2;eleNodeCoor(iP2,:);eleNodeCoor(iP4,:);eleNodeCoor(iP3,:);point3];
    elseif sum(eleNodeLS>0)==2 && abs(indexLS(end)-indexLS(end-1))==2 && sum(eleNodeLS<0)==2
        iP1 = min(indexLS(end),indexLS(end-1)); iP11 = max(indexLS(end),indexLS(end-1));
        if iP1==1&&iP11==3
            iP2 = 2; iP3 = 4; iP12 = 4; iP13 = 2;
        elseif iP1==2&&iP11==4
            iP2 = 3; iP3 = 1; iP12 = 1; iP13 = 3;
        end
        if eleNodeLS(iP1)>eleNodeLS(iP2)
            point2 = eleNodeCoor(iP1,:)-(eleNodeCoor(iP1,:)-eleNodeCoor(iP2,:))*eleNodeLS(iP1)/(eleNodeLS(iP1)-eleNodeLS(iP2));
        else
            point2 = eleNodeCoor(iP2,:)-(eleNodeCoor(iP2,:)-eleNodeCoor(iP1,:))*eleNodeLS(iP2)/(eleNodeLS(iP2)-eleNodeLS(iP1));
        end
        if eleNodeLS(iP1)>eleNodeLS(iP3)
            point3 = eleNodeCoor(iP1,:)-(eleNodeCoor(iP1,:)-eleNodeCoor(iP3,:))*eleNodeLS(iP1)/(eleNodeLS(iP1)-eleNodeLS(iP3));
        else
            point3 = eleNodeCoor(iP3,:)-(eleNodeCoor(iP3,:)-eleNodeCoor(iP1,:))*eleNodeLS(iP3)/(eleNodeLS(iP3)-eleNodeLS(iP1));
        end
        if eleNodeLS(iP11)>eleNodeLS(iP12)
            point12 = eleNodeCoor(iP11,:)-(eleNodeCoor(iP11,:)-eleNodeCoor(iP12,:))*eleNodeLS(iP11)/(eleNodeLS(iP11)-eleNodeLS(iP12));
        else
            point12 = eleNodeCoor(iP12,:)-(eleNodeCoor(iP12,:)-eleNodeCoor(iP11,:))*eleNodeLS(iP12)/(eleNodeLS(iP12)-eleNodeLS(iP11));
        end
        if eleNodeLS(iP11)>eleNodeLS(iP13)
            point13 = eleNodeCoor(iP11,:)-(eleNodeCoor(iP11,:)-eleNodeCoor(iP13,:))*eleNodeLS(iP11)/(eleNodeLS(iP11)-eleNodeLS(iP13));
        else
            point13 = eleNodeCoor(iP13,:)-(eleNodeCoor(iP13,:)-eleNodeCoor(iP11,:))*eleNodeLS(iP13)/(eleNodeLS(iP13)-eleNodeLS(iP11));
        end
        if sum(eleNodeLS([iP1;iP11]))>=sum(abs(eleNodeLS))/2
            %Hexagon
            %Global coordinate of polygon
            polygonPoint_solid = [eleNodeCoor(iP1,:);point2;point13;eleNodeCoor(iP11,:);point12;point3];
        else
            %Two Triangles
            %Local coordinate of polygon
            polygonPoint_solid = [eleNodeCoor(iP1,:);point2;point3];
            polygonPoint_solid(:,:,2) = [eleNodeCoor(iP11,:);point12;point13];
        end
    end
    eleVol = 0;
    for iSolid = 1:size(polygonPoint_solid,3)
        [~,eleVol_iSolid] = PolygonBaryCenter(polygonPoint_solid(:,:,iSolid));
        eleVol = eleVol+eleVol_iSolid;%Area of solid under global coordinate
    end
    %Quad4
    elementNode = [1,2,3,4];
    a = eleNodeCoor(elementNode(:,2),:)-eleNodeCoor(elementNode(:,1),:);
    b = eleNodeCoor(elementNode(:,4),:)-eleNodeCoor(elementNode(:,1),:);
    c = eleNodeCoor(elementNode(:,2),:)-eleNodeCoor(elementNode(:,3),:);
    d = eleNodeCoor(elementNode(:,4),:)-eleNodeCoor(elementNode(:,3),:);
    vol = abs(a(:,1).*b(:,2)-a(:,2).*b(:,1))/2+abs(c(:,1).*d(:,2)-c(:,2).*d(:,1))/2;
    eleDen = eleVol/vol;
    eleMat.E = eleDen.^penal*(eleMat.youngE-eleMat.youngEmin)+eleMat.youngEmin;
    sen = penal*eleDen.^(penal-1)*(eleMat.youngE-eleMat.youngEmin);
end
eleNodeCoor = eleNodeCoor0;
[eleComp,eleSen] = ElementResponse(eleNodeDisp,eleNodeCoor,eleMat,sen,eleProfile);
end