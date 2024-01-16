function Plot3DSurf(varargin)
% An ODE-driven Level-Set Density Method for Topology Optimization of Shell Structures
% Zuyu Li, Email:lizuyu0091@163.com; Yang Liu, Email:285548029@qq.com
%% Contour filling of 3D surface
nodeCoor = varargin{1};
val = varargin{2};
level = varargin{3};
elementNode = varargin{4};
axisLim = varargin{5};
nEle = size(elementNode,1);
nNode = size(nodeCoor,1);
nPatch = 0;
isoEleNode = zeros(nEle*4,3);
isoNodeCoor = [nodeCoor;zeros(nEle*4,size(nodeCoor,2))];
hold on;
set(gcf,'Renderer','zbuffer','Color',[1,1,1],'Position',[60,60,1200,600]);
view(3); axis equal; axis manual; axis vis3d; axis(axisLim); axis off;
for iEle = 1:nEle
    eleNode = elementNode(iEle,:);
    nodeC = nodeCoor(eleNode,:);
    eleNodeVal =  val(eleNode,1);
    if min(eleNodeVal)>=0
        isoEleNode(nPatch+1,:) = [eleNode(1),eleNode(2),eleNode(3)];
        isoEleNode(nPatch+2,:) = [eleNode(1),eleNode(3),eleNode(4)];
        nPatch = nPatch+2;
    else
        %Cross Element
        [~,indexVal] = sort(eleNodeVal);
        if sum(eleNodeVal>level(1))==1 && any(eleNodeVal<level(1))
            %Triangle
            iP1 = indexVal(end);
            if iP1==1
                iP2 = 2; iP3 = 4;
            elseif iP1==2
                iP2 = 3; iP3 = 1;
            elseif iP1==3
                iP2 = 4; iP3 = 2;
            elseif iP1==4
                iP2 = 1; iP3 = 3;
            end
            if eleNodeVal(iP1)>eleNodeVal(iP2)
                point2 = nodeC(iP1,:)-(nodeC(iP1,:)-nodeC(iP2,:))*eleNodeVal(iP1)/(eleNodeVal(iP1)-eleNodeVal(iP2));
            else
                point2 = nodeC(iP2,:)-(nodeC(iP2,:)-nodeC(iP1,:))*eleNodeVal(iP2)/(eleNodeVal(iP2)-eleNodeVal(iP1));
            end
            if eleNodeVal(iP1)>eleNodeVal(iP3)
                point3 = nodeC(iP1,:)-(nodeC(iP1,:)-nodeC(iP3,:))*eleNodeVal(iP1)/(eleNodeVal(iP1)-eleNodeVal(iP3));
            else
                point3 = nodeC(iP3,:)-(nodeC(iP3,:)-nodeC(iP1,:))*eleNodeVal(iP3)/(eleNodeVal(iP3)-eleNodeVal(iP1));
            end
            isoEleNode(nPatch+1,:) = [eleNode(iP1),nNode+1,nNode+2];
            nPatch = nPatch+1;
            isoNodeCoor(nNode+1:nNode+2,:) = [point2;point3];
            nNode = nNode+2;
        elseif sum(eleNodeVal>level(1))==2 && abs(indexVal(end)-indexVal(end-1))~=2 && any(eleNodeVal<level(1))
            %Quadrangle
            iP1 = min(indexVal(end),indexVal(end-1)); iP2 = max(indexVal(end),indexVal(end-1));
            if iP1==1&&iP2==2
                iP3 = 3; iP4 = 4;
            elseif iP1==2&&iP2==3
                iP3 = 4; iP4 = 1;
            elseif iP1==3&&iP2==4
                iP3 = 1; iP4 = 2;
            elseif iP1==1&&iP2==4
                iP1 = 4; iP2 = 1; iP3 = 2; iP4 = 3;
            end
            if eleNodeVal(iP2)>eleNodeVal(iP3)
                point3 = nodeC(iP2,:)-(nodeC(iP2,:)-nodeC(iP3,:))*eleNodeVal(iP2)/(eleNodeVal(iP2)-eleNodeVal(iP3));
            else
                point3 = nodeC(iP3,:)-(nodeC(iP3,:)-nodeC(iP2,:))*eleNodeVal(iP3)/(eleNodeVal(iP3)-eleNodeVal(iP2));
            end
            if eleNodeVal(iP1)>eleNodeVal(iP4)
                point4 = nodeC(iP1,:)-(nodeC(iP1,:)-nodeC(iP4,:))*eleNodeVal(iP1)/(eleNodeVal(iP1)-eleNodeVal(iP4));
            else
                point4 = nodeC(iP4,:)-(nodeC(iP4,:)-nodeC(iP1,:))*eleNodeVal(iP4)/(eleNodeVal(iP4)-eleNodeVal(iP1));
            end
            isoEleNode(nPatch+1,:) = [eleNode(iP1),eleNode(iP2),nNode+2];
            isoEleNode(nPatch+2,:) = [eleNode(iP2),nNode+1,nNode+2];
            nPatch = nPatch+2;
            isoNodeCoor(nNode+1:nNode+2,:) = [point3;point4];
            nNode = nNode+2;
        elseif (sum(eleNodeVal>level(1))==3 || (sum(eleNodeVal>level(1))==2 && abs(indexVal(end)-indexVal(end-1))==2)) && sum(eleNodeVal<level(1))==1
            %Pentagon
            iP1 = indexVal(1);
            if iP1==1
                iP2 = 2; iP3 = 4; iP4 = 3;
            elseif iP1==2
                iP2 = 3; iP3 = 1; iP4 = 4;
            elseif iP1==3
                iP2 = 4; iP3 = 2; iP4 = 1;
            elseif iP1==4
                iP2 = 1; iP3 = 3; iP4 = 2;
            end
            if eleNodeVal(iP1)>eleNodeVal(iP2)
                point2 = nodeC(iP1,:)-(nodeC(iP1,:)-nodeC(iP2,:))*eleNodeVal(iP1)/(eleNodeVal(iP1)-eleNodeVal(iP2));
            else
                point2 = nodeC(iP2,:)-(nodeC(iP2,:)-nodeC(iP1,:))*eleNodeVal(iP2)/(eleNodeVal(iP2)-eleNodeVal(iP1));
            end
            if eleNodeVal(iP1)>eleNodeVal(iP3)
                point3 = nodeC(iP1,:)-(nodeC(iP1,:)-nodeC(iP3,:))*eleNodeVal(iP1)/(eleNodeVal(iP1)-eleNodeVal(iP3));
            else
                point3 = nodeC(iP3,:)-(nodeC(iP3,:)-nodeC(iP1,:))*eleNodeVal(iP3)/(eleNodeVal(iP3)-eleNodeVal(iP1));
            end
            isoEleNode(nPatch+1,:) = [eleNode(iP2),eleNode(iP4),nNode+1];
            isoEleNode(nPatch+2,:) = [nNode+1,eleNode(iP4),nNode+2];
            isoEleNode(nPatch+3,:) = [nNode+2,eleNode(iP4),eleNode(iP3)];
            nPatch = nPatch+3;
            isoNodeCoor(nNode+1:nNode+2,:) = [point2;point3];
            nNode = nNode+2;
        elseif sum(eleNodeVal>level(1))==2 && abs(indexVal(end)-indexVal(end-1))==2 && sum(eleNodeVal<level(1))==2
            iP1 = min(indexVal(end),indexVal(end-1)); iP11 = max(indexVal(end),indexVal(end-1));
            if iP1==1&&iP11==3
                iP2 = 2; iP3 = 4; iP12 = 4; iP13 = 2;
            elseif iP1==2&&iP11==4
                iP2 = 3; iP3 = 1; iP12 = 1; iP13 = 3;
            end
            if eleNodeVal(iP1)>eleNodeVal(iP2)
                point2 = nodeC(iP1,:)-(nodeC(iP1,:)-nodeC(iP2,:))*eleNodeVal(iP1)/(eleNodeVal(iP1)-eleNodeVal(iP2));
            else
                point2 = nodeC(iP2,:)-(nodeC(iP2,:)-nodeC(iP1,:))*eleNodeVal(iP2)/(eleNodeVal(iP2)-eleNodeVal(iP1));
            end
            if eleNodeVal(iP1)>eleNodeVal(iP3)
                point3 = nodeC(iP1,:)-(nodeC(iP1,:)-nodeC(iP3,:))*eleNodeVal(iP1)/(eleNodeVal(iP1)-eleNodeVal(iP3));
            else
                point3 = nodeC(iP3,:)-(nodeC(iP3,:)-nodeC(iP1,:))*eleNodeVal(iP3)/(eleNodeVal(iP3)-eleNodeVal(iP1));
            end
            if eleNodeVal(iP11)>eleNodeVal(iP12)
                point12 = nodeC(iP11,:)-(nodeC(iP11,:)-nodeC(iP12,:))*eleNodeVal(iP11)/(eleNodeVal(iP11)-eleNodeVal(iP12));
            else
                point12 = nodeC(iP12,:)-(nodeC(iP12,:)-nodeC(iP11,:))*eleNodeVal(iP12)/(eleNodeVal(iP12)-eleNodeVal(iP11));
            end
            if eleNodeVal(iP11)>eleNodeVal(iP13)
                point13 = nodeC(iP11,:)-(nodeC(iP11,:)-nodeC(iP13,:))*eleNodeVal(iP11)/(eleNodeVal(iP11)-eleNodeVal(iP13));
            else
                point13 = nodeC(iP13,:)-(nodeC(iP13,:)-nodeC(iP11,:))*eleNodeVal(iP13)/(eleNodeVal(iP13)-eleNodeVal(iP11));
            end
            if sum(eleNodeVal([iP1;iP11]))>=sum(abs(eleNodeVal))/2
                %Hexagon
                isoEleNode(nPatch+1,:) = [nNode+3,nNode+2,eleNode(iP11)];
                isoEleNode(nPatch+2,:) = [nNode+2,eleNode(iP1),eleNode(iP11)];
                isoEleNode(nPatch+3,:) = [eleNode(iP1),nNode+1,eleNode(iP11)];
                isoEleNode(nPatch+4,:) = [nNode+1,nNode+4,eleNode(iP11)];
                nPatch = nPatch+4;
                isoNodeCoor(nNode+1:nNode+4,:) = [point2;point3;point12;point13];
                nNode = nNode+4;
            else
                %Two Triangles
                isoEleNode(nPatch+1,:) = [eleNode(iP1),nNode+1,nNode+2];
                isoEleNode(nPatch+2,:) = [eleNode(iP11),nNode+3,nNode+4];
                nPatch = nPatch+2;
                isoNodeCoor(nNode+1:nNode+4,:) = [point2;point3;point12;point13];
                nNode = nNode+4;
            end
        end
    end
end
isoEleNode = isoEleNode(1:nPatch,:);
isoNodeCoor = isoNodeCoor(1:nNode,:);
[isoNodeCoor,~,indexn] = unique(isoNodeCoor,'rows');
isoEleNode = indexn(isoEleNode);
if size(isoEleNode,2)==1
    isoEleNode = isoEleNode';
end
isoEleNode = isoEleNode(:,[1,3,2]);
trisurf(isoEleNode,isoNodeCoor(:,1),isoNodeCoor(:,2),isoNodeCoor(:,3),'EdgeColor','none','FaceColor',[0.5,0.5,0.5],'Facelighting','phong');
camlight('headlight');
set(findobj(gca,'type','patch'),...
    'FaceLighting','phong','AmbientStrength',0.3,'DiffuseStrength',0.8,...
    'SpecularStrength',0.5,'SpecularColorReflectance',1,'SpecularExponent',15,'BackFaceLighting','unlit');
hold off;
end