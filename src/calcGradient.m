function G = calcGradient(p, phi, ovData, A, S, edgeWeights)
global eps_;

NPoints=size(phi,1);
% 2D
G=zeros(NPoints, 2);

if size(S,2) > 1
    isUseSxyAsMainAxis=true;
    Sx=S(:,1);
    Sy=S(:,2);
else
    isUseSxyAsMainAxis=false;
    Sx=S;
    Sy=S;
end

ovData5=ovData(:,5);


for idxP=1:NPoints
    
    crtPhi=phi(idxP);
    crtCoordinate=p(idxP,:);
    ovEndPhi=phi(ovData5{idxP}{2});
    isOnZeroContour=(crtPhi > eps_ && any(ovEndPhi < -eps_)) ...
        || (crtPhi < -eps_ && any(ovEndPhi > eps_));
    % set isOnZeroContour = false here to disable explicit searching branch
    %isOnZeroContour = false;
    
    myOVEdges=ovData5{idxP}{1};
    myOVEdgesMid=ovData5{idxP}{3};
    myOVEdgesMidIndices=ovData5{idxP}{4};

    
    if isOnZeroContour
        [GExplicit,OK]=calcGradientViaExplicitGeometry(p, idxP, myOVEdgesMidIndices, crtCoordinate, crtPhi, phi, isUseSxyAsMainAxis);
        if OK
            G(idxP,:)=GExplicit;
        else
            error('calcGradientViaExplicitGeometry Error')
        end

    else
        
        changingAmount=ovEndPhi-crtPhi;
        changingRate=(changingAmount)./(myOVEdges(:,5));

        
        % find the vertex aiming at zero or aiming to upwind of S
        if isUseSxyAsMainAxis
            aimingDirection=-normr([Sx(idxP), Sy(idxP)]);
        else
            if crtPhi > 0
                [~,im]=min(changingRate);
            else
                [~,im]=max(changingRate);
            end
            aimingDirection=myOVEdges(im,3:4);
        end
        
        [GAiming1,OK1] = estimateGradientOnOneDirection(aimingDirection, myOVEdges, changingAmount,...
            myOVEdgesMid, myOVEdgesMidIndices, crtPhi, phi);
        [GAiming2,OK2] = estimateGradientOnOneDirection(-aimingDirection, myOVEdges, changingAmount,...
            myOVEdgesMid, myOVEdgesMidIndices, crtPhi, phi);
        
        if OK1 && OK2
            [~, tag]=Godunov(dot(GAiming2, aimingDirection, 2),dot(GAiming1, aimingDirection, 2), S(idxP));
            if tag == 1
                G(idxP,:)=GAiming1;
            elseif tag == -1
                G(idxP,:)=GAiming2;
            else
                G(idxP,:)=[0, 0];
            end
        elseif OK1 && (~OK2)
            G(idxP,:)=GAiming1;
        elseif (~OK1) && OK2
            G(idxP,:)=GAiming2;
        else
            error('aimingDirection unavailable');
        end
    end
end

end

function [G,OK]=calcGradientViaExplicitGeometry(p, idxP, myOVEdgesMidIndices, crtCoordinate, crtPhi, phi, isUseSxyAsMainAxis)
global eps_;
% zeroContourEdges: [x1, y1, x2, y2];
zeroContourEdges1=zeros(size(myOVEdgesMidIndices));
zeroContourEdges2=zeros(size(myOVEdgesMidIndices));
nZeroContourEdges=0;

for i = 1:size(myOVEdgesMidIndices,1)
    triangleOtherVertices=myOVEdgesMidIndices(i,:);
    phiArray=[crtPhi;phi(triangleOtherVertices)];
    
    % see if this triangle is on zero contour
    if sum(phiArray>eps_)>0 && sum(phiArray<-eps_)>0
        
        % see if we need to filter the triangle
        if isUseSxyAsMainAxis
            divergence=dot([Sx(idxP) Sy(idxP)], myOVEdgesMid(i, 3:4));
            
            % check for safe
            if abs(divergence) > 1
                error('divergence cannot > 1')
            end
            
            if divergence < divergenceThreshold + 0.1
                continue;
            end
        end
        
        nZeroContourEdges=nZeroContourEdges+1;
        
        % find the edge and add to zeroContourEdges
        [lineSeg1, lineSeg2] = findContourInTriangle(p([idxP, triangleOtherVertices],:), phiArray);
        zeroContourEdges1(nZeroContourEdges,:)=lineSeg1;
        zeroContourEdges2(nZeroContourEdges,:)=lineSeg2;
    end
end

zeroContourEdges1(nZeroContourEdges+1:end,:)=[];
zeroContourEdges2(nZeroContourEdges+1:end,:)=[];

% find the nearest point
nearestPoint=[];
nearestDistance=Inf;
for i=1:nZeroContourEdges
    [xy,dis,~]=distance2curve([zeroContourEdges1(i,:);zeroContourEdges2(i,:)],crtCoordinate);
    if dis < nearestDistance
        nearestPoint=xy;
        nearestDistance=dis;
    end
end

% calc real gradient
gradientDirection=normr(nearestPoint-crtCoordinate);
G=gradientDirection*(-crtPhi)/nearestDistance;

OK=true;
end

function [G] = calcGradient1Axis1Vertex(phi,crtPhi,OVTags,OVMidTags,dA,weight)
changingAmount=(phi(OVTags)-crtPhi);

phiMid=(phi(OVMidTags(:,1))+phi(OVMidTags(:,2)))./2;
changingAmountMid=(phiMid-crtPhi);

b=[changingAmount;changingAmountMid] .* weight;
G=dA\b;
end

function [G,OK] = estimateGradientOnOneDirection(aimingDirection, myOVEdges, changingAmount,...
    myOVEdgesMid, myOVEdgesMidIndices, crtPhi, phi)

global divergenceThreshold;
crtDivergenceThreshold=divergenceThreshold;

% find the region where gradient vector may exist
% construct A
realEdgeDivergence=dot(myOVEdges(:,3:4),repmat(aimingDirection,size(myOVEdges,1),1),2);
midEdgeDivergence=dot(myOVEdgesMid(:,3:4),repmat(aimingDirection,size(myOVEdgesMid,1),1),2);
isRealEdgeDivergenceAccepted=realEdgeDivergence>crtDivergenceThreshold;
isMidEdgeDivergenceAccepted=midEdgeDivergence>crtDivergenceThreshold;

OK = (sum(isRealEdgeDivergenceAccepted) + sum(isMidEdgeDivergenceAccepted)) >= 2;

if ~OK
    G=[0,0];
    return;
end

weights=[realEdgeDivergence(isRealEdgeDivergenceAccepted);...
    midEdgeDivergence(isMidEdgeDivergenceAccepted)];
A=[myOVEdges(isRealEdgeDivergenceAccepted,1:2);myOVEdgesMid(isMidEdgeDivergenceAccepted,1:2)]...
    .*weights;
b=[changingAmount(isRealEdgeDivergenceAccepted);...
    (phi(myOVEdgesMidIndices(isMidEdgeDivergenceAccepted,1)) +...
    phi(myOVEdgesMidIndices(isMidEdgeDivergenceAccepted,2)))./2 - crtPhi]...
    .*weights;
gradientEstimatedVia1stStep=A\b;
G=gradientEstimatedVia1stStep';
end

function [R, returnTag] = Godunov(negativeDiff,positiveDiff,S)
if S>0
    if (negativeDiff <= 0 && positiveDiff <= 0)
        returnTag=1;
        R=positiveDiff;
        return;
    end
    if (negativeDiff >= 0 && positiveDiff >= 0)
        returnTag=-1;
        R=negativeDiff;
        return;
    end
    if (negativeDiff > 0 && positiveDiff < 0)
        if (abs(negativeDiff) > abs(positiveDiff))
            returnTag=-1;
            R=negativeDiff;
            return;
        else
            returnTag=1;
            R=positiveDiff;
            return;
        end
    else
        returnTag=0;
        R=0;
        return;
    end
else
    if (negativeDiff <= 0 && positiveDiff <= 0)
        returnTag=-1;
        R=negativeDiff;
        return;
    end
    if (negativeDiff >= 0 && positiveDiff >= 0)
        returnTag=1;
        R=positiveDiff;
        return;
    end
    if (negativeDiff < 0 && positiveDiff > 0)
        if (abs(negativeDiff) > abs(positiveDiff))
            returnTag=-1;
            R=negativeDiff;
            return;
        else
            returnTag=1;
            R=positiveDiff;
            return;
        end
    else
        returnTag=0;
        R=0;
        return;
    end
end
end
