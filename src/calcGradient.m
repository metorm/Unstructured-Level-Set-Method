function G = calcGradient(phi, ovData, A, S)

NPoints=size(phi,1);
% 2D
G=zeros(NPoints, 2);
eps=1e-3;

if size(S,2) > 1
    Sx=S(:,1);
    Sy=S(:,2);
else
    Sx=S;
    Sy=S;
end

A1=A(:,1);
A2=A(:,2);
A3=A(:,3);
A4=A(:,4);

ovData1=ovData(:,1);
ovData2=ovData(:,2);
ovData3=ovData(:,3);
ovData4=ovData(:,4);

parfor idxP=1:NPoints
    crtPhi=phi(idxP);
    GXPOK=A1(idxP).OK;
    if GXPOK
        GXP=calcGradient1Axis1Vertex(phi,crtPhi,...
            ovData1{idxP}{1,2},ovData1{idxP}{1,4},A1(idxP).A);
    end
    
    GXNOK=A2(idxP).OK;
    if GXNOK
        [GXN,]=calcGradient1Axis1Vertex(phi,crtPhi,...
            ovData2{idxP}{1,2},ovData2{idxP}{1,4},A2(idxP).A);
    end
    
    GYPOK=A3(idxP).OK;
    if GYPOK
        [GYP,]=calcGradient1Axis1Vertex(phi,crtPhi,...
            ovData3{idxP}{1,2},ovData3{idxP}{1,4},A3(idxP).A);
    end
    
    GYNOK=A4(idxP).OK;
    if GYNOK
        [GYN,]=calcGradient1Axis1Vertex(phi,crtPhi,...
            ovData4{idxP}{1,2},ovData4{idxP}{1,4},A4(idxP).A);
    end
    
    GXOK=GXPOK && GXNOK;
    GXPMissing=(~GXPOK) && GXNOK;
    GXNMissing=GXPOK && (~GXNOK);
    if (~GXPOK) && (~GXNOK)
        error('Unavailabel!')
    end
    
    GYOK=GYPOK && GYNOK;
    GYPMissing=(~GYPOK) && GYNOK;
    GYNMissing=GYPOK && (~GYNOK);
    if (~GYPOK) && (~GYNOK)
        error('Unavailabel!')
    end
    
    if GXOK
        GX=Godunov(GXN(1),GXP(1),Sx(idxP));
    else
        if GXNMissing
            GX=GXP(1);
        end
        if GXPMissing
            GX=GXN(1);
        end
    end
    
    if GYOK
        GY=Godunov(GYN(2),GYP(2),Sy(idxP));
    else
        if GYNMissing
            GY=GYP(2);
        end
        if GYPMissing
            GY=GYN(2);
        end
    end
    
    if GXPMissing || GXNMissing || GYPMissing || GYNMissing
        if (GXNMissing) && (abs(GY) < eps) && GX > 0
            GX=max(0,2-GX);
        elseif (GYNMissing) && (abs(GX) < eps) && GY > 0
            GY=max(0,2-GY);
        elseif (GXPMissing) && (abs(GY) < eps) && GX < 0
            GX=min(0,-2-GX);
        elseif (GYPMissing) && (abs(GX) < eps) && GY < 0
            GY=max(0,-2-GY);
        end
    end
    
    G(idxP,:)=[GX,GY];
end

end

function [G] = calcGradient1Axis1Vertex(phi,crtPhi,OVTags,OVMidTags,dA)
changingAmount=(phi(OVTags)-crtPhi);

phiMid=(phi(OVMidTags(:,1))+phi(OVMidTags(:,2)))./2;
changingAmountMid=(phiMid-crtPhi);

b=[changingAmount;changingAmountMid];
G=dA\b;
end

function R = Godunov(negativeDiff,positiveDiff,S)
if S>0
    if (negativeDiff <= 0 && positiveDiff <= 0)
        R=positiveDiff;
        return;
    end
    if (negativeDiff >= 0 && positiveDiff >= 0)
        R=negativeDiff;
        return;
    end
    if (negativeDiff > 0 && positiveDiff < 0)
        if (abs(negativeDiff) > abs(positiveDiff))
            R=negativeDiff;
            return;
        else
            R=positiveDiff;
            return;
        end
    else
        R=0;
        return;
    end
else
    if (negativeDiff <= 0 && positiveDiff <= 0)
        R=negativeDiff;
        return;
    end
    if (negativeDiff >= 0 && positiveDiff >= 0)
        R=positiveDiff;
        return;
    end
    if (negativeDiff < 0 && positiveDiff > 0)
        if (abs(negativeDiff) > abs(positiveDiff))
            R=negativeDiff;
            return;
        else
            R=positiveDiff;
            return;
        end
    else
        R=0;
        return;
    end
end
end
