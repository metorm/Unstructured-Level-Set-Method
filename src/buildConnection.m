function [C, NC, CMid, NCMid, NeighbourCells] = buildConnection(t)

MAX_CONNECTION_PER_VERTEX=10;
MAX_NEIGHBOR_CELLS_PER_VERTEX=8;

NPoints=max(max(t));
C=zeros(NPoints, MAX_CONNECTION_PER_VERTEX);
CMid=zeros(NPoints, MAX_NEIGHBOR_CELLS_PER_VERTEX, 2);
NeighbourCells=zeros(NPoints, MAX_NEIGHBOR_CELLS_PER_VERTEX);

for i=1:NPoints
    
    logicalInvolvedTag=sum(t==i, 2)>0;
    [r,~,~]=find(logicalInvolvedTag);
    NeighbourCells(i,1:numel(r))=r';
    
    sourceMat=t(logicalInvolvedTag,:);
    writeCounter=1;
    
    for p=1:size(sourceMat,1)
        triTags=sourceMat(p,:);
        b=ismember(triTags,i);
        for bi=1:size(b,2)
            if((~b(bi)) && (~sum(ismember(C(i,:),triTags(bi)))))
                if(writeCounter>MAX_CONNECTION_PER_VERTEX)
                    error(['writeCounter > ' num2str(MAX_CONNECTION_PER_VERTEX)]);
                else
                    C(i,writeCounter)=triTags(bi);
                    writeCounter=writeCounter+1;
                end
            end
        end
        CMid(i,p,:)=triTags(~b);
    end
end

NC=sum(C~=0,2);
NCMid=sum(sum(CMid~=0,3)/2,2);
end