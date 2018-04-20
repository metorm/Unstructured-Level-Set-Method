%% load result data
%load resultData/U_0.025_Extra.mat
load resultData/2Circles_0.025_Explicit.mat
rExtra=Result;
%load resultData/U_0.025_Explicit.mat
rExplicit=Result;

%% config
gap=3;

%% draw
nPoints=size(p,1);

hold on;
for r=1:gap:numel(rExtra)
    if numel(rExtra{r}) == nPoints && numel(rExplicit{r}) == nPoints
        %% draw
        [~, h1]=tricontour(p,t,rExtra{r},[0 0]);
        [~, h2]=tricontour(p,t,rExplicit{r},[0 0]);
        
        drawnow;
        
        %% set color
        set(h1,'EdgeColor',[0 0 1])
        set(h2,'EdgeColor',[1 0 0])
    else
        break;
    end
end