function [stats,areas,polygons] = voronoiBounded(P,domain)
% Bounded planar Voronoi: exact half-plane clipping in a rectangle.
% P: N x 2 centroids (m); domain=[xmin xmax ymin ymax].
% No periodic boundary conditions; normalized per-frame cell areas.
assert(size(P,2)==2&&size(P,1)>=2&&all(isfinite(P(:))));
assert(all(isfinite(domain)) && numel(domain)==4 && ...
   domain(2)>domain(1)&&domain(4)>domain(3));
assert(size(unique(P,'rows'),1)==size(P,1),'Duplicate vortex centroids.');
assert(all(P(:,1)>=domain(1)&P(:,1)<=domain(2)& ...
    P(:,2)>=domain(3)&P(:,2)<=domain(4)), 'Centroids outside bounded domain.');
n=size(P,1);areas=zeros(n,1);polygons=cell(n,1);
rect=[domain(1) domain(3);domain(2) domain(3); ...
      domain(2) domain(4);domain(1) domain(4)];
for i=1:n
    polygon=rect;
    for j=1:n
        if i==j,continue;end
        a=P(j,:)-P(i,:);
        b=(dot(P(j,:),P(j,:))-dot(P(i,:),P(i,:)))/2;
        polygon=clipHalfplane(polygon,a,b);
        if isempty(polygon),break;end
    end
    if size(polygon,1)>=3,areas(i)=polyarea(polygon(:,1),polygon(:,2));end
    polygons{i}=polygon;
end
areaDomain=(domain(2)-domain(1))*(domain(4)-domain(3));
partitionError=abs(sum(areas)-areaDomain)/areaDomain;
assert(all(areas>0) && partitionError<1e-6, ...
   'Bounded Voronoi failed area-partition check.');
a=areas/mean(areas);
stats=struct('NVortices',n,'VoronoiCV',std(a,0)/mean(a),...
    'SmallCellFraction',mean(a<0.5),'LargeCellFraction',mean(a>2),...
    'NormalizedAreaQ95',pof.percentileInclusive(a,95), ...
    'PartitionAreaError',partitionError,'NormalizedAreas',a);
end
function Q=clipHalfplane(P,a,b)
if isempty(P),Q=P;return;end
s=P*a(:)-b;
Q=zeros(2*size(P,1),2);m=0;
for k=1:size(P,1)
    j=mod(k,size(P,1))+1;si=s(k);sj=s(j);
    in1=(si<=1e-12);in2=(sj<=1e-12);
    if in1 && in2
        m=m+1;Q(m,:)=P(j,:);
    elseif in1 && ~in2
        t=si/(si-sj);m=m+1;Q(m,:)=P(k,:)+t*(P(j,:)-P(k,:));
    elseif ~in1 && in2
        t=si/(si-sj);m=m+1;Q(m,:)=P(k,:)+t*(P(j,:)-P(k,:));
        m=m+1;Q(m,:)=P(j,:);
    end
end
Q=Q(1:m,:);
if size(Q,1)>1
    dd=hypot(diff(Q(:,1)),diff(Q(:,2)));
    Q=Q([true;dd>1e-12],:);
end
end
