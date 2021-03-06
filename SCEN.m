function AE = SCEN(data,boxsize)

%The function performs the transformation from gene expression matrix to association-entropy matrix (AE)
%data: Gene expression matrix (normalized count), rows = genes, columns = cells
%boxsize: Size of neighborhood, Default = 0.1 (nx(k) = ny(k) = 0.1*n)

if nargin < 2 || isempty(boxsize)
    boxsize = 0.1;
end

%Define the neighborhood of each plot
[m,n] = size(data);
upper = zeros(m,n);
lower = zeros(m,n);
for i = 1 : m
    [s1,s2] = sort(data(i,:));
    n0 = n - sum(sign(s1));
    h = round(boxsize*(n - n0)/2);
    k = 1;
    while k <= n
        s = sum(s1 == s1(k)) - 1;
        if s >= h
            upper(i,s2(k : k + s)) = s1(k);
            lower(i,s2(k : k + s)) = s1(k);
        else
            upper(i,s2(k : k + s)) = s1(min(n, k + s + h));
            lower(i,s2(k : k + s)) = s1(max(n0*(n0 > h) + 1, k - h));
        end
        k = k + s + 1;
    end
end

%Construct single-cell entropy network
ent = zeros(sum(sum(sign(data))),2);
mi = zeros(m,m);
f = [1; zeros(n,1)];
for k = 1 : n
    tic
    b = bsxfun(@and,(bsxfun(@le,data,upper(:,k)) & ...
        bsxfun(@ge,data,lower(:,k))),data(:,k) > 0);
    a = sum(b,2);
    b = sparse(double(b));
    b = b*b';
    scen = b/n.*log(n*b./(a*a' + eps) + eps);
    mi = mi + scen;
    scen(scen < 0) = 0;
    e = full(sum(scen,2) - diag(scen));
    f(k + 1) = f(k) + sum(e > 0);
    ent(f(k) : f(k + 1) - 1, :) = [find(e), e(e > 0)];
    disp(['Cell ' num2str(k) ' is completed']);
    toc
end
clear data upper lower a b csn

%Construct association-entropy matrix
AE = zeros(m,n);
k = 0;
for i = 1 : length(ent)
    k = k + (i == f(k + 1));
    AE(ent(i,1),k) = ent(i,2);
end
AE = bsxfun(@rdivide,AE,sum(AE))*sum(sum(mi,2) - diag(mi))/n;