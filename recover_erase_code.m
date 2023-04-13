% Recover set of point from erased code
%% 1. Generate random points
silent = true; % do not draw pictures
if silent
    Scale = 20000000; % generate coords in range [-Scale;Scale]
    N = 19; % number of points
    seed = 51; rng(seed); % fix RNG seed for repeatability 
    %XY_orig = randi(Scale*2+1, [N,2]) - Scale - 1; % generate Nx2 array
    XY_orig = rand([N,2]);
else
    Scale = 1000; % generate coords in range [-Scale;Scale]
    N = 5; % number of points
    seed = 42; rng(seed); % fix RNG seed for repeatability 
    XY_orig = randi(Scale*2+1, [N,2]) - Scale - 1; % generate Nx2 array
    
end
close all;

if ~silent 
    figure;
    subplot(1,3,1);
    plot(XY_orig(:,1), XY_orig(:,2), 'rx');
    hold on; 
    text(XY_orig(:,1), XY_orig(:,2),...
        arrayfun(@(x) {['a_' num2str(x)]}, 1:N));
    title('Original points');
end
%% 2. Generate erased code and find min(S/S)
[CodeE, Sorig] = erased_code(XY_orig);   
% Sorig - areas of all triangle - for correctness check only!
%TODO: check the conditions of Obshee polozhenie!
CodeE = CodeE(CodeE>=1); % select only S_1/S_2 > 1 (drop S_2/S_1)
CodeE = CodeE(:); % make column vector
%minS = min(CodeE);
%CodeEnorm = CodeE ./ minS;
%% 3. Find set of areas (two variants)
XY_recover_2way = cell(1,2);
tic
for way = 1:2
    Sset = find_set_areas(CodeE, way, N); 
    Sset = sort(Sset);
    if isempty(Sset), continue; end
    if all(near(sort(Sorig ) , Sset*min(Sorig)))
        fprintf('Recovering original S\n');
    else
        fprintf('Recovering ghost S\n');
    end
    %% 4. Find all respective triples (allign at central triangle)
    % each triple is (i1,i2,i3) - indices in Sset;
    % each triple type is 1 or 2
    try
        [Triples, TripleTypes] = find_triples(Sset, N); 
    catch
        fprintf('Recovery not succeed!\n');
        continue;
    end
    assert(size(Triples,1) == N-3 && size(Triples,2) == 3);
    TrianglesInTriplesInd = unique(Triples(:));% find triangles used in the triples
    SsetTriangleInTriples = Sset(TrianglesInTriplesInd); % 

    %% 5. Factorize in 3 class of adjacency with the central triangle
    % Determine one of 6 types
    ClassInd = factorize_adj_class(Sset, Triples);
    % return array of 1,2,3 or 0 - adjacency class, 0 is for not in triples
    %assert(numel(ClassInd)== numel(SsetTriangleInTriples));
    ClassMap = containers.Map(Sset, ClassInd);% TODO: check syntax!
    %% 6. Recover point coordinates from each triples
    % 6.1 Let the initial triangle has coordinates (0,0)-(0,1)-(1,0)
    % From the Class determine the point position with respect to axis
    [x_sign, y_sign] = getSignFromClass(Triples, TripleTypes, ClassInd);    
    % TODO: need other params???, Implement me!
    % Let the initial triangle has coordinates (0,0)-(0,1)-(1,0)
    % then S'=S/2 (as it has area=0.5)
    % Suppose S1 - prilegaet k (0,0)-(0,1) - triple of class 1.
    % then geometric set is a line such that 1/2*|y|*1 = 1/2S1 =>
    % y=S1*y_sign(i), ditto for x = S2*x_sign(i)
    XY_recover = zeros(N,2);
    XY_recover(2,:) = [1,0]; 
    XY_recover(3,:) = [0,1];
    for i = 4:N
        % recover i-th point from (i-3)th triple
        i_t = i-3;
        triple = Triples(i_t,:);
        x_s =x_sign(i_t);
        y_s = y_sign(i_t);
        assert(x_s~=0 && y_s ~=0);
        for tr = triple 
            cls_tr = ClassMap(Sset(tr));
            switch cls_tr
                case 1
                    assert(y_s~=0);
                    XY_recover(i,2) =  Sset(tr)*y_s;
                case 2
                    XY_recover(i,1) = - Sset(tr)*x_s;%% TODO: why minus here???
                case 3 % do nothing
                otherwise
                    error('Wrong class ID = %d', cls_tr);
            end
        end
    end
    
    XY_recover_2way{way} =  XY_recover;
    if ~silent
        subplot(1,3,2); 
        plot(XY_recover(:,1),   XY_recover(:,2), 'ro');
        title('Recovered points');
    end
    [CErec, Srec] = erased_code(XY_recover);
    Srec = sort(Srec); Sset = sort(Sset);
    for k = 1:10
        if find(near(Srec / Srec(1), Sset(k)))
            %fprintf('%d found\n', k);
        else
            fprintf('%d not found\n', k);
        end
    end
    
    CErec=CErec(CErec>1);
    CErec = sort(CErec);CodeE = sort(CodeE);
    assert(all(near(Srec / min(Srec) , Sset/min(Sset) )));
    fprintf('Areas match original ones!\n');
    assert(all(near(CErec, CodeE, 1.0e-6)));
    fprintf('Erased code matches the original one!\n')
    %% try to align 
    pos_orig = [4,3,1];% Map restored point to the respective original points
    % TODO: make it automatic
    e_1 = XY_orig(pos_orig(2), :) - XY_orig(pos_orig(1), :) ;
    e_2 = XY_orig(pos_orig(3), :) - XY_orig(pos_orig(1), :) ;
    A_aff  = [e_1', e_2'];
    b_aff = XY_orig(pos_orig(1));
    XY_recover_aff = (A_aff*XY_recover' +b_aff)';
    if ~silent
        subplot(1,3,3); 
        plot(XY_orig(:,1), XY_orig(:,2), 'ro');
        hold on;
        plot(XY_recover_aff(:,1), XY_recover_aff(:,2), 'bp');
        title('Aligned recovered to original');
    end
    s=0;
    
end
fprintf('%d points ', N); toc
%% 7. Display the original and the recovered sets
if 0
    figure; 
    subplot(1,3,1)
    plot(XY_orig(:,1), XY_orig(:,2), 'ro');
    title('Original');
    for w = 1:2
        subplot(1,3,w+1)
        XY_r = XY_recover_2way{w};
        if isempty(XY_r) , continue; end
        plot(XY_r(:,1), XY_r(:,2), 'ro');
        title(['Recoverd#', num2str(w)]);
    end
end
%% 8. TODO: make affine transform to align recovered with the original
% automatic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [CE, Sorig] = erased_code(XY)
%% input - set of points [ x_1 y_1; x_2 y_2;...;x_n y_n]
% output - the sorted array  of erased codes
n = size(XY,1) ; assert(n>3);
assert(size(XY,2)==2);
n_triangles = nchoosek(n,3);
triangle_ind = nchoosek(1:n,3);
S = zeros(1, n_triangles);
for i = 1:n_triangles
    triple = triangle_ind(i,:);
    S(i) = polyarea( XY(triple, 1), XY(triple, 2));
end
S_n = S ./ min(S);

Sdebug=[];
for i =1:n_triangles
    %fprintf('%d %d %d %g\n', triangle_ind(i,1), triangle_ind(i,2),...
        %triangle_ind(i,3), S_n(i));
    Sdebug = [Sdebug; S_n(i), triangle_ind(i,:)];
end
Sdebug = sortrows(Sdebug);
%disp(Sdebug);
for i = 1:size(Sdebug,1)
    %fprintf('%d) %d %d %d %g\n', i,  Sdebug(i,2:end),Sdebug(i));
end
n_code = nchoosek(n_triangles, 2)*2;
pairs =  nchoosek(1:n_triangles, 2);
CE = zeros(n_code/2,2);
for i = 1:n_code/2
    pp=pairs(i,:);
    CE(i,1) = S(pp(1))/S(pp(2));
    CE(i,2) = S(pp(2))/S(pp(1));
end
CE = sort(CE(:));
CE = unique(CE);
if numel(CE) < n_code
    error('Not general position triangles');
end
if nargout >=2, Sorig=S; end
end

%function b = near(x,y)
%b = abs(x-y) < 1.0e-5;
%end


function Sset = find_set_areas(CodeE, way, n)
% n - number of points, to check only!!
CodeE=sort(CodeE);
Smax = CodeE(end);
M=length(CodeE);
CodeE_rev = Smax ./ CodeE(M:-1:1);
% must be sorted too
assert(all(CodeE_rev(1:end-1) < CodeE_rev(2:end)));
%assert(CodeEnorm(1)==1);

% find pairs
i=1;j=1;
pairs0=[];
while i <= M && j <= M
    if near(CodeE(i), CodeE_rev(j))
        if near(CodeE(i)*CodeE(M+1-j), Smax)
            pairs0 = [pairs0; sort([i, M+1-j])];
        end
        i=i+1;j=j+1;
    else
        if CodeE(i) < CodeE_rev(j)
            i=i+1;
        else
            j=j+1;
        end
    end
end
pairs0=unique(pairs0, 'rows');

compare_with_old = n <= 7;
if compare_with_old 
    pairs=[];
    for i = 1:M-1
        for j = i+1:M
            if near(CodeE(i)*CodeE(j), Smax)
                pairs = [pairs; i j];
            end
        end
    end
    P1 = sortrows(pairs0);
    P2 = sortrows(pairs);
    assert(all(P1(:)==P2(:)));
else
    pairs = pairs0;
end
assert(size(pairs, 1) == nchoosek(n,3)-2); % - 2, because Smin,Smax excluded
% then test for Sbase - either CodeEnorm(i) or CodeEnorm(j)
Sbase=CodeE(pairs(1,3-way)); %TODO: remove "3- " after Debug
assert(Sbase>1);
Sset = [1, Sbase, Smax];
for p = 2:size(pairs,1)
    i = pairs(p,1); j = pairs(p,2);
    S_i = max(CodeE(i)/Sbase, Sbase/CodeE(i));
    S_j = max(CodeE(j)/Sbase, Sbase/CodeE(j));
    b_i = any(near(S_i, CodeE));
    b_j = any(near(S_j, CodeE));
    assert(b_i+b_j < 2);
    if b_i+b_j == 0% recover failes
        Sset=[];
        return;
    end
    if b_i, Sset=[Sset, CodeE(i)];end
    if b_j, Sset=[Sset, CodeE(j)];end
end
assert(length(Sset) == nchoosek(n,3));
end


function [Triples, TripleTypes] = find_triples(Sset, Npt)
% Npt - number of points
n = numel(Sset); % number of triangles, should be C_N^3
assert(n==nchoosek(Npt,3));
Triples =[]; TripleTypes=[];

%Sset = sort(Sset);
assert(Sset(1)==1);
for main = 2:n % main triangle
    subs = nchoosek([2:main-1, main+1:n], 2); % additional triangle
    for sub = subs'
        % type 1 Main = S1+S2+1, tyhpe 2 Main = S1+S2-1
        Smain = Sset(main); S12 = sum(Sset(sub));
        if near(Smain, S12+1) || near(Smain, S12-1)
            Triples = [Triples; main, sub'];
            if  near(Smain, S12+1)
                TripleTypes=[TripleTypes,1];
            else
                TripleTypes=[TripleTypes,2];
            end
        end
    end
end
% each triple correspond to point (exclude central triangle)
assert(size(Triples,1) == Npt-3, 'Not general position triangles!');
assert(length(TripleTypes) == Npt-3);
end


function com_side = have_common_side(i1,i2, Sset)
% checks if Sset(i1) and Sset(i2) have a common side
% firs sign always "+", others - 7 combinations
SignsCombine = [...
    1,  1,  1, -1;...
    1,  1, -1,  1;...
    1, -1,  1,  1;... % 1"-"
    1,  1, -1, -1;...
    1, -1, -1,  1;...
    1, -1,  1, -1;... % 2"-"
    1, -1, -1, -1;... % 3"-"
    ];
nTr = numel(Sset);
ind = setdiff(2:nTr, [i1,i2]);
pairs = nchoosek(ind, 2);
S = zeros(4, 1);
S(1) = Sset(i1); S(2) = Sset(i2);
for pair = pairs'
    S(3:4) = Sset(pair);
    if any(near(SignsCombine*S, 0))
        com_side = true;
        return;
    end
end
com_side = false;
end


function ClassIndExt = factorize_adj_class(Sset, Triples)
%% triangles fall in same class iff
    % 1) They have common side
    % 2) They are NOT in the same triple (maybe this is enough? - NO!)
n_Triples = size(Triples,1);
pos_max = length(Sset);%max(Triples(:));
conn = zeros(pos_max);
% set -1 for pos if two triangles in the same triple
for i = 1:n_Triples 
    for j1 = 1:3
        for j2 = 1:3
            if j1==j2,continue; end
            conn(Triples(i, j1), Triples(i, j2)) = -1;
        end
    end
end

TriplesInd  =unique(Triples(:));

for i = 1:pos_max
    if all(i~= TriplesInd), continue; end % only check for triangles in triples
    for j = 1:pos_max
        if i==j, continue; end
        if all(j~= TriplesInd), continue; end
        if conn(i,j) == -1, continue; end
        if have_common_side(i, j, Sset) % TODO: Implement me!
            conn(i,j) = 1;
        end
    end
end
conn(conn==-1)=0;
G = graph(conn(TriplesInd,TriplesInd));
ClassInd = conncomp(G);
ClassIndExt = zeros(1,pos_max);
ClassIndExt(TriplesInd) = ClassInd;
%% check that there are 3 types in each triple
for i = 1:n_Triples
    tri_bin = ClassIndExt(Triples(i,: ));
    assert( all(sort(tri_bin) == 1:3), 'Not 1,2,3 class in triple');
end

end

function [x_s, y_s] = type_class2sign(triangle_type, triangle_class)
switch triangle_class
    case 0
        x_s=0;y_s=0;
    case 1
        if triangle_type == 1
            x_s=-1;y_s=1;
        else
            x_s=1;y_s=-1;
        end
    case 2
        if triangle_type == 1
            x_s=-1;y_s=-1;
        else
            x_s=1;y_s=1;
        end
    case 3
        if triangle_type == 1
            x_s=1;y_s=-1;
        else
            x_s=-1;y_s=1;
        end
end

end

function [x_sign, y_sign] = getSignFromClassOld(Triples, TripleTypes, ClassInd)
n_triples = size(Triples,1);
assert(length(TripleTypes)==n_triples);
n_triangle = length(ClassInd);
assert(max(Triples(:))<=n_triangle);
triangle_types = zeros(1,n_triangle);
for i = 1:length(TripleTypes)
    triangle_types(Triples(i,:)) = TripleTypes(i);
end
[x_sign, y_sign] = arrayfun(@type_class2sign, triangle_types, ClassInd);

end


%%% signs are assigned to the triples! not triangles
function [x_sign, y_sign] = getSignFromClass(Triples, TripleTypes, ClassInd)
n_triples = size(Triples,1);
assert(length(TripleTypes)==n_triples);
n_triangle = length(ClassInd);
assert(max(Triples(:))<=n_triangle);
[x_sign, y_sign] = arrayfun(@type_class2sign, TripleTypes, ClassInd(Triples(:,1)));

end
