clear;
data = load('hand2.txt');


x_data = data(:,1);
y_data = data(:,2);
z_data = data(:,3);

rbf = @(e,r) sqrt((e*r).^2+1);
c = 700;
N = size(data, 1);
neval = 50;
ctrs = [x_data y_data z_data];
ctrs = [ctrs;
    [0.133894,0.280227,0.111215];
    [0.226034,0.305576,0.16492];
    [0.404979,0.31825,0.119331];
    [0.501432,-0.0239546,0.0346656];
    [0.114395,0.252269,0.112196];
    [0.318766,0.296531,0.16621];
    [0.538522,0.0752226,0.146095]];
rhs = [zeros(N,1); ones(7,1)];
DM_data = pdist2(ctrs, ctrs);
IM = sqrt((c * DM_data).^2 + 1);
coe = IM \ rhs;
bmin = min(ctrs, [], 1);
bmax = max(ctrs, [], 1);
xgrid = linspace(bmin(1)+0.02, bmax(1)+0.05, neval);
ygrid = linspace(bmin(2)-0.1, bmax(2)+0.1, neval);
zgrid = linspace(bmin(3)-0.05, bmax(3)+0.05, neval);
[xe, ye, ze] = meshgrid(xgrid, ygrid, zgrid);
epoints = [xe(:) ye(:) ze(:)];
block_size = 500; 
% Determine the number of blocks
num_blocks = ceil(size(epoints, 1) / block_size);
pf = zeros(neval^3, 1);
for block = 1:num_blocks
    start_idx = (block - 1) * block_size + 1;
    end_idx = min(block * block_size, size(epoints, 1));
    current_epoints = epoints(start_idx:end_idx, :);
    DM_eval = pdist2(current_epoints, ctrs);
    EM = rbf(c, DM_eval);
    pf_block = EM * coe;
    pf(start_idx:end_idx) = pf(start_idx:end_idx) + pf_block;
end
figure(1);
hold on;
plot3(ctrs(:,1), ctrs(:,2), ctrs(:,3), 'ro', 'Marker', 'o', 'MarkerFaceColor', 'blue');
view(45, 30);
axis off;
pf = reshape(pf, neval, neval, neval);
% Generate isosurface
f = figure;
p = patch(isosurface(xe, ye, ze, pf, 0));
% p.FaceColor = 'green';
% p.EdgeColor = 'none';
set(p,'FaceLighting','gouraud','FaceColor','g','EdgeColor','none');

% Compute normals
isonormals(xe, ye, ze, pf, p);
isocaps(xe,ye,ze,pf,0,'below')
light('Position',[1 0 1],'Style','infinite');
light('Position',[-1 0 -1],'Style','infinite');

axis equal; colormap jet; colorbar;
view(3);
axis off;