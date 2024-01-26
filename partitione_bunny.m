% Clear workspace
clear;
data = load('Data3D_Bunny2.mat');
dsites = data.dsites;
x_data = dsites(:,1);
y_data = dsites(:,2);
z_data = dsites(:,3);
rbf = @(e,r) sqrt((e*r).^2+1);
c = 700;
N = size(dsites, 1);
neval = 50;
% Exterior Point
data1 = load('external_buny.txt'); 
x0 = data1(:,1);
y0 = data1(:,2); 
z0 = data1(:,3);
ctrs = [x_data y_data z_data];
ctrs = [ctrs;[x0,y0,z0]];
rhs = [zeros(N,1); ones(18,1)];  
DM_data = pdist2(ctrs,ctrs);
IM = sqrt((c*DM_data).^2+1);
coe = IM \ rhs; %RBF coefficient
bmin = min(ctrs,[],1);
bmax = max(ctrs,[],1);
xgrid = linspace(bmin(1)-0.001, bmax(1)+0.001, neval);
ygrid = linspace(bmin(2)-0.001, bmax(2)+0.001, neval);
zgrid = linspace(bmin(3)-0.005, bmax(3)+0.005, neval);
[xe, ye, ze] = meshgrid(xgrid, ygrid, zgrid);
epoints = [xe(:) ye(:) ze(:)];   
number = 100;
% Define block size
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
plot3(ctrs(1:N,1), ctrs(1:N,2), ctrs(1:N,3), 'bo', 'Marker', 'o', 'MarkerFaceColor', 'blue');
plot3(data1(:,1), data1(:,2), data1(:,3), 'ro', 'Marker', 'o', 'MarkerFaceColor', 'red');
view(45, 30);
axis off;
pf = reshape(pf, neval, neval, neval);
%disp(pf)
% Generate isosurface
f = figure;
p = patch(isosurface(xe, ye, ze, pf, 0));
p.FaceColor = 'yellow';
p.EdgeColor = 'none';
% Compute normals
isonormals(xe, ye, ze, pf, p);
set(p,'FaceLighting','gouraud','EdgeColor','none');
light('Position',[1 0 1],'Style','infinite');
light('Position',[-1 0 -1],'Style','infinite');
%axis equal; colormap jet; colorbar;
view(3);
axis equal;
axis off;