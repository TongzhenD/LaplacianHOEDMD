% This code is for generating the simulated dynamic manifold system that
% diffusion along the x_axis, and revising the diffusion factor functions
% could generate the data diffusing along the y-axis.
clc
clear
n=1000;
[u, v] = meshgrid(linspace(-1, 1, 50));
x = u;
y = v;
z = u.^2 - v.^2;

t_start = 0;
t_end = pi;
T=42;
T_lst=linspace(t_start,t_end,T);
T_lst=T_lst(2:T-1);
T=length(T_lst); % T=40

omega = 1;  % angular velocity

% define the diffusion factor functions
kx = @(t) 2*sin(t);   % along the x-axis
ky = @(t) 1;   % along the y-axis
kz = @(t) 1; % along the z-axis

% generate 1000 random points within and on the unit circle
rng('default');
theta = 2 * pi * rand(1, n); 
r = sqrt(rand(1, n));
x_points = r .* cos(theta);
y_points = r .* sin(theta);

% generate a colormap based on coordinates
c=zeros(n,3);
for i =1:n
    if x_points(i)<=y_points(i) && y_points(i)<=-x_points(i)
        c(i,:)=[1 0 0];
    elseif x_points(i)<=y_points(i) && y_points(i)>-x_points(i)
        c(i,:)=[0 1 1];
    elseif x_points(i)>y_points(i) && y_points(i)>-x_points(i)
        c(i,:)=[1 0 1];
    else 
        c(i,:)=[1 1 0];
    end
end
% display D
figure
scatter(x_points,y_points,10,c,'filled')
xlim([-1.2,1.2]);
ylim([-1.2,1.2]);
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
axis square

data=zeros(n,3,T);

videoWriter = VideoWriter('n_1000_dynamic_xdim_manifold_pi_rotated.avi');
videoWriter.FrameRate = 5;
open(videoWriter);

% draw dynamic saddle surfaces and embedded points
figure;
for i = 1:T
    
    t=T_lst(i);
    theta = omega * t;
    
    % stretch or shrink along x-axis and y-axis
    x_scale = x .* kx(t);
    y_scale = y .* ky(t);

    % make a rotation
    xt = x_scale * cos(theta) - y_scale * sin(theta);
    yt = x_scale * sin(theta) + y_scale * cos(theta);

    zt = ((x .* kx(t)).^2 - (y .* ky(t)).^2) .* kz(t);
    
    % calculate the coordinates of the point on the saddle surface
    x_points_scale=x_points .* kx(t);
    y_points_scale=y_points .* ky(t);

    x_points_t = x_points_scale * cos(theta) - y_points_scale * sin(theta);
    y_points_t = x_points_scale * sin(theta) + y_points_scale * cos(theta);

    z_points_t = ((x_points .* kx(t)).^2 - (y_points .* ky(t)).^2) .* kz(t);

    data(:,1,i)=x_points_t';
    data(:,2,i)=y_points_t';
    data(:,3,i)=z_points_t';
    
    % draw the saddle surface
    surf(xt, yt, zt);
    hold on;
    
    % draw the point cloud
    scatter3(x_points_t, y_points_t, z_points_t, 10, c, 'filled');
    
    hold off;

    xlim([-2 2]);
    ylim([-2 2]);
    zlim([-1 3]);
    title(['t = ' num2str(i)],'FontSize',14);
    xlabel('x','FontSize',14);
    ylabel('y','FontSize',14);
    zlabel('z','FontSize',14);
   
    frame = getframe(gcf);
    writeVideo(videoWriter, frame);
    pause(0.1);
end
close(videoWriter);