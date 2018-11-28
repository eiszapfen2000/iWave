clear all;
close all;

deltaTime = 0.03;
deltaAlpha = 0.03;

gravity = 9.81;
gravitydtdt = 9.81 * deltaTime * deltaTime;
onealphat = 1 + deltaAlpha * deltaTime;

g_n = 10000;
g_deltaQ = 0.001;
kernel_radius = 10;
resolution = 512;
border_mode = 'symmetric';

heights     = zeros(resolution);
prevHeights = zeros(resolution);
derivative  = zeros(resolution);
sources     = zeros(resolution);
obstruction = ones(resolution);
% depth       = zeros(resolution);

sources(20:40, 80:100) = 0.25;
sources(25:35, 85:95) = 0.5;
sources(128, 128) = 1;

% obstruction(60:100, 80:100) = 0;
% obstruction = obstruction - 1;

depthr = 1:10/(resolution - 1):11;
depth = repmat(depthr, resolution, 1);
gkernel = G(kernel_radius, g_n, g_deltaQ);

% subplot(2,2,1);
% imshow(sources,[]);
% subplot(2,2,2);
% imshow(obstruction,[]);
% subplot(2,2,3);
% imshow(heights,[]);

hf = figure();
ha = axes('Parent',hf);
hi = imshow(heights, 'DisplayRange', [], 'Parent', ha);
% hd = surf(-depth, 'Parent', ha, 'FaceColor', [1 0 0]);
% hold on
% hs = surf(heights, 'Parent', ha, 'FaceColor', [0 1 0]);
% shading(ha,'flat');
% alpha(hs, 0.75);
% axis(ha, 'equal');
% zlim(ha, [-10 5]);

while true
    
    heights = heights + sources;
    heights = heights .* obstruction;
    
%     depthDerivative = conv2(depth, gkernel, 'zero');
%     shallowheights = tanh(depthDerivative) .* heights;
%     derivative = conv2(shallowheights, gkernel, 'zero');

    depth_p = padarraymirror(depth,kernel_radius, kernel_radius);
    depthDerivative = conv2(depth_p, gkernel, 'valid');
    shallowheights = tanh(depthDerivative) .* heights;
    shallowheights_p = padarraymirror(shallowheights,kernel_radius,kernel_radius);
    derivative = conv2(shallowheights_p, gkernel, 'valid');

    
    temp = heights;
    
    heights = heights .* ((2 - deltaAlpha * deltaTime) / onealphat);
    heights = heights - prevHeights .* (1 / onealphat);
    heights = heights - derivative .* (gravitydtdt / onealphat);
    
    prevHeights = temp;
    
    set(hi, 'CData', heights);
%     waitforbuttonpress;
%     set(hs, 'ZData', heights);
    pause(eps);
    drawnow();
    
end

