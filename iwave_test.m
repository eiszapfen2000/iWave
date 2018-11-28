clear all;
close all;

g_n = 10000;
g_deltaQ = 0.001;
kernel_radius = 10;
resolution = 512;

heights     = zeros(resolution);
prevHeights = zeros(resolution);
derivative  = zeros(resolution);
sources     = zeros(resolution);
obstruction = ones(resolution);
depth       = ones(resolution);

forceIn = fspecial('gauss',15,1)*5;

% sources(20:40, 80:100) = 0.25;
% sources(25:35, 85:95) = 0.5;
% sources(128, 128) = 1;

sources(253:267,253:267) = forceIn;

obstruction(200-kernel_radius/2:200+kernel_radius/2,:) = 0;
obstruction(200-kernel_radius/2:200+kernel_radius/2,200:250) = 1;

% depthr = 1:10/(resolution - 1):11;
% depth = repmat(depthr, resolution, 1);
gkernel = G(kernel_radius, g_n, g_deltaQ);

hf = figure();
% ha = axes('Parent',hf,'Units','normalized');
has = subplot(2,2,1,'Parent',hf,'Units','normalized');
hao = subplot(2,2,2,'Parent',hf,'Units','normalized');
hah = subplot(2,2,3,'Parent',hf,'Units','normalized');
hs = imshow(sources, 'DisplayRange', [0 1], 'Parent', has,'InitialMagnification','fit','Border','tight');
ho = imshow(1 - obstruction, 'DisplayRange', [0 1], 'Parent', hao,'InitialMagnification','fit','Border','tight');
hh = imshow(heights, 'DisplayRange', [], 'Parent', hah,'InitialMagnification','fit','Border','tight');
% hd = surf(-depth, 'Parent', ha, 'FaceColor', [1 0 0]);
% hold on
% hs = surf(heights, 'Parent', ha, 'FaceColor', [0 1 0]);
% shading(ha,'flat');
% alpha(hs, 0.75);
% axis(ha, 'equal');
% zlim(ha, [-10 5]);

depth_p = padarraymirror(depth, kernel_radius, kernel_radius);
depthDerivative = conv2(depth_p, gkernel, 'valid');

deltaAlpha = 0.03;
startTime = rem(now(),1);

while true

    endTime = rem(now(),1);
    deltaTime = (endTime - startTime) * 1e5;
    startTime = endTime;
    
    [heights, prevHeights] = iWave(gkernel,heights,prevHeights,sources,obstruction,depthDerivative,deltaTime,deltaAlpha);

    set(hh, 'CData', heights);
    pause(eps);
    drawnow();
    
end

