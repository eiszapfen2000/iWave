function [heights, prev_heights] = iWave(gkernel, heights, prev_heights, sources, obstructions, depth, dt, da)
%%
gravity = 9.81;
gravitydtdt = gravity * dt * dt;
onealphat = 1 + da * dt;

heights = heights + sources;
heights = heights .* obstructions;

borderSize = floor(size(gkernel,1)/2);
depth_p = padarraymirror(depth,borderSize, borderSize);
depthDerivative = conv2(depth_p, gkernel, 'valid');
shallowheights = tanh(depthDerivative) .* heights;
shallowheights_p = padarraymirror(shallowheights,borderSize,borderSize);
derivative = conv2(shallowheights_p, gkernel, 'valid');

temp = heights;

heights = heights .* ((2 - da * dt) / onealphat);
heights = heights - prev_heights .* (1 / onealphat);
heights = heights - derivative .* (gravitydtdt / onealphat);

prev_heights = temp;
end