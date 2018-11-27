function o = padarraymirror(x,r,c)
%%
o = padarray(x,[r+1 c+1],'both','symmetric');
o(end-r,:) = [];
o(:,end-c) = [];
o(r+1,:) = [];
o(:,c+1) = [];
end