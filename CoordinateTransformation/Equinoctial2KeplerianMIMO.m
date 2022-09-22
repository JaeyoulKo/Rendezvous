function [a,e,i,o,w,nu] = Equinoctial2KeplerianMIMO(p,f,g,h,k,L)
% every element = column or row vector

a = p./(1-f.^2-g.^2);
e = sqrt(f.^2+g.^2);
i = 2*atan(sqrt(h.^2+k.^2));
o = atan2(k,h);
w = atan2(g,f)-o; %?
nu = L-o-w;
% if o<0
%     o=o+2*pi;
% end
% if w<0
%     w=w+2*pi;
% end
% if nu<0
%     nu=nu+2*pi;
% end

end
