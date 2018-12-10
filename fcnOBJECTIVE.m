function [out, foil] = fcnOBJECTIVE(z, points, LE, chord)

truncation = z(end - 1);
eps = z(end);

p = z(1:end-2);
foil = fcnPARSEC2COORD(p,100).*(chord.*(1 + truncation));

u_x = 0;
u_y = 0;
u_z = 0;
theta = eps;
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];

for kk = 1:size(foil)
    foil(kk,:) = R*foil(kk,:)';
end

foil = foil + LE(1:2:3);


tmp1 = repmat(foil, 1, 1, size(points,1));
tmp2 = repmat(reshape(points', 1, 2, size(points,1)), size(foil,1),1,1);
lens = sqrt(sum((tmp2 - tmp1).^2, 2));
out = sum(min(lens,[],1));

end