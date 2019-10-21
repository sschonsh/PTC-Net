function T = Rotation(v,u,n)
%%%%   rotation from v to u, preserve n, which is prependicular to v and u

%n = cross(v,u); n = n/norm(n);
v2 = cross(n,v); v2 = v2/norm(v2);
Q = [v, v2, n];
T = Q*[dot(u,v) -dot(u,v2) 0;dot(u,v2),dot(u,v),0;0 0 1]*Q';
end