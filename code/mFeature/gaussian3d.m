function g = gaussian3d(mu,sigma)
    gx = @(X) 1/(sigma(1)*sqrt(2*pi)) * exp(-(X-mu(1)).^2/(2*sigma(1)^2));
    gy = @(Y) 1/(sigma(2)*sqrt(2*pi)) * exp(-(Y-mu(2)).^2/(2*sigma(2)^2));
    gz = @(Z) 1/(sigma(3)*sqrt(2*pi)) * exp(-(Z-mu(3)).^2/(2*sigma(3)^2));
    g = @(X,Y,Z) gx(X).*gy(Y).*gz(Z);
end