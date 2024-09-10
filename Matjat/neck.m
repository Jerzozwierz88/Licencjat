function [res] = neck(p, r, n0, n1, u)
if (n0 == n1) && (n0 == 1)
    res = (-2).*p.*((-1)+r).*(p+(-1).*p.*r.^2+2.*r.*cosh(p)+(1+r.^2).*sinh( ...
  p)).^(-1);
else
    res = (-2).*n0.*((-1)+n1).^(1/2).*n1.*((-1)+r).*((-1).*((-1)+n0).^(1/2) ...
  .*r.*sin((((-1)+n1).*n1.^(-1)).^(1/2).*p.*u).*sinh(n0.^(-1/2).*p)+ ...
  sin((((-1)+n0).*n0.^(-1)).^(1/2).*p).*(((-1)+r).*sin((((-1)+n1).* ...
  n1.^(-1)).^(1/2).*p.*u)+((-1)+n1).^(1/2).*sinh(n1.^(-1/2).*p.*u))) ...
  .*(2.*((((-1)+n0).*n0.*n1).^(1/2)+(-1).*(((-1)+n0).*n0.*n1.^3).^( ...
  1/2)).*r.*cos((((-1)+n0).*n0.^(-1)).^(1/2).*p).*cos((((-1)+n1).* ...
  n1.^(-1)).^(1/2).*p.*u)+(-2).*((((-1)+n0).*n0.*n1).^(1/2)+(-1).*(( ...
  (-1)+n0).*n0.*n1.^3).^(1/2)).*r.*cosh(n0.^(-1/2).*p).*cosh(n1.^( ...
  -1/2).*p.*u)+((-1)+n1).^(1/2).*((-1).*n0.*r.^2+n1.*((-1)+2.*n0.*r) ...
  ).*sin((((-1)+n0).*n0.^(-1)).^(1/2).*p).*sin((((-1)+n1).*n1.^(-1)) ...
  .^(1/2).*p.*u)+((-1)+n1).*(n1+(-1).*n0.*r.^2).*sin((((-1)+n0).* ...
  n0.^(-1)).^(1/2).*p).*sinh(n1.^(-1/2).*p.*u)+sinh(n0.^(-1/2).*p).* ...
  ((((-1)+n0).*((-1)+n1)).^(1/2).*((-1).*n1+n0.*r.^2).*sin((((-1)+ ...
  n1).*n1.^(-1)).^(1/2).*p.*u)+((-1)+n0).^(1/2).*((-1)+n1).*(n1+n0.* ...
  r.^2).*sinh(n1.^(-1/2).*p.*u))).^(-1);
end
end