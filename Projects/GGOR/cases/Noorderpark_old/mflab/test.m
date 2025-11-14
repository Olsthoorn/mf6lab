% test van de limit van vergeijking q1 vs qt

w   = 2;
c   = logspace(-10,10,21);
b   = 50;
D1  = 10;
k1  = 10;
kD1 = k1*D1;
N   = 1;
gamma  = 1;
lambda = sqrt(kD1 .* c);
qt = 1;

G = w./c * b/D1 + b./lambda .* coth(b./lambda) - 1;


q1 = (qt * gamma - N .* c .* G)./(c.*(G+1)+gamma)

(qt * gamma - N * w * b/D1) ./ (gamma + w * b /D1)
