h = 5 %height OA
x = pi / 6 % angle


w = 2 * pi %angular velocity of P
a = -5 * pi %angular acceleration of P


c = tan(x)* h % position


u = w * (c^2+h^2)/ h % velocity


v = (((a * ((c ^ 2 + h ^ 2) ^ 2) / h) + 2 * c * (u ^ 2)) / (c ^ 2 + h ^ 2)) % acceleration


%compare
