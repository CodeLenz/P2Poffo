function Quadrado(a)
       Iz = a^4 / 12;
       Jeq = 0.1406*a^4
       A  = a^2
       return A, Iz, Jeq
end

function Circulo(r)
       Iz = pi*r^4 / 4;
       Jeq = pi*r^4 / 2;
       A  = pi*r^2
       return A, Iz, Jeq
end
