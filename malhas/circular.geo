//
// Seção circular 
//

// Raio da seção 
R = 1E-2;

// tamanho do elemento
lc = R/20;

// Centro, direita e acima esquerda e abaixo
Point(1) = {   R ,   R,   0, lc};
Point(2) = { 2*R ,   R,   0, lc};
Point(3) = {   R , 2*R,   0, lc};
Point(4) = {   0,    R,   0, lc};
Point(5) = {   R,    0,   0, lc};


// Circulo
Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};


// Superfície
Curve Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};

// Material
Physical Surface("Material,aço,1,210E9,0.3,7850.0") = {1};

// Prende todos os nós do contornos
Physical Curve("U,1,0.0") = {1:4};

// Até aqui, informações para gerar a malha

// Converte os triângulos para retângulos
Recombine Surface{:};

// Algoritmo para geração de malha
Mesh.Algorithm = 8;

// Cria a malha
Mesh 2;

// Grava a malha
Save "circular.msh";
