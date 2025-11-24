//
// Seção em L
//

// tamanho da Arestas
a = 1E-2;
b = 1e-3;

// tamanho do elemento
lc = a/400;

// Pontos
Point(1) = {   0 ,   0,   0, lc};
Point(2) = {   a ,   0,   0, lc};
Point(3) = {   a ,   b,   0, lc};
Point(4) = {   b,    b,   0, lc};
Point(5) = {   b,    a,   0, lc};
Point(6) = {   0,    a,   0, lc};

// Arestas
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,1};


// Superfície
Curve Loop(1) = {1,2,3,4,5,6};
Plane Surface(1) = {1};

// Material
Physical Surface("Material,aço,1,210E9,0.3,7850.0") = {1};

// Prende todos os nós do contorno
Physical Curve("U,1,0.0") = {1:6};

// Até aqui, informações para gerar a malha

// Converte os triângulos para retângulos
Recombine Surface{:};

// Algoritmo para geração de malha
Mesh.Algorithm = 8;

// Cria a malha
Mesh 2;

// Grava a malha
Save "L.msh";