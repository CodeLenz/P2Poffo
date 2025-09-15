//
// Seção quadrada
//

// Lado da seção 
a = 1E-2;

// tamanho do elemento
lc = a/20;

// 
Point(1) = {   0 ,   0,   0, lc};
Point(2) = {   a ,   0,   0, lc};
Point(3) = {   a ,   a,   0, lc};
Point(4) = {   0,    a,   0, lc};


// Arestas
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};


// Superfície
Curve Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};

// Material
Physical Surface("") = {1};

// Prende todos os nós do contorno
Physical Curve("U,1,0.0") = {1:4};

// Até aqui, informações para gerar a malha

// Converte os triângulos para retângulos
Recombine Surface{:};

// Algoritmo para geração de malha
Mesh.Algorithm = 8;

// Cria a malha
Mesh 2;

// Grava a malha
Save "quadrado.msh";