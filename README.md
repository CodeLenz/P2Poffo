# P2Poffo
Programa Pré e Pós processamento de seções transversais utilizando MEF.

## Installation:

```bash
]add https://github.com/CodeLenz/P2Poffo
```

## Pré-Processamento
Esta etapa do programa tem como objetivo utilizar o software livre Gmsh para a geração da malha, de modo a obter as propriedades da seção transversal no sistema principal.

### Formato do arquivo de entrada
Os arquivos de entrada podem estar nos formatos .geo ou .msh. No caso de arquivos .geo, o programa utiliza a biblioteca interna do Gmsh para a geração da malha, enquanto os arquivos .msh já contêm diretamente as informações da malha.

'''bash
//
// Arquivo .geo para a leitura no GMSH e geração de malha, exemplo com seção circular
//

// Raio da seção transversal
R = 1E-2;

// tamanho do elemento
lc = R/30;

// Pontos da circuferencia com o sistema de referêcia na esquerda baixo
//Point(ID) = {coordenada em x, coordenada em y, coordenada em z, tamanho do lemento ao redor do nó}
Point(1) = {   R ,   R,   0, lc};
Point(2) = { 2*R ,   R,   0, lc};
Point(3) = {   R , 2*R,   0, lc};
Point(4) = {   0,    R,   0, lc};
Point(5) = {   R,    0,   0, lc};


// Circulo
//Circle(ID) = {nó inicial, nó central da circunferencia, nó final};
Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};


// Superfície
// Curve Loop(ID) = {curva1, curva2, curva3, curva4};
Curve Loop(1) = {1,2,3,4};

//Plane Surface(ID) = {ID};
Plane Surface(1) = {1};

// Material Precisa para o conversor do Lgmsh
//Physical Surface("Material,nome do paterial,ID,E,Nu,G") = {ID};
Physical Surface("Material,aço,1,210E9,0.3,7850.0") = {1};

// Prende todos os nós do contornos
// Physical Curve("")= {conectivadaes(1 ao 4)}; 
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

'''