# P2Poffo
Programa Pré e Pós processamento de seções transversais utilizando MEF.

## Installation:

```bash
]add https://github.com/CodeLenz/P2Poffo
```

## Pré-Processamento
Esta etapa do programa tem como objetivo utilizar o software livre Gmsh para a geração da malha, de modo a obter as propriedades da seção transversal no sistema principal.

Retorno do Programa:
```bash
    (cx,cy) - Centroides da seção
    area    - Área da seção
    Izl     - Segundo momento de área em relação ao eixo z 
    Iyl     - Segundo momento de área em relação ao eixo y 
    Jeq     - Momento de inércia polar 
    α       - Ângulo formado entre o eixo de referência original e o eixo de referência principal
    ∇Φ      - Vetor com Gradiente da função de Airy de todos os elementos da malha
```

### Formato do arquivo de entrada
Os arquivos de entrada podem estar nos formatos .geo ou .msh. No caso de arquivos .geo, o programa utiliza a biblioteca interna do Gmsh para a geração da malha, enquanto os arquivos .msh já contêm diretamente as informações da malha.


Arquivo .geo para a leitura no GMSH e geração de malha, exemplo com seção circular

Raio da seção transversal
```bash
R = 1E-2;
```

Tamanho do elemento
```bash
lc = R/30;
```

Pontos da circuferencia com o sistema de referêcia na esquerda baixo.

Point(ID) = {coordenada em x, coordenada em y, coordenada em z, tamanho do lemento ao redor do nó}
```bash
Point(1) = {   R,    R,   0, lc};
Point(2) = { 2*R,    R,   0, lc};
Point(3) = {   R,  2*R,   0, lc};
Point(4) = {   0,    R,   0, lc};
Point(5) = {   R,    0,   0, lc};
```


Circle(ID) = {nó inicial, nó central da circunferencia, nó final};
```bash
Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};
```


Curve Loop(ID) = {curva1, curva2, curva3, curva4};
```bash
Curve Loop(1) = {1,2,3,4};
```

Plane Surface(ID) = {ID};
```bash
Plane Surface(1) = {1};
```

Material: Precisa para o conversor do Lgmsh

Physical Surface("Material,nome do paterial,ID,E,ν,ρ") = {ID};
```bash
Physical Surface("Material,aço,1,210E9,0.3,7850.0") = {1};
```

Prende todos os nós do contornor, ordem para tem que ser do conversor.

Physical Curve("U,1,0,0")= {conectivadaes(1 ao 4)}; 
```bash
Physical Curve("U,1,0.0") = {1:4};
```

Até aqui, informações para gerar a malha, agora a malha será manipulada.

Converte os triângulos para retângulos
```bash
Recombine Surface{:};
```
Algoritmo para geração de malha
```bash
Mesh.Algorithm = 8;
```
Cria a malha
```bash
Mesh 2;
```
Grava a malha com .msh
```bash
Save "circular.msh";
```

### Exemplo
Dada uma seção em L com as seguintes dimensões
```bash
a = 1 cm
b = 1 mm
```

<p align="center">
  <img src="Imagens/Seção em L.png" alt="Figura 1" width="50%">
</p>

Criando o arquivo L.geo
#### L.geo
```bash
a = 1E-2;
b = 1e-3;

lc = a/20;

Point(1) = {   0 ,   0,   0, lc};
Point(2) = {   a ,   0,   0, lc};
Point(3) = {   a ,   b,   0, lc};
Point(4) = {   b,    b,   0, lc};
Point(5) = {   b,    a,   0, lc};
Point(6) = {   0,    a,   0, lc};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,1};

Curve Loop(1) = {1,2,3,4,5,6};
Plane Surface(1) = {1};

Physical Surface("Material,aço,1,210E9,0.3,7850.0") = {1};

Physical Curve("U,1,0.0") = {1:6};

Recombine Surface{:};

Mesh.Algorithm = 8;

Mesh 2;

Save "L.msh";
```
#### Rodando o exemplo
No prompt de comando
```bash
cd Local arquivo
```
Inicie o julia
```bash
julia
```
Utilize o pacote e verifica as dependência
```bash
using P2Poffo
```
Rode o arquivo do exemplo
```bash
Pre_processamento("malhas/L.geo") 
```
#### Resultados
Para os resultados relacionados à seção transversal, observam-se diferenças quanto à dependência da discretização por elementos finitos:

1) Propriedades independentes de elementos finitos:

a) Centroide 

b) Área 

c) Segundos momento de área(Izl,Iyl)

d) Ângulo 𝛼

Essas propriedades são obtidas diretamente a partir dos primeiros momentos de área, ou seja, não dependem do tamanho ou da densidade da malha de elementos finitos.

2) Propriedades dependente de elementos finitos:

a) Função de Airy

b) Jeq da seção 

Com isso, notamos que a convergência do Jeq 
<p align="center">
  <img src="Imagens/Convergencia seção L.png" alt="Figura 1" width="80%">
</p>

A distribuição da função de airy na seção
<p align="center">
  <img src="Imagens/Phi L.png" alt="Figura 1" width="80%">
</p>

Resultados das propriedade da seção tranvesal:
```bash
    cx    = 0.002868421052631579  [m]
    cy    = 0.0028684210526315787 [m]
    area  = 1.8999999999999998e-5 [m]
    Izl   = 2.865833333333282e-10 [m] 
    Iyl   = 7.342543859648254e-11 [m]
    Jeq   = 6.178847989325534e-12 [m] 
    α     = 45.0                  [°]
```