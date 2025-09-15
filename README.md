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