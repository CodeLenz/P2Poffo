#
# Lê um arquivo .msh do gmsh e converte para a nossa entrada de dados
#
function ConversorFEM1(arquivo,esp=0.01)

    # Evita que o usuário seja tanso
    contains(arquivo,".geo") && error("Usar o .msh e não o .geo")

    # Processa o arquivo usando o lgmsh
    malha = Lgmsh.Parsemsh_FEM_Solid(arquivo)

    # Como temos mais informações do que estamos utilizando agora,
    # temos que dar uma "processada"

    # Número de nós
    nn = malha.nn

    # Número de elementos 
    ne = malha.ne

    # Coordenadas
    XY = malha.coord

    # Conectividades
    IJ = malha.connect[:,3:end]

    # Tipos de elementos (para exportar para o gmsh depois)
    etypes = malha.connect[:,1]

    # Baseado nos tipo, podemos descobrir os elementos do tipo 2
    # (triângulos) e repetir o terceiro nó
    for ele=1:ne

        # Se for um triângulo
        if etypes[ele]==2
            #println("Elemento $ele é um triângulo")
            #@show IJ[ele,:]
            IJ[ele,end] = IJ[ele,end-1]
        end 

    end


    # Informações sobre o material - Por simplicidade, vamos assumir 
    # que temos somente um material 
    MAT       = zeros(ne,2)

    # E
    MAT[:,1] .= malha.materials[1,1]

    # ν
    MAT[:,2] .= malha.materials[1,2]

    # Espessura 
    ESP = esp*ones(ne)

    # Forças concentradas
    nf = malha.nfc
    FC = malha.FC

    # Forças distribuídas
    np = malha.nft
    P  = malha.FT

    # Forças de corpo
    nfb = malha.nfb
    FB  = malha.FB

    # Apoios
    na = malha.nap
    AP = malha.AP

    # centroides
    centroides = malha.centroids

    # Retorna como estava em Dados()
    return nn,XY,ne,IJ,MAT,ESP,nf,FC,np,P,na,AP,nfb,FB,etypes,centroides

end