#
# Programa principal
#

function AnaliseTorcao(arquivo)

    # Entrada de dados
    nn,XY,ne,IJ,MAT,ESP,nf,FC,np,P,na,AP,nfb,FB,etypes,centroides = ConversorFEM1(arquivo)

    #
    #                      Processamento
    #

    # Matriz de rigidez global
    K = RigidezGlobal(nn,ne,XY,IJ)

    # Vetor de forças de corpo
    VFB = ForcabGlobal(nn,ne,XY,IJ)

    # Aplica as condições de contorno essenciais homogêneas
    K,F = AplicaCCH(nn,na,AP,K,VFB)

    # Solução do sistema linear de Equações 
    Φ = K\FFB

    # Calcula o J_eq para a seção transversal e também 
    # devolve um vetor com a área de cada elemento da malha
    Jeq, A = Jequivalente(Φ, ne, IJ,XY)

    # Área total da seção 
    area = sum(A)
    
    # Agora podemos calcular as outras propriedades da seção transversal
    # como posição do centróide, Ix, Iy, Ixy e α
    α, Izl , Iyl = MomentosInercia(A,centroides)

    # Podemos calcular as tensões tangenciais no centro de cada elemento da 
    # malha. 
    #τ =  Tensoes(Φ,ne,IJ,XY)

    # Visualização dos resultados
    Lgmsh_export_init("saida.pos",nn,ne,XY,etypes,IJ) 

    # Exporta o campo Φ
    Lgmsh_export_nodal_scalar("saida.pos",Φ,"Φ")

    # Exporta as tensões τ_zx
    #Lgmsh_export_element_scalar("saida.pos",τ[:,1],"τzx")

    # Exporta o campo τ_zy
    #Lgmsh_export_element_scalar("saida.pos",τ[:,2],"τzy")

    # Exporta o campo τ
    # Lgmsh_export_element_scalar("saida.pos",sqrt.(τ[:,1].^2 + τ[:,2].^2),"τ")
    
    # Retorna os valores calculados para a seção
    return Jeq, area, Izl, Iyl, α

end
