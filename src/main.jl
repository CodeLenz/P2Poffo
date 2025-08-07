#
# Programa principal
#

function AnaliseTorcao(arquivo)

    # Entrada de dados
    nn,XY,ne,IJ,MAT,ESP,nf,FC,np,P,na,AP,nfb,FB,etypes = ConversorFEM1(arquivo)

    #
    #                      Processamento
    #

    # Matriz de rigidez global
    K = RigidezGlobal(nn,ne,XY,IJ)

    # Vetor de forças de corpo
    VFB = ForcabGlobal(nn,ne,XY,IJ)

    # Soma as contribuições das forças 
    F = VFB

    # Aplica as condições de contorno essenciais homogêneas
    K,F = AplicaCCH(nn,na,AP,K,F)

    # Solução do sistema linear de Equações 
    Φ = K\F

    # Calcula o J_eq para a seção transversal
    Jeq, area = Jequivalente(Φ, ne, IJ,XY)

    println("Momento de inércia polar equivalente para a seção:  ",Jeq)
    println("Área da seção ", area)
        
    # Podemos calcular as tensões tangenciais no centro de cada elemento da 
    # malha. 
    τ =  Tensoes(Φ,ne,IJ,XY)

    # Visualização dos resultados
    Lgmsh_export_init("saida.pos",nn,ne,XY,etypes,IJ) 

    # Exporta o campo Φ
    Lgmsh_export_nodal_scalar("saida.pos",Φ,"Φ")

    # Exporta as tensões τ_zx
    Lgmsh_export_element_scalar("saida.pos",τ[:,1],"τzx")

    # Exporta o campo τ_zy
    Lgmsh_export_element_scalar("saida.pos",τ[:,2],"τzy")

    # Exporta o campo τ
    Lgmsh_export_element_scalar("saida.pos",sqrt.(τ[:,1].^2 + τ[:,2].^2),"τ")
    
    # Retorna os valores calculados
    return Jeq, area

end
