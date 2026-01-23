#
# Rotinas de integração
#
#
# V(x) <= V_sup = vf * V0
# 

# Parametrização SIMP 
simp(x,p=3) = x^p

# Derivada da parametrização 
dsimp(x,p=3) = p*x^(p-1)

#
# Rotina principal
#
function Main_Otim_OC(arquivo::AbstractString, fkparam::Function, fdkparam::Function, posfile=true; verbose=false,vf = 0.5, niter=100)

    # Chama o Analise3D com o nome do arquivo, para receber a estrutura de malha
    U0, malha = Analise3D(arquivo,posfile,verbose=verbose)

    # Numero de elementos 
    ne = malha.ne

    # Estimativa inicial das densidades relativas
    x0 = vf*ones(ne)

    # Calcula o volume de cada elemento sem considerar a 
    # parametrização 
    V = Volumes(malha)

    # Volume total ta estrutura
    V0 = sum(V)

    # Volume limite
    V_sup = vf * V0
    
    # Derivada do volume é fixa 
    dV = V
    
    # Loop externo de otimização 
    for iter=1:niter

        # Calcula os deslocamentos
        U, _ = Analise3D(malha,posfile,x0=x0,kparam=[fkparam])
       
        # @show iter, sum(x0.*dV), U[8]

        # Deriva da compliance
        dC = dCompliance(malha,U,x0,fdkparam)    

        # Atualiza as densidades relativas utilizando o OC
        x0 .= OC(x0,dC,dV,V_sup,ne)
        
    end # loop externo
    
    # Adiciona uma vista com as densidades relativas ao arquivo de saída
    if posfile

        # Nome do arquivo com .pos 
        arquivo_pos =  replace(arquivo,".yaml"=>".pos")
        
        # Exporta a vista escalar 
        Lgmsh.Lgmsh_export_element_scalar(arquivo_pos,x0,"x")

    end

   # Retorna as densidades relativas 
   return x0

end






