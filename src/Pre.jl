#
# Programa principal
#

function Pre_processamento(arquivo, gera_pos=true)


    # Se o arquivo for um .geo, geramos um .msh utilizando a biblioteca
    # do gmsh

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
    K,VFB = AplicaCCH(nn,na,AP,K,VFB)

    # Solução do sistema linear de Equações 
    Φ = K\VFB

    # Calcula as propriedaes da seção em relação a origem do 
    # sistema de referência e também o vetor A com a área de cada elemento da malha
    A, Qx, Qy, Ix, Iy, Ixy, Jeq = Propriedades_secao(Φ, ne, IJ, XY)

    # Área total da seção 
    area = sum(A)

    # Centróide da seção 
    An = A./area
    cx = sum(centroides[:,1].*An)
    cy = sum(centroides[:,2].*An)

   
    #
    # Por fim, podemos converter os valores para o centróide da seção
    #
    Ix  = Ix - Qx^2 / area
    Iy  = Iy - Qy^2 / area
    Ixy = Ixy - Qx*Qy / area
    
    # Mudamos a notação para zy
    Iz  = Ix
    Izy = -Ixy 

    println("área",area)
    println("Iz",Iz)
    println("Iy",Iy)
    println("Izy",Izy)
    println("Jeq",Jeq)
    println("")
   

    # Se o produto de inércia for nulo, não precisamos calcular a rotação do sistema 
    # de referência
    if isapprox(Izy, 0.0) #, atol=1E-16)
 
       α = 0.0
       Izl = Iz
       Iyl = Iy
 
    else
 
       # Podemos calcular o α da direção principal;
       # Evitamos divisão por zero
       if isapprox(Iz,Iy) #,atol=1E-15)
          α = sign(Izy)*45.0
       else
          α = 0.5*atand( 2*Izy/(Iz-Iy) )
       end
 
       # Com isso, temos os valores extremos dados por
       Im = (Iz + Iy) / 2
       R = sqrt( ((Iz-Iy)/2)^2 + Ixy^2 )
       Izl = Im + R
       Iyl = Im - R
    end

    # Podemos calcular as tensões tangenciais no centro de cada elemento da 
    # malha. 
    #τ =  Tensoes(Φ,ne,IJ,XY)

    # Visualização dos resultados
    if gera_pos

       # Inicializa o arquivo de saída
       Lgmsh_export_init("saida.pos",nn,ne,XY,etypes,IJ) 

       # Exporta o campo Φ
       Lgmsh_export_nodal_scalar("saida.pos",Φ,"Φ")

      # Exporta as tensões τ_zx
      #Lgmsh_export_element_scalar("saida.pos",τ[:,1],"τzx")

      # Exporta o campo τ_zy
      #Lgmsh_export_element_scalar("saida.pos",τ[:,2],"τzy")

      # Exporta o campo τ
      # Lgmsh_export_element_scalar("saida.pos",sqrt.(τ[:,1].^2 + τ[:,2].^2),"τ")

   end

    
    # Retorna os valores calculados para a seção
    return (cx,cy), area, Izl, Iyl, Jeq, α

end
