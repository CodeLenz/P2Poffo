#
# Programa principal
#

function Pre_processamento(arquivo, gera_pos=true)


    # Se o arquivo for um .geo, geramos um .msh utilizando a biblioteca
    # do gmsh
    if occursin(".geo",arquivo)
       
       # Gera a malha
       gmsh.initialize()
       gmsh.open(arquivo)
       gmsh.model.mesh.generate(2)
       
       # Cria o mesmo nome, mas com .msh
       mshfile = replace(arquivo,".geo"=>".msh")

       # Cria o .msh
       gmsh.write(mshfile)
      
    else 

       # Assumimos que já passaram o .msh (seria bom testar...)
       mshfile = arquivo

    end

    # Entrada de dados, lendo do arquivo .msh
    nn,XY,ne,IJ,MAT,na,AP,etypes,centroides = ConversorFEM(mshfile)

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

    println("área ",area)
    println("Iz   ",Iz)
    println("Iy   ",Iy)
    println("Izy  ",Izy)
    println("Jeq  ",Jeq)
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

    #
    # CUIDADO...esses são só as derivadas de Φ em relação a x e a y
    #
    ∇Φ =  GradΦ(Φ,ne,IJ,XY)

    # Visualização dos resultados
    if gera_pos
      
      # Nome do arquivo 
      pos_file = mshfile[1:end-4]*".pos"

      # Caminho para a pasta POS
      #caminho2 = pathof(P2Poffo)[1:end-14]*"\\Pos"

      # Retira os caminhos do nome do arquivo
      #mshfile2 = "∇Φ" * basename(mshfile)
      
      # Gera o nome do arquivo .pos 
      #posfile = replace(mshfile2,".msh"=>".pos")
      
      # Cria o arquivo completo do .pos com o nome do yaml
      #nome_pos = joinpath(caminho2,posfile)

      # Inicializa o arquivo de saída
      Lgmsh_export_init(pos_file,nn,ne,XY,etypes,IJ) 

      # Exporta o campo Φ
      Lgmsh_export_nodal_scalar(pos_file,Φ,"Φ")

      # Exporta o gradiente
      Lgmsh_export_element_scalar(pos_file,∇Φ[:,1],"∇Φx")

      # Exporta o gradiente
      Lgmsh_export_element_scalar(pos_file,∇Φ[:,2],"∇yΦ")

      # Exporta o módulo do gradiente
      Lgmsh_export_element_scalar(pos_file,sqrt.(∇Φ[:,1].^2 + ∇Φ[:,2].^2),"∇Φ")

   end

    # Retorna os valores calculados para a seção
    return (cx,cy), area, Izl, Iyl, Jeq, α, ∇Φ

end
