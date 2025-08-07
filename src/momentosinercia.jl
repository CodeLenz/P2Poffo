function MomentosInercia(nn,XY,ne,IJ,A,centroides)


    # Incializa os valores
    area = sum(A)
    
    # Centróide da seção 
    centroide_x = sum(centroides[:,1].*A)/area;
    centroide_y = sum(centroides[:,2].*A)/area;

    # Podemos calcular os momentos de inércia
    Ix   = 0.0 
    Iy   = 0.0
    Ixy  = 0.0
    for ele in LinearIndices(A)
        Ix  = Ix  + A[ele]*(centroide_y - centroides[ele,2])^2
        Iy  = Iy  + A[ele]*(centroide_x - centroides[ele,1])^2
        Ixy = Ixy + A[ele]*(centroide_x - centroides[ele,1])*(centroide_y - centroides[ele,2])
    end

    # Mudamos a nomenclatura para zy
    Iz  = Ix
    Izy = -Ixy
     
    # Se o produto de inércia for nulo, não precisamos calcular a rotação do sistema 
    # de referência
    if isapprox(Izy,0.0, atol=1E-12)
 
       α = 0.0
       Izl = Iz
       Iyl = Iy
 
    else
 
       # Podemos calcular o α da direção principal;
       # Evitamos divisão por zero
       if isapprox(Iz,Iy,atol=1E-12)
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
 
    
    # Retorna as propriedades da seção 
    return α, Izl , Iyl

end