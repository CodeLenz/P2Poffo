#
# Monta a matriz C^v_EPT
#
function MontaCEPT(E,ν)

    # Termo comum
    c = E/(1-ν^2)

    # Retorna a matriz
    Cv = [c    c*ν      0 ; 
         c*ν   c       0 ;
          0    0  c*(1-ν)/2 ]

    # Retorna a matriz
    return Cv      
end


