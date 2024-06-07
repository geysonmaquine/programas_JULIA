using Roots,DelimitedFiles

function f0(gd) # autovalores nao degenerados
    f(x)=(-sqrt(2)*sin((gd+1)*x))-(sin(gd*x))
    raizes=find_zeros(f,0,pi)
    raizes=raizes[findall(raizes.>0)]
    autovalores=zeros(gd)

    for i in 1:(gd)
        autovalores[i]=3-(2*sqrt(2)*cos(raizes[i])) 
    end
    return autovalores
end

function f4(gd) # complemento para gd-n-1 (nao esta em uso por enquanto...)
   f(x)=(sinh((gd+1)*x))-(sqrt(2)*sinh(gd*x))
    raizes=find_zeros(f,0,pi)
    println(raizes)
    #raizes=raizes[findall(raizes.>0)]
    println(length(raizes))

    autovalores=[]
    autovalores=vcat(autovalores,3-(2*sqrt(2)*cosh(raizes[length(raizes)]))) 
    return autovalores
end

function f3(gd) # autovalores degenerados
    r=[]
    um=[]

    for m in 0:gd-1
        f(x)=(sin((gd+1-m)*x))-(sqrt(2)*sin((gd-m)*x))
        raizes=find_zeros(f,0,pi)
        raizes=raizes[findall(raizes.>0)]
        autovalores=[]

        for i in 1:gd-m
            autovalores=vcat(autovalores,3-(2*sqrt(2)*cos(raizes[i]))) 
        end  
    
        if m==0
            r=vcat(r,autovalores,autovalores)
        elseif m==gd-1
            um=vcat(um,ones(3*2^(gd-2)))
        else
            r=vcat(r,repeat(autovalores,3*(2^(m-1))))
        end
    end
    r=vcat(um,r)

    return r
end

# raizes para o dual sierpinki

function DSGRSD(gs,gd)
    ags=[]
    for i in 1:gs
        if i==1
            agd=vcat(f0(gd),f3(gd))
            lambda1=([5].+sqrt.([25].-([4].*agd)))./2
            lambda2=([5].-sqrt.([25].-([4].*agd)))./2
            ags=vcat(lambda1,lambda2)
            #ags=vcat(ags,repeat([3],(2^(i-1)*(3^i+3)-3^(i-1)))) # adiciona os autovalores 3                  
        else
            lambda1=([5].+sqrt.([25].-([4].*ags)))./2
            lambda2=([5].-sqrt.([25].-([4].*ags)))./2
            ags=vcat(lambda1,lambda2)
            #ags=vcat(ags,repeat([5],(2^(gd-1)*(3^i-3)-3^(i-1)+1))) # adiciona os autovalores 5
            #ags=vcat(ags,repeat([3],(2^(gd-1)*(3^i+3)-3^(i-1)))) # adiciona os autovalores 3
        end
            ags=vcat(ags,repeat([5],(2^(gd-1)*(3^i-3)-3^(i-1)+1))) # adiciona os autovalores 5
            ags=vcat(ags,repeat([3],(2^(gd-1)*(3^i+3)-3^(i-1)))) # adiciona os autovalores 3
    end
    ags=vcat([0],ags)

    
    return ags
end

function multicamada_analitico(gs,gd,L)
    lambdaM=[]
    mono=DSGRSD(gs,gd)
    for j  in 0:L-1
        r=[(2-2*cos((j*pi)/L))].+mono
        lambdaM=vcat(lambdaM,r)
    end
    
    autovalores=sort(lambdaM) # retorna os autovalores em ordenados

    #writedlm("/content/drive/MyDrive/transport_DSG/autovalores_analiticos/autovalores_gs_$(string(gs))_gd_$(string(gd))_L_$(string(L)).csv",autovalores,";")

    return autovalores
    
end