using LinearAlgebra

function gerar_sequencia(n)
    a=zeros(n)
    a[1]=1
    for i in 2:n
        a[i]=a[i-1]+(3^(i-1))
    end
    return a
end

function MATRIZ_ADJACENCIA_gs(gs)
    N=3^gs
    A=zeros(N,N)
for i in 3:3:N
    A[i,i-2]=1
    A[i-2,i]=1
    A[i,i-1]=1
    A[i-1,i]=1
    A[i-1,i-2]=1
    A[i-2,i-1]=1        
end

if(gs>1)
   seq=gerar_sequencia(gs-1)
   for i in 1:(gs-1)
       aux2=collect(3^(i):3^(i):N)
        for j in 1:3:length(aux2)
            A[aux2[j]+1,aux2[j]-Int(seq[i])]=1
            A[aux2[j]-Int(seq[i]),aux2[j]+1]=1 
            A[aux2[j],aux2[j+1]+1]=1
            A[aux2[j+1]+1,aux2[j]]=1
            A[aux2[j+2]-Int(seq[i]),aux2[j+1]]=1
            A[aux2[j+1],aux2[j+2]-Int(seq[i])]=1
        end
    end

end
    return A
end

function MATRIZ_ADJACENCIA_gd(gs,gd)
Nd = 3^gs * ((3 * 2^gd) - 2)  # numero de nos dendrimero
A=MATRIZ_ADJACENCIA_gs(gs)
N=length(A[1,:])
Ad=kron(Matrix{Float64}(I,Int(Nd/N),Int(Nd/N)),A)


result_d = findall(sum.(eachrow(Ad)) .== 2)  # procura os indices com valores 1 na matriz adjacencia
laterais_grau_dois = transpose(reshape(result_d,3,Int(length(result_d)/3)))

Ad[laterais_grau_dois[1, 1], laterais_grau_dois[2, 1]] = 1
Ad[laterais_grau_dois[2, 1], laterais_grau_dois[1, 1]] = 1
Ad[laterais_grau_dois[1, 2], laterais_grau_dois[3, 1]] = 1
Ad[laterais_grau_dois[3, 1], laterais_grau_dois[1, 2]] = 1
Ad[laterais_grau_dois[1, 3], laterais_grau_dois[4, 1]] = 1
Ad[laterais_grau_dois[4, 1], laterais_grau_dois[1, 3]] = 1



    if gd>1
        cont10 = 2
        for i in 5:length(laterais_grau_dois[:,1])
            if i % 2 == 0
                Ad[laterais_grau_dois[i, 1], laterais_grau_dois[cont10, 2]] = 1
                Ad[laterais_grau_dois[cont10, 2], laterais_grau_dois[i, 1]] = 1
                cont10 += 1        
            else
                Ad[laterais_grau_dois[i, 1], laterais_grau_dois[cont10, 3]] = 1
                Ad[laterais_grau_dois[cont10, 3], laterais_grau_dois[i, 1]] = 1
        
            end
        end
    end
return Ad
end

function MULTICAMADAS_FRACTAL(gs,gd,L)
    Ad=MATRIZ_ADJACENCIA_gd(gs,gd)
    Nd=length(Ad[1,:])
    NL=Nd*L
    if L>1
        Ad_multi=kron(Matrix{Float64}(I,L,L),Ad)
       indices=findall(Ad.==1)
       cont1=0
       cont2=Nd
       for i in 2:L
            for j in 1:length(indices)
                Ad_multi[cont1+Int(indices[j][1]),cont2+Int(indices[j][1])]=1
                Ad_multi[cont1+Int(indices[j][2]),cont2+Int(indices[j][2])]=1
                Ad_multi[cont2+Int(indices[j][1]),cont1+Int(indices[j][1])]=1
                Ad_multi[cont2+Int(indices[j][2]),cont1+Int(indices[j][2])]=1
            end
            cont1=cont1+Nd
            cont2=cont2+Nd
       end
      
    else
        Ad_multi=Ad
    end
    return Ad_multi
end

function LAPLACIANA_MULTICAMADA(gs, gd, L)
    Ad_multi = MULTICAMADAS_FRACTAL(gs, gd, L)
    Ld_multi = -Ad_multi  # Laplaciana do dual sierpinski forma dendrimero multicamadas
    Nd_multi = size(Ld_multi, 2)
    sum_Ad_multi = abs.(sum(Ad_multi, dims=1))

    Ld_multi[diagind(Ld_multi)] = sum_Ad_multi

    return Ld_multi
end