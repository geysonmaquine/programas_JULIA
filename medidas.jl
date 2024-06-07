using LinearAlgebra,DataFrames,DelimitedFiles

function zera_listas()
    global gs_gd_L=[]
    global lista_autovalores=[]
    global lista_autovetores=[]
    global lista_base=[]
end

function meshgrid(x, y)
    X = [i for i in x, j in 1:length(y)]
    Y = [j for i in 1:length(x), j in y]
    return X, Y
end

function delta(x,y)
   return x == y ? 1 : 0
end


function prob_med_classico(t, autovalores)
    p = zeros(length(t))
    for i in 1:length(t)
       p[i] = sum(exp.(-autovalores.*t[i]))
    end
    return p./length(autovalores)
end
function prob_transicao_classica(t,base,autovalores,autovetores)
    N=length(autovalores)
    W=zeros(N,N)
    P=zeros(N,N)
    for a in 1:N
        for i in 1:N
            W[a,i]=sum(base[:,a].*autovetores[:,i])
        end
    end
    for a in 1:N
        for b in 1:N
            P[a,b]=sum(exp.(-autovalores.*t).*W[a,:].*W[b,:])
        end
    end
    return P
end

function prob_media_retorno_quantico(t,base,autovalores,autovetores)
    pi=zeros(length(t))
    N=length(autovalores)
    W=zeros(N,N)
    for a in 1:N
        for i in 1:N
            W[a,i]=sum(base[:,a].*autovetores[:,i])^2
        end
    end

    for i in 1:length(t)
        aux=0
        for a in 1:N
            aux=aux+(abs(sum(exp.(-im.*autovalores.*t[i]).*W[a,:]))^2)
        end
        pi[i]=(1/N)*aux
    end
    return pi
end
function limite_inferior_prob_media_retorno_quantico(t,autovalores)
    PI = zeros(length(t))
    for i in 1:length(t)
       PI[i] = abs((1/length(autovalores))*sum(exp.(-im.*autovalores.*t[i])))^2
    end
    return PI
end

function prob_transicao_quantica(t,base,autovalores,autovetores)
    N=length(autovalores)
    W=zeros(N,N)
    PI=zeros(N,N)
    for a in 1:N
        for i in 1:N
            W[a,i]=sum(base[:,a].*autovetores[:,i])
        end
    end
    for a in 1:N
        for b in 1:N
            PI[a,b]=abs(sum(exp.(-im.*autovalores.*t).*W[a,:].*W[b,:]))^2
        end
    end
    return PI
end
function valor_med_transicao(base,autovalores,autovetores)
    N=length(autovalores)
    X=zeros(N,N)
    for a in 1:N
        aux=0
        for b in 1:N
            for m in 1:N
                for n in 1:N
                    aux=aux+(delta(autovalores[n],autovalores[m])*sum(base[:,b].*autovetores[:,n])*sum(base[:,a].*autovetores[:,n])*sum(base[:,a].*autovetores[:,m])*sum(base[:,b].*autovetores[:,m]))
                end    
            end
        end
        X[a,b]=aux
    end
    return X
end

function valor_med_retorno(base,autovalores,autovetores)
    N=length(autovalores)
    X=zeros(N)
    for a in 1:N
        aux=0
        for n in 1:N
            for m in 1:N
                aux=aux+(delta(autovalores[n],autovalores[m])*(sum(abs.(base[:,a].*autovetores[:,n]))^2)*(sum(abs.(base[:,a].*autovetores[:,m]))^2))
            end
        end
        X[a]=aux
    end
    return sum(X)/N
end

function valor_med_retorno_limite_inferior(autovalores) # vesao original
   N=length(autovalores)
    X=0
    for n in 1:N
        for m in 1:N
            X=X+(delta(autovalores[n],autovalores[m]))
       end
    end
    return X/(N^2)
end

function prob_med_retorno_limite_inferior(t,autovalores)
    PI=zeros(length(t))
    N=length(autovalores)
    for i in 1:length(t)
        PI[i]=abs((1/N)*sum(exp.(-im.*autovalores.*t[i])))^2
    end
    return PI
end

function quantidade_autovalores(autovalores)
    valores=unique(autovalores)
    quantidade=zeros(length(valores))
    for i in 1:length(valores)
        quantidade[i]=count(autovalores.==valores[i])
    end
    return valores,quantidade
end

#function valor_med_retorno_limite_inferior(autovalores) # versao otimizada
 #   N=length(autovalores)
  #  valores,quantidade=quantidade_autovalores(autovalores)
   # X=0
    #for i in 1:length(quantidade)
     # X=X+((quantidade[i]^2))
    #end
    #return X/(N^2)
# #end
function diferenca_relativa(autovalores) 
    N=length(autovalores)
    valores,quantidade=quantidade_autovalores(autovalores)
    max=findall(x -> x == maximum(quantidade), quantidade)
    w=0
  
    for i in 1:length(max)
      w=w+((quantidade[max[i]]^2))
    end
    aux=(1/(N))*w
    X=valor_med_retorno_limite_inferior(autovalores)
    X1=(aux^2)+((1/N)*(1-aux))
            
    return((X-X1)/X)
end

function listas(gs,gd,L) # gera a lista de todos os autovalores e autovetores variando a geracao
    zera_listas() # zera todas as listas de resultados (use caso sua ide salve as variavies)
    global gs_gd_L
    global lista_autovalores
    global lista_autovetores
    global lista_base
    for k in 1:length(L)
        for j in 1:length(gd)
            for i in 1:length(gs)
                push!(gs_gd_L,[gs[i],gd[j],L[k]])
                autovalores, autovetores=eigen(LAPLACIANA_MULTICAMADA(gs[i],gd[j],L[k]))
                autovalores=round.(autovalores,digits=4)
                autovetores=round.(autovetores,digits=4)
                base=Matrix{Int64}(I,length(autovalores),length(autovalores))
                push!(lista_autovalores,autovalores)
                push!(lista_autovetores,autovetores)
                push!(lista_base,base)
                #writedlm("/content/drive/MyDrive/transport_DSG/autovalores/autovalores_gs_$(string(gs[i]))_gd_$(string(gd[j]))_L_$(string(L[k])).csv",autovalores,";")
                #writedlm("/content/drive/MyDrive/transport_DSG/autovetores/autovetores_gs_$(string(gs[i]))_gd_$(string(gd[j]))_L_$(string(L[k])).csv",autovetores,";")
            end
        end
    end
  
end

