using LinearAlgebra, PyPlot,LaTeXStrings, DelimitedFiles

#================================================================================================#
######### ########    ########  #######   #   #########    ##########    #########
#         #       #   #      #  #         #   #            #        #    #
#   ####  ########    ########  #####     #   #            #        #    #########
#      #  #       #   #      #  #         #   #            #        #            #
########  #        #  #      #  #         #   #########    ##########    #########
#================================================================================================#
function Grafico_espectro_autovalores()
    PyPlot.clf()
    global gs_gd_L
    global lista_autovalores
    
   for i in 1:length(gs_gd_L)
        valores,quantidade=quantidade_autovalores(lista_autovalores[i])
        gs="$(string(gs_gd_L[i][1]))"
        gd="$(string(gs_gd_L[i][2]))"
        L="$(string(gs_gd_L[i][3]))"
        PyPlot.plot(1:length(valores),valores,label="\$g_s=\$ $gs, \$g_d=\$ $gd, \$L=\$ $L")
    end
    PyPlot.xlabel("Eigenvalue order",fontsize=20)
    PyPlot.ylabel(L"$\lambda$",fontsize=20)
    PyPlot.grid(linestyle=":")
    PyPlot.legend(loc="upper center", bbox_to_anchor=(1.3, .9 ), ncol=1)
    PyPlot.savefig("autovalores_espectro.jpg",dpi=300, bbox_inches = "tight")
end

function Grafico_prob_media_retorno_quantico(t)
    PyPlot.clf()
    global gs_gd_L
    global lista_autovalores
    global lista_autovetores
    global lista_base
    for i in 1:length(gs_gd_L)
        PyPlot.clf()
        pi=prob_media_retorno_quantico(t,lista_base[i],lista_autovalores[i],lista_autovetores[i])
        #pi=limite_inferior_prob_media_retorno_quantico(t,lista_autovalores[i])
        gs="$(string(gs_gd_L[i][1]))"
        gd="$(string(gs_gd_L[i][2]))"
        L="$(string(gs_gd_L[i][3]))"
        PyPlot.loglog(t,pi,label="\$g_s=\$ $gs, \$g_d=\$ $gd, \$L=\$ $L")
        PyPlot.xlabel("t",fontsize=20)
        PyPlot.ylabel(L"$\bar{\pi}(t)$",fontsize=20)
        PyPlot.grid(linestyle=":")
        PyPlot.legend(loc="upper center", bbox_to_anchor=(1.3, .9 ), ncol=1)
        PyPlot.savefig("Prob_medio_retorno_Quantico_gs_$(gs)_gd_$(gd)_L_$(L).jpg",dpi=300, bbox_inches = "tight")
        #writedlm("/content/drive/MyDrive/transport_DSG/Probabilidade_media_retorno_Quantico/dados/Prob_medio_retorno_Quantico_gs_$(gs)_gd_$(gd)_L_$(L).csv", hcat(t,pi), ";")

    end
    #PyPlot.xlabel("t",fontsize=20)
    #PyPlot.ylabel(L"$\bar{\pi}(t)$",fontsize=20)
    #PyPlot.grid(linestyle=":")
    #PyPlot.legend(loc="upper center", bbox_to_anchor=(1.3, .9 ), ncol=1)
    #PyPlot.savefig("autovalores_espectro.jpg",dpi=300, bbox_inches = "tight")
end

function Grafico_prob_media_retorno_classico(t)
    PyPlot.clf()
    global gs_gd_L
    global lista_autovalores
    for i in 1:length(gs_gd_L)
        PyPlot.clf()
        P=prob_med_classico(t,lista_autovalores[i])
        gs="$(string(gs_gd_L[i][1]))"
        gd="$(string(gs_gd_L[i][2]))"
        L="$(string(gs_gd_L[i][3]))"
        PyPlot.loglog(t,P,label="\$g_s=\$ $gs, \$g_d=\$ $gd, \$L=\$ $L")
        PyPlot.xlabel("t",fontsize=20)
        PyPlot.ylabel(L"$\bar{p}(t)$",fontsize=20)
        PyPlot.grid(linestyle=":")
        PyPlot.legend(loc="upper center", bbox_to_anchor=(1.3, .9 ), ncol=1)
        PyPlot.savefig("/content/drive/MyDrive/transport_DSG/Probabilidade_media_de_retorno_classico/Prob_medio_retorno_Classico_gs_$(gs)_gd_$(gd)_L_$(L).jpg",dpi=300, bbox_inches = "tight")
        writedlm("/content/drive/MyDrive/transport_DSG/Probabilidade_media_de_retorno_classico/dados/Prob_medio_retorno_Clasico_gs_$(gs)_gd_$(gd)_L_$(L).csv", hcat(t,P), ";")

    end
    #PyPlot.xlabel("t",fontsize=20)
    #PyPlot.ylabel(L"$\bar{p}(t)$",fontsize=20)
    #PyPlot.grid(linestyle=":")
    #PyPlot.legend(loc="upper center", bbox_to_anchor=(1.3, .9 ), ncol=1)
    #PyPlot.savefig("autovalores_espectro.jpg",dpi=300, bbox_inches = "tight")
end

function Graficos_valor_medio()
    PyPlot.clf()
    global gs_gd_L
    global lista_autovalores
    global lista_autovetores
    global lista_base
    X=valor_med_retorno(lista_base[1],lista_autovalores[1],lista_autovetores[1])
    XA=valor_med_retorno_limite_inferior(lista_autovalores[1])
    gs="$(string(gs_gd_L[1][1]))"
    gd="$(string(gs_gd_L[1][2]))"
    L="$(string(gs_gd_L[1][3]))"
    PyPlot.bar("\$g_s=\$$gs, \n \$g_d=\$$gd, \n \$L=\$$L",X,color="darkblue",label=L"$\bar{\chi}$")
    PyPlot.text("\$g_s=\$$gs, \n \$g_d=\$$gd, \n \$L=\$$L",(X+0.0005),"$(string(round(X,digits=3)))",horizontalalignment="center")
    PyPlot.bar("\$g_s=\$$gs, \n \$g_d=\$$gd, \n \$L=\$$L",XA,color="gray",label=L"$\bar{\chi}^*$")
    PyPlot.text("\$g_s=\$$gs, \n \$g_d=\$$gd, \n \$L=\$$L",(XA-0.01),"$(string(round(XA,digits=3)))",horizontalalignment="center")
    if length(gs_gd_L)>1
        for i in 2:length(gs_gd_L)
            X=valor_med_retorno(lista_base[i],lista_autovalores[i],lista_autovetores[i])
            XA=valor_med_retorno_limite_inferior(lista_autovalores[i])
            gs="$(string(gs_gd_L[i][1]))"
            gd="$(string(gs_gd_L[i][2]))"
            L="$(string(gs_gd_L[i][3]))"
            PyPlot.bar("\$g_s=\$$gs, \n \$g_d=\$$gd, \n \$L=\$$L",X,color="darkblue")
            PyPlot.text("\$g_s=\$$gs, \n \$g_d=\$$gd, \n \$L=\$$L",(X+0.0005),"$(string(round(X,digits=3)))",horizontalalignment="center")
            PyPlot.bar("\$g_s=\$$gs, \n \$g_d=\$$gd, \n \$L=\$$L",XA,color="gray")
            PyPlot.text("\$g_s=\$$gs, \n \$g_d=\$$gd, \n \$L=\$$L",(XA-0.01),"$(string(round(XA,digits=3)))",horizontalalignment="center")
        end
    end
    PyPlot.grid(linestyle=":")
    PyPlot.xlabel("Structures",fontsize=20)
    PyPlot.ylabel(L"$\bar{\chi},\bar{\chi}^*$",fontsize=20)
    PyPlot.legend(loc="upper right")

    PyPlot.savefig("/content/drive/MyDrive/transport_DSG/Valor_medio_de_retorno/Valor_medio_retorno_quantico_fixo_L_$(L).jpg",dpi=300, bbox_inches = "tight")

end

function Grafico_probabilidade_transicao_quantica(gs,gd,L,t)
    
    global gs_gd_L
    global lista_autovalores
    global lista_autovetores
    global lista_base
    posicao = findall(x -> x == [gs, gd, L], gs_gd_L)[1]
    x = collect(1:length(lista_autovalores[posicao]))
    y = collect(1:length(lista_autovalores[posicao]))
    X,Y=meshgrid(x,y)
    for i in 1:length(t)
      aux=prob_transicao_quantica(t[i],lista_base[posicao],lista_autovalores[posicao],lista_autovetores[posicao])
      PyPlot.clf()
      p=PyPlot.pcolormesh(X,Y,aux,norm=matplotlib[:colors][:LogNorm](vmin=10^(-3),vmax=10^0),cmap="plasma")
      bar=PyPlot.colorbar(p)
      bar.set_label(L"$\pi_{j,k}$",fontsize=20)
      PyPlot.title("\$g_s=\$$(string(gs)),\$g_d=\$$(string(gd)), \$L=\$$(string(L)) \n t=$(string(t[i]))")
      PyPlot.xlabel("j",fontsize=20)
      PyPlot.ylabel("k",fontsize=20)
      PyPlot.savefig("/content/drive/MyDrive/transport_DSG/Probabilidade_transicao_quantica/Prob_transicao_quantica_gs_$(gs)_gd_$(gd)_L_$(L)_t_$(string(t[i])).jpg",dpi=300, bbox_inches = "tight")
      writedlm("/content/drive/MyDrive/transport_DSG/Probabilidade_transicao_quantica/dados/Prob_transicao_quantica_gs_$(gs)_gd_$(gd)_L_$(L)_t_$(string(t[i])).csv", aux, ";")
  end

magalu notebook
    
    function Grafico_valor_medio_geral(gs,gd,L)
    global gs_gd_L
    global lista_autovalores
    global lista_autovetores
    global lista_base
    if length(gs)==1
        x=1:length(gd)
        y=1:length(L)
        Z=zeros(length(x),length(y))
        for i in 1:length(x)
            for j in 1:length(y)
                posicao = findall(x -> x == [gs[1], gd[i], L[j]], gs_gd_L)[1]
                Z[i,j]=valor_med_retorno(lista_base[posicao],lista_autovalores[posicao],lista_autovetores[posicao])
                
            end
        end
    end
    if length(L)==1
        x=1:length(gs)
        y=1:length(gd)
        Z=zeros(length(x),length(y))
        for i in 1:length(x)
            for j in 1:length(y)
                posicao = findall(x -> x == [gs[i], gd[j], L[1]], gs_gd_L)[1]
                Z[i,j]=valor_med_retorno(lista_base[posicao],lista_autovalores[posicao],lista_autovetores[posicao])
                
            end
        end
    end
    if length(gd)==1
        x=1:length(gs)
        y=1:length(L)
        Z=zeros(length(x),length(y))
        for i in 1:length(x)
            for j in 1:length(y)
                posicao = findall(x -> x == [gs[i], gd[1], L[j]], gs_gd_L)[1]
                Z[i,j]=valor_med_retorno(lista_base[posicao],lista_autovalores[posicao],lista_autovetores[posicao])
                
            end
        end
    end
    if gs==gd && length(L)>1
        x=1:length(gs)
        y=1:length(L)
        Z=zeros(length(x),length(y))
        for i in 1:length(x)
            for j in 1:length(y)
                posicao = findall(x -> x == [gs[i], gd[i], L[j]], gs_gd_L)[1]
                Z[i,j]=valor_med_retorno(lista_base[posicao],lista_autovalores[posicao],lista_autovetores[posicao])
                
            end
        end
    end
    if gs==L && length(gd)>1
        x=1:length(gs)
        y=1:length(gd)
        Z=zeros(length(x),length(y))
        for i in 1:length(x)
            for j in 1:length(y)
                posicao = findall(x -> x == [gs[i], gd[j], L[i]], gs_gd_L)[1]
                Z[i,j]=valor_med_retorno(lista_base[posicao],lista_autovalores[posicao],lista_autovetores[posicao])
                
            end
        end
    end
    if gd==L && length(gs)>1
        x=1:length(gd)
        y=1:length(gs)
        Z=zeros(length(x),length(y))
        for i in 1:length(x)
            for j in 1:length(y)
                posicao = findall(x -> x == [gs[j], gd[i], L[i]], gs_gd_L)[1]
                Z[i,j]=valor_med_retorno(lista_base[posicao],lista_autovalores[posicao],lista_autovetores[posicao])
                
            end
        end
    end
    X,Y=meshgrid(x,y)
    PyPlot.clf()
    p=PyPlot.pcolormesh(X,Y,Z,norm=matplotlib[:colors][:LogNorm](vmin=10^(-3),vmax=10^0),cmap="plasma")
    bar=PyPlot.colorbar(p)
    bar.set_label(L"$\bar{\chi}$",fontsize=20)
    if length(gs)==1
      PyPlot.title("\$g_s=\$$(string(gs[1]))")
      PyPlot.xlabel(L"$g_d$",fontsize=20)
      PyPlot.ylabel("L",fontsize=20)
      PyPlot.savefig("/content/drive/MyDrive/transport_DSG/Valor_medio_de_retorno/valor_medio_geral_fixo_gs_$(gs[1]).jpg",dpi=300, bbox_inches = "tight")
      writedlm("/content/drive/MyDrive/transport_DSG/Valor_medio_de_retorno/dados/valor_medio_geral_fixo_gs_$(gs[1]).csv", Z, ";")
  
    end
    if length(gd)==1
      PyPlot.title("\$g_d=\$$(string(gd[1]))")
      PyPlot.xlabel(L"$g_s$",fontsize=20)
      PyPlot.ylabel("L",fontsize=20)
      PyPlot.savefig("/content/drive/MyDrive/transport_DSG/Valor_medio_de_retorno/valor_medio_geral_fixo_gd_$(gd[1]).jpg",dpi=300, bbox_inches = "tight")
      writedlm("/content/drive/MyDrive/transport_DSG/Valor_medio_de_retorno/dados/valor_medio_geral_fixo_gd_$(gd[1]).csv", Z, ";")
  
    end
    if length(L)==1
      PyPlot.title("\$L=\$$(string(L[1]))")
      PyPlot.xlabel(L"$g_s$",fontsize=20)
      PyPlot.ylabel(L"$g_d$",fontsize=20)
      PyPlot.savefig("/content/drive/MyDrive/transport_DSG/Valor_medio_de_retorno/valor_medio_geral_fixo_L_$(L[1]).jpg",dpi=300, bbox_inches = "tight")
      writedlm("/content/drive/MyDrive/transport_DSG/Valor_medio_de_retorno/dados/valor_medio_geral_fixo_L_$(L[1]).csv", Z, ";")
  
    end
    if gd==gs && length(L)>1
      PyPlot.xlabel(L"$g_s,g_d$",fontsize=20)
      PyPlot.ylabel("L",fontsize=20)
      PyPlot.savefig("/content/drive/MyDrive/transport_DSG/Valor_medio_de_retorno/valor_medio_geral_fixo_gs_gd.jpg",dpi=300, bbox_inches = "tight")
      writedlm("/content/drive/MyDrive/transport_DSG/Valor_medio_de_retorno/dados/valor_medio_geral_fixo_gs_gd.csv",Z, ";")
  
    end
    if gd==L && length(gs)>1
      PyPlot.xlabel(L"$g_d,L$",fontsize=20)
      PyPlot.ylabel(L"$g_s$",fontsize=20)
      PyPlot.savefig("/content/drive/MyDrive/transport_DSG/Valor_medio_de_retorno/valor_medio_geral_fixo_gd_L.jpg",dpi=300, bbox_inches = "tight")
      writedlm("/content/drive/MyDrive/transport_DSG/Valor_medio_de_retorno/dados/valor_medio_geral_fixo_gd_L.csv",Z, ";")
  
    end
    if gs==L && length(gd)>1
      PyPlot.xlabel(L"$g_s,L$",fontsize=20)
      PyPlot.ylabel(L"$g_d$",fontsize=20)
      PyPlot.savefig("/content/drive/MyDrive/transport_DSG/Valor_medio_de_retorno/valor_medio_geral_fixo_gs_L.jpg",dpi=300, bbox_inches = "tight")
      writedlm("/content/drive/MyDrive/transport_DSG/Valor_medio_de_retorno/dados/valor_medio_geral_fixo_gs_L.csv",Z, ";")
  
    end
   
end

function Grafico_diferenca_relativa(gs,gd,L)
    global gs_gd_L
    global lista_autovalores
    if length(gs)==1
        x=1:length(gd)
        y=1:length(L)
        Z=zeros(length(x),length(y))
        for i in 1:length(x)
            for j in 1:length(y)
                posicao = findall(x -> x == [gs[1], gd[i], L[j]], gs_gd_L)[1]
                Z[i,j]=diferenca_relativa(lista_autovalores[posicao])
                
            end
        end
    end
    if length(L)==1
        x=1:length(gs)
        y=1:length(gd)
        Z=zeros(length(x),length(y))
        for i in 1:length(x)
            for j in 1:length(y)
                posicao = findall(x -> x == [gs[i], gd[j], L[1]], gs_gd_L)[1]
                Z[i,j]=diferenca_relativa(lista_autovalores[posicao])
                
            end
        end
    end
    if length(gd)==1
        x=1:length(gs)
        y=1:length(L)
        Z=zeros(length(x),length(y))
        for i in 1:length(x)
            for j in 1:length(y)
                posicao = findall(x -> x == [gs[i], gd[1], L[j]], gs_gd_L)[1]
                Z[i,j]=diferenca_relativa(lista_autovalores[posicao])
                
            end
        end
    end
    if gs==gd && length(L)>1
        x=1:length(gs)
        y=1:length(L)
        Z=zeros(length(x),length(y))
        for i in 1:length(x)
            for j in 1:length(y)
                posicao = findall(x -> x == [gs[i], gd[i], L[j]], gs_gd_L)[1]
                Z[i,j]=diferenca_relativa(lista_autovalores[posicao])
                
            end
        end
    end
    if gs==L && length(gd)>1
        x=1:length(gs)
        y=1:length(gd)
        Z=zeros(length(x),length(y))
        for i in 1:length(x)
            for j in 1:length(y)
                posicao = findall(x -> x == [gs[i], gd[j], L[i]], gs_gd_L)[1]
                Z[i,j]=diferenca_relativa(lista_autovalores[posicao])
                
            end
        end
    end
    if gd==L && length(gs)>1
        x=1:length(gd)
        y=1:length(gs)
        Z=zeros(length(x),length(y))
        for i in 1:length(x)
            for j in 1:length(y)
                posicao = findall(x -> x == [gs[j], gd[i], L[i]], gs_gd_L)[1]
                Z[i,j]=diferenca_relativa(lista_autovalores[posicao])
                
            end
        end
    end
    X,Y=meshgrid(x,y)
    PyPlot.clf()
    p=PyPlot.pcolormesh(X,Y,Z,cmap="plasma")
    bar=PyPlot.colorbar(p)
    bar.set_label(L"$(\bar{\chi}-\bar{\chi}^*)/\bar{\chi}$",fontsize=20)
    if length(gs)==1
      PyPlot.title("\$g_s=\$$(string(gs[1]))")
      PyPlot.xlabel(L"$g_d$",fontsize=20)
      PyPlot.ylabel("L",fontsize=20)
      PyPlot.savefig("/content/drive/MyDrive/transport_DSG/Eficiencia_relativa/ef_relativa_fixo_gs_$(gs[1]).jpg",dpi=300, bbox_inches = "tight")
      writedlm("/content/drive/MyDrive/transport_DSG/Eficiencia_relativa/dados/ef_relativa_fixo_gs_$(gs[1]).csv",Z,";")
    end
    if length(gd)==1
      PyPlot.title("\$g_d=\$$(string(gd[1]))")
      PyPlot.xlabel(L"$g_s$",fontsize=20)
      PyPlot.ylabel("L",fontsize=20)
      PyPlot.savefig("/content/drive/MyDrive/transport_DSG/Eficiencia_relativa/ef_relativa_fixo_gd_$(gd[1]).jpg",dpi=300, bbox_inches = "tight")
      writedlm("/content/drive/MyDrive/transport_DSG/Eficiencia_relativa/dados/ef_relativa_fixo_gd_$(gd[1]).csv",Z,";")
  
    end
    if length(L)==1
      PyPlot.title("L=$(string(L[1]))")
      PyPlot.xlabel(L"$g_s$",fontsize=20)
      PyPlot.ylabel(L"$g_d$",fontsize=20)
      PyPlot.savefig("/content/drive/MyDrive/transport_DSG/Eficiencia_relativa/ef_relativa_fixo_L_$(L[1]).jpg",dpi=300, bbox_inches = "tight")
      writedlm("/content/drive/MyDrive/transport_DSG/Eficiencia_relativa/dados/ef_relativa_fixo_L_$(L[1]).csv",Z,";")
  
    end
    if gd==gs && length(L)>1
      PyPlot.xlabel(L"$g_s,g_d$",fontsize=20)
      PyPlot.ylabel("L",fontsize=20)
      PyPlot.savefig("/content/drive/MyDrive/transport_DSG/Eficiencia_relativa/ef_relativa_fixo_gs_gd.jpg",dpi=300, bbox_inches = "tight")
      writedlm("/content/drive/MyDrive/transport_DSG/Eficiencia_relativa/dados/ef_relativa_fixo_gs_gd.csv",Z,";")
    end
    if gd==L && length(gs)>1
      PyPlot.xlabel(L"$g_d,L$",fontsize=20)
      PyPlot.ylabel(L"$g_s$",fontsize=20)
      PyPlot.savefig("/content/drive/MyDrive/transport_DSG/Eficiencia_relativa/ef_relativa_fixo_gd_L.jpg",dpi=300, bbox_inches = "tight")
      writedlm("/content/drive/MyDrive/transport_DSG/Eficiencia_relativa/dados/ef_relativa_fixo_gd_L.csv", Z, ";")
    end
    if gs==L && length(gd)>1
      PyPlot.xlabel(L"$g_s,L$",fontsize=20)
      PyPlot.ylabel(L"$g_d$",fontsize=20)
      PyPlot.savefig("/content/drive/MyDrive/transport_DSG/Eficiencia_relativa/ef_relativa_fixo_gs_L.jpg",dpi=300, bbox_inches = "tight")
      writedlm("/content/drive/MyDrive/transport_DSG/Eficiencia_relativa/dados/ef_relativa_gs_L.csv", Z,";")
  
    end
end

