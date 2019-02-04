# -*- coding: utf-8 -*-
"""
Metodo = Metodo iterativo de Gauss Jacobi
Professora = Marciana Lima Goes
"""
import numpy as np
import matplotlib.pyplot as plt

def jacobi(n, A, b, iterMax, tol, x):

    
    aray.append(x.copy())
    print(x)
    k = 0
    erro = 2
    
    while(k < iterMax and erro >= tol):
        
        erro = 0
        
        xant = x.copy()
        
        for i in range(n):
            soma = 0
            for j in range(n):
                if(j != i):
                    soma = soma + float(A[i,j])*xant[j]
            x[i] = (b[i] - soma)/float(A[i,i])
            if(abs(x[i]-xant[i])>erro):
                erro = abs(x[i] - xant[i])
        k +=1

        aray.append(x.copy())
       
        
    if(erro < tol):
        print('iteracoes = ', k)
        print('erro = ', erro)
        print("x = " , x)
        sist_array = np.array(aray, dtype=float)
        return k,sist_array

    else:
        print('Nao houve convergencia')
        return 0, 0
        
        

    
def crit_convergencia(n,m):

    #Simplificando o sistema
    for i in range(m):
        b[i] = b[i]/A[i,i]
        A[i] = A[i]/A[i,i]
        

#######Criterio de Linha       
    max_somatoria = 0
    for i in range(m):
        somatoria = 0
        for j in range(n):
            somatoria += abs(A[i,j])
            
        if(somatoria > max_somatoria): max_somatoria = somatoria
    
    if(max_somatoria < 2): 
        print("Criterio de Linha foi Aprovado!")
        return True
    
#######Criterio de Diagonal dominante
    diagonais = 0
    for i in range(m):
        if(A[i,i]>A[:,i].sum()-A[i,i] and A[i,i]>A[i,:].sum()-A[i,i]):
            diagonais +=1
    if diagonais == n:
        print('Criterio da Diagonal Dominante Aprovado!')
	return True
        
#######Criterio de Sassenfeld
    max_sassenfeld = 0
    bj = [0]
    for i in range(m):
        bi = 0
        if i > 0:
            for j in range(m - 1):

                bi += abs(A[i,j]) * bj[j]
        
        for j in range(i + 1, n):
            bi += abs(A[i,j])
            
        bj.append(bi)
            
        if(bi >  max_sassenfeld):
            max_sassenfeld = bi
            
    if(max_sassenfeld < 1):
        print('Criterio de Sassenfeld foi Aprovado!')
        return True    
    print("Nenhum criterio foi atendido!")

A = np.array([[1,2,3],
              [2,5,6],
              [2,3,4]], dtype=float)
    
    
chute_inicial = [0,0,0]
b = [23.5,50,36]
aray = []
iterMax = 1000
tol = 1E-3
n = A[0,:].size
m = A[:,0].size

### Plotagem do Grafico e Execucao do Codigo
if(crit_convergencia(n,m)):
    total_iter,sstm = jacobi(n,A,b,iterMax,tol,chute_inicial)
    
    for i in range(n):
        plt.plot(range(tot_iter+1),sstm[:,i],'-o',label='X{}'.format(i+1))
    plt.xlabel("Iteracoes")
    plt.ylabel("Aproximacao")
    plt.title("Gauss Jacob")
    plt.grid(True, linestyle='-.')
    plt.legend()
    plt.savefig("gaus.png")



