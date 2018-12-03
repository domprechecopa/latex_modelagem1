# -*- coding: utf-8 -*-
"""
Nome = Pedro Pacheco Mendes Filho
Turma = CTEC2017
Matricula = 201701029
Professora = Marciana Lima Goes
Metodo = Metodo iterativo de Gauss Jacobi
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
                    soma = soma + float(A[i,j])*x[j]
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
        
        

    
def crit_convergencia(n):

    #Simplificando o sistema
    for i in range(A[:,0].size):
        b[i] = b[i]/A[i,i]
        A[i] = A[i]/A[i,i]
        

#######Criterio de Linha       
    max_somatoria = 0
    for i in range(A[:,0].size):
        somatoria = 0
        for j in range(n):
            somatoria += abs(A[i,j])
            
        if(somatoria > max_somatoria): max_somatoria = somatoria
    
    if(max_somatoria < 2): 
        print("Criterio de Linha foi Aprovado!")
        return True
    
#######Criterio de Diagonal dominante
    diagonais = 0
    for i in range(A[:,0].size):
        if(A[i,i]>A[:,i].sum()-A[i,i] and A[i,i]>A[i,:].sum()-A[i,i]):
            diagonais +=1
    if diagonais == n:
        print('Criterio da Diagonal Dominante Aprovado!')
        
#######Criterio de Sassenfeld
    max_sassenfeld = 0
    bj = [0]
    for i in range(A[:,0].size):
        bi = 0
        if i > 0:
            for j in range(A[:,0].size - 1):

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

A = np.array([[-4,1,0,1,0,0],
              [1,-4,1,0,1,0],
              [0,1,-4,0,0,1],
              [1,0,0,-4,1,0],
              [0,1,0,1,-4,1],
              [0,0,1,0,1,-4]], dtype=float)
    
    
chute_inicial = [0,0,0,0,0,0]
b = [-30,-20,-60,-30,-20,-60]
aray = []
iterMax = 1000
tol = 1E-3
n = A[0,:].size
if(crit_convergencia(n)):
    total_iter,sstm = jacobi(n,A,b,iterMax,tol,chute_inicial)
    
    for i in range(n):
        plt.plot(range(total_iter+1),sstm[:,i],'-o',label='X{}'.format(i+1))
    plt.xlabel("Iteracoes")
    plt.ylabel("Aproximacao")
    plt.title("Gauss Jacob")
    plt.grid(True, linestyle='-.')
    plt.legend()
    plt.savefig("gaus.png")



