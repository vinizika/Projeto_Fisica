import numpy as np
from math import * 

def menu():
    print("Simulador Crocodilo de Perdizes\n")
    print("[1] Determinação da função de onda quântica e outros parâmetros\n")
    print("[2] Cálculo dos parâmetros da caixa e partícula, dada a função de onda\n")
    print("[3] Encerrar Programa\n")
    n = int(input("Escolha a sua opcao: "))
    return n

while(True):
    n = menu()
    if n == 1:
        print("Opcao [1] selecionada")
    elif n == 2:
        A = float(input("Digite o valor da amplitude (A): "))
        k = float(input("Digite o valor do número de onda (k): "))
        Xp = float(input("Digite a posição específica (Xp) onde deseja calcular a probabilidade: "))

        L = 2/A**2
        n = round(k * L / np.pi)  

        probabilidade_xp = L*Xp
        probabilidade = A**2 * np.sin(k * probabilidade_xp)**2

        print(f"Largura da caixa (L): {L} metros")
        print(f"Nível quântico da partícula (n): {n}")
        print(f"Probabilidade de encontrar a partícula na posição {probabilidade_xp} é: {probabilidade:.5f}")
    elif n == 3:
        break
    else:
        print("Entrada invalida")
