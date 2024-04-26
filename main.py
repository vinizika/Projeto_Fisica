
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
        print("Opcao [2] selecionada")
    elif n == 3:
        break
    else:
        print("Entrada invalida")
