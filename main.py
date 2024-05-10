from math import * 
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.patches as patches
import numpy as np
from scipy.integrate import quad

def menu():
    print("Simulador Física\n")
    print("[1] Determinação da função de onda quântica e outros parâmetros\n")
    print("[2] Cálculo dos parâmetros da caixa e partícula, dada a função de onda\n")
    print("[3] Simulacao\n")
    print("[4] Encerrar Programa\n")
    n = int(input("Escolha a sua opcao: "))
    return n

def menu_posicao():
    print("Você deseja entrar com a posicao em:\n")
    print("[1] Probabilidade (0 a 1)\n")
    print("[2] Metros\n")
    n = int(input())
    return n

def psi(n, L):
    amplitude = sqrt(2 / L)
    k = n*pi / L
    return amplitude, k

def energy_n(n, L):
    hbar = 1.0545718e-34  
    print("Voce deseja calcular a energia de que particula?\n")
    print("[1] Proton")
    print("[2] Eletron")
    i = int(input())
    if i == 1:
        m = 1.6726219e-27  
    else:
        m = 9.11e-31 
    E_joules = ((n ** 2) * (pi ** 2) * (hbar ** 2)) / (2 * m * (L ** 2))
    E_eV = E_joules / 1.60218e-19 
    return E_joules, E_eV

def menu_energia():
    print("Voce deseja calcular a energia de que particula?\n")
    print("[1] Proton")
    print("[2] Eletron")
    n = int(input())
    return n

def wave_function(x, n, L):
    return np.sqrt(2 / L) * np.sin(n * np.pi * x / L)

def probability(a, b, n, L):
    integrand = lambda x: wave_function(x, n, L)**2
    result, _ = quad(integrand, a, b)
    return result

def velocidade_electron(n):
    hbar = 1.0545718e-34  
    print(f"Voce deseja calcular a velocidade no nivel {n} de que particula?\n")
    print("[1] Proton")
    print("[2] Eletron")
    i = int(input())
    if i == 1:
        m = 1.6726219e-27  
    else:
        m = 9.11e-31 
    E_joules = ((n ** 2) * (pi ** 2) * (hbar ** 2)) / (2 * m * (L ** 2))
    E_eV = E_joules / 1.60218e-19  
    p = sqrt(2 * m * E_joules)
    v_n = p / m
    return v_n

def grafico_funcao_de_onda(L, n_inicial, n_final):
    num_x = np.linspace(0, L, 1000)
    amplitude_inicial, k_inicial = psi(n_inicial, L)
    amplitude_final, k_final = psi(n_final, L)
    psi_inicial = amplitude_inicial * np.sin(k_inicial * num_x)
    psi_final = amplitude_final * np.sin(k_final * num_x)

    plt.figure(figsize=(10, 6))
    plt.plot(num_x, psi_inicial, label='Nível Inicial')
    plt.plot(num_x, psi_final, label='Nível Final')
    plt.xlabel('Posição x (m)')
    plt.ylabel('ψ(x)')
    plt.title('Gráfico da função de onda inicial e final')
    plt.legend()
    plt.grid(True)
    plt.show()

def grafico_distribuicao_probabilidades(L, n_inicial, n_final):
    num_x = np.linspace(0, L, 1000)
    amplitude_inicial, k_inicial = psi(n_inicial, L)
    amplitude_final, k_final = psi(n_final, L)
    
    psi_inicial = amplitude_inicial * np.sin(k_inicial * num_x)
    psi_final = amplitude_final * np.sin(k_final * num_x)
    
    probabilidade_inicial_vals = psi_inicial**2
    probabilidade_final_vals = psi_final**2
    
    plt.figure(figsize=(10, 6))
    plt.plot(num_x, probabilidade_inicial_vals, label='Distribuição Inicial')
    plt.plot(num_x, probabilidade_final_vals, label='Distribuição Final')
    plt.xlabel('Posição x (m)')
    plt.ylabel('|ψ(x)|²')
    plt.title('Gráfico da Distribuição de Probabilidade')
    plt.legend()
    plt.grid(True)
    plt.ylim(0, np.max([probabilidade_inicial_vals, probabilidade_final_vals]) * 1.1)  # Definindo o limite superior para melhor visualização
    plt.show()

largura = 1.0

while(True):
    print(' ')
    print('Gabriel Koiama - RA: 24.123.051-5')
    print('João Pedro Peterutto - RA: 24.123.045-7')
    print('João Pedro Lopes - RA: 24.123.071-3')
    print('Vinicius Duarte - RA: 24.123.073-9\n')
    print('O projeto "Simulador Física" é uma ferramenta educativa em Python desenvolvida para explorar conceitos de mecânica quântica através da simulação de uma partícula em uma caixa. O software utiliza a biblioteca Matplotlib para visualizações gráficas interativas, permitindo aos usuários visualizar funções de onda, distribuições de probabilidade e dinâmicas de transição de estado quântico. O código inclui um menu interativo para facilitar a navegação entre diferentes cálculos e visualizações, como determinação de função de onda, cálculo de energia, e simulações de movimento de partículas em níveis quânticos. Essencialmente, o design foca na usabilidade e na integração de elementos interativos para promover um aprendizado mais intuitivo e visual dos princípios da física quântica.\n')
    print('É importante que você veja atentamente o menu e durante todo o código será perguntado se você deseja calcular as informações de um Próton (1) ou Elétron (2), para que não haja nenhum erro de interpretação do código')
    n = menu()
    if n == 1:
        L = float(input("Digite a largura da caixa (L): "))
        largura = L
        while(True):
            n_inicial = int(input("Digite o nível quântico inicial (n inicial): "))
            n_final = int(input("Digite o nível quântico final (n final): "))
            if n_final < 0 or n_inicial < 0:
                print("O valor de n precisa ser maior que 0")
            else:
                break
        if n_final > n_inicial:
            print("Foton absorvido para o salto quantico")
        else:
            print("Foton emitido para o salto quantico")
        while True:
            a = float(input("Digite o limite inferior (a): "))
            b = float(input("Digite o limite superior (b): "))
            if 0 <= a <= b <= L:
                break
            else:
                print(f"Por favor, garanta que 0 <= a <= b <= L ({L}). Tente novamente.")
        amplitude_inicial, k_inicial = psi(n_inicial, L)
        amplitude_final, k_final =psi(n_final,L)
        print(f"Função de onda para n inicial: 𝜓(x) = {amplitude_inicial:.3e} sin({k_inicial:.3e} x)\n")
        print(f"Função de onda para n final: 𝜓(x) = {amplitude_final:.3e} sin({k_final:.3e} x)\n")

        E_inicial_joules, E_inicial_eV = energy_n(n_inicial, L)
        E_final_joules, E_final_eV = energy_n(n_final, L)

        E_inicial_joules, E_inicial_eV = energy_n(n_inicial, L)
        E_final_joules, E_final_eV = energy_n(n_final, L)
        
        delta_E_joules = abs(E_final_joules - E_inicial_joules)
        delta_E_eV = abs(E_final_eV - E_inicial_eV)
        
        
        h = 6.62607015e-34  
        foton_frequencia = delta_E_joules / h

      
        c = 3e8  
        foton_lambda = c / foton_frequencia if foton_frequencia != 0 else 0
        
        print(f"Energia do fóton (𝐸_𝑓ó𝑡𝑜𝑛): {delta_E_eV:.4e} eV\n")
        print(f"Frequência do fóton (𝑓): {foton_frequencia:.4e} Hz\n")
        print(f"Comprimento de onda do fóton (𝜆): {foton_lambda:.4e} m\n")

        print(f"Energia do nível quântico inicial (E_i): {E_inicial_joules:.4e} Joules ({E_inicial_eV:.4e} eV)\n")
        print(f"Energia do nível quântico final (E_f): {E_final_joules:.4e} Joules ({E_final_eV:.4e} eV)\n")

        print(f"Velocidade no nível {n_inicial}: {velocidade_electron(n_inicial):.2f} m/s\n")
        print(f"Velocidade no nível {n_final}: {velocidade_electron(n_final):.2f} m/s\n")

        print(f"Comprimento de onda de De Broglie no nível inicial ({n_inicial}): {2*L/n_inicial} m\n")
        print(f"Comprimento de onda de De Broglie no nível final ({n_final}): {2*L/n_final} m\n")

        prob = probability(a, b, n_inicial, L)
        probf = probability(a, b, n_final, L)
        print(f"A probabilidade de encontrar a partícula entre {a} m e {b} m no nível n = {n_inicial} é {prob:.4f}.\n")
        print(f"A probabilidade de encontrar a partícula entre {a} m e {b} m no nível n = {n_final} é {probf:.4f}.\n")

        grafico_funcao_de_onda(L, n_inicial, n_final)
        grafico_distribuicao_probabilidades(L, n_inicial, n_final)

    elif n == 2:
        if menu_posicao() == 1:
            A = float(input("Digite o valor da amplitude (A): \n"))
            k = float(input("Digite o valor do número de onda (k): \n"))
            Xp = float(input("Digite a posição específica (Xp) onde deseja calcular a probabilidade: \n"))

            L = 2/A**2
            largura = L
            n = round((k*L)/pi)  

            probabilidade_xp = L*Xp
            probabilidade = A**2 * sin(k * probabilidade_xp)**2

            print(f"Largura da caixa (L): {L} metros\n")
            print(f"Nível quântico da partícula (n): {n}\n")
            print(f"Probabilidade de encontrar a partícula na posição {probabilidade_xp} é: {probabilidade:.5f}\n")
        elif menu_posicao() == 2:
            A = float(input("Digite o valor da amplitude (A): \n"))
            k = float(input("Digite o valor do número de onda (k): \n"))
            Xp = float(input("Digite a posição específica (Xp) onde deseja calcular a probabilidade: \n"))

            L = 2/A**2
            largura = L
            n = round((k*L)/pi)  
            probabilidade = A**2 * sin(k * Xp)**2

            print(f"Largura da caixa (L): {L:.2f} metros\n")
            print(f"Nível quântico da partícula (n): {n}\n")
            print(f"Probabilidade de encontrar a partícula na posição {Xp} é: {probabilidade:.5f}\n")
            
        else:
            print("Opcao invalida")
    elif n == 3:
        print("Simulacao")

        n_levels = 5
        energy_values = np.linspace(0.1, 0.4, n_levels)

       
        movement_sequence = [1, 1, 4, 2, 5, 3, 1, 3, 5, 2, 4, 1, 5, 3, 1, 4, 1]
        current_index = 0  
        current_level = movement_sequence[current_index] 

        def init():
            particle.set_data([], [])
            glow.set_radius(0)
            for i, line in enumerate(lines):
                line.set_data([0, 5], [energy_values[i]] * 2)
            return [particle, glow] + lines

        def animate(i):
            global current_index, current_level
            x = i % 500 / 100
            next_level_index = (current_index + 1) % len(movement_sequence)
            next_level = movement_sequence[next_level_index]

            if x == 4.9:  
                if next_level > current_level:
                    glow.set_color('blue')  
                elif next_level < current_level:
                    glow.set_color('green')  
                else:
                    glow.set_color('none')
                glow.set_radius(0.05)
            elif x == 0:  
                current_level = next_level
                current_index = next_level_index
            else:
                glow.set_radius(max(0, glow.get_radius() - 0.001))  
            particle.set_data(x, energy_values[current_level - 1])
            glow.center = (x, energy_values[current_level - 1])
            return [particle, glow] + lines

        fig, ax = plt.subplots()
        ax.set_xlim(0, 5)
        ax.set_ylim(0, 0.5)
        ax.set_ylabel("E (eV)")
        ax.set_xticks([])
        ax.set_yticks(energy_values)
        ax.set_yticklabels([f'E{i+1}' for i in range(n_levels)])
        ax.grid(True)

        lines = [ax.plot([], [], 'b-')[0] for _ in range(n_levels)]
        particle, = ax.plot([], [], 'ro', ms=8)
        glow = patches.Circle((0, 0), 0.05, color='none', transform=ax.transData)
        ax.add_patch(glow)

        ani = animation.FuncAnimation(fig, animate, init_func=init, frames=500, interval=3, blit=True, repeat=True)

        colors = {'blue': 'Absorve: Subir', 'green': 'Emite: Descer'}
        patches = [patches.Patch(color=color, label=label) for color, label in colors.items()]

        #PROTON
        E5 = (((5 ** 2) * (pi ** 2) * (1.0545718e-34 ** 2)) / (2 * 1.6726219e-27 * (largura ** 2)))
        E5f = "{:.4e}".format(E5)
        E4 = (((4 ** 2) * (pi ** 2) * (1.0545718e-34 ** 2)) / (2 * 1.6726219e-27 * (largura ** 2)))
        E4f = "{:.4e}".format(E4)
        E3 = (((3 ** 2) * (pi ** 2) * (1.0545718e-34 ** 2)) / (2 * 1.6726219e-27 * (largura ** 2)))
        E3f = "{:.4e}".format(E3)
        E2 = (((2 ** 2) * (pi ** 2) * (1.0545718e-34 ** 2)) / (2 * 1.6726219e-27 * (largura ** 2)))
        E2f = "{:.4e}".format(E2)
        E1 = (((1 ** 2) * (pi ** 2) * (1.0545718e-34 ** 2)) / (2 * 1.6726219e-27 * (largura ** 2)))
        E1f = "{:.4e}".format(E1)

        #ELETRON
        E5e = (((5 ** 2) * (pi ** 2) * (1.0545718e-34 ** 2)) / (2 * 9.11e-31 * (largura ** 2)))
        E5fe = "{:.4e}".format(E5e)
        E4e = (((4 ** 2) * (pi ** 2) * (1.0545718e-34 ** 2)) / (2 * 9.11e-31 * (largura ** 2)))
        E4fe = "{:.4e}".format(E4e)
        E3e = (((3 ** 2) * (pi ** 2) * (1.0545718e-34 ** 2)) / (2 * 9.11e-31 * (largura ** 2)))
        E3fe = "{:.4e}".format(E3e)
        E2e = (((2 ** 2) * (pi ** 2) * (1.0545718e-34 ** 2)) / (2 * 9.11e-31 * (largura ** 2)))
        E2fe = "{:.4e}".format(E2e)
        E1e = (((1 ** 2) * (pi ** 2) * (1.0545718e-34 ** 2)) / (2 * 9.11e-31 * (largura ** 2)))
        E1fe = "{:.4e}".format(E1e)

        legend_info = [f'E5: proton={str(E5f)} eletron={str(E5fe)}', f'E4: proton={str(E4f)} eletron={str(E4fe)}', f'E3: proton={str(E3f)} eletron={str(E3fe)}', f'E2: proton={str(E2f)} eletron={str(E2fe)}', f'E1: proton={str(E1f)} eletron={str(E1fe)}', f'Largura: {largura}']
        ax.legend(handles=patches, title="\n".join(legend_info), loc='upper right')

        plt.show()
    elif n == 4:
        print("Obrigado")
        break
    else:
        print("Entrada invalida")
