from math import * 
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.patches as patches
import numpy as np

def menu():
    print("Simulador Crocodilo de Perdizes\n")
    print("[1] Determinação da função de onda quântica e outros parâmetros\n")
    print("[2] Cálculo dos parâmetros da caixa e partícula, dada a função de onda\n")
    print("[3] Grafico\n")
    print("[4] Encerrar Programa\n")
    n = int(input("Escolha a sua opcao: "))
    return n

def menu_posicao():
    print("Você deseja entrar com a posicao em:\n")
    print("[1] Probabilidade (0 a 1)\n")
    print("[2] Metros\n")
    n = int(input())
    return n

def psi(x, n, L):
    amplitude = sqrt(2 / L)
    k = n*pi / L
    return amplitude, k

def energy_n(n, L):
    hbar = 1.0545718e-34  # Constante de Planck reduzida, em J*s
    print("Voce deseja calcular a energia de que particula?\n")
    print("[1] Proton")
    print("[2] Eletron")
    i = int(input())
    if i == 1:
        m = 1.6726219e-27  # Massa do proton, em kg
    else:
        m = 9.11e-31 # Massa do eletron, em kg
    E_joules = ((n ** 2) * (pi ** 2) * (hbar ** 2)) / (2 * m * (L ** 2))
    E_eV = E_joules / 1.60218e-19  # Conversão de Joule para eV
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
    hbar = 1.0545718e-34  # Constante de Planck reduzida, em J*s
    print(f"Voce deseja calcular a velocidade no nivel {n} de que particula?\n")
    print("[1] Proton")
    print("[2] Eletron")
    i = int(input())
    if i == 1:
        m = 1.6726219e-27  # Massa do proton, em kg
    else:
        m = 9.11e-31 # Massa do eletron, em kg
    E_joules = ((n ** 2) * (pi ** 2) * (hbar ** 2)) / (2 * m * (L ** 2))
    E_eV = E_joules / 1.60218e-19  # Conversão de Joule para eV
    p = sqrt(2 * m * E_joules)
    v_n = p / m
    return v_n

def grafico_funcao_de_onda(L, A, ki, kf):
    # Valores de x para o grafico
    x_vals = np.linspace(0, L * 1e-9, 1000)
    psi_inicial_vals = A * np.sin(ki * x_vals)
    psi_final_vals = A * np.sin(kf * x_vals)

    # Normalização das funções de onda para evitar distorções no gráfico
    psi_inicial_vals_normalizado = psi_inicial_vals / np.max(np.abs(psi_inicial_vals))
    psi_final_vals_normalizado = psi_final_vals / np.max(np.abs(psi_final_vals))

    # Plotar as funções de onda
    plt.figure(figsize=(10, 6))
    plt.plot(x_vals, psi_inicial_vals_normalizado, label='Nível Inicial')
    plt.plot(x_vals, psi_final_vals_normalizado, label='Nível Final')
    plt.xlabel('Posição x (m)')
    plt.ylabel('𝜓(x)')
    plt.title('Gráfico da função de onda inicial e final')
    plt.legend()
    plt.grid(True)
    plt.ylim(-1.1, 1.1) # Definir os limites ideias para o eixo y
    plt.show()

def grafico_distribuicao_probabilidades(L, A, ki, kf):
    # Valores de x para o gráfico
    x_vals = np.linspace(0, L * 1e-9, 1000)
    psi_inicial_vals = (A * np.sin(ki * x_vals))**2
    psi_final_vals = (A * np.sin(kf * x_vals))**2

    plt.figure(figsize=(10, 6))
    plt.plot(x_vals, psi_inicial_vals, label='Nível Inicial')
    plt.plot(x_vals, psi_final_vals, label='Nível Final')
    plt.xlabel('Posição x (m)')
    plt.ylabel('𝜓(x)^2')
    plt.title('Gráfico da distribuição de probabilidade inicial e final')
    plt.legend()
    plt.grid(True)
    plt.show()

while(True):
    n = menu()
    if n == 1:
        L = float(input("Digite a largura da caixa (L): "))
        n_inicial = int(input("Digite o nível quântico inicial (n inicial): "))
        n_final = int(input("Digite o nível quântico final (n final): "))
        while True:
            a = float(input("Digite o limite inferior (a): "))
            b = float(input("Digite o limite superior (b): "))
            if 0 <= a <= b <= L:
                break
            else:
                print(f"Por favor, garanta que 0 <= a <= b <= L ({L}). Tente novamente.")
        amplitude_inicial, k_inicial = psi(a, n_inicial, L)
        amplitude_final, k_final = psi(b, n_final, L)

        print(f"Função de onda para n inicial: 𝜓(x) = {amplitude_inicial:.3e} sin({k_inicial:.3e} x)")
        print(f"Função de onda para n final: 𝜓(x) = {amplitude_final:.3e} sin({k_final:.3e} x)")

        E_inicial_joules, E_inicial_eV = energy_n(n_inicial, L)
        E_final_joules, E_final_eV = energy_n(n_final, L)

        # Cálculos das energias iniciais e finais
        E_inicial_joules, E_inicial_eV = energy_n(n_inicial, L)
        E_final_joules, E_final_eV = energy_n(n_final, L)
        
        # Diferença de energia entre os níveis
        delta_E_joules = abs(E_final_joules - E_inicial_joules)
        delta_E_eV = abs(E_final_eV - E_inicial_eV)
        
        # Cálculo da frequência do fóton
        h = 6.62607015e-34  # Constante de Planck em J*s
        foton_frequencia = delta_E_joules / h

        # Cálculo do comprimento de onda do fóton
        c = 3e8  # Velocidade da luz em m/s
        foton_lambda = c / foton_frequencia if foton_frequencia != 0 else 0
        
        print(f"Energia do fóton (𝐸_𝑓ó𝑡𝑜𝑛): {delta_E_eV:.4e} eV")
        print(f"Frequência do fóton (𝑓): {foton_frequencia:.4e} Hz")
        print(f"Comprimento de onda do fóton (𝜆): {foton_lambda:.4e} m")

        print(f"Energia do nível quântico inicial (E_i): {E_inicial_joules:.4e} Joules ({E_inicial_eV:.4e} eV)")
        print(f"Energia do nível quântico final (E_f): {E_final_joules:.4e} Joules ({E_final_eV:.4e} eV)")

        print(f"Velocidade no nível {n_inicial}: {velocidade_electron(n_inicial):.2f} m/s")
        print(f"Velocidade no nível {n_final}: {velocidade_electron(n_final):.2f} m/s")

        print(f"Comprimento de onda de De Broglie no nível inicial ({n_inicial}): {2*L/n_inicial} unidades")
        print(f"Comprimento de onda de De Broglie no nível final ({n_final}): {2*L/n_final} unidades")

        prob = probability(a, b, n_inicial, L)
        print(f"A probabilidade de encontrar a partícula entre {a} m e {b} m no nível n = {n_inicial} é {prob:.4f}.")

    elif n == 2:
        if menu_posicao() == 1:
            A = float(input("Digite o valor da amplitude (A): "))
            k = float(input("Digite o valor do número de onda (k): "))
            Xp = float(input("Digite a posição específica (Xp) onde deseja calcular a probabilidade: "))

            L = 2/A**2
            n = round((k*L)/pi)  

            probabilidade_xp = L*Xp
            probabilidade = A**2 * sin(k * probabilidade_xp)**2

            print(f"Largura da caixa (L): {L} metros")
            print(f"Nível quântico da partícula (n): {n}")
            print(f"Probabilidade de encontrar a partícula na posição {probabilidade_xp} é: {probabilidade:.5f}")
        elif menu_posicao() == 2:
            A = float(input("Digite o valor da amplitude (A): "))
            k = float(input("Digite o valor do número de onda (k): "))
            Xp = float(input("Digite a posição específica (Xp) onde deseja calcular a probabilidade: "))

            L = 2/A**2
            n = round((k*L)/pi)  
            probabilidade = A**2 * sin(k * Xp)**2

            print(f"Largura da caixa (L): {L} metros")
            print(f"Nível quântico da partícula (n): {n}")
            print(f"Probabilidade de encontrar a partícula na posição {Xp} é: {probabilidade:.5f}")
        else:
            print("Opcao invalida")
    elif n == 3:
        print("Simulacao")

        n_levels = 5
        energy_values = np.linspace(0.1, 0.4, n_levels)

        # Sequência específica de níveis para o movimento da partícula
        movement_sequence = [1, 1, 4, 2, 5, 3, 1, 3, 5, 2, 4, 1, 5, 3, 1, 4, 1]
        current_index = 0  # Índice inicial na sequência de movimento
        current_level = movement_sequence[current_index]  # Nível inicial baseado no primeiro índice

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

            if x == 4.9:  # Atualiza o brilho para o próximo movimento
                if next_level > current_level:
                    glow.set_color('blue')  # Azul para subir
                elif next_level < current_level:
                    glow.set_color('green')  # Verde para descer
                else:
                    glow.set_color('none')
                glow.set_radius(0.05)
            elif x == 0:  # Atualiza o nível da partícula no começo da linha
                current_level = next_level
                current_index = next_level_index
            else:
                glow.set_radius(max(0, glow.get_radius() - 0.001))  # Faz o brilho diminuir gradualmente
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

        colors = {'blue': 'Próximo: Subir', 'green': 'Próximo: Descer'}
        patches = [patches.Patch(color=color, label=label) for color, label in colors.items()]
        legend_info = ['E5', 'E4', 'E3', 'E2', 'E1', 'L=Length']
        ax.legend(handles=patches, title="\n".join(legend_info), loc='upper right')

        plt.show()
    elif n == 4:
        break
    else:
        print("Entrada invalida")
