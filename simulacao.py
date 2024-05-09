import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.patches as patches
import numpy as np

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
ax.set_yticklabels([f'{i}' for i in range(n_levels)])
ax.grid(True)

lines = [ax.plot([], [], 'b-')[0] for _ in range(n_levels)]
particle, = ax.plot([], [], 'ro', ms=8)
glow = patches.Circle((0, 0), 0.05, color='none', transform=ax.transData)
ax.add_patch(glow)

ani = animation.FuncAnimation(fig, animate, init_func=init, frames=500, interval=5, blit=True, repeat=True)

colors = {'blue': 'Absorver: Subir', 'green': 'Emitir: Descer'}
patches = [patches.Patch(color=color, label=label) for color, label in colors.items()]
legend_info = ['L=Length']
ax.legend(handles=patches, title="\n".join(legend_info), loc='upper right')

plt.show()
