import numpy as np
import matplotlib.pyplot as plt

__all__ = [
    "diagram_transitions",
    "diagram_energies",
]

def diagram_transitions(multiplet, Bmax, EL=""):
    max_strength = max(gf_ for _, gf_, *_ in multiplet.transitions(0))
    B = np.array([0, Bmax])
    for tr in multiplet.transitions(B):
        x_line, gf_, *_, dmJ, _ = tr
        alpha = gf_/max_strength
        c = 'C4' if dmJ == 0 else 'C0'
        plt.plot(x_line, B, f'{c}-', lw=1.5, alpha=alpha)
    plt.xlabel(r"Wavelength [\AA]")
    plt.ylabel("Field Strength [kG]")
    plt.title(EL)
    plt.show()

def diagram_energies(multiplet, Bmax, EL=""):
    B = np.array([0, Bmax])
    for i, level in enumerate(multiplet.Upper_and_Lower_states):
        plt.subplot(2, 1, i+1)
        for state in level:
            for substate in state.energies(B):
                plt.plot(B, substate.k, f'C{i}-')
        plt.ylabel("Energy [1/cm]")
    plt.xlabel("Field Strength [kG]")
    plt.title(EL)
    plt.show()
