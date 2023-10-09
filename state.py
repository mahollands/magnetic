import functools
from dataclasses import dataclass
from collections import namedtuple
import numpy as np

SubState = namedtuple("SubState", "mJ k")

@dataclass(frozen=True)
class State:
    """
    This class is an energy level requiring 4 quantities:

    k0 = unperturbed energy of the state (cm-1)
    J  = total angular momenutm qunatum number
    L  = orbintal angular momenutm qunatum number
    S  = spin angular momenutm qunatum number

    Example:
    >>> LVL = State(15000., 1.5, 2, 0.5)

    On intitialisation, the class variable 'g', the Lande-g factor
    is also calculated and can get be accessed as LVL.g

    Various methods are included to calculate perturbations of
    energy levels in a 'weak' magnetic field:

    LVL.splitting(mJ, B)
    LVL.splittings(B)
    LVL.energy(mJ, B)
    LVL.energies(B)

    'splittings' and 'energies' do not require an mJ value,
    but return a list of (mJ, dk) or (mJ, k) tuples.
    """
    k0: float
    J: float
    L: int
    S: float

    def __post_init__(self):
        """
        Check J L S values are okay
        """
        if any(X < 0 for X in self.JLS):
            raise ValueError("JLS must be non-negative")
        if self.J not in np.arange(abs(self.L-self.S), self.L+self.S+0.5, 1.0):
            raise ValueError("Invalid JLS combination")


    def __str__(self):
        """
        Returns a formatted Ket string as "|(2S+1)L_J>"
        """
        llist = 'SPDFGHIJ'
        if self.L >= len(llist):
            raise ValueError("L too large")
        suplist, sublist = '⁰¹²³⁴⁵⁶⁷⁸⁹', '₀₁₂₃₄₅₆₇₈₉'
        ss = suplist[self.S_multiplicity]
        ll = llist[self.L]
        jj = f"{sublist[int(self.J)]}" if (self.J % 1.0 == 0.0) else \
            f"{sublist[int(2*self.J)]}⸝₂"
        return f"|{ss}{ll}{jj}⟩"

    @property
    def J_multiplicity(self):
        return int(2*self.J)+1

    @property
    def S_multiplicity(self):
        return int(2*self.S)+1

    @functools.cached_property
    def g(self):
        if self.S == 0:
            return 1
        if self.L == 0:
            return 2
        if self.L == self.S:
            return 3/2
        if self.J == 0:
            return 0
        JJ1, LL1, SS1 = (X*(X+1) for X in self.JLS)
        return 1 + (JJ1 + SS1 - LL1) / (2*JJ1)

    @property
    def JLS(self):
        """
        Property to unpack J, L, S angular momenta simultaneously
        """
        return self.J, self.L, self.S

    @property
    def mJs(self):
        """
        Property to get all values of mJ between -mJ and +mJ
        """
        return np.arange(-self.J, self.J+0.5, 1.0)

    def splitting(self, mJ, B):
        """
        For a specific mJ and field strength B in kG, calculate the energy shift in 1/cm
        """
        return 0.046686 * B * mJ * self.g

    def splittings(self, B):
        """
        For a field strength B in kG, calculate the energy shifts in 1/cm
        for all substates mJ
        """
        return [SubState(mJ, self.splitting(mJ, B)) for mJ in self.mJs]

    def energy(self, mJ, B):
        """
        For a specific mJ and field strength B, calculate the total energy in 1/cm
        """
        return self.k0 + self.splitting(mJ, B)

    def energies(self, B):
        """
        For a field strength B, calculate the total energies in 1/cm
        for all substates mJ
        """
        return [SubState(mJ, self.energy(mJ, B)) for mJ in self.mJs]

