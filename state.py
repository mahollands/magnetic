"""
Classes for dealing with atomic states
"""

import functools
from collections import namedtuple
from fractions import Fraction
import numpy as np

SubState = namedtuple("SubState", "mJ k")

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
    def __init__(self, k0, J, L, S):
        self.k0 = float(k0)
        self.J = Fraction(J)
        self.L = int(L)
        self.S = Fraction(S)

        if any(X < 0 for X in self.JLS):
            raise ValueError("JLS must be non-negative")
        if self.J not in np.arange(abs(self.L-self.S), self.L+self.S+1):
            raise ValueError("Invalid JLS combination")


    def __repr__(self):
        Jdiv2 = "/2" if self.J.denominator == 2 else ""
        Sdiv2 = "/2" if self.S.denominator == 2 else ""
        Jstr = f"{self.J.numerator}{Jdiv2}"
        Sstr = f"{self.S.numerator}{Sdiv2}"
        return f"State({self.k0}, {Jstr}, {self.L}, {Sstr})"

    def __str__(self):
        """
        String conversion is given as a formatted Ket string, i.e. "|(2S+1)L_J>"
        """
        llist = 'SPDFGHIJ'
        if self.L >= len(llist):
            raise ValueError("L too large")
        suplist, sublist = '⁰¹²³⁴⁵⁶⁷⁸⁹', '₀₁₂₃₄₅₆₇₈₉'
        ss = suplist[self.S_multiplicity]
        ll = llist[self.L]
        div2 = "⸝₂" if self.J.denominator == 2 else ""
        jj = f"{sublist[self.J.numerator]}{div2}"
        return f"|{ss}{ll}{jj}⟩"

    @property
    def J_multiplicity(self):
        return int(2*self.J+1)

    @property
    def S_multiplicity(self):
        return int(2*self.S+1)

    @functools.cached_property
    def g(self):
        """
        Lande g-factor. Calculated as a fraction.
        """
        if self.S == 0:
            return Fraction(1, 1)
        if self.L == 0:
            return Fraction(2, 1)
        if self.L == self.S:
            return Fraction(3, 2)
        if self.J == 0:
            return Fraction(0, 1)
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
        return np.arange(-self.J, self.J+1)

    def splitting(self, mJ, B):
        """
        For a specific mJ and field strength B in kG, calculate the energy shift in 1/cm
        """
        return 0.046686 * B * float(mJ * self.g)

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
