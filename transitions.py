"""
Contains a class for dealing with multiplets between an levels, determining the
Zeeman splitting, and constructing synthetic line profiles.
"""
from itertools import product
from collections import defaultdict
from typing import List, Dict
import numpy as np
from scipy.special import voigt_profile
from .state import State

__all__ = ["Multiplet"]

class Multiplet:
    """
    Multiplet class, containing lists of lower and upper states, and oscillator
    strengths (log[gf]) between them.
    """

    def __init__(self,
        Lower_states: List[State],
        Upper_states: List[State],
        log_gf: Dict[tuple,float],
    ):
        self.Lower_states = Lower_states
        self.Upper_states = Upper_states
        self.log_gf = log_gf

    @property
    def Upper_and_Lower_states(self):
        return self.Upper_states, self.Lower_states

    @property
    def Lower_and_Upper_states(self):
        return self.Lower_states, self.Upper_states

    def states_outer_product(self):
        return product(self.Upper_states, self.Lower_states)

    @property
    def log_gf(self):
        return self._log_gf

    @log_gf.setter
    def log_gf(self, value):
        self._log_gf = defaultdict(lambda: -10.0, value)

    @staticmethod
    def relative_strength(Ji, Jf, mi, mf):
        """
        Calculate the relative strengths of Zeeman components from Höhn et al
        (1925). See also Table B1 from Dorsch et al. (2022)
        """
        dJ = int(Jf - Ji)
        dm = int(mf - mi)

        if abs(dm) > 1:
            return 0

        m_ = dm * mi

        if dJ == 0:
            return mi**2 if dm == 0 else (Ji-m_)*(Ji+m_+1)/4
        if dJ == 1:
            return (Ji+1)**2 - mi**2 if dm == 0 else (Ji+m_+1)*(Ji+m_+2)/4
        if dJ == -1:
            return Ji**2 - mi**2 if dm == 0 else (Ji-m_)*(Ji-m_-1)/4
        raise ValueError

    def transitions(self, B, T=6000.):
        """
        Generator for components in a Zeeman multiplet for an input field strength.
        Temperature can also be used to obtain a Boltzman factor from the lower
        energy.
        """
        for uS, lS in self.states_outer_product():
            if abs(uS.J - lS.J) > 1: #dJ selection rule
                continue

            if (uS.J, lS.J) == (0, 0): #J=0 to J=0 transitions forbidden
                continue

            l_substates, u_substates = lS.energies(B), uS.energies(B)

            total_strength = sum(self.relative_strength(lS.J, uS.J, ls.mJ, us.mJ) \
                for ls, us in product(l_substates, u_substates))

            gf0 = 10**self.log_gf[lS.J, uS.J]

            for ls, us in product(l_substates, u_substates):
                if abs(dmJ := int(us.mJ - ls.mJ)) > 1: #dmJ selection rule
                    continue
                boltz = np.exp(-ls.k/(0.695*T))
                w_line = 1e8/(us.k - ls.k)
                rel_strength = self.relative_strength(lS.J, uS.J, ls.mJ, us.mJ)
                gf_ = gf0 * rel_strength/total_strength
                yield w_line, gf_, lS, uS, dmJ, boltz

    def line_profile(self, B, x, res_l, res_g, psi, T, rv):
        """
        Generator for line profiles of components in the Zeeman multiplet.
        """
        z = 1 + rv/2.998e5
        cos2psi, sin2psi = np.cos(psi)**2, np.sin(psi)**2
        for tr in self.transitions(B, T=T):
            x_line, gf_, *_, dmJ, boltz = tr
            rot_factor = sin2psi if dmJ == 0 else 1+cos2psi
            V = voigt_profile(x-x_line*z, res_g/2.355, res_l/2)
            yield boltz*gf_*rot_factor * V

    def profile(self, B, x, strength, res_l, res_g, psi=1, T=6000., rv=0):
        """
        Caclulates a synthetic spectrum of a multiplet in a magnetic field
        B : field strength in [kG]
        x: vacuum wavelengths [AA]
        strength: dimensionless strength scaling parameter
        res_l: Lorentzian FWHM [AA]
        res_g: Gaussian FWHM [AA]
        psi: viewing angle [radians]
        T: temperature [K]
        rv: radial velocity [km/s]
        """
        ylines = sum(self.line_profile(B, x, res_l, res_g, psi, T, rv))
        return np.exp(-strength*ylines)
