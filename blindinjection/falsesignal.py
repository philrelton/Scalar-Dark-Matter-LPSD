#!/usr/bin/env python3

import numpy as np
import scipy.constants as consts

class FalseSignal(dict):
    def __init__(self, frequency, amplitude, Nfreqs=100, phase_seed=None, FWHM=1e-6, day=None, leap_year=False):
        self.frequency = frequency
        self.amplitude = amplitude
        self.Nfreqs = Nfreqs
        if phase_seed is None:
            self.phase_seed = int(np.random.randint(0, 2**32, 1)[0])
        else:
            self.phase_seed = phase_seed
        self.FWHM = FWHM
        self.day = day
        self.leap_year = leap_year
        _properties = ["frequencies", "phases", "amplitudes"]
        for _prop in _properties:
            self[_prop] = dict()
        self.doppler_shift()
        self.freq_array()
        self.phase_array()
        self.amplitude_array()

    @classmethod
    def random(cls):
        frequency = float(np.random.uniform(1000, 6000, 1))
        amplitude = float(np.random.uniform(5e-22, 1e-19, 1))
        Nfreqs = 100
        phase_seed = np.random.randint(0, 2**32, 1)
        FWHM = 1e-6
        day = np.random.randint(0, 365, 1)
        leap_year = False
        return cls(frequency, amplitude, Nfreqs, phase_seed, FWHM, day, leap_year)

    def doppler_shift(self):
        if self.day is not None:
            v_sun  = 2.4e5
            v_earth = 3.0e4
            if self.leap_year:
                t0 = 152
                ndays = 366
            else:
                t0 = 151
                ndays = 365
            v_obs = (v_sun
                    + (0.5
                       * v_earth
                       * np.cos(2 * np.pi * (self.day - t0) / ndays))
                    )
            # This mass equation is taken from Derevianko eq6
            m0 = (self.frequency * 1e-9 * consts.e
                 / (2 * np.pi * 2.42e5 * consts.c ** 2))
            delta_f = m0 * v_obs ** 2 / (2 * consts.hbar)
            self.frequency += delta_f
            return f"The change in frequency is {delta_f}, \
                     the new peak frequency is {self.frequency}"

    def freq_array(self):
        """Finds the position of frequencies with respect to the Full
           Width Half Maximum (FWHM).
           
           This depends on the representation in Figure 2 of Derevianko
           (2018). The curve here extends from -0.5 to 8 when the axis
           is centered on 0. Therefore one seventeenth of the curve is 
           before the true frequency.
           
           The FWHM of this curve is calculated to take up roughly 
           29.24274% of this range of the curve. Therefore from -0.5 
           to 8 the width is 3.419652194014651 times the FWHM. See
           curve_shape_calculation.ipynb for the calculation
        """

        signal_width = 3.419652194014651 * self.FWHM
        r = (self.frequency * signal_width) / self.Nfreqs
        start = self.frequency - np.ceil(self.Nfreqs / 17) * r
        end = start + self.Nfreqs * r
        self["frequencies"] = np.linspace(start, end, self.Nfreqs).reshape(self.Nfreqs)

    def phase_array(self):
        np.random.seed(self.phase_seed)
        self["phases"] = np.random.uniform(0, 2*np.pi, self.Nfreqs)

    def line_shape(self, eta=1):
        """This creates the shape of the curve in frequency space
           the equation is taken from Derevianko (2018) equation 10
        
           We have taken eta to be 1 as in the paper.
           x0 is the 'true' angular frequency of the signal
           colen is the coherence length of the signal
           x is the array of observed angular frequencies in the curve
           The values here are arbitrary as this is just to get the
           shape, the signal is scaled correctly elsewhere.
        """
        x0 = 1
        x = np.linspace(x0/2, 9*x0, self.Nfreqs)
        colen = 1/x0
        xt = colen * (x - x0)
        coeffs = colen * np.exp(-(eta ** 2)) / (np.sqrt(2 * np.pi) * eta)
        exponential = np.exp(-xt)
        hyperbole = np.sinh(eta * np.sqrt(eta ** 2 + 2 * xt))
        curve = coeffs * exponential * hyperbole
        norm = max(curve) * np.sqrt(self.Nfreqs)
        return curve / norm
    
    def amplitude_array(self):
        """This matches the curve to the desired signal amplitude"""
        curve = self.line_shape()
        const = float(self.amplitude)
        # The 1 / 0.317052 constant is derived from experimentation
        self["amplitudes"] = const * curve / 0.317052
