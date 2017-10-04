#!/usr/bin/env python3
"""
Destexhe 2009 in NEST.

File: destexhe2009.py

Copyright 2017 Ankur Sinha
Author: Ankur Sinha <sanjay DOT ankur AT gmail DOT com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import nest
import nest.voltage_trace
import pylab
import matplotlib
import matplotlib.pyplot as plt


class Destexhe2009:

    """Replicate Destexhe2009."""

    def __init__(self):
        """Init params."""
        nest.ResetKernel()
        # not setting a, b here - they differ for different neuron sets so
        # we'll set them there. These are the common properties.
        self.S = 20000e-12  # m^2
        self.C_m = 1e6/1e-4  # pF/m^2
        self.g_L = 0.05e6/1e-4  # nS/m^2
        self.E_L = -60.  # mV
        self.V_reset = self.E_L
        self.Delta_T = 2.5  # mV
        self.V_th = -50.  # mV
        self.tau_w = 600.  # ms
        self.t_ref = 2.5  # ms
        self.I_e = 0.  # pA
        self.V_peak = 0.  # mV

        self.neuron_dict_common = {
            'C_m': self.C_m * self.S,
            't_ref': self.t_ref,
            'V_reset': self.V_reset,
            'E_L': self.E_L,
            'g_L': self.g_L * self.S,
            'I_e': self.I_e,
            'Delta_T': self.Delta_T,
            'tau_w': self.tau_w,
            'V_th': self.V_th,
            'V_peak': self.V_peak
        }

        # cell set specific dicts
        # remember that assignment does not copy, only refers
        self.dict_RS_strong = dict(self.neuron_dict_common)
        self.dict_RS_strong.update(
            {'a': 0.001e3,  # nS
             'b': 0.04e3  # pA
             }
        )
        nest.CopyModel('aeif_cond_exp', 'RS_strong_cell')
        nest.SetDefaults('RS_strong_cell', self.dict_RS_strong)

        self.dict_RS_weak = dict(self.neuron_dict_common)
        self.dict_RS_weak.update(
            {'a': 0.001e3,  # nS
             'b': 0.005e3  # pA
             }
        )
        nest.CopyModel('aeif_cond_exp', 'RS_weak_cell')
        nest.SetDefaults('RS_weak_cell', self.dict_RS_weak)

        self.dict_FS = dict(self.neuron_dict_common)
        self.dict_FS.update(
            {'a': 0.001e3,  # nS
             'b': 0.  # pA
             }
        )
        nest.CopyModel('aeif_cond_exp', 'FS_cell')
        nest.SetDefaults('FS_cell', self.dict_FS)

        self.dict_LTS = dict(self.neuron_dict_common)
        self.dict_LTS.update(
            {'a': 0.02e3,  # nS
             'b': 0.  # pA
             }
        )
        nest.CopyModel('aeif_cond_exp', 'LTS_cell')
        nest.SetDefaults('LTS_cell', self.dict_LTS)

        self.dict_TC = dict(self.neuron_dict_common)
        self.dict_TC.update(
            {'a': 0.04e3,  # nS
             'b': 0.  # pA
             }
        )
        nest.CopyModel('aeif_cond_exp', 'TC_cell')
        nest.SetDefaults('TC_cell', self.dict_TC)

        self.dict_RE = dict(self.neuron_dict_common)
        self.dict_RE.update(
            {'a': 0.08e3,  # nS
             'b': 0.03e3  # pA
             }
        )
        nest.CopyModel('aeif_cond_exp', 'RE_cell')
        nest.SetDefaults('RE_cell', self.dict_RE)

    def figure1(self):
        """Figure 1."""
        dc_stim = nest.Create('dc_generator')
        dc_properties = [
            {'amplitude': 0.25e3,  # pA
             'start': 200.,
             'stop': 400.
             }]
        voltmeter_properties = {'withgid': True,
                                'withtime': True,
                                'interval': 0.1}

        nest.SetStatus(dc_stim, dc_properties)

        RS_strong_cell = nest.Create('RS_strong_cell', 1)
        nest.Connect(dc_stim, RS_strong_cell, 'all_to_all')
        RS_strong_voltmeter = nest.Create('voltmeter')
        nest.SetStatus(RS_strong_voltmeter, voltmeter_properties)
        nest.Connect(RS_strong_voltmeter, RS_strong_cell)

        RS_weak_cell = nest.Create('RS_weak_cell', 1)
        nest.Connect(dc_stim, RS_weak_cell, 'all_to_all')
        RS_weak_voltmeter = nest.Create('voltmeter')
        nest.SetStatus(RS_weak_voltmeter, {'withgid': True,
                                           'withtime': True,
                                           'interval': 0.1})
        nest.Connect(RS_weak_voltmeter, RS_weak_cell)

        FS_cell = nest.Create('FS_cell', 1)
        nest.Connect(dc_stim, FS_cell, 'all_to_all')
        FS_voltmeter = nest.Create('voltmeter')
        nest.SetStatus(FS_voltmeter, {'withgid': True,
                                      'withtime': True,
                                      'interval': 0.1})
        nest.Connect(FS_voltmeter, FS_cell)

        LTS_cell = nest.Create('LTS_cell', 1)
        nest.Connect(dc_stim, LTS_cell, 'all_to_all')
        LTS_voltmeter = nest.Create('voltmeter')
        nest.SetStatus(LTS_voltmeter, {'withgid': True,
                                       'withtime': True,
                                       'interval': 0.1})
        nest.Connect(LTS_voltmeter, LTS_cell)

        nest.Simulate(600)

        events = nest.GetStatus(voltmeter, 'events')[0]
        times = events['times']
        V_m = events['V_m']
        print(times)
        print(V_m)


if __name__ == "__main__":
    sim = Destexhe2009()
    sim.figure1()
