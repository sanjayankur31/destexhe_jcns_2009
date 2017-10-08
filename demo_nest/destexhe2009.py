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
import subprocess
import numpy


class Destexhe2009:

    """Replicate Destexhe2009."""

    def __init__(self):
        """Init params."""
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
        self.V_peak = 40.  # mV
        self.E_ex = 0.  # mV
        self.E_in = -80.  # mV
        self.tau_syn_ex = 5.  # ms
        self.tau_syn_in = 10.  # ms

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
            'V_peak': self.V_peak,
            'E_ex': self.E_ex,
            'E_in': self.E_in,
            'tau_syn_ex': self.tau_syn_ex,
            'tau_syn_in': self.tau_syn_in

        }

        # cell set specific dicts
        # remember that assignment does not copy, only refers
        self.dict_RS_strong = dict(self.neuron_dict_common)
        self.dict_RS_strong.update(
            {'a': 0.001e3,  # nS
             'b': 0.04e3  # pA
             }
        )
        self.dict_RS_weak = dict(self.neuron_dict_common)
        self.dict_RS_weak.update(
            {'a': 0.001e3,  # nS
             'b': 0.005e3  # pA
             }
        )
        self.dict_FS = dict(self.neuron_dict_common)
        self.dict_FS.update(
            {'a': 0.001e3,  # nS
             'b': 0.  # pA
             }
        )
        self.dict_LTS = dict(self.neuron_dict_common)
        self.dict_LTS.update(
            {'a': 0.02e3,  # nS
             'b': 0.  # pA
             }
        )
        self.dict_TC = dict(self.neuron_dict_common)
        self.dict_TC.update(
            {'a': 0.04e3,  # nS
             'b': 0.  # pA
             }
        )
        self.dict_RE = dict(self.neuron_dict_common)
        self.dict_RE.update(
            {'a': 0.08e3,  # nS
             'b': 0.03e3  # pA
             }
        )

    def __setup(self):
        """Setup neuron models."""
        nest.ResetKernel()
        nest.SetKernelStatus(
            {
                'resolution': 0.1,
                'overwrite_files': True,
            })
        nest.CopyModel('aeif_cond_exp', 'RS_strong_cell')
        nest.SetDefaults('RS_strong_cell', self.dict_RS_strong)

        nest.CopyModel('aeif_cond_exp', 'RS_weak_cell')
        nest.SetDefaults('RS_weak_cell', self.dict_RS_weak)

        nest.CopyModel('aeif_cond_exp', 'FS_cell')
        nest.SetDefaults('FS_cell', self.dict_FS)

        nest.CopyModel('aeif_cond_exp', 'LTS_cell')
        nest.SetDefaults('LTS_cell', self.dict_LTS)

        nest.CopyModel('aeif_cond_exp', 'TC_cell')
        nest.SetDefaults('TC_cell', self.dict_TC)

        nest.CopyModel('aeif_cond_exp', 'RE_cell')
        nest.SetDefaults('RE_cell', self.dict_RE)

    def figure1(self):
        """Figure 1."""
        self.__setup()
        # set up depolarizing step current
        dc_stim_depol = nest.Create('dc_generator')
        dc_depol_properties = [
            {'amplitude': 0.25e3,  # pA
             'start': 200.,
             'stop': 600.
             }]
        nest.SetStatus(dc_stim_depol, dc_depol_properties)

        # set up hyperpolarizing step current
        dc_stim_hyperpol = nest.Create('dc_generator')
        dc_hyperpol_properties = [
            {'amplitude': -0.25e3,  # pA
             'start': 200.,
             'stop': 600.
             }]
        nest.SetStatus(dc_stim_hyperpol, dc_hyperpol_properties)

        # Set up voltmeters - one each for each cell
        voltmeter_properties = {'withgid': True,
                                'withtime': True,
                                'interval': 0.1,
                                'to_file': True,
                                }
        voltmeter1 = nest.Create('voltmeter')
        voltmeter2 = nest.Create('voltmeter')
        voltmeter3 = nest.Create('voltmeter')
        voltmeter4 = nest.Create('voltmeter')
        voltmeter5 = nest.Create('voltmeter')
        voltmeter6 = nest.Create('voltmeter')
        voltmeter7 = nest.Create('voltmeter')
        voltmeter8 = nest.Create('voltmeter')
        voltmeter9 = nest.Create('voltmeter')
        nest.SetStatus(voltmeter1, voltmeter_properties)
        nest.SetStatus(voltmeter2, voltmeter_properties)
        nest.SetStatus(voltmeter3, voltmeter_properties)
        nest.SetStatus(voltmeter4, voltmeter_properties)
        nest.SetStatus(voltmeter5, voltmeter_properties)
        nest.SetStatus(voltmeter6, voltmeter_properties)
        nest.SetStatus(voltmeter7, voltmeter_properties)
        nest.SetStatus(voltmeter8, voltmeter_properties)
        nest.SetStatus(voltmeter9, voltmeter_properties)

        # individual neurons
        RS_strong_cell = nest.Create('RS_strong_cell', 1)
        print("GID of RS_strong_cell is: {}".format(RS_strong_cell))
        nest.Connect(dc_stim_depol, RS_strong_cell, 'all_to_all')
        nest.Connect(voltmeter1, RS_strong_cell)

        RS_weak_cell = nest.Create('RS_weak_cell', 1)
        print("GID of RS_weak_cell is: {}".format(RS_weak_cell))
        nest.Connect(dc_stim_depol, RS_weak_cell, 'all_to_all')
        nest.Connect(voltmeter2, RS_weak_cell)

        FS_cell = nest.Create('FS_cell', 1)
        print("GID of FS_cell is: {}".format(FS_cell))
        nest.Connect(dc_stim_depol, FS_cell, 'all_to_all')
        nest.Connect(voltmeter3, FS_cell)

        LTS_cell = nest.Create('LTS_cell', 1)
        print("GID of LTS_cell is: {}".format(LTS_cell))
        nest.Connect(dc_stim_depol, LTS_cell, 'all_to_all')
        nest.Connect(voltmeter4, LTS_cell)

        TC_cell = nest.Create('TC_cell', 1)
        print("GID of TC_cell is: {}".format(TC_cell))
        nest.Connect(dc_stim_depol, TC_cell, 'all_to_all')
        nest.Connect(voltmeter5, TC_cell)

        RE_cell = nest.Create('RE_cell', 1)
        print("GID of RE_cell is: {}".format(RE_cell))
        nest.Connect(dc_stim_depol, RE_cell, 'all_to_all')
        nest.Connect(voltmeter6, RE_cell)

        LTS_cell_2 = nest.Create('LTS_cell', 1)
        print("GID of LTS_cell_2 is: {}".format(LTS_cell_2))
        nest.Connect(dc_stim_hyperpol, LTS_cell_2, 'all_to_all')
        nest.Connect(voltmeter7, LTS_cell_2)

        TC_cell_2 = nest.Create('TC_cell', 1)
        print("GID of TC_cell_2 is: {}".format(TC_cell_2))
        nest.Connect(dc_stim_hyperpol, TC_cell_2, 'all_to_all')
        nest.Connect(voltmeter8, TC_cell_2)

        RE_cell_2 = nest.Create('RE_cell', 1)
        print("GID of RE_cell_2 is: {}".format(RE_cell_2))
        nest.Connect(dc_stim_hyperpol, RE_cell_2, 'all_to_all')
        nest.Connect(voltmeter9, RE_cell_2)

        nest.Simulate(1000)

        # plot the graph
        args = ['gnuplot', 'figure1.plt']
        subprocess.call(args)

    def figure2(self):
        """
        Figure 2.


        Unable to replicate - stimulus unclear.
        """
        self.__setup()
        TC_cells = nest.Create('TC_cell', 2)
        RE_cells = nest.Create('RE_cell', 2)

        nest.Connect(TC_cells, RE_cells,
                     conn_spec={'rule': 'one_to_one'},
                     syn_spec={'model': 'static_synapse',
                               'weight': 30.}
                     )
        nest.Connect(TC_cells, list(reversed(RE_cells)),
                     conn_spec={'rule': 'one_to_one'},
                     syn_spec={'model': 'static_synapse',
                               'weight': 30.}
                     )
        for i in range(0, 4):
            nest.Connect(RE_cells, TC_cells,
                         conn_spec={'rule': 'one_to_one'},
                         syn_spec={'model': 'static_synapse',
                                   'weight': -30.}
                         )
            nest.Connect(RE_cells, list(reversed(TC_cells)),
                         conn_spec={'rule': 'one_to_one'},
                         syn_spec={'model': 'static_synapse',
                                   'weight': -30.}
                         )
        for i in range(0, 8):
            nest.Connect(RE_cells, list(reversed(RE_cells)),
                         conn_spec={'rule': 'one_to_one'},
                         syn_spec={'model': 'static_synapse',
                                   'weight': -30.}
                         )
        print("From RE: {} ".format(nest.GetConnections(source=RE_cells)))
        print("From TC: {}".format(nest.GetConnections(source=TC_cells)))

        voltmeter_properties = {'withgid': True,
                                'withtime': True,
                                'interval': 0.1,
                                'to_file': True,
                                }
        voltmeter1 = nest.Create('voltmeter')
        voltmeter2 = nest.Create('voltmeter')
        voltmeter3 = nest.Create('voltmeter')
        voltmeter4 = nest.Create('voltmeter')
        nest.SetStatus(voltmeter1, voltmeter_properties)
        nest.SetStatus(voltmeter2, voltmeter_properties)
        nest.SetStatus(voltmeter3, voltmeter_properties)
        nest.SetStatus(voltmeter4, voltmeter_properties)

        nest.Connect(voltmeter1, [RE_cells[0]])
        nest.Connect(voltmeter2, [RE_cells[1]])
        nest.Connect(voltmeter3, [TC_cells[0]])
        nest.Connect(voltmeter4, [TC_cells[1]])

        mystim = nest.Create('poisson_generator', 1,
                             {'rate': 200.,
                              'start': 150., 'stop': 200.})
#
        nest.Connect(mystim, [TC_cells[0]])
        nest.Connect(mystim, [TC_cells[1]])
        nest.Simulate(1500.)

        # plot the graph
        args = ['gnuplot', '-e',
                'outputfile="Figure2.png"', 'figure2.plt']
        subprocess.call(args)

    def figure3(self):
        """
        Figure 3.

        Not replicated yet.
        """
        self.__setup()
        TC_cells = nest.Create('TC_cell', 10)
        RE_cells = nest.Create('RE_cell', 10)

        # fixed indegree - to ensure that scaling doesn't affect network much
        nest.Connect(RE_cells, TC_cells,
                     {'rule': 'pairwise_bernoulli',
                      'p': 0.08},
                     syn_spec={'model': 'static_synapse',
                               'weight': -67.}
                     )
        nest.Connect(RE_cells, RE_cells,
                     {'rule': 'pairwise_bernoulli',
                      'p': 0.08},
                     syn_spec={'model': 'static_synapse',
                               'weight': -67.}
                     )
        nest.Connect(TC_cells, RE_cells,
                     {'rule': 'pairwise_bernoulli',
                      'p': 0.02},
                     syn_spec={'model': 'static_synapse',
                               'weight': 6.}
                     )

        stim = nest.Create('poisson_generator', 1, {'rate': 300., 'start': 0.,
                                                    'stop': 50.})
        stim_TC_cells = TC_cells
        nest.Connect(stim, stim_TC_cells,
                     {'rule': 'fixed_indegree', 'indegree': 20},
                     {'model': 'static_synapse', 'weight': 6.}
                     )

        detector = nest.Create('spike_detector', params={
                                   'to_file': True,
                                   'to_memory': False,
                                }
                               )
        nest.Connect(RE_cells, detector)
        nest.Connect(TC_cells, detector)

        nest.Simulate(51000.)

        # plot the graph
        args = ['gnuplot', 'figure3.plt']
        subprocess.call(args)

    def figure4(self):
        """
        Figure 4.

        Not replicated yet.
        """
        self.__setup()
        TC_cells = nest.Create('TC_cell', 50)
        RE_cells = nest.Create('RE_cell', 50)

        # fixed indegree - to ensure that scaling doesn't affect network much
        nest.Connect(RE_cells, TC_cells,
                     {'rule': 'pairwise_bernoulli',
                      'p': 0.08},
                     syn_spec={'model': 'static_synapse',
                               'weight': -67.}
                     )
        nest.Connect(RE_cells, RE_cells,
                     {'rule': 'pairwise_bernoulli',
                      'p': 0.08},
                     syn_spec={'model': 'static_synapse',
                               'weight': -67.}
                     )
        nest.Connect(TC_cells, RE_cells,
                     {'rule': 'pairwise_bernoulli',
                      'p': 0.02},
                     syn_spec={'model': 'static_synapse',
                               'weight': 6.}
                     )

        stim = nest.Create('poisson_generator', 1, {'rate': 400., 'origin': 0.,
                                                    'start': 0.,
                                                    'stop': 50.})
        stim_TC_cells = TC_cells
        nest.Connect(stim, stim_TC_cells,
                     {'rule': 'fixed_indegree', 'indegree': 20},
                     {'model': 'static_synapse', 'weight': 6.}
                     )

        detector = nest.Create('spike_detector', params={
                                   'to_file': True,
                                   'to_memory': False,
                                }
                               )
        nest.Connect(RE_cells, detector)
        nest.Connect(TC_cells, detector)

        nest.Simulate(51000.)

        # plot the graph
        args = ['gnuplot', 'figure4.plt']
        subprocess.call(args)


if __name__ == "__main__":
    sim = Destexhe2009()
    #  sim.figure1()
    #  sim.figure2(isi=i)
    #  sim.figure3()
    #  sim.figure4()
