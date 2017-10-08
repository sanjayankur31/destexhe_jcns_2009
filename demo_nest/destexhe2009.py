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
import subprocess
import random
import logging


class Destexhe2009:

    """Replicate Destexhe2009."""

    def __init__(self):
        """Init params."""
        # not setting a, b here - they differ for different neuron sets so
        self.dt = 0.1

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
             'b': 0.01e3  # pA
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
        logging.info("Setting up NEST and neuron models")
        nest.ResetKernel()
        nest.SetKernelStatus(
            {
                'resolution': self.dt,
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
        logging.info("FIGURE 1: Different spiking behaviour traces.")
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
        nest.Connect(dc_stim_depol, RS_strong_cell, 'all_to_all')
        nest.Connect(voltmeter1, RS_strong_cell)

        RS_weak_cell = nest.Create('RS_weak_cell', 1)
        nest.Connect(dc_stim_depol, RS_weak_cell, 'all_to_all')
        nest.Connect(voltmeter2, RS_weak_cell)

        FS_cell = nest.Create('FS_cell', 1)
        nest.Connect(dc_stim_depol, FS_cell, 'all_to_all')
        nest.Connect(voltmeter3, FS_cell)

        LTS_cell = nest.Create('LTS_cell', 1)
        nest.Connect(dc_stim_depol, LTS_cell, 'all_to_all')
        nest.Connect(voltmeter4, LTS_cell)

        TC_cell = nest.Create('TC_cell', 1)
        nest.Connect(dc_stim_depol, TC_cell, 'all_to_all')
        nest.Connect(voltmeter5, TC_cell)

        RE_cell = nest.Create('RE_cell', 1)
        nest.Connect(dc_stim_depol, RE_cell, 'all_to_all')
        nest.Connect(voltmeter6, RE_cell)

        LTS_cell_2 = nest.Create('LTS_cell', 1)
        nest.Connect(dc_stim_hyperpol, LTS_cell_2, 'all_to_all')
        nest.Connect(voltmeter7, LTS_cell_2)

        TC_cell_2 = nest.Create('TC_cell', 1)
        nest.Connect(dc_stim_hyperpol, TC_cell_2, 'all_to_all')
        nest.Connect(voltmeter8, TC_cell_2)

        RE_cell_2 = nest.Create('RE_cell', 1)
        nest.Connect(dc_stim_hyperpol, RE_cell_2, 'all_to_all')
        nest.Connect(voltmeter9, RE_cell_2)

        nest.Simulate(1000)

        # plot the graph
        args = ['gnuplot', 'figure1.plt']
        subprocess.call(args)

    def figure2(self, num_TC=100, num_RE=100):
        """
        Figure 2.

        """
        self.__setup()
        TC_cells = nest.Create('TC_cell', num_TC)
        RE_cells = nest.Create('RE_cell', num_RE)

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
                             {'rate': 1000./70.,
                              'start': 0., 'stop': 50.})
#
        nest.Connect(mystim, [TC_cells[0]])
        nest.Connect(mystim, [TC_cells[1]])
        nest.Simulate(1500.)

        # plot the graph
        args = ['gnuplot', '-e',
                'outputfile="Figure2.png"', 'figure2.plt']
        subprocess.call(args)

    def thalamic(self, num_neurons=200, outputfile="Figure3"):
        """
        Figure 3 and 4 - oscillations in thalamic neuronas.

        Not replicated yet.
        """
        self.__setup()
        num_TC = num_neurons/2
        num_RE = num_neurons/2
        scale = 200/num_neurons
        TC_cells = nest.Create('TC_cell', num_TC)
        RE_cells = nest.Create('RE_cell', num_RE)

        nest.Connect(RE_cells, TC_cells,
                     {'rule': 'pairwise_bernoulli',
                      'autapses': False,
                      'multapses': False,
                      'p': 0.08*scale},
                     syn_spec={'model': 'static_synapse',
                               'delay': self.dt,
                               'weight': -67.}
                     )
        nest.Connect(RE_cells, RE_cells,
                     {'rule': 'pairwise_bernoulli',
                      'autapses': False,
                      'multapses': False,
                      'p': 0.08*scale},
                     syn_spec={'model': 'static_synapse',
                               'delay': self.dt,
                               'weight': -67.}
                     )
        nest.Connect(TC_cells, RE_cells,
                     {'rule': 'pairwise_bernoulli',
                      'autapses': False,
                      'multapses': False,
                      'p': 0.02*scale},
                     syn_spec={'model': 'static_synapse',
                               'delay': self.dt,
                               'weight': 6.}
                     )

        stim_cells = random.sample(TC_cells, int(len(TC_cells)/5))
        for stim_cell in stim_cells:
            stim = nest.Create('poisson_generator', 20,
                               {'rate': 1000./70.,
                                'start': 0., 'stop': 50.})
            nest.Connect(stim, [stim_cell],
                         syn_spec={'model': 'static_synapse',
                                   'delay': self.dt, 'weight': 6.*scale}
                         )

        detector = nest.Create('spike_detector', params={
                                   'to_file': True,
                                   'to_memory': False,
                                }
                               )
        nest.Connect(RE_cells, detector)
        nest.Connect(TC_cells, detector)

        print("TC -> RE: {} synapses per neuron".format(
            len(nest.GetConnections(source=TC_cells,
                                    target=RE_cells))/num_RE))
        print("From RE -> RE: {} synapses per neuron".format(
            len(nest.GetConnections(
                source=RE_cells, target=RE_cells))/num_RE))
        print("From RE -> TC: {} synapses per neuron".format(
            len(nest.GetConnections(
                source=RE_cells, target=TC_cells))/num_TC))

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

        nest.Simulate(500.)

        # plot the graph
        args = ['gnuplot',
                '-e', 'outputfile="{}"'.format(outputfile),
                '-e', 'sd="{}"'.format(detector[0]),
                '-e', 'v1="{}"'.format(voltmeter1[0]),
                '-e', 'v2="{}"'.format(voltmeter2[0]),
                '-e', 'v3="{}'.format(voltmeter3[0]),
                '-e', 'v4="{}"'.format(voltmeter4[0]),
                'thalamic.plt'
                ]
        subprocess.call(args)

    def corticals6(self, num_neurons=2000, outputfile="Figure6"):
        """
        Both cortical diagrams.
        """
        self.cortical6('RS_strong_cell', num_neurons, "Figure6a")
        self.cortical6('RS_weak_cell', num_neurons, "Figure6b")

    def cortical6(self, PY_cell_model, num_neurons, outputfile):
        """
        Figure 6.
        """
        self.__setup()
        scale = 2000/num_neurons
        num_PY = int(0.8*num_neurons)
        num_IN = int(0.2*num_neurons)
        PY_cells = nest.Create(PY_cell_model, num_PY)
        IN_cells = nest.Create('FS_cell', num_IN)

        # the indegree must remain the same as for a complete circuit
        # this can either be done by increasing the probability of connections
        # by multiplying it by the scale factor, or by using a constant
        # indegree
        nest.Connect(IN_cells, PY_cells,
                     {'rule': 'pairwise_bernoulli',
                      'multapses': False,
                      'autapses': False,
                      'p': 0.02*scale},
                     syn_spec={'model': 'static_synapse',
                               'delay': self.dt,
                               'weight': -67.}
                     )
        nest.Connect(IN_cells, IN_cells,
                     {'rule': 'pairwise_bernoulli',
                      'autapses': False,
                      'multapses': False,
                      'p': 0.02*scale},
                     syn_spec={'model': 'static_synapse',
                               'delay': self.dt,
                               'weight': -67.}
                     )
        nest.Connect(PY_cells, IN_cells,
                     {'rule': 'pairwise_bernoulli',
                      'multapses': False,
                      'autapses': False,
                      'p': 0.02*scale},
                     syn_spec={'model': 'static_synapse',
                               'delay': self.dt,
                               'weight': 6.}
                     )
        nest.Connect(PY_cells, PY_cells,
                     {'rule': 'pairwise_bernoulli',
                      'multapses': False,
                      'autapses': False,
                      'p': 0.02*scale},
                     syn_spec={'model': 'static_synapse',
                               'delay': self.dt,
                               'weight': 6.}
                     )

        stim_cells = random.sample(PY_cells, int(len(PY_cells)/5))
        for stim_cell in stim_cells:
            stim = nest.Create('poisson_generator', 20,
                               {'rate': 1000./70.,
                                'start': 0., 'stop': 50.})
            nest.Connect(stim, [stim_cell],
                         syn_spec={'model': 'static_synapse',
                                   'delay': self.dt, 'weight': 6.*scale}
                         )

        detector = nest.Create('spike_detector', params={
                                   'to_file': True,
                                   'to_memory': False,
                                }
                               )
        nest.Connect(IN_cells, detector)
        nest.Connect(PY_cells, detector)

        print("PY -> IN: {} synapses per neuron".format(
            len(nest.GetConnections(source=PY_cells,
                                    target=IN_cells))/num_IN))
        print("PY -> PY: {} synapses per neuron".format(
            len(nest.GetConnections(source=PY_cells,
                                    target=PY_cells))/num_PY))
        print("From IN -> IN: {} synapses per neuron".format(
            len(nest.GetConnections(
                source=IN_cells, target=IN_cells))/num_IN))
        print("From IN -> PY: {} synapses per neuron".format(
            len(nest.GetConnections(
                source=IN_cells, target=PY_cells))/num_PY))

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

        nest.Connect(voltmeter1, [IN_cells[0]])
        nest.Connect(voltmeter2, [IN_cells[1]])
        nest.Connect(voltmeter3, [PY_cells[0]])
        nest.Connect(voltmeter4, [PY_cells[1]])

        nest.Simulate(500.)

        # plot the graph
        args = ['gnuplot',
                '-e', 'outputfile="{}"'.format(outputfile),
                '-e', 'sd="{}"'.format(detector[0]),
                '-e', 'v1="{}"'.format(voltmeter1[0]),
                '-e', 'v2="{}"'.format(voltmeter2[0]),
                '-e', 'v3="{}"'.format(voltmeter3[0]),
                '-e', 'v4="{}"'.format(voltmeter4[0]),
                'cortical6.plt'
                ]
        subprocess.call(args)

    def cortical7(self, num_neurons=2000):
        """
        Figure 7.
        """
        self.__setup()
        num_PY = int(num_neurons * 0.8)
        num_IN = int(num_neurons * 0.2)
        num_PY_RS = int(0.95 * num_PY)
        num_PY_LTS = int(0.05 * num_PY)
        PY_cells_RS = nest.Create('RS_weak_cell', num_PY_RS)
        PY_cells_LTS = nest.Create('LTS_cell', num_PY_LTS)
        PY_cells = PY_cells_RS + PY_cells_LTS
        IN_cells = nest.Create('FS_cell', num_IN)

        scale = 2000/num_neurons

        nest.Connect(IN_cells, PY_cells,
                     {'rule': 'pairwise_bernoulli',
                      'multapses': False,
                      'autapses': False,
                      'p': 0.02*scale},
                     syn_spec={'model': 'static_synapse',
                               'delay': self.dt,
                               'weight': -67.}
                     )
        nest.Connect(IN_cells, IN_cells,
                     {'rule': 'pairwise_bernoulli',
                      'multapses': False,
                      'autapses': False,
                      'p': 0.02*scale},
                     syn_spec={'model': 'static_synapse',
                               'delay': self.dt,
                               'weight': -67.}
                     )
        nest.Connect(PY_cells, IN_cells,
                     {'rule': 'pairwise_bernoulli',
                      'multapses': False,
                      'autapses': False,
                      'p': 0.02*scale},
                     syn_spec={'model': 'static_synapse',
                               'delay': self.dt,
                               'weight': 6.}
                     )
        nest.Connect(PY_cells, PY_cells,
                     {'rule': 'pairwise_bernoulli',
                      'multapses': False,
                      'autapses': False,
                      'p': 0.02*scale},
                     syn_spec={'model': 'static_synapse',
                               'delay': self.dt,
                               'weight': 6.}
                     )

        stim_cells = random.sample(PY_cells, int(len(PY_cells)/5))
        for stim_cell in stim_cells:
            stim = nest.Create('poisson_generator', 20,
                               {'rate': 1000./70.,
                                'start': 0., 'stop': 50.})
            nest.Connect(stim, [stim_cell],
                         syn_spec={'model': 'static_synapse',
                                   'delay': self.dt, 'weight': 6.*scale}
                         )

        detector = nest.Create('spike_detector', params={
                                   'to_file': True,
                                   'to_memory': False,
                                }
                               )
        nest.Connect(IN_cells, detector)
        nest.Connect(PY_cells, detector)

        print("PY -> IN: {} synapses per neuron".format(
            len(nest.GetConnections(source=PY_cells,
                                    target=IN_cells))/num_IN))
        print("PY -> PY: {} synapses per neuron".format(
            len(nest.GetConnections(source=PY_cells,
                                    target=PY_cells))/num_PY))
        print("From IN -> IN: {} synapses per neuron".format(
            len(nest.GetConnections(
                source=IN_cells, target=IN_cells))/num_IN))
        print("From IN -> PY: {} synapses per neuron".format(
            len(nest.GetConnections(
                source=IN_cells, target=PY_cells))/num_PY))

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

        nest.Connect(voltmeter1, [PY_cells_RS[0]])
        nest.Connect(voltmeter2, [PY_cells_RS[20]])
        nest.Connect(voltmeter3, [PY_cells_LTS[0]])
        nest.Connect(voltmeter4, [PY_cells_LTS[10]])

        nest.Simulate(5000.)

        # plot the graph
        args = ['gnuplot',
                '-e', 'outputfile="{}"'.format("Figure7b"),
                '-e', 'sd="{}"'.format(detector[0]),
                '-e', 'v1="{}"'.format(voltmeter1[0]),
                '-e', 'v2="{}"'.format(voltmeter2[0]),
                '-e', 'v3="{}'.format(voltmeter3[0]),
                '-e', 'v4="{}"'.format(voltmeter4[0]),
                'cortical7.plt'
                ]
        subprocess.call(args)


if __name__ == "__main__":
    sim = Destexhe2009()
    sim.figure1()
    #  sim.figure2(isi=i)
    #  sim.thalamic(num_TC=2, num_RE=2, outputfile="Figure3")
    #  sim.thalamic(num_TC=20, num_RE=20, outputfile="Figure3")
    #  sim.thalamic(num_TC=30, num_RE=30, outputfile="Figure4b")
    #  sim.thalamic(num_TC=50, num_RE=50, outputfile="Figure4c")
    sim.corticals6(num_neurons=2000)
    #  sim.cortical7(num_neurons=500)
