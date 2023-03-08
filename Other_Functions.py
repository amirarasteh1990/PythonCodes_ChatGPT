import numpy as np
import math
import cmath
from components import AcCable


class CableSegmental:
    def __init__(self, project_data, optimisation_data, df_components, operating_point):
        self.project_data = project_data
        self.df_c0 = df_components['df c0']
        self.df_c1 = df_components['df c1']
        self.df_c2 = df_components['df c2']
        self.df_c3 = df_components['df c3']
        self.temp_c0 = operating_point['T_C0 (oC)']
        self.temp_c1 = operating_point['T_C1 (oC)']
        self.temp_c2 = operating_point['T_C2 (oC)']
        self.temp_c3 = operating_point['T_C3 (oC)']
        self.n_parallel = optimisation_data['n_parallel']
        self.d_start_c0 = 0
        self.d_c0 = project_data['d_C0 (km)']
        self.d_start_c1 = self.d_c0
        self.d_c1 = project_data['d_C1 (km)']
        self.d_start_c2 = self.d_c0 + self.d_c1
        self.d_c2 = project_data['d_C2 (km)']
        self.d_start_c3 = self.d_c0 + self.d_c1 + self.d_c2
        self.d_c3 = project_data['d_C3 (km)']
        self.d_comp_mid = project_data['d_Comp (km) from offshore']

    def vi(self, v_onshore, i_onshore):
        project_data = self.project_data
        df_c0 = self.df_c0
        df_c1 = self.df_c1
        df_c2 = self.df_c2
        df_c3 = self.df_c3
        temp_c0 = self.temp_c0
        temp_c1 = self.temp_c1
        temp_c2 = self.temp_c2
        temp_c3 = self.temp_c3
        n_parallel = self.n_parallel
        d_start_c0 = self.d_start_c0
        d_c0 = self.d_c0
        d_start_c1 = self.d_start_c1
        d_c1 = self.d_c1
        d_start_c2 = self.d_start_c2
        d_c2 = self.d_c2
        d_start_c3 = self.d_start_c3
        d_c3 = self.d_c3
        d_comp_mid = self.d_comp_mid

        v_onshore = v_onshore / math.sqrt(3)
        i_onshore = i_onshore

        c0_analyse = AcCable(project_data, n_parallel, df_c0, d_start_c0, d_c0, d_comp_mid, temp_c0)

        c1_analyse = AcCable(project_data, n_parallel, df_c1, d_start_c1, d_c1, d_comp_mid, temp_c1)

        c2_analyse = AcCable(project_data, n_parallel, df_c2, d_start_c2, d_c2, d_comp_mid, temp_c2)

        c3_analyse = AcCable(project_data, n_parallel, df_c3, d_start_c3, d_c3, d_comp_mid, temp_c3)

        vi_offshore_c3 = c3_analyse.inputs_vionshore_outputs_vioffshore(v_onshore, i_onshore)

        vi_offshore_c2 = c2_analyse.inputs_vionshore_outputs_vioffshore(vi_offshore_c3['v_offshore'],
                                                                        vi_offshore_c3['i_offshore'])
        vi_offshore_c1 = c1_analyse.inputs_vionshore_outputs_vioffshore(vi_offshore_c2['v_offshore'],
                                                                        vi_offshore_c2['i_offshore'])
        vi_offshore_c0 = c0_analyse.inputs_vionshore_outputs_vioffshore(vi_offshore_c1['v_offshore'],
                                                                        vi_offshore_c1['i_offshore'])
        i_offshore_c3 = vi_offshore_c3['i_offshore']
        v_offshore_c3 = vi_offshore_c3['v_offshore']
        i_offshore_c2 = vi_offshore_c2['i_offshore']
        v_offshore_c2 = vi_offshore_c2['v_offshore']
        i_offshore_c1 = vi_offshore_c1['i_offshore']
        v_offshore_c1 = vi_offshore_c1['v_offshore']
        i_offshore_c0 = vi_offshore_c0['i_offshore']
        v_offshore_c0 = vi_offshore_c0['v_offshore']

        ########################################################
        # onshore compensation (check if i_onshore of the cables is in the limit)
        if ((d_c0 > 0 and abs(vi_offshore_c0['i_offshore']) > df_c0['Ampacity 20DegC (kA)'] * n_parallel) or
                (d_c0 > 0 and abs(vi_offshore_c1['i_offshore']) > df_c0['Ampacity 20DegC (kA)'] * n_parallel) or
                (d_c1 > 0 and abs(vi_offshore_c1['i_offshore']) > df_c1['Ampacity 20DegC (kA)'] * n_parallel) or
                (d_c1 > 0 and abs(vi_offshore_c2['i_offshore']) > df_c1['Ampacity 20DegC (kA)'] * n_parallel) or
                (d_c2 > 0 and abs(vi_offshore_c2['i_offshore']) > df_c2['Ampacity 20DegC (kA)'] * n_parallel) or
                (d_c2 > 0 and abs(vi_offshore_c3['i_offshore']) > df_c2['Ampacity 20DegC (kA)'] * n_parallel) or
                (d_c3 > 0 and abs(vi_offshore_c3['i_offshore']) > df_c3['Ampacity 20DegC (kA)'] * n_parallel) or
                (d_c3 > 0 and abs(i_onshore) > df_c3['Ampacity 20DegC (kA)'] * n_parallel)):
            required_compensation = 'Onshore'
            # Shifting and scaling input voltage and current
            v_onshore_phase = cmath.phase(v_onshore)
            v_onshore = abs(v_onshore)
            i_onshore = i_onshore * cmath.exp(-v_onshore_phase * 1j)
            # ip_onshore = i_onshore.real
            # iq_onshore = i_onshore.imag
            #
            # i_onshore = ip_onshore - iq_onshore * 1j
            vi_offshore_c3 = c3_analyse.inputs_vionshore_outputs_vioffshore(v_onshore, i_onshore)
            vi_offshore_c2 = c2_analyse.inputs_vionshore_outputs_vioffshore(vi_offshore_c3['v_offshore'], vi_offshore_c3['i_offshore'])
            vi_offshore_c1 = c1_analyse.inputs_vionshore_outputs_vioffshore(vi_offshore_c2['v_offshore'], vi_offshore_c2['i_offshore'])
            vi_offshore_c0 = c0_analyse.inputs_vionshore_outputs_vioffshore(vi_offshore_c1['v_offshore'], vi_offshore_c1['i_offshore'])
            i_offshore = vi_offshore_c0['i_offshore']

            i_max_nc3 = df_c3['Ampacity 20DegC (kA)'] * n_parallel

            if d_c3 > 0 and abs(i_onshore.real) > i_max_nc3:
                return {'comments': 'Error: ip > ampacity_c3 due to voltage drop (possible solution: tap changer of onshore transformer)'}
            elif d_c3 > 0 and abs(v_onshore) > df_c3['Insulation Voltage (kV)']:
                return {'comments': 'Error: V > V_c3 due to voltage drops (possible solution: tap changer of onshore transformer)'}
            elif d_c3 > 0 and abs(i_onshore) > i_max_nc3:
                i_onshore = i_onshore.real + 1j * np.sign(i_onshore.imag) * math.sqrt(i_max_nc3 ** 2 - i_onshore.real ** 2)
                vi_offshore_c3 = c3_analyse.inputs_vionshore_outputs_vioffshore(v_onshore, i_onshore)
                vi_offshore_c2 = c2_analyse.inputs_vionshore_outputs_vioffshore(vi_offshore_c3['v_offshore'], vi_offshore_c3['i_offshore'])
                vi_offshore_c1 = c1_analyse.inputs_vionshore_outputs_vioffshore(vi_offshore_c2['v_offshore'], vi_offshore_c2['i_offshore'])
                vi_offshore_c0 = c0_analyse.inputs_vionshore_outputs_vioffshore(vi_offshore_c1['v_offshore'], vi_offshore_c1['i_offshore'])
                i_offshore = vi_offshore_c0['i_offshore']

            while (((d_c0 > 0 and abs(vi_offshore_c0['i_offshore']) > df_c0['Ampacity 20DegC (kA)'] * n_parallel) or
                   (d_c1 > 0 and abs(vi_offshore_c1['i_offshore']) > df_c1['Ampacity 20DegC (kA)'] * n_parallel) or
                   (d_c2 > 0 and abs(vi_offshore_c2['i_offshore']) > df_c2['Ampacity 20DegC (kA)'] * n_parallel) or
                   (d_c3 > 0 and abs(vi_offshore_c3['i_offshore']) > df_c3['Ampacity 20DegC (kA)'] * n_parallel))): # and (d_c3 > 0 and abs(i_onshore) + 0.001 < df_c3['Ampacity 20DegC (kA)'] * n_parallel)):
                i_onshore -= 0.001j
                vi_offshore_c3 = c3_analyse.inputs_vionshore_outputs_vioffshore(v_onshore, i_onshore)
                vi_offshore_c2 = c2_analyse.inputs_vionshore_outputs_vioffshore(vi_offshore_c3['v_offshore'], vi_offshore_c3['i_offshore'])
                vi_offshore_c1 = c1_analyse.inputs_vionshore_outputs_vioffshore(vi_offshore_c2['v_offshore'], vi_offshore_c2['i_offshore'])
                vi_offshore_c0 = c0_analyse.inputs_vionshore_outputs_vioffshore(vi_offshore_c1['v_offshore'], vi_offshore_c1['i_offshore'])
                i_offshore = vi_offshore_c0['i_offshore']

            if (((d_c0 > 0 and abs(vi_offshore_c0['i_offshore']) > df_c0['Ampacity 20DegC (kA)'] * n_parallel) or
                 (d_c0 > 0 and abs(vi_offshore_c1['i_offshore']) > df_c0['Ampacity 20DegC (kA)'] * n_parallel)) or
                 ((d_c1 > 0 and abs(vi_offshore_c1['i_offshore']) > df_c1['Ampacity 20DegC (kA)'] * n_parallel) or
                 (d_c1 > 0 and abs(vi_offshore_c2['i_offshore']) > df_c1['Ampacity 20DegC (kA)'] * n_parallel)) or
                 ((d_c2 > 0 and abs(vi_offshore_c2['i_offshore']) > df_c2['Ampacity 20DegC (kA)'] * n_parallel) or
                 (d_c2 > 0 and abs(vi_offshore_c3['i_offshore']) > df_c2['Ampacity 20DegC (kA)'] * n_parallel)) or
                 ((d_c3 > 0 and abs(vi_offshore_c3['i_offshore']) > df_c3['Ampacity 20DegC (kA)'] * n_parallel) or
                 (d_c3 > 0 and abs(i_onshore) > df_c3['Ampacity 20DegC (kA)'] * n_parallel))):
                required_compensation = 'Mid point'

            if required_compensation == 'Mid point' and d_comp_mid == 0:
                return {'comments': 'Error: midpoint compensation is required'}
            elif required_compensation == 'Mid point':
                return {'comments': 'Error: midpoint compensation is required'}

            # outputs = {'i_onshore_cable': i_onshore * cmath.exp(v_onshore_phase * 1j),
            #            'comments': 'onshore_compensation is required'}
            # return outputs

            i_offshore_c3 = vi_offshore_c3['i_offshore']
            v_offshore_c3 = vi_offshore_c3['v_offshore']
            i_offshore_c2 = vi_offshore_c2['i_offshore']
            v_offshore_c2 = vi_offshore_c2['v_offshore']
            i_offshore_c1 = vi_offshore_c1['i_offshore']
            v_offshore_c1 = vi_offshore_c1['v_offshore']
            i_offshore_c0 = vi_offshore_c0['i_offshore']
            v_offshore_c0 = vi_offshore_c0['v_offshore']

            # Shifting back the input voltage and current

            v_onshore = v_onshore * cmath.exp(v_onshore_phase * 1j)
            i_onshore = i_onshore * cmath.exp(v_onshore_phase * 1j)

            v_offshore_c3 = v_offshore_c3 * cmath.exp(v_onshore_phase * 1j)
            i_offshore_c3 = i_offshore_c3 * cmath.exp(v_onshore_phase * 1j)

            v_offshore_c2 = v_offshore_c2 * cmath.exp(v_onshore_phase * 1j)
            i_offshore_c2 = i_offshore_c2 * cmath.exp(v_onshore_phase * 1j)

            v_offshore_c1 = v_offshore_c1 * cmath.exp(v_onshore_phase * 1j)
            i_offshore_c1 = i_offshore_c1 * cmath.exp(v_onshore_phase * 1j)

            v_offshore_c0 = v_offshore_c0 * cmath.exp(v_onshore_phase * 1j)
            i_offshore_c0 = i_offshore_c0 * cmath.exp(v_onshore_phase * 1j)

        #####################################################################

        # No load analysis (with maximum voltage offshore and zero current offshore)
        c0_vionshore_nl = c0_analyse.inputs_vioffshore_outputs_vionshore(df_c0['Insulation Voltage (kV)'] / math.sqrt(3), 0)
        if abs(c0_vionshore_nl['v_onshore']) > df_c0['Insulation Voltage (kV)'] / math.sqrt(3):
            return {'comments': 'Error: c0 Insulation voltage is exceeded'}

        if d_c0 == 0:
            c1_vioffshore_nl = {'v_offshore': df_c1['Insulation Voltage (kV)'] / math.sqrt(3),
                                'i_offshore': 0
                                }
        else:
            c1_vioffshore_nl = {'v_offshore': c0_vionshore_nl['v_onshore'],
                             'i_offshore': c0_vionshore_nl['i_onshore']
                             }

        c1_vionshore_nl = c1_analyse.inputs_vioffshore_outputs_vionshore(c1_vioffshore_nl['v_offshore'], c1_vioffshore_nl['i_offshore'])

        if abs(c1_vionshore_nl['v_onshore']) > df_c1['Insulation Voltage (kV)'] / math.sqrt(3):
            return {'comments': 'Error: c1 Insulation voltage is exceeded'}

        if d_c0 == 0 and d_c1 == 0:
            c2_vioffshore_nl = {'v_offshore': df_c2['Insulation Voltage (kV)'] / math.sqrt(3),
                                'i_offshore': 0}
        else:
            c2_vioffshore_nl = {'v_offshore': c1_vionshore_nl['v_onshore'],
                                'i_offshore': c1_vionshore_nl['i_onshore']
                                }

        c2_vionshore_nl = c2_analyse.inputs_vioffshore_outputs_vionshore(c2_vioffshore_nl['v_offshore'], c2_vioffshore_nl['i_offshore'])

        if abs(c2_vionshore_nl['v_onshore']) > df_c2['Insulation Voltage (kV)'] / math.sqrt(3):
            return {'comments': 'Error: c2 Insulation voltage is exceeded'}

        if d_c0 == 0 and d_c1 == 0 and d_c2 == 0:
            c3_vioffshore_nl = {'v_offshore': df_c3['Insulation Voltage (kV)'] / math.sqrt(3),
                                'i_offshore': 0}
        else:
            c3_vioffshore_nl = {'v_offshore': c2_vionshore_nl['v_onshore'],
                                'i_offshore': c2_vionshore_nl['i_onshore']}

        c3_vionshore_nl = c3_analyse.inputs_vioffshore_outputs_vionshore(c3_vioffshore_nl['v_offshore'], c3_vioffshore_nl['i_offshore'])

        if abs(c3_vionshore_nl['v_onshore']) > df_c3['Insulation Voltage (kV)'] / math.sqrt(3):
            return {'comments': 'Error: c3 Insulation voltage is exceeded'}

        v_nl_onshore = c3_vionshore_nl['v_onshore']
        i_nl_onshore = c3_vionshore_nl['i_onshore']

        # scaling voltages and currents
        v_onshore = v_onshore * math.sqrt(3)
        v_nl_onshore = v_nl_onshore * math.sqrt(3)
        v_offshore_c3 = v_offshore_c3 * math.sqrt(3)
        v_offshore_c2 = v_offshore_c2 * math.sqrt(3)
        v_offshore_c1 = v_offshore_c1 * math.sqrt(3)
        v_offshore_c0 = v_offshore_c0 * math.sqrt(3)
        # v_nl_onshore = abs(v_nl_onshore) # already converter to ll

        outputs = {'v_onshore': v_onshore,
                   'i_onshore': i_onshore,
                   'v_offshore_c3': v_offshore_c3,
                   'i_offshore_c3': i_offshore_c3,
                   'v_offshore_c2': v_offshore_c2,
                   'i_offshore_c2': i_offshore_c2,
                   'v_offshore_c1': v_offshore_c1,
                   'i_offshore_c1': i_offshore_c1,
                   'v_offshore': v_offshore_c0,
                   'i_offshore': i_offshore_c0,
                   'v_nl_onshore': v_nl_onshore,
                   'i_nl_onshore': i_nl_onshore,
                   'q_comp_mid': 0,
                   'comments': ''}

        return outputs

    def p_limit(self, outputs_vi):
        project_data = self.project_data
        df_c0 = self.df_c0
        df_c1 = self.df_c1
        df_c2 = self.df_c2
        df_c3 = self.df_c3
        temp_c0 = self.temp_c0
        temp_c1 = self.temp_c1
        temp_c2 = self.temp_c2
        temp_c3 = self.temp_c3
        n_parallel = self.n_parallel
        d_start_c0 = self.d_start_c0
        d_c0 = self.d_c0
        d_start_c1 = self.d_start_c1
        d_c1 = self.d_c1
        d_start_c2 = self.d_start_c2
        d_c2 = self.d_c2
        d_start_c3 = self.d_start_c3
        d_c3 = self.d_c3
        d_comp_mid = self.d_comp_mid
        v_onshore_c3 = outputs_vi['v_onshore']
        v_offshore_c3 = outputs_vi['v_offshore_c3']
        i_offshore_c3 = outputs_vi['i_offshore_c3']
        v_onshore_c2 = v_offshore_c3
        v_offshore_c2 = outputs_vi['v_offshore_c2']
        i_offshore_c2 = outputs_vi['i_offshore_c2']
        v_onshore_c1 = v_offshore_c2
        v_offshore_c1 = outputs_vi['v_offshore_c1']
        i_offshore_c1 = outputs_vi['i_offshore_c1']
        v_onshore_c0 = v_offshore_c1
        v_offshore_c0 = outputs_vi['v_offshore']
        i_offshore_c0 = outputs_vi['i_offshore']

        s_c0 = math.sqrt(3) * v_offshore_c0 * np.conj(i_offshore_c0)
        p_c0 = s_c0.real
        s_c1 = math.sqrt(3) * v_offshore_c1 * np.conj(i_offshore_c1)
        p_c1 = s_c1.real
        s_c2 = math.sqrt(3) * v_offshore_c2 * np.conj(i_offshore_c2)
        p_c2 = s_c2.real
        s_c3 = math.sqrt(3) * v_offshore_c3 * np.conj(i_offshore_c3)
        p_c3 = s_c3.real

        c0 = AcCable(project_data, n_parallel, df_c0, d_start_c0, d_c0, d_comp_mid, temp_c0)
        c1 = AcCable(project_data, n_parallel, df_c1, d_start_c1, d_c1, d_comp_mid, temp_c1)
        c2 = AcCable(project_data, n_parallel, df_c2, d_start_c2, d_c2, d_comp_mid, temp_c2)
        c3 = AcCable(project_data, n_parallel, df_c3, d_start_c3, d_c3, d_comp_mid, temp_c3)

        p_max_c0 = c0.p_max(v_onshore_c0, v_offshore_c0)
        p_max_c1 = c1.p_max(v_onshore_c1, v_offshore_c1)
        p_max_c2 = c2.p_max(v_onshore_c2, v_offshore_c2)
        p_max_c3 = c3.p_max(v_onshore_c3, v_offshore_c3)

        if p_max_c0 < p_c0 or p_max_c1 < p_c1 or p_max_c2 < p_c2 or p_max_c3 < p_c3:
            return True
        else:
            return False


def pq_cal(vi_nodes):
    #
    v_offshore = vi_nodes['v_offshore']
    i_offshore_wf = vi_nodes['i_offshore_wf']
    i_offshore_comp_wfside = vi_nodes['i_offshore_comp_wfside']
    i_offshore_tr_wfside = vi_nodes['i_offshore_tr_wfside']

    i_offshore_tr_cableside = vi_nodes['i_offshore_tr_cableside']

    i_offshore_comp_cableside = vi_nodes['i_offshore_comp_cableside']
    v_offshore_cable = vi_nodes['v_offshore_cable']
    i_offshore_cable = vi_nodes['i_offshore_cable']
    v_offshore_c3 = vi_nodes['v_offshore_c3']
    i_offshore_c3 = vi_nodes['i_offshore_c3']
    v_offshore_c2 = vi_nodes['v_offshore_c2']
    i_offshore_c2 = vi_nodes['i_offshore_c2']
    v_offshore_c1 = vi_nodes['v_offshore_c1']
    i_offshore_c1 = vi_nodes['i_offshore_c1']
    v_onshore_cable = vi_nodes['v_onshore_cable']
    i_onshore_cable = vi_nodes['i_onshore_cable']
    i_onshore_comp_cableside = vi_nodes['i_onshore_comp_cableside']
    i_onshore_tr_cableside = vi_nodes['i_onshore_tr_cableside']

    i_onshore_tr_gridside = vi_nodes['i_onshore_tr_gridside']
    i_onshore_comp_gridside = vi_nodes['i_onshore_comp_gridside']

    v_poc = vi_nodes['v_poc']
    i_poc = vi_nodes['i_poc']

    v_grid = vi_nodes['v_grid']
    i_grid = vi_nodes['i_grid']
    ########################
    # P and Q of the nodes #
    ########################
    # offshore_wf
    s_offshore_wf = math.sqrt(3) * v_offshore * np.conj(i_offshore_wf)

    # offshore_comp_wfside
    s_offshore_comp_wfside = math.sqrt(3) * v_offshore * np.conj(i_offshore_comp_wfside)

    # offshore_tr_wfside
    s_offshore_tr_wfside = math.sqrt(3) * v_offshore * np.conj(i_offshore_tr_wfside)

    # offshore_tr_cableside
    s_offshore_tr_cableside = math.sqrt(3) * v_offshore_cable * np.conj(i_offshore_tr_cableside)

    # offshore_comp_cableside
    s_offshore_comp_cableside = math.sqrt(3) * v_offshore_cable * np.conj(i_offshore_comp_cableside)

    # offshore_cable
    s_offshore_cable = math.sqrt(3) * v_offshore_cable * np.conj(i_offshore_cable)

    # c2_c3
    s_c2c3 = math.sqrt(3) * v_offshore_c3 * np.conj(i_offshore_c3)

    # c1_c2
    s_c1c2 = math.sqrt(3) * v_offshore_c2 * np.conj(i_offshore_c2)

    # c0_c1
    s_c0c1 = math.sqrt(3) * v_offshore_c1 * np.conj(i_offshore_c1)

    # onshore_cable
    s_onshore_cable = math.sqrt(3) * v_onshore_cable * np.conj(i_onshore_cable)

    # onshore_comp_cableside
    s_onshore_comp_cableside = math.sqrt(3) * v_onshore_cable * np.conj(i_onshore_comp_cableside)

    # onshore_tr_cableside
    s_onshore_tr_cableside = math.sqrt(3) * v_onshore_cable * np.conj(i_onshore_tr_cableside)

    # onshore_tr_gridside
    s_onshore_tr_gridside = math.sqrt(3) * v_poc * np.conj(i_onshore_tr_gridside)

    # onshore_comp_gridside
    s_onshore_comp_gridside = math.sqrt(3) * v_poc * np.conj(i_onshore_comp_gridside)

    # poc
    s_poc = math.sqrt(3) * v_poc * np.conj(i_poc)

    # grid
    s_grid = math.sqrt(3) * v_grid * np.conj(i_grid)

    # outputs
    outputs = {'s_offshore_wf': s_offshore_wf,
               's_offshore_comp_wfside': s_offshore_comp_wfside,
               's_offshore_tr_wfside': s_offshore_tr_wfside,
               's_offshore_tr_cableside': s_offshore_tr_cableside,
               's_offshore_comp_cableside': s_offshore_comp_cableside,
               's_offshore_cable': s_offshore_cable,
               's_c0c1': s_c0c1,
               's_c1c2': s_c1c2,
               's_c2c3': s_c2c3,
               's_onshore_cable': s_onshore_cable,
               's_onshore_comp_cableside': s_onshore_comp_cableside,
               's_onshore_tr_cableside': s_onshore_tr_cableside,
               's_onshore_tr_gridside': s_onshore_tr_gridside,
               's_onshore_comp_gridside': s_onshore_comp_gridside,
               's_poc': s_poc,
               's_grid': s_grid
               }

    return outputs
