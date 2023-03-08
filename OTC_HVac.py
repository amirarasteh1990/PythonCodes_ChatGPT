import numpy as np
import cmath
import math
from components import AcTransformer, AcGrid
from Other_Functions import CableSegmental


def hvac(project_data, optimisation_data, df_components, operating_point):

    ################
    # project data #
    ################
    # q_project = project_data['s']
    # p_project = project_data['p']

    ##########################
    # data of the components #
    ##########################
    df_offshore_tr = df_components['df offshore transformer']
    df_onshore_tr = df_components['df onshore transformer']

    tap_offshore_tr_pos = 0
    tap_offshore_tr_neg = 0
    tap_onshore_tr_pos = 0
    tap_onshore_tr_neg = 0

    #####################
    # optimisation data #
    #####################
    p_poc = operating_point['P_poc (MW)']
    q_poc = operating_point['Q_poc (MVAR)']
    v_poc_pu = operating_point['V_poc (pu)']
    v_offshore = optimisation_data['V_wf (kV)']
    v_export = optimisation_data['V_export (kV)']
    v_onshore = optimisation_data['V_grid (kV)']
    comp_offshore_location = optimisation_data['Offshore comp']
    comp_midpoint_inout = optimisation_data['Midpoint comp']
    comp_onshore_location = optimisation_data['Onshore comp']
    dv_poc_pos = optimisation_data['DVmaxOnshore positive (%)']
    dv_poc_neg = optimisation_data['DVmaxOnshore negative (%)']
    q_poc_min = optimisation_data['Q_grid min (MVAR)']
    q_poc_max = optimisation_data['Q_grid max (MVAR)']
    n_tr_offshore = optimisation_data['n_parallel']
    n_tr_onshore = optimisation_data['n_parallel']

    ###########
    # AC Grid #
    ###########
    #v_grid = v_onshore *
    ac_grid = AcGrid(project_data, optimisation_data)
    # p_poc = p_project
    # q_poc = ac_grid.q_poc(q_project)

    #vi_poc = ac_grid.vi_poc2(p_poc, q_poc, v_poc_pu)
    v_grid = optimisation_data['V_grid (kV)'] #vi_poc['v_grid']
    #i_grid = vi_poc['i_grid']
    v_poc = v_grid * v_poc_pu#vi_poc['v_poc']
    #i_poc = vi_poc['i_poc']
    i_poc = np.conj((p_poc + 1j * q_poc) / (math.sqrt(3) * v_poc))
    i_grid = i_poc
    #i_grid = i_poc
    i_onshore_tr_gridside = i_poc

    tr_onshore = AcTransformer(project_data, df_onshore_tr, n_tr_onshore)

    if v_onshore > v_export:
        inputside_tr_onshore = 'high voltage'
    else:
        inputside_tr_onshore = 'low voltage'

    if project_data['Tap changer'] == 'Automatic Tap Changing':
        tap_tr_onshore = tr_onshore.tap_changer(v_poc, -i_onshore_tr_gridside, inputside_tr_onshore, v_export)
    else:
        tap_tr_onshore = 0

    vi_onshore_tr_cableside = tr_onshore.vi_output(v_poc, -i_onshore_tr_gridside, inputside_tr_onshore, tap_tr_onshore)
    v_onshore_tr_cableside = vi_onshore_tr_cableside['v_output']
    i_onshore_tr_cableside = -vi_onshore_tr_cableside['i_output']

    v_onshore_cable = v_onshore_tr_cableside
    i_onshore_cable = i_onshore_tr_cableside

    ##########################################
    # Voltages and currents across the cable #
    ##########################################
    export_cable = CableSegmental(project_data, optimisation_data, df_components, operating_point)
    outputs_exportcable = export_cable.vi(v_onshore_cable, i_onshore_cable)

    if 'Error' in outputs_exportcable['comments']:
        return {'comments': outputs_exportcable['comments']}

    i_onshore_cable = outputs_exportcable['i_onshore']

    #if comp_onshore_location == 'Export cable side':
    i_onshore_comp_cableside = i_onshore_tr_cableside - i_onshore_cable
    i_onshore_comp_gridside = 0

    if comp_onshore_location == 'Grid side' and i_onshore_comp_cableside != 0:
        q_comp_onshore = math.sqrt(3) * v_onshore_cable * np.conj(i_onshore_comp_cableside)
        i_onshore_comp_gridside = np.conj(q_comp_onshore / (math.sqrt(3) * v_poc))
        i_onshore_tr_gridside_new = i_poc - i_onshore_comp_gridside
        vi_onshore_tr_cableside_new = tr_onshore.vi_output(v_poc, -i_onshore_tr_gridside_new, inputside_tr_onshore, 0)
        i_onshore_tr_cableside_new = -vi_onshore_tr_cableside_new['i_output']
        d_q_onshore = q_comp_onshore
        d_i_cableside = abs(i_onshore_tr_cableside_new - i_onshore_tr_cableside)
        d_i_cableside_2 = abs(i_onshore_tr_cableside_new - i_onshore_cable)
        k_delta = d_i_cableside_2 / d_i_cableside
        d_q_onshore_2 = k_delta * d_q_onshore
        q_comp_onshore = q_comp_onshore - d_q_onshore_2

        i_onshore_comp_gridside = np.conj(q_comp_onshore / (math.sqrt(3) * v_poc))
        i_onshore_tr_gridside = i_poc - i_onshore_comp_gridside
        i_onshore_comp_cableside = 0

        if v_onshore > v_export:
            inputside_tr_onshore = 'high voltage'
        else:
            inputside_tr_onshore = 'low voltage'

        vi_onshore_tr_cableside = tr_onshore.vi_output(v_poc, -i_onshore_tr_gridside, inputside_tr_onshore, 0)
        v_tr_onshore_cableside = vi_onshore_tr_cableside['v_output']
        i_onshore_tr_cableside = -vi_onshore_tr_cableside['i_output']

        # while abs(abs(i_onshore_tr_cableside) - abs(i_onshore_cable)) > 0.001:
        #     q_comp_onshore -= np.sign(abs(i_onshore_tr_cableside) - abs(i_onshore_cable)) * 0.01j
        #     i_onshore_comp_gridside = np.conj(q_comp_onshore / (math.sqrt(3) * v_poc))
        #     i_onshore_tr_gridside = i_poc - i_onshore_comp_gridside
        #     i_onshore_comp_cableside = 0
        #
        #     if v_onshore > v_export:
        #         inputside_tr_onshore = 'high voltage'
        #     else:
        #         inputside_tr_onshore = 'low voltage'
        #
        #     vi_onshore_tr_cableside = tr_onshore.vi_output(v_poc, -i_onshore_tr_gridside, inputside_tr_onshore, 0)
        #     v_tr_onshore_cableside = vi_onshore_tr_cableside['v_output']
        #     i_onshore_tr_cableside = -vi_onshore_tr_cableside['i_output']

        v_onshore_cable = v_tr_onshore_cableside
        i_onshore_cable = i_onshore_tr_cableside
        outputs_exportcable = export_cable.vi(v_onshore_cable, i_onshore_cable)

    try:
        v_offshore_cable = outputs_exportcable['v_offshore']
        i_offshore_cable = outputs_exportcable['i_offshore']
        v_onshore_cable = outputs_exportcable['v_onshore']
        i_onshore_cable = outputs_exportcable['i_onshore']
        q_comp_mid = outputs_exportcable['q_comp_mid']
    except TypeError and KeyError:
        return outputs_exportcable

    s_onshore_tr_gridside = math.sqrt(3) * v_poc * np.conj(i_onshore_tr_gridside)
    s_onshore_tr_cableside = math.sqrt(3) * v_onshore_tr_cableside * np.conj(i_onshore_tr_cableside)
    s_onshore_tr_abs = max(abs(s_onshore_tr_gridside), abs(s_onshore_tr_cableside))
    if s_onshore_tr_abs > df_onshore_tr['Capacity (MVA)'] * n_tr_onshore + 10:
        return {'comments': "Error: s_onshore_tr > df_onshore_tr[Capacity (MVA)] * n_tr_onshore"}

    ####
    if export_cable.p_limit(outputs_exportcable):
        return {'comments': 'Error: Maximum active power transmission limit is exceeded'}
    ############
    # OFFSHORE #
    ############
    tr_offshore = AcTransformer(project_data, df_offshore_tr, n_tr_offshore)

    if v_export > v_offshore:
        inputside_tr_offshore = 'high voltage'
    else:
        inputside_tr_offshore = 'low voltage'

    v_offshore_tr_cableside = v_offshore_cable
    i_offshore_tr_cableside = i_offshore_cable
    vi_tr_offshore_wfside = tr_offshore.vi_output(v_offshore_cable, -i_offshore_tr_cableside, inputside_tr_offshore, 0)
    v_offshore_tr_wfside = vi_tr_offshore_wfside['v_output']
    i_offshore_tr_wfside = -vi_tr_offshore_wfside['i_output']

    v_offshore = v_offshore_tr_wfside
    i_offshore_wf = i_offshore_tr_wfside

    s_offshore_wf = math.sqrt(3) * v_offshore_tr_wfside * np.conj(i_offshore_wf)
    if s_offshore_wf.real > project_data['P (MW)'] + 20:
        return {'comments': 'Error: P_wf > P_project'}

    try:
        q_available_wf = math.sqrt(project_data['S (MVA)'] ** 2 - s_offshore_wf.real ** 2) * abs(v_offshore / optimisation_data['V_wf (kV)']) ** 2
    except ValueError:
        q_available_wf = 0

    if s_offshore_wf.imag > q_available_wf:
        q_comp_offshore = (s_offshore_wf.imag - q_available_wf) * 1j
    elif s_offshore_wf.imag < -q_available_wf:
        q_comp_offshore = (s_offshore_wf.imag + q_available_wf) * 1j
    else:
        q_comp_offshore = 0

    #if comp_offshore_location == 'Wind Farm side':
    i_offshore_comp_wfside = np.conj(q_comp_offshore / (math.sqrt(3) * v_offshore_tr_wfside))
    i_offshore_comp_cableside = 0

    if comp_offshore_location == 'Export cable side' and i_offshore_comp_wfside != 0:

        s_offshore_tr_wfside = math.sqrt(3) * v_offshore_tr_wfside * np.conj(i_offshore_wf)
        q_offshore_tr_wfside = s_offshore_tr_wfside.imag
        while abs(q_offshore_tr_wfside - np.sign(q_offshore_tr_wfside) * q_available_wf) > 1:
            q_comp_offshore += np.sign(q_offshore_tr_wfside - np.sign(q_offshore_tr_wfside) * q_available_wf) * 0.1j
            i_offshore_comp_cableside = np.conj(q_comp_offshore / (math.sqrt(3) * v_offshore_cable))
            i_offshore_tr_cableside = i_offshore_cable - i_offshore_comp_cableside
            i_offshore_comp_wfside = 0
            if v_export > v_offshore:
                inputside_tr_offshore = 'high voltage'
            else:
                inputside_tr_offshore = 'low voltage'

            vi_offshore_tr_wfside = tr_offshore.vi_output(v_offshore_cable, -i_offshore_tr_cableside, inputside_tr_offshore, 0)
            v_offshore_tr_wfside = vi_offshore_tr_wfside['v_output']
            i_offshore_tr_wfside = -vi_offshore_tr_wfside['i_output']
            v_offshore = v_offshore_tr_wfside

            s_offshore_tr_wfside = math.sqrt(3) * v_offshore_tr_wfside * np.conj(i_offshore_tr_wfside)
            q_offshore_tr_wfside = s_offshore_tr_wfside.imag

    i_offshore_wf = i_offshore_tr_wfside - i_offshore_comp_wfside

    s_offshore_tr_wfside = math.sqrt(3) * v_offshore_tr_wfside * np.conj(i_offshore_tr_wfside)
    s_offshore_tr_cableside = math.sqrt(3) * v_offshore_tr_cableside * np.conj(i_offshore_tr_cableside)
    s_offshore_tr_abs = max(abs(s_offshore_tr_wfside), abs(s_offshore_tr_cableside))
    if s_offshore_tr_abs > df_offshore_tr['Capacity (MVA)'] * n_tr_offshore + 10:
        return {'comments': 'Error: s_offshore_tr > df_offshore_tr[Capacity (MVA)] * n_tr_offshore'}


    ############################################
    # tap changer of the transformer (no load) #
    ############################################
    vi_poc_nl = tr_onshore.vi_output(outputs_exportcable['v_nl_onshore'], outputs_exportcable['i_nl_onshore'], 'high voltage', 0)
    i_poc_nl = vi_poc_nl['i_output']
    v_poc_nl = ac_grid.v_poc(i_poc_nl)
    vi_nl_onshore_2 = tr_onshore.vi_output(v_poc_nl, i_poc_nl, 'low voltage', 0)
    v_nl_onshore_2 = vi_nl_onshore_2['v_output']

    if abs(outputs_exportcable['v_nl_onshore']) < abs(v_nl_onshore_2):
        tap_tr_onshore_nl = tr_onshore.tap_changer(outputs_exportcable['v_nl_onshore'], outputs_exportcable['i_nl_onshore'], 'low voltage', v_poc_nl)
    else:
        tap_tr_onshore_nl = 0

    #########################################
    # Return outputs (voltage and currents) #
    #########################################
    outputs = {'v_offshore': v_offshore_tr_wfside,
               'i_offshore_wf': i_offshore_wf,
               'i_offshore_comp_wfside': i_offshore_comp_wfside,
               'i_offshore_tr_wfside': i_offshore_tr_wfside,
               'i_offshore_tr_cableside': i_offshore_tr_cableside,
               'i_offshore_comp_cableside': i_offshore_comp_cableside,
               'v_offshore_cable': v_offshore_cable,
               'i_offshore_cable': i_offshore_cable,
               'v_onshore_cable': v_onshore_cable,
               'i_onshore_cable': i_onshore_cable,
               'v_offshore_c3': outputs_exportcable['v_offshore_c3'],
               'i_offshore_c3': outputs_exportcable['i_offshore_c3'],
               'v_offshore_c2': outputs_exportcable['v_offshore_c2'],
               'i_offshore_c2': outputs_exportcable['i_offshore_c2'],
               'v_offshore_c1': outputs_exportcable['v_offshore_c1'],
               'i_offshore_c1': outputs_exportcable['i_offshore_c1'],
               'i_onshore_comp_cableside': i_onshore_comp_cableside,
               'i_onshore_tr_cableside': i_onshore_tr_cableside,
               'i_onshore_tr_gridside': i_onshore_tr_gridside,
               'i_onshore_comp_gridside': i_onshore_comp_gridside,
               'v_poc': v_poc,
               'i_poc': i_poc,
               'v_grid': v_grid,
               'i_grid': i_grid,
               'v_nl_onshore': outputs_exportcable['v_nl_onshore'],
               'tap_tr_onshore_nl': tap_tr_onshore_nl,
               'tap_tr_onshore': tap_tr_onshore,
               'q_comp_mid': q_comp_mid,
               #'q_comp_onshore': q_comp_onshore,
               'comments': outputs_exportcable['comments']
               }

    return outputs


def hvac_initialized(project_data, optimisation_data, df_components, initial_values):

    ################
    # project data #
    ################
    # q_project = project_data['s']
    # p_project = project_data['p']

    ##########################
    # data of the components #
    ##########################
    df_offshore_tr = df_components['df offshore transformer']
    df_onshore_tr = df_components['df onshore transformer']

    tap_offshore_tr_pos = 0
    tap_offshore_tr_neg = 0
    tap_onshore_tr_pos = 0
    tap_onshore_tr_neg = 0

    #####################
    # optimisation data #
    #####################
    p_poc = operating_point['P_poc (MW)']
    q_poc = operating_point['Q_poc (MVAR)']
    v_poc_pu = operating_point['V_poc (pu)']
    v_offshore = optimisation_data['V_wf (kV)']
    v_export = optimisation_data['V_export (kV)']
    v_onshore = optimisation_data['V_grid (kV)']
    comp_offshore_location = optimisation_data['Offshore comp']
    comp_midpoint_inout = optimisation_data['Midpoint comp']
    comp_onshore_location = optimisation_data['Onshore comp']
    dv_poc_pos = optimisation_data['DVmaxOnshore positive (%)']
    dv_poc_neg = optimisation_data['DVmaxOnshore negative (%)']
    q_poc_min = optimisation_data['Q_grid min (MVAR)']
    q_poc_max = optimisation_data['Q_grid max (MVAR)']
    n_tr_offshore = optimisation_data['n_parallel']
    n_tr_onshore = optimisation_data['n_parallel']

    ###########
    # AC Grid #
    ###########
    #v_grid = v_onshore *
    ac_grid = AcGrid(project_data, optimisation_data)
    # p_poc = p_project
    # q_poc = ac_grid.q_poc(q_project)

    #vi_poc = ac_grid.vi_poc2(p_poc, q_poc, v_poc_pu)
    v_grid = optimisation_data['V_grid (kV)'] #vi_poc['v_grid']
    #i_grid = vi_poc['i_grid']
    v_poc = v_grid * v_poc_pu#vi_poc['v_poc']
    #i_poc = vi_poc['i_poc']
    i_poc = np.conj((p_poc + 1j * q_poc) / (math.sqrt(3) * v_poc))
    i_grid = i_poc
    #i_grid = i_poc
    i_onshore_tr_gridside = i_poc

    tr_onshore = AcTransformer(project_data, df_onshore_tr, n_tr_onshore)

    if v_onshore > v_export:
        inputside_tr_onshore = 'high voltage'
    else:
        inputside_tr_onshore = 'low voltage'

    if project_data['Tap changer'] == 'Automatic Tap Changing':
        tap_tr_onshore = tr_onshore.tap_changer(v_poc, -i_onshore_tr_gridside, inputside_tr_onshore, v_export)
    else:
        tap_tr_onshore = 0

    vi_onshore_tr_cableside = tr_onshore.vi_output(v_poc, -i_onshore_tr_gridside, inputside_tr_onshore, tap_tr_onshore)
    v_onshore_tr_cableside = vi_onshore_tr_cableside['v_output']
    i_onshore_tr_cableside = -vi_onshore_tr_cableside['i_output']

    v_onshore_cable = v_onshore_tr_cableside
    i_onshore_cable = i_onshore_tr_cableside

    ##########################################
    # Voltages and currents across the cable #
    ##########################################
    export_cable = CableSegmental(project_data, optimisation_data, df_components, operating_point)
    outputs_exportcable = export_cable.vi(v_onshore_cable, i_onshore_cable)

    if 'Error' in outputs_exportcable['comments']:
        return {'comments': outputs_exportcable['comments']}

    i_onshore_cable = outputs_exportcable['i_onshore']

    #if comp_onshore_location == 'Export cable side':
    i_onshore_comp_cableside = i_onshore_tr_cableside - i_onshore_cable
    i_onshore_comp_gridside = 0

    if comp_onshore_location == 'Grid side' and i_onshore_comp_cableside != 0:
        q_comp_onshore = math.sqrt(3) * v_onshore_cable * np.conj(i_onshore_comp_cableside)
        i_onshore_comp_gridside = np.conj(q_comp_onshore / (math.sqrt(3) * v_poc))
        i_onshore_tr_gridside_new = i_poc - i_onshore_comp_gridside
        vi_onshore_tr_cableside_new = tr_onshore.vi_output(v_poc, -i_onshore_tr_gridside_new, inputside_tr_onshore, 0)
        i_onshore_tr_cableside_new = -vi_onshore_tr_cableside_new['i_output']
        d_q_onshore = q_comp_onshore
        d_i_cableside = abs(i_onshore_tr_cableside_new - i_onshore_tr_cableside)
        d_i_cableside_2 = abs(i_onshore_tr_cableside_new - i_onshore_cable)
        k_delta = d_i_cableside_2 / d_i_cableside
        d_q_onshore_2 = k_delta * d_q_onshore
        q_comp_onshore = q_comp_onshore - d_q_onshore_2

        i_onshore_comp_gridside = np.conj(q_comp_onshore / (math.sqrt(3) * v_poc))
        i_onshore_tr_gridside = i_poc - i_onshore_comp_gridside
        i_onshore_comp_cableside = 0

        if v_onshore > v_export:
            inputside_tr_onshore = 'high voltage'
        else:
            inputside_tr_onshore = 'low voltage'

        vi_onshore_tr_cableside = tr_onshore.vi_output(v_poc, -i_onshore_tr_gridside, inputside_tr_onshore, 0)
        v_tr_onshore_cableside = vi_onshore_tr_cableside['v_output']
        i_onshore_tr_cableside = -vi_onshore_tr_cableside['i_output']

        # while abs(abs(i_onshore_tr_cableside) - abs(i_onshore_cable)) > 0.001:
        #     q_comp_onshore -= np.sign(abs(i_onshore_tr_cableside) - abs(i_onshore_cable)) * 0.01j
        #     i_onshore_comp_gridside = np.conj(q_comp_onshore / (math.sqrt(3) * v_poc))
        #     i_onshore_tr_gridside = i_poc - i_onshore_comp_gridside
        #     i_onshore_comp_cableside = 0
        #
        #     if v_onshore > v_export:
        #         inputside_tr_onshore = 'high voltage'
        #     else:
        #         inputside_tr_onshore = 'low voltage'
        #
        #     vi_onshore_tr_cableside = tr_onshore.vi_output(v_poc, -i_onshore_tr_gridside, inputside_tr_onshore, 0)
        #     v_tr_onshore_cableside = vi_onshore_tr_cableside['v_output']
        #     i_onshore_tr_cableside = -vi_onshore_tr_cableside['i_output']

        v_onshore_cable = v_tr_onshore_cableside
        i_onshore_cable = i_onshore_tr_cableside
        outputs_exportcable = export_cable.vi(v_onshore_cable, i_onshore_cable)

    try:
        v_offshore_cable = outputs_exportcable['v_offshore']
        i_offshore_cable = outputs_exportcable['i_offshore']
        v_onshore_cable = outputs_exportcable['v_onshore']
        i_onshore_cable = outputs_exportcable['i_onshore']
        q_comp_mid = outputs_exportcable['q_comp_mid']
    except TypeError and KeyError:
        return outputs_exportcable

    s_onshore_tr_gridside = math.sqrt(3) * v_poc * np.conj(i_onshore_tr_gridside)
    s_onshore_tr_cableside = math.sqrt(3) * v_onshore_tr_cableside * np.conj(i_onshore_tr_cableside)
    s_onshore_tr_abs = max(abs(s_onshore_tr_gridside), abs(s_onshore_tr_cableside))
    if s_onshore_tr_abs > df_onshore_tr['Capacity (MVA)'] * n_tr_onshore + 10:
        return {'comments': "Error: s_onshore_tr > df_onshore_tr[Capacity (MVA)] * n_tr_onshore"}

    ####
    if export_cable.p_limit(outputs_exportcable):
        return {'comments': 'Error: Maximum active power transmission limit is exceeded'}
    ############
    # OFFSHORE #
    ############
    tr_offshore = AcTransformer(project_data, df_offshore_tr, n_tr_offshore)

    if v_export > v_offshore:
        inputside_tr_offshore = 'high voltage'
    else:
        inputside_tr_offshore = 'low voltage'

    v_offshore_tr_cableside = v_offshore_cable
    i_offshore_tr_cableside = i_offshore_cable
    vi_tr_offshore_wfside = tr_offshore.vi_output(v_offshore_cable, -i_offshore_tr_cableside, inputside_tr_offshore, 0)
    v_offshore_tr_wfside = vi_tr_offshore_wfside['v_output']
    i_offshore_tr_wfside = -vi_tr_offshore_wfside['i_output']

    v_offshore = v_offshore_tr_wfside
    i_offshore_wf = i_offshore_tr_wfside

    s_offshore_wf = math.sqrt(3) * v_offshore_tr_wfside * np.conj(i_offshore_wf)
    if s_offshore_wf.real > project_data['P (MW)'] + 20:
        return {'comments': 'Error: P_wf > P_project'}

    try:
        q_available_wf = math.sqrt(project_data['S (MVA)'] ** 2 - s_offshore_wf.real ** 2) * abs(v_offshore / optimisation_data['V_wf (kV)']) ** 2
    except ValueError:
        q_available_wf = 0

    if s_offshore_wf.imag > q_available_wf:
        q_comp_offshore = (s_offshore_wf.imag - q_available_wf) * 1j
    elif s_offshore_wf.imag < -q_available_wf:
        q_comp_offshore = (s_offshore_wf.imag + q_available_wf) * 1j
    else:
        q_comp_offshore = 0

    #if comp_offshore_location == 'Wind Farm side':
    i_offshore_comp_wfside = np.conj(q_comp_offshore / (math.sqrt(3) * v_offshore_tr_wfside))
    i_offshore_comp_cableside = 0

    if comp_offshore_location == 'Export cable side' and i_offshore_comp_wfside != 0:

        s_offshore_tr_wfside = math.sqrt(3) * v_offshore_tr_wfside * np.conj(i_offshore_wf)
        q_offshore_tr_wfside = s_offshore_tr_wfside.imag
        while abs(q_offshore_tr_wfside - np.sign(q_offshore_tr_wfside) * q_available_wf) > 1:
            q_comp_offshore += np.sign(q_offshore_tr_wfside - np.sign(q_offshore_tr_wfside) * q_available_wf) * 0.1j
            i_offshore_comp_cableside = np.conj(q_comp_offshore / (math.sqrt(3) * v_offshore_cable))
            i_offshore_tr_cableside = i_offshore_cable - i_offshore_comp_cableside
            i_offshore_comp_wfside = 0
            if v_export > v_offshore:
                inputside_tr_offshore = 'high voltage'
            else:
                inputside_tr_offshore = 'low voltage'

            vi_offshore_tr_wfside = tr_offshore.vi_output(v_offshore_cable, -i_offshore_tr_cableside, inputside_tr_offshore, 0)
            v_offshore_tr_wfside = vi_offshore_tr_wfside['v_output']
            i_offshore_tr_wfside = -vi_offshore_tr_wfside['i_output']
            v_offshore = v_offshore_tr_wfside

            s_offshore_tr_wfside = math.sqrt(3) * v_offshore_tr_wfside * np.conj(i_offshore_tr_wfside)
            q_offshore_tr_wfside = s_offshore_tr_wfside.imag

    i_offshore_wf = i_offshore_tr_wfside - i_offshore_comp_wfside

    s_offshore_tr_wfside = math.sqrt(3) * v_offshore_tr_wfside * np.conj(i_offshore_tr_wfside)
    s_offshore_tr_cableside = math.sqrt(3) * v_offshore_tr_cableside * np.conj(i_offshore_tr_cableside)
    s_offshore_tr_abs = max(abs(s_offshore_tr_wfside), abs(s_offshore_tr_cableside))
    if s_offshore_tr_abs > df_offshore_tr['Capacity (MVA)'] * n_tr_offshore + 10:
        return {'comments': 'Error: s_offshore_tr > df_offshore_tr[Capacity (MVA)] * n_tr_offshore'}


    ############################################
    # tap changer of the transformer (no load) #
    ############################################
    vi_poc_nl = tr_onshore.vi_output(outputs_exportcable['v_nl_onshore'], outputs_exportcable['i_nl_onshore'], 'high voltage', 0)
    i_poc_nl = vi_poc_nl['i_output']
    v_poc_nl = ac_grid.v_poc(i_poc_nl)
    vi_nl_onshore_2 = tr_onshore.vi_output(v_poc_nl, i_poc_nl, 'low voltage', 0)
    v_nl_onshore_2 = vi_nl_onshore_2['v_output']

    if abs(outputs_exportcable['v_nl_onshore']) < abs(v_nl_onshore_2):
        tap_tr_onshore_nl = tr_onshore.tap_changer(outputs_exportcable['v_nl_onshore'], outputs_exportcable['i_nl_onshore'], 'low voltage', v_poc_nl)
    else:
        tap_tr_onshore_nl = 0

    #########################################
    # Return outputs (voltage and currents) #
    #########################################
    outputs = {'v_offshore': v_offshore_tr_wfside,
               'i_offshore_wf': i_offshore_wf,
               'i_offshore_comp_wfside': i_offshore_comp_wfside,
               'i_offshore_tr_wfside': i_offshore_tr_wfside,
               'i_offshore_tr_cableside': i_offshore_tr_cableside,
               'i_offshore_comp_cableside': i_offshore_comp_cableside,
               'v_offshore_cable': v_offshore_cable,
               'i_offshore_cable': i_offshore_cable,
               'v_onshore_cable': v_onshore_cable,
               'i_onshore_cable': i_onshore_cable,
               'v_offshore_c3': outputs_exportcable['v_offshore_c3'],
               'i_offshore_c3': outputs_exportcable['i_offshore_c3'],
               'v_offshore_c2': outputs_exportcable['v_offshore_c2'],
               'i_offshore_c2': outputs_exportcable['i_offshore_c2'],
               'v_offshore_c1': outputs_exportcable['v_offshore_c1'],
               'i_offshore_c1': outputs_exportcable['i_offshore_c1'],
               'i_onshore_comp_cableside': i_onshore_comp_cableside,
               'i_onshore_tr_cableside': i_onshore_tr_cableside,
               'i_onshore_tr_gridside': i_onshore_tr_gridside,
               'i_onshore_comp_gridside': i_onshore_comp_gridside,
               'v_poc': v_poc,
               'i_poc': i_poc,
               'v_grid': v_grid,
               'i_grid': i_grid,
               'v_nl_onshore': outputs_exportcable['v_nl_onshore'],
               'tap_tr_onshore_nl': tap_tr_onshore_nl,
               'tap_tr_onshore': tap_tr_onshore,
               'q_comp_mid': q_comp_mid,
               #'q_comp_onshore': q_comp_onshore,
               'comments': outputs_exportcable['comments']
               }

    return outputs
