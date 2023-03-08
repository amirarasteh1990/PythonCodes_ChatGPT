import cmath
import math
import numpy as np
import pandas as pd


class AcCable:
    def __init__(self, project_data, n_cable, df_ac_cable, d_start_cable, d_cable, d_comp_mid, temp):
        self.s = project_data['S (MVA)'] / n_cable
        self.p = int(project_data['P (MW)']) / n_cable
        self.d_start = d_start_cable
        self.d = d_cable
        self.d_comp_mid = d_comp_mid
        self.f = project_data['F (Hz)']
        self.w = 2 * math.pi * self.f
        self.n_cable = int(n_cable)

        # Extracting the cable parameters
        if df_ac_cable['Conductor Material'] == 'Cu':
            temp_coef = 3.93e-3
        else:
            temp_coef = 4.03e-3

        self.r_cable_km = df_ac_cable['Resistance 20°C (Ohms/km)'] * (1 + temp_coef * (temp - 20))
        self.l_cable_km = df_ac_cable['Inductance (mH/km)'] * 1e-3
        self.c_cable_km = df_ac_cable['Capacitance (µF/km)'] * 1e-6
        self.zc = cmath.sqrt((self.r_cable_km + self.l_cable_km * self.w * 1j) / (self.c_cable_km * self.w * 1j))
        self.yp = cmath.sqrt((self.r_cable_km + self.l_cable_km * self.w * 1j) * (self.c_cable_km * self.w * 1j))
        self.v_cable = df_ac_cable['Voltage (kV)'] / math.sqrt(3)
        self.v_cable_insulation = df_ac_cable['Insulation Voltage (kV)'] / math.sqrt(3)
        self.i_max_cable = df_ac_cable['Ampacity 20DegC (kA)'] / 1000
        self.key_cable = df_ac_cable['Key']
        self.capacity_mva = df_ac_cable['Capacity (MVA)']
        self.df_ac_cable = df_ac_cable

    def inputs_vioffshore_outputs_vionshore(self, v_offshore, i_offshore):
        length_cable = self.d
        n_cable = self.n_cable
        zc = self.zc
        yp = self.yp
        v_offshore = v_offshore
        i_offshore = i_offshore / n_cable

        i_onshore = i_offshore * cmath.cosh(yp * length_cable) - v_offshore * cmath.sinh(yp * length_cable) / zc
        v_onshore = v_offshore * cmath.cosh(yp * length_cable) - i_offshore * cmath.sinh(yp * length_cable) * zc

        v_onshore = v_onshore
        i_onshore = n_cable * i_onshore

        outputs = {'v_onshore': v_onshore,
                   'i_onshore': i_onshore
                   }

        return outputs

    def inputs_vionshore_outputs_vioffshore(self, v_onshore, i_onshore):
        length_cable = self.d
        n_cable = self.n_cable
        zc = self.zc
        yp = self.yp
        v_onshore = v_onshore
        i_onshore = i_onshore / n_cable

        i_offshore = i_onshore * cmath.cosh(yp * length_cable) + v_onshore * cmath.sinh(yp * length_cable) / zc
        v_offshore = v_onshore * cmath.cosh(yp * length_cable) + i_onshore * cmath.sinh(yp * length_cable) * zc

        i_offshore = n_cable * i_offshore

        outputs = {'v_offshore': v_offshore,
                   'i_offshore': i_offshore
                   }

        return outputs

    def inputs_ioffshore_vonshore_outputs_voffshore_ionshore(self, i_offshore, v_onshore):
        length_cable = self.d
        n_cable = self.n_cable
        zc = self.zc
        yp = self.yp
        i_offshore = i_offshore / n_cable

        v_offshore = (v_onshore + i_offshore * zc * cmath.sinh(yp * length_cable)) / cmath.cosh(yp * length_cable)
        i_onshore = i_offshore * cmath.cosh(yp * length_cable)  - v_onshore * cmath.sinh(yp * length_cable) / zc

        i_onshore = n_cable * i_onshore

        outputs = {'v_offshore': v_offshore,
                   'i_onshore': i_onshore
                   }

        return outputs

    def inputs_ionshore_voffshore_outputs_vonshore_ioffshore(self, i_onshore, v_offshore):
        length_cable = self.d
        n_cable = self.n_cable
        zc = self.zc
        yp = self.yp
        i_onshore = i_onshore / n_cable

        v_onshore = (v_offshore - i_onshore * zc * cmath.sinh(yp * length_cable)) / cmath.cosh(yp * length_cable)
        i_offshore = i_onshore * cmath.cosh(yp * length_cable) + v_onshore * cmath.sinh(yp * length_cable) / zc

        i_offshore = n_cable * i_offshore

        outputs = {'v_onshore': v_onshore,
                   'i_offshore': i_offshore
                   }

        return outputs

    def p_max(self, v_onshore, v_offshore):
        length_cable = self.d
        zc = self.zc
        yp = self.yp
        alpha = cmath.cosh(yp * length_cable)
        beta = zc * cmath.sinh(yp * length_cable)
        if length_cable > 0:
            p_max = (abs(v_onshore) * abs(v_offshore) - abs(alpha) * abs(v_onshore) ** 2 * math.cos(cmath.phase(beta)-cmath.phase(alpha))) / abs(beta)
        else:
            p_max = float('inf')
        return p_max

    def resistance(self, cond_max_temp, Depth, SoilResist, al_cu):
        df_ac_cables = self.df_ac_cable
        cond_section_mm2 = df_ac_cables['Section (mm2)']
        cond_diam_mm = df_ac_cables['Diam (mm)']
        space_mm = df_ac_cables['Diam over insulation (mm)']
        rho_elec_cu = 1.7241e-8
        rho_elec_al = 2.8264e-8
        temp_coef_cu = 3.93e-3
        temp_coef_al = 4.03e-3
        ks_cu = 1
        ks_al = 1
        kp_cu = 1
        kp_al = 1
        # Calculation of the Steady - State permissible current(SSCurrent) for a three - phase group of single-core cables according to IEC - 60287 - 1 - 1 Page 25
        # T4 = ((SoilResist) / (2 * pi)) * ( (log((4 * Depth) / (Cable.CableD * 10 ** (-3)))) + (log(1 + ((2 * Depth) / (Space)) ^ 2))) # T4 Soil Thermal Resistance(K.m / W)

        r_dc_20 = 1e3 * (rho_elec_cu) / (cond_section_mm2 * 10 ** (-6))
        r_dc_max = r_dc_20 * (
                1 + temp_coef_cu * (cond_max_temp - 20))  # Maximum DC resistance of the conductor(Ohm / m)

        xs2 = ((8 * math.pi * self.f) / (r_dc_max)) * 10 ** (-7) * ks_cu
        ys = (xs2 ** 2) / (192 + 0.8 * xs2 ** 2)

        xp2 = ((8 * math.pi * self.f) / r_dc_max) * 10 ** (-7) * kp_cu
        expression_xp2 = (xp2 ** 2) / (192 + 0.8 * xp2 ** 2)
        yp = expression_xp2 * (cond_diam_mm / space_mm) ** 2 * ((0.312 * (cond_diam_mm / space_mm) ** 2) + (
                1.18 / (expression_xp2 + 0.27)))

        r_ac_max = r_dc_max * (1 + ys + yp)  # Maximum AC resistance of the conductor(Ohm / m)

        # [CondDiam, InsuDiam, ~] = Cable.Diameters()

        # ElectCap = (Cable.XLPEe) * 10 ^ (-9) / (18 * log((InsuDiam) / (CondDiam))) # Insulation equivalent capacitance(F / m)

        # wd = 2 * pi * f * ElectCap * (Cable.Voltagelevel * 1000 / (sqrt(3))) ^ 2 * Cable.XLPETanDelta # Insulation Dielectric losses(W / m)

        # [T1, T2, T3] = Cable.EqThermalResis()

        # current = sqrt((CondMaxTemp - (20) - wd * (0.5 * T1 + (T2 + T3 + T4))) / ( r_ac_max * T1 + r_ac_max * T2 * (1 + Cable.Lambda1) + r_ac_max * (1 + Cable.Lambda1 + Cable.Lambda2) * ( T3 + T4)))
        return r_dc_max, r_ac_max

    def max_vi(self):
        length_cable = self.d
        n_cable = self.n_cable
        p = self.p
        q = self.q_wf
        v_start = self.v_offshore
        length_segment = 1
        i_start = (p + q * 1j) / (math.sqrt(3) * v_start * n_cable)
        v_start = v_start / math.sqrt(3)
        abs_v_max = 0
        abs_i_max = 0
        #
        r_cable = self.r_cable_km
        l_cable = self.l_cable_km
        c_cable = self.c_cable_km
        w = self.w
        zc = cmath.sqrt((r_cable + l_cable * w * 1j) / (c_cable * w * 1j))
        yc = cmath.sqrt((r_cable + l_cable * w * 1j) * (c_cable * w * 1j))

        for n in range(length_cable):
            v_end = v_start * cmath.cosh(yc * length_segment) - i_start * zc * cmath.sinh(yc * length_segment)
            i_end = i_start * cmath.cosh(yc * length_segment) - (v_start / zc) * cmath.sinh(yc * length_segment)
            v_start = v_end
            i_start = i_end

            if abs(v_end) > abs_v_max:
                abs_v_max = abs(v_end)

            if abs(i_end) > abs_i_max:
                abs_i_max = abs(i_end)

        abs_v_max = abs_v_max * math.sqrt(3)
        return abs_v_max, abs_i_max

    # def losses(self, cable_list, n_cable, p, q, v_start):
    #     length_cable = self.d
    #     i_start = (p + q * 1j) / (math.sqrt(3) * v_start * n_cable)
    #     v_start = v_start / math.sqrt(3)
    #     # r_cable = cable_list['Resistance 90°C (Ohms/km)']
    #     # l_cable = cable_list['Inductance (mH/km)']
    #     # c_cable = cable_list['Capacitance (µF/km)']
    #     # z_cable = math.sqrt(r_cable ** 2 + l_cable ** 2)
    #     # zc = math.sqrt(z_cable / c_cable)
    #     # yc = math.sqrt(z_cable * c_cable)
    #
    #     # sam
    #     list_res = []
    #     for n in [1]:  # ind, data in cable_list.iterrows():
    #         r_cable = 0.01  # data['Resistance 90°C (Ohms/km)']
    #         l_cable = 0.4 * 1e-3  # data['Inductance (mH/km)'] * 1e-3
    #         c_cable = 0.2 * 1e-6  # data['Capacitance (µF/km)'] * 1e-6
    #         z_cable = math.sqrt(r_cable ** 2 + l_cable ** 2)
    #         w = 2 * math.pi * self.f
    #         zc = cmath.sqrt((r_cable + l_cable * w * 1j) / (c_cable * w * 1j))
    #         yc = cmath.sqrt((r_cable + l_cable * w * 1j) * (c_cable * w * 1j))
    #         v_end = v_start * cmath.cosh(yc * length_cable) - i_start * zc * cmath.sinh(yc * length_cable)
    #         i_end = i_start * cmath.cosh(yc * length_cable) - (v_start / zc) * cmath.sinh(yc * length_cable)
    #         v_end_oc = v_start / cmath.cosh(yc * length_cable)
    #         v_end_nl = v_start / cmath.cosh(yc * length_cable)
    #         v_end = v_end * math.sqrt(3)
    #         abs_v_end = abs(v_end)
    #         phase_v_end = cmath.phase(v_end) * 180 / math.pi
    #         dv_cable_percent = 100 * (abs(v_start) - abs(v_end)) / abs(v_start)
    #         abs_i_end = abs(i_end)
    #         phase_i_end = cmath.phase(i_end)
    #         s_end = math.sqrt(3) * v_end * np.conj(i_end) * n_cable
    #         p_end = s_end.real
    #         q_end = s_end.imag
    #         dp_cable = p - p_end
    #         dq_cable = q - q_end
    #         list_res.append([abs_v_end, phase_v_end, dv_cable_percent, abs_i_end,
    #                          phase_i_end, p_end, q_end, dp_cable, dq_cable])
    #         # df_results = pd.concat([df_results, pd.DataFrame([abs_v_end, phase_v_end, dv_cable_percent, abs_i_end,
    #         #                                                   phase_i_end, p_end, q_end, dp_cable, dq_cable])])
    #     df_results = pd.DataFrame(list_res, columns=['abs_v_end', 'phase_v_end', 'dv_cable_percent', 'abs_i_end',
    #                                                  'phase_i_end', 'p_end', 'q_end', 'dp_cable', 'dq_cable'])
    #     df_results.insert(0, 'Key', cable_list.Key)
    #     return df_results
    #
    # def losses_piecewise(self, cable_list, n_cable, p, q, v_start):
    #     length_cable = self.d
    #     i_start = (p + q * 1j) / (math.sqrt(3) * v_start * n_cable)
    #     list_res = []
    #     for ind, data in cable_list.iterrows():
    #         r_cable = 0.08  # data['Resistance 90°C (Ohms/km)']
    #         l_cable = 0.44 * 1e-3  # data['Inductance (mH/km)'] * 1e-3
    #         c_cable = 0.14 * 1e-6  # data['Capacitance (µF/km)'] * 1e-6
    #         z_cable = math.sqrt(r_cable ** 2 + l_cable ** 2)
    #         v_end = v_start
    #         i_end = i_start
    #         for n in range(length_cable):
    #             v_end += (r_cable + l_cable * 1e-3 * 2 * math.pi * self.f * 1j) * i_end
    #             i_cap_cable = v_end * (c_cable * 1e-6 * 2 * math.pi * self.f) * 1j
    #             i_end += i_cap_cable
    #
    #         abs_v_end = abs(v_end)
    #         phase_v_end = cmath.phase(v_end)
    #         dv_cable_percent = 100 * (abs(v_start) - abs(v_end)) / abs(v_start)
    #         abs_i_end = abs(i_end)
    #         phase_i_end = cmath.phase(i_end)
    #         s_end = math.sqrt(3) * v_end * numpy.conj(i_end) * n_cable
    #         p_end = s_end.real
    #         q_end = s_end.imag
    #         dp_cable = p - p_end
    #         dq_cable = q - q_end
    #         list_res.append([abs_v_end, phase_v_end, dv_cable_percent, abs_i_end,
    #                          phase_i_end, p_end, q_end, dp_cable, dq_cable])
    #         # df_results = pd.concat([df_results, pd.DataFrame([abs_v_end, phase_v_end, dv_cable_percent, abs_i_end,
    #         #                                                   phase_i_end, p_end, q_end, dp_cable, dq_cable])])
    #     df_results = pd.DataFrame(list_res, columns=['abs_v_end', 'phase_v_end', 'dv_cable_percent', 'abs_i_end',
    #                                                  'phase_i_end', 'p_end', 'q_end', 'dp_cable', 'dq_cable'])
    #     df_results.insert(0, 'Key', cable_list.Key)
    #     return df_results
    #     # i_cable = (p + q * 1j) / (v * n_cable)
    #     #
    #     # abs_v_end = abs(v_end)
    #     # phase_v_end = cmath.phase(v_end)
    #     # dv_cable_percent = 100 * (abs(v) - abs(v_end)) / abs(v)
    #     # abs_i_end = abs(i_end)
    #     # phase_i_end = cmath.phase(i_end)
    #     # s_end = math.sqrt(3) * v_end * numpy.conj(i_end) * n_cable
    #     # p_end = s_end.real
    #     # q_end = s_end.imag
    #     # dp_cable = p - p_end
    #     # dq_cable = q - q_end
    #     # return abs_v_end, phase_v_end, dv_cable_percent, abs_i_end, phase_i_end, p_end, q_end, dp_cable, dq_cable
    #
    # def sizing0(self):
    #     df_cables = self.df_ac_cable
    #     cable_v = df_cables['Voltage (kV)']
    #     i_cable = self.s / (self.n_cable * self.v_offshore)
    #     s_cable = self.s / self.n_cable
    #     df_ac_cable = self.df_ac_cable
    #     df_ac_cable_v = df_ac_cable['Voltage (kV)']
    #     df_ac_cable_a = df_ac_cable['Ampacity 20DegC']
    #     df_ac_cable_s = df_ac_cable['Capacity (MVA)']
    #     df_ac_cable_losses = self.losses(df_cables, self.n_cable, self.p, 0, self.v_offshore)
    #     df_ac_cable_insulationvoltage = df_ac_cable['Insulation Voltage (kV)']
    #     v_min = 0  # df_ac_cable_losses[0]
    #     amp_min = 0  # df_ac_cable_losses[3]
    #     s_min = 0  # abs(df_ac_cable_losses[5] + df_ac_cable_losses[6])
    #     ac_cable_valid_list = df_ac_cable[(df_ac_cable_v > v_min) & (df_ac_cable_a > amp_min) & (df_ac_cable_s > s_min)]
    #     ac_cables_valid_min = min(ac_cable_valid_list['Capacity (MVA)'])
    #     ac_cable_opt = ac_cable_valid_list[ac_cable_valid_list['Capacity (MVA)'] == ac_cables_valid_min]
    #     return ac_cable_opt, df_ac_cable_losses
    # def q_max_onshore(self, s_onshore):
    #     p_onshore = s_onshore.real
    #     capacity_mva = self.capacity_mva * self.n_cable
    #     try:
    #         q_max_onshore = math.sqrt(capacity_mva ** 2 - p_onshore ** 2)
    #     except ValueError:
    #         q_max_onshore = 0
    #     return q_max_onshore
    #
    # def sizing(self):
    #     v_onshore_shifted = self.v_onshore_shifted
    #     v_onshore_phase = self.v_onshore_phase
    #     i_onshore_shifted = self.i_onshore_shifted
    #     length_cable = self.d
    #     n_cable = self.n_cable
    #     ip_onshore = i_onshore_shifted.real
    #     iq_onshore_shifted = i_onshore_shifted.imag
    #     i_max_cable = self.i_max_cable
    #     iq_max_cable_onshore = self.iq_max_cable_onshore
    #     zc = self.zc
    #     yp = self.yp
    #     v_cable = self.v_cable
    #     v_cable_insulation = self.v_cable_insulation
    #     key_cable = self.key_cable
    #
    #     if v_cable < v_onshore_shifted:
    #         return 'Error: the cable is not fitting (rated voltage)'
    #
    #     if ip_onshore < i_max_cable:
    #         iq_offshore0 = 0
    #         ip_onshore_shifted = ip_onshore
    #     else:
    #         return f"Error: the current limit is exceeded (p), ip_offshore = {round(ip_offshore, 3)} and i_max_cable {round(i_max_cable, 3)}."
    #
    #     # Initialisation
    #     q_svc_mid = 0
    #     y_svc_mid = 0
    #     cable_comments = 'None'
    #
    #     ###########################################################
    #     # q_wf (Compensation by the wind farm remaining capacity) #
    #     ###########################################################
    #
    #     required_compensation = 'Wind Farm'
    #
    #     if abs(iq_max_cable_onshore) < abs(iq_onshore_shifted):
    #         cable_comments = 'Notice: iq_max_cable < iq_offshore'
    #
    #     i_onshore_shifted = ip_onshore + min(iq_onshore_shifted, iq_onshore_shifted) * 1j
    #
    #     i_onshore_shifted = i_onshore_shifted * cmath.cosh(yp * length_cable) - v_offshore_shifted * cmath.sinh(
    #         yp * length_cable) / zc
    #     v_onshore_shifted = v_offshore_shifted * cmath.cosh(yp * length_cable) - i_offshore_shifted * cmath.sinh(
    #         yp * length_cable) * zc
    #
    #     v_onshore_shifted_min = v_onshore_shifted
    #     i_onshore_shifted_min = i_onshore_shifted
    #     v_onshore_shifted_max = v_onshore_shifted
    #     i_onshore_shifted_max = i_onshore_shifted
    #
    #     # ------------------------------#
    #     # range of vi_onshore with q_wf #
    #     # ------------------------------#
    #     i_onshore_shifted2 = i_onshore_shifted
    #     v_onshore_shifted2 = v_onshore_shifted
    #
    #     if abs(i_onshore_shifted) < i_max_cable:
    #         cable_comments = 'Notice: i_onshore < i_max_cable'
    #         i_offshore_shifted2 = ip_offshore - iq_offshore_shifted * 1j
    #         i_onshore_shifted2 = i_offshore_shifted2 * cmath.cosh(yp * length_cable) - v_offshore_shifted * cmath.sinh(
    #             yp * length_cable) / zc
    #
    #         while abs(i_onshore_shifted2) < i_max_cable:
    #             i_offshore_shifted2 += 0.01j
    #             i_onshore_shifted2 = i_offshore_shifted2 * cmath.cosh(
    #                 yp * length_cable) - v_offshore_shifted * cmath.sinh(yp * length_cable) / zc
    #
    #         i_onshore_shifted2 = i_offshore_shifted2 * cmath.cosh(
    #             yp * length_cable) - v_offshore_shifted * cmath.sinh(yp * length_cable) / zc
    #         v_onshore_shifted2 = v_offshore_shifted * cmath.cosh(yp * length_cable) - i_offshore_shifted2 * cmath.sinh(
    #             yp * length_cable) * zc
    #
    #     ##########################################
    #     # q_svc_offshore (offshore compensation) #
    #     ##########################################
    #
    #     if abs(i_onshore_shifted) > i_max_cable:
    #         required_compensation = 'Offshore'
    #         if abs(i_offshore_shifted) < i_max_cable:
    #             while abs(i_onshore_shifted) > i_max_cable:
    #                 i_offshore_shifted += 0.01j
    #                 i_onshore_shifted = i_offshore_shifted * cmath.cosh(
    #                     yp * length_cable) - v_offshore_shifted * cmath.sinh(yp * length_cable) / zc
    #                 if abs(i_offshore_shifted) > i_max_cable:
    #                     required_compensation = 'Mid point'
    #                     break
    #
    #     # ------------------------------#
    #     # range of q_onshore2 with q_wf #
    #     # ------------------------------#
    #
    #     # if abs(i_offshore_shifted.imag) < min(iq_max_wf, iq_max_offshore):
    #     #     i_offshore_shifted2 = ip - iq_max_wf * 1j
    #     #     i_onshore_shifted2 = i_offshore_shifted2 * cmath.cosh(yp * length_cable) - v_offshore_shifted * cmath.sinh(
    #     #         yp * length_cable) / zc
    #     #     v_onshore_shifted2 = v_offshore_shifted * cmath.cosh(yp * length_cable) - i_offshore_shifted2 * cmath.sinh(
    #     #         yp * length_cable) * zc
    #     #
    #     #     s_onshore2 = 3 * v_onshore_shifted2 * np.conj(i_onshore_shifted2)
    #     #     q_onshore2 = -s_onshore2.imag
    #
    #     ######################################
    #     # q_svc_mid (Mid point compensation) #
    #     ######################################
    #
    #     if required_compensation == 'Mid point':
    #
    #         return 'Error: Midpoint compensation is required'
    #
    #         i_offshore_shifted = complex(ip_offshore, iq_offshore_shifted)
    #
    #         v_mid_shifted = v_offshore_shifted * cmath.cosh(
    #             yp * length_cable / 2) - i_offshore_shifted * zc * cmath.sinh(yp * length_cable / 2)
    #         i_offshore_mid = i_offshore_shifted * cmath.cosh(
    #             yp * length_cable / 2) - v_offshore_shifted * cmath.sinh(yp * length_cable / 2) / zc
    #
    #         if abs(i_offshore_mid) > i_max_cable:
    #             return 'Error: the length limit is exceeded (offshore)'
    #
    #         v_mid_shifted_phase = cmath.phase(v_mid_shifted)
    #         v_mid_shifted2 = abs(v_mid_shifted)
    #         s_offshore_mid = 3 * v_mid_shifted * np.conj(i_offshore_mid)
    #         p_offshore_mid = s_offshore_mid.real
    #         ip_offshore_mid = p_offshore_mid / (3 * v_mid_shifted2)
    #
    #         ip_mid_onshore = ip_offshore_mid
    #         iq_mid_onshore = 0
    #
    #         i_mid_onshore = ip_mid_onshore + iq_mid_onshore * 1j
    #         i_onshore_shifted2 = i_mid_onshore * cmath.cosh(yp * length_cable / 2) - v_mid_shifted2 * cmath.sinh(
    #             yp * length_cable / 2) / zc
    #
    #         while abs(i_onshore_shifted2) > i_max_cable:
    #             i_mid_onshore += 0.01j
    #             i_onshore_shifted2 = i_mid_onshore * cmath.cosh(
    #                 yp * length_cable / 2) - v_mid_shifted2 * cmath.sinh(yp * length_cable / 2) / zc
    #             if abs(i_mid_onshore) > i_max_cable:
    #                 return 'Error: the length limit is exceeded (onshore)'
    #
    #         s_offshore_mid = 3 * v_mid_shifted * np.conj(i_offshore_mid)
    #         q_svc_offshore_mid = s_offshore_mid.imag
    #
    #         s_mid_onshore = 3 * v_mid_shifted * np.conj(i_mid_onshore)
    #         q_svc_mid_onshore = s_mid_onshore.imag
    #
    #         q_svc_mid = q_svc_offshore_mid + q_svc_mid_onshore
    #
    #         y_svc_mid = ((3 * v_mid_shifted ** 2) / q_svc_mid) * 1j
    #
    #         i_onshore_shifted = i_onshore_shifted2 * cmath.exp(
    #             1j * v_mid_shifted_phase)
    #         v_onshore_shifted2 = v_mid_shifted * cmath.cosh(
    #             yp * length_cable / 2) - i_mid_onshore * zc * cmath.sinh(
    #             yp * length_cable / 2)
    #         v_onshore_shifted = v_onshore_shifted2 * cmath.exp(
    #             1j * v_mid_shifted_phase)
    #
    #     #############################################################
    #     # No load voltage and current + Maximum voltage and current #
    #     #############################################################
    #
    #     if required_compensation != 'Mid point':
    #         v_onshore_noload = abs(v_offshore_shifted / cmath.cosh(yp * length_cable))
    #         i_offshore_noload = v_onshore_noload * cmath.sinh(yp * length_cable) / zc
    #         i_max = abs(i_onshore_shifted)
    #         v_max = abs(v_onshore_shifted)
    #         v_onshore_noload_limit = abs(v_offshore_shifted * cmath.cosh(yp * length_cable))
    #     else:
    #         matrix_a = np.array([[cmath.cosh(yp * length_cable / 2), zc * cmath.sinh(yp * length_cable / 2)],
    #                              [cmath.sinh(yp * length_cable / 2) / zc, cmath.cosh(yp * length_cable / 2)]])
    #         matrix_b = np.array([[1, 0], [y_svc_mid, 1]])
    #         matrix_equivalent = np.matmul(matrix_a, matrix_b, matrix_a)
    #         v_onshore_noload = v_offshore_shifted / matrix_equivalent[0][0]
    #         i_offshore_noload = v_onshore_noload * matrix_equivalent[0][1]
    #         [v_max, i_max] = [0, 0] #self.max_vi()
    #
    #     if v_cable_insulation < v_max:
    #         return 'Error: the cable is not fitting (insulation voltage)'
    #
    #
    #
    #     # if v_cable_insulation < v_onshore_noload:
    #     #     return 'The cable is not fitting (insulation voltage - no load)'
    #
    #     # ##########################################
    #     # # Offshore and Onshore Power calculation #
    #     # ##########################################
    #     # s_offshore = 3 * v_offshore_shifted * np.conj(i_offshore_shifted)
    #     # q_svc_offshore = -s_offshore.imag - self.q_wf
    #     #
    #     # s_onshore = 3 * v_onshore_shifted * np.conj(i_onshore_shifted)
    #     # q_onshore = -s_onshore.imag
    #     #
    #     # s_loss = s_offshore - s_onshore
    #     # p_loss = s_loss.real
    #
    #     ###########################################
    #     # Shifting back the voltages and currents #
    #     ###########################################
    #
    #     v_offshore = v_offshore_shifted * cmath.exp(v_offshore_phase * 1j)
    #     v_onshore = v_onshore_shifted * cmath.exp(v_offshore_phase * 1j)
    #     i_offshore = i_offshore_shifted * cmath.exp(v_offshore_phase * 1j)
    #     i_onshore = i_onshore_shifted * cmath.exp(v_offshore_phase * 1j)
    #
    #     v_onshore_min = v_onshore_shifted_min * cmath.exp(v_offshore_phase * 1j)
    #     i_onshore_min = i_onshore_shifted_min * cmath.exp(v_offshore_phase * 1j)
    #     v_onshore_max = v_onshore_shifted_max * cmath.exp(v_offshore_phase * 1j)
    #     i_onshore_max = i_onshore_shifted_max * cmath.exp(v_offshore_phase * 1j)
    #
    #     v_offshore = math.sqrt(3) * v_offshore
    #     v_onshore = math.sqrt(3) * v_onshore
    #     i_offshore = n_cable * i_offshore
    #     i_onshore = n_cable * i_onshore
    #
    #     v_onshore_min = math.sqrt(3) * v_onshore_min
    #     i_onshore_min = n_cable * np.conj(i_onshore_min)
    #     v_onshore_max = math.sqrt(3) * v_onshore_max
    #     i_onshore_max = n_cable * np.conj(i_onshore_max)
    #
    #     v_onshore_noload = math.sqrt(3) * v_onshore_noload
    #     v_onshore_noload_limit = math.sqrt(3) * v_onshore_noload_limit
    #     v_max = 0 # math.sqrt(3) * v_max
    #
    #     # outputs0 = {'p_loss': round(p_loss, 2),
    #     #             'q_svc_offshore': round(q_svc_offshore, 2),
    #     #             'q_svc_mid': round(q_svc_mid, 2),
    #     #             'q_onshore': round(q_onshore, 2),
    #     #             'q_onshore2': round(q_onshore2, 2),
    #     #             'v_offshore': complex(round(v_offshore.real, 2), round(v_offshore.imag, 2)),
    #     #             'v_onshore': complex(round(v_onshore.real, 2), round(v_onshore.imag, 2)),
    #     #             'i_offshore': complex(round(i_offshore.real, 2), -round(i_offshore.imag, 2)),
    #     #             'i_onshore': complex(round(i_onshore.real, 2), -round(i_onshore.imag, 2)),
    #     #             'i_max_cable': i_max_cable,
    #     #             'v_onshore_noload': round(abs(v_onshore_noload), 2),
    #     #             'i_offshore_noload': round(abs(i_offshore_noload), 2),
    #     #             'v_max': round(v_max, 2),
    #     #             'i_max': round(i_max, 2),
    #     #             'comments': cable_comments}
    #
    #     outputs = {'v_offshore': v_offshore,
    #                'v_onshore': v_onshore,
    #                'i_offshore': i_offshore,
    #                'i_onshore': i_onshore,
    #                'v_onshore_min': v_onshore_min,
    #                'v_onshore_max': v_onshore_max,
    #                'i_onshore_min': abs(i_onshore_min),
    #                'i_onshore_max': abs(i_onshore_max),
    #                'v_onshore_noload': abs(v_onshore_noload),
    #                'i_offshore_noload': abs(i_offshore_noload),
    #                'v_onshore_noload_limit': v_onshore_noload_limit,
    #                'v_max': v_max,
    #                'i_max': i_max,
    #                'q_svc_mid': q_svc_mid,
    #                'y_svc_mid': y_svc_mid,
    #                'comments': cable_comments}
    #
    #     return outputs
    #
    #     #######################
    #     # Power calcaultation #
    #     #######################
    #
    #     # s_offshore = 3 * v_offshore_abs * np.conj(i_offshore_shifted)
    #     # s_onshore = 3 * v_onshore_shifted * np.conj(i_onshore_shifted)
    #     # q_offshore = s_offshore.imag
    #     # q_onshore = s_onshore.imag
    #     # q_onshore2 = 0
    #     # q_svc_offshore = 0
    #     # q_svc_midpoint = 0
    #     #
    #     #
    #     #
    #     # if abs(i_onshore_shifted2) > i_max_cable:
    #     #
    #     #     # q_svc_onshore (Onshore compensation)
    #     #     i_onshore_phase2 = -math.pi / 2
    #     #     i_offshore_shifted2 = 0
    #     #
    #     #     while abs(i_offshore_shifted2.real - ip) > 0.01:
    #     #         i_onshore_shifted2 = i_max_cable * cmath.exp(1j * i_onshore_phase2)
    #     #         i_offshore_shifted2 = (zc * i_onshore_shifted2 + v_offshore_shifted * cmath.sinh(yp * length_cable)) / (
    #     #                 zc * cmath.cosh(yp * length_cable))
    #     #         i_onshore_phase2 += 0.01
    #     #         if i_onshore_phase2 > math.pi / 2:
    #     #             return 'i_offshore is not found'
    #     #
    #     #     i_onshore_shifted2 = i_offshore_shifted2 * cmath.cosh(yp * length_cable) - v_offshore_shifted * cmath.sinh(
    #     #         yp * length_cable) / zc
    #     #     v_onshore_shifted2 = v_offshore_shifted * cmath.cosh(yp * length_cable) - i_offshore_shifted2 * cmath.sinh(
    #     #         yp * length_cable) * zc
    #     #
    #     #     s_offshore2 = 3 * v_offshore_shifted * np.conj(i_offshore_shifted2)
    #     #     s_onshore2 = 3 * v_onshore_shifted2 * np.conj(i_onshore_shifted2)
    #     #     q_onshore2 = s_onshore2.imag
    #     #
    #     #     if abs(i_onshore_shifted) > i_max_cable:
    #     #
    #     #         i_onshore_phase = -math.pi / 2
    #     #         i_offshore_shifted = 0
    #     #
    #     #         while abs(i_offshore_shifted.real - ip) > 0.01:
    #     #             i_onshore_shifted = i_max_cable * cmath.exp(1j * i_onshore_phase)
    #     #             i_offshore_shifted = (zc * i_onshore_shifted + v_offshore_abs * cmath.sinh(yp * length_cable)) / (
    #     #                     zc * cmath.cosh(yp * length_cable))
    #     #             i_onshore_phase += 0.01
    #     #             if i_onshore_phase > math.pi / 2:
    #     #                 return 'i_offshore is not found'
    #     #
    #     #         while abs(i_offshore_shifted.real - ip) > 0.001:
    #     #             i_onshore_shifted = i_max_cable * cmath.exp(1j * i_onshore_phase)
    #     #             i_offshore_shifted = (zc * i_onshore_shifted + v_offshore_abs * cmath.sinh(yp * length_cable)) / (
    #     #                     zc * cmath.cosh(yp * length_cable))
    #     #             i_onshore_phase += 0.001
    #     #
    #     #         i_offshore_shifted = ip + i_offshore_shifted.imag * 1j
    #     #         i_onshore_shifted = i_offshore_shifted * cmath.cosh(yp * length_cable) - v_offshore_abs * cmath.sinh(
    #     #             yp * length_cable) / zc
    #     #         v_onshore_shifted = v_offshore_abs * cmath.cosh(
    #     #             yp * length_cable) - i_offshore_shifted * zc * cmath.sinh(
    #     #             yp * length_cable)
    #     #         abs_i_onshore = abs(i_onshore_shifted)
    #     #         abs_i_offshore = abs(i_offshore_shifted)
    #     #
    #     #         s_offshore = 3 * v_offshore_abs * np.conj(i_offshore_shifted)
    #     #         s_onshore = 3 * v_onshore_shifted * np.conj(i_onshore_shifted)
    #     #         q_offshore = -s_offshore.imag
    #     #         q_onshore = s_onshore.imag
    #     #
    #     #         v_onshore_noload = v_offshore_abs / cmath.cosh(yp * length_cable)
    #     #         i_offshore_noload = v_onshore_noload * cmath.sinh(yp * length_cable) / zc
    #     #         i_max = abs(i_onshore_shifted)
    #     #         v_max = abs(v_onshore_shifted)
    #     #         iq_offshore = i_offshore_shifted.imag
    #     #
    #     #         q_svc_offshore = 0
    #     #
    #     #         if abs(i_offshore_shifted) < i_max_cable:
    #     #             if q_offshore > 0:
    #     #                 q_svc_offshore = q_offshore - q_max_wf
    #     #             else:
    #     #                 q_svc_offshore = q_offshore + q_max_wf
    #     #             q_svc_midpoint = 0
    #     #
    #     #         elif abs_i_offshore > i_max_cable:
    #     #             i_offshore_shifted = complex(ip, iq_max_offshore)
    #     #             s_offshore = 3 * v_offshore_shifted * np.conj(i_offshore_shifted)
    #     #             q_svc_offshore = -s_offshore.imag - self.q_wf
    #     #
    #     #             v_mid_shifted = v_offshore_shifted * cmath.cosh(
    #     #                 yp * length_cable / 2) - i_offshore_shifted * zc * cmath.sinh(yp * length_cable / 2)
    #     #             i_mid_offshore_shifted = i_offshore_shifted * cmath.cosh(
    #     #                 yp * length_cable / 2) - v_offshore_shifted * cmath.sinh(yp * length_cable / 2) / zc
    #     #
    #     #             v_mid_phase = cmath.phase(v_mid_shifted)
    #     #             v_mid_shifted = abs(v_mid_shifted)
    #     #             s_mid_offshore = 3 * v_mid_shifted * np.conj(i_mid_offshore_shifted)
    #     #             p_mid_offshore = s_mid_offshore.real
    #     #             ip_mid_offshore = p_mid_offshore / (3 * v_offshore_shifted)
    #     #
    #     #             i_onshore_phase = -math.pi / 2
    #     #             i_mid_onshore_shifted = 0
    #     #
    #     #             while abs(i_mid_onshore_shifted.real - ip_mid_offshore) > 0.01:
    #     #                 i_onshore_shifted = i_max_cable * cmath.exp(1j * i_onshore_phase)
    #     #                 i_mid_onshore_shifted = (zc * i_onshore_shifted + v_offshore_abs * cmath.sinh(
    #     #                     yp * length_cable / 2)) / (
    #     #                                                 zc * cmath.cosh(yp * length_cable / 2))
    #     #                 i_onshore_phase += 0.01
    #     #                 if i_onshore_phase > math.pi / 2:
    #     #                     return 'i_offshore is not found in the second loop'
    #     #
    #     #             while abs(i_mid_onshore_shifted.real - ip_mid_offshore) > 0.001:
    #     #                 i_onshore_shifted = i_max_cable * cmath.exp(1j * i_onshore_phase)
    #     #                 i_mid_onshore_shifted = (zc * i_onshore_shifted + v_offshore_abs * cmath.sinh(
    #     #                     yp * length_cable / 2)) / (
    #     #                                                 zc * cmath.cosh(yp * length_cable / 2))
    #     #                 i_onshore_phase += 0.001
    #     #
    #     #             i_mid_onshore = i_mid_onshore_shifted * cmath.exp(
    #     #                 1j * v_mid_phase)  # (zc * i_onshore_shifted + v_mid * cmath.sinh(yp * length_cable / 2)) / (zc * cmath.cosh(yp * length_cable / 2))
    #     #             v_onshore_shifted = v_mid_shifted * cmath.cosh(
    #     #                 yp * length_cable / 2) - i_mid_onshore * zc * cmath.sinh(
    #     #                 yp * length_cable / 2)
    #     #
    #     #             if abs(i_mid_offshore_shifted) > i_max_cable or abs(i_mid_onshore) > i_max_cable:
    #     #                 return 'The length limit is exceeded'
    #     #
    #     #             s_mid_onshore = 3 * v_mid_shifted * np.conj(i_mid_onshore)
    #     #             q_svc_mid_onshore = s_mid_onshore.imag
    #     #
    #     #             q_svc_mid_offshore = s_mid_offshore.imag
    #     #             q_svc_midpoint = q_svc_mid_onshore + q_svc_mid_offshore
    #     #
    #     #             s_onshore = 3 * v_onshore_shifted * np.conj(i_onshore_shifted)
    #     #             q_onshore = s_onshore.imag
    #     #
    #     #             y_midsvc = q_svc_midpoint / (math.sqrt(3) * v_mid_shifted)
    #     #             matrix_a = np.array([[cmath.cosh(yp * length_cable / 2), zc * cmath.sinh(yp * length_cable / 2)],
    #     #                                  [cmath.sinh(yp * length_cable / 2) / zc, cmath.cosh(yp * length_cable / 2)]])
    #     #             matrix_b = np.array([[1, 0],
    #     #                                  [y_midsvc, 1]])
    #     #
    #     #             matrix_equivalent = np.matmul(matrix_a, matrix_b, matrix_a)
    #     #
    #     #             v_onshore_noload = v_offshore_abs / matrix_equivalent[0][0]
    #     #             i_offshore_noload = v_onshore_noload * matrix_equivalent[0][1]
    #     #             [v_max, i_max] = self.max_vi()
    #     #
    #     #         q_onshore2 = q_onshore
    #     #
    #     #         if v_cable_insulation < v_max:
    #     #             return 'The cable is not fitting (insulation voltage)'
    #     #
    #     #         if v_cable_insulation < v_onshore_noload:
    #     #             return 'The cable is not fitting (insulation voltage - no load)'
    #     #
    #     # i_offshore = i_offshore_shifted * cmath.exp(v_offshore_phase * 1j)
    #     # v_offshore = v_offshore_shifted * cmath.exp(v_offshore_phase * 1j)
    #     # i_onshore = i_onshore_shifted * cmath.exp(v_offshore_phase * 1j)
    #     # v_onshore = v_onshore_shifted * cmath.exp(v_offshore_phase * 1j)
    #     #
    #     # s_loss = 3 * v_offshore * np.conj(i_offshore) - 3 * v_onshore * np.conj(i_onshore)
    #     # p_loss = s_loss.real
    #     #
    #     # v_offshore = math.sqrt(3) * v_offshore
    #     # v_onshore = math.sqrt(3) * v_onshore
    #     # v_onshore_noload = math.sqrt(3) * v_onshore_noload
    #     # v_max = math.sqrt(3) * v_max
    #     #
    #     # outputs = {'p_loss': round(p_loss, 2),
    #     #            'q_svc_offshore': round(q_svc_offshore, 2),
    #     #            'q_svc_midpoint': round(q_svc_midpoint, 2),
    #     #            'q_onshore': round(q_onshore, 2),
    #     #            'q_onshore2': round(q_onshore2, 2),
    #     #            'v_offshore': complex(round(v_offshore.real, 2), round(v_offshore.imag, 2)),
    #     #            'v_onshore': complex(round(v_onshore.real, 2), round(v_onshore.imag, 2)),
    #     #            'i_offshore': complex(round(i_offshore.real, 2), round(i_offshore.imag, 2)),
    #     #            'i_onshore': complex(round(i_onshore.real, 2), round(i_onshore.imag, 2)),
    #     #            'i_max_cable': i_max_cable,
    #     #            'v_onshore_noload': round(abs(v_onshore_noload), 2),
    #     #            'i_offshore_noload': round(abs(i_offshore_noload), 2),
    #     #            'v_max': round(v_max, 2),
    #     #            'i_max': round(i_max, 2)}
    #     #
    #     # return outputs
    #
    # def comp_onshore(self, v_onshore, i_onshore):
    #     v_offshore_shifted = self.v_offshore_shifted
    #     v_offshore_phase = self.v_offshore_phase
    #     i_offshore_shifted = self.i_offshore_shifted
    #     length_cable = self.d
    #     n_cable = self.n_cable
    #     ip_offshore = i_offshore_shifted.real
    #     iq_offshore_shifted = i_offshore_shifted.imag
    #     i_max_cable = self.i_max_cable
    #     iq_max_cable_offshore = self.iq_max_cable_offshore
    #     zc = self.zc
    #     yp = self.yp
    #     v_cable = self.v_cable
    #     v_cable_insulation = self.v_cable_insulation
    #     key_cable = self.key_cable
    #
    #     if v_cable < v_offshore_shifted:
    #         return 'Error: the cable is not fitting (rated voltage)'
    #
    #     if ip_offshore < i_max_cable:
    #         iq_offshore0 = 0
    #         ip_onshore_shifted = ip_offshore
    #     else:
    #         return f"Error: the current limit is exceeded (p), ip_offshore = {round(ip_offshore, 3)} and i_max_cable {round(i_max_cable, 3)}."
    #
    #     # Initialisation
    #     q_svc_mid = 0
    #     y_svc_mid = 0
    #     cable_comments = 'None'
    #
    #     ###########################################################
    #     # q_wf (Compensation by the wind farm remaining capacity) #
    #     ###########################################################
    #
    #     required_compensation = 'Wind Farm'
    #
    #     if abs(iq_max_cable_offshore) < abs(iq_offshore_shifted):
    #         cable_comments = 'Notice: iq_max_cable < iq_offshore'
    #
    #     i_offshore_shifted = ip_offshore + min(iq_offshore_shifted, iq_offshore_shifted) * 1j
    #
    #     i_onshore_shifted = i_offshore_shifted * cmath.cosh(yp * length_cable) - v_offshore_shifted * cmath.sinh(
    #         yp * length_cable) / zc
    #     v_onshore_shifted = v_offshore_shifted * cmath.cosh(yp * length_cable) - i_offshore_shifted * cmath.sinh(
    #         yp * length_cable) * zc
    #
    #     v_onshore_shifted_min = v_onshore_shifted
    #     i_onshore_shifted_min = i_onshore_shifted
    #     v_onshore_shifted_max = v_onshore_shifted
    #     i_onshore_shifted_max = i_onshore_shifted
    #
    #     # ------------------------------#
    #     # range of vi_onshore with q_wf #
    #     # ------------------------------#
    #     i_onshore_shifted2 = i_onshore_shifted
    #     v_onshore_shifted2 = v_onshore_shifted
    #
    #     if abs(i_onshore_shifted) < i_max_cable:
    #         cable_comments = 'Notice: i_onshore < i_max_cable'
    #         i_offshore_shifted2 = ip_offshore - iq_offshore_shifted * 1j
    #         i_onshore_shifted2 = i_offshore_shifted2 * cmath.cosh(yp * length_cable) - v_offshore_shifted * cmath.sinh(
    #             yp * length_cable) / zc
    #
    #         while abs(i_onshore_shifted2) < i_max_cable:
    #             i_offshore_shifted2 += 0.01j
    #             i_onshore_shifted2 = i_offshore_shifted2 * cmath.cosh(
    #                 yp * length_cable) - v_offshore_shifted * cmath.sinh(yp * length_cable) / zc
    #
    #         i_onshore_shifted2 = i_offshore_shifted2 * cmath.cosh(
    #             yp * length_cable) - v_offshore_shifted * cmath.sinh(yp * length_cable) / zc
    #         v_onshore_shifted2 = v_offshore_shifted * cmath.cosh(yp * length_cable) - i_offshore_shifted2 * cmath.sinh(
    #             yp * length_cable) * zc
    #
    #     ##########################################
    #     # q_svc_offshore (offshore compensation) #
    #     ##########################################
    #
    #     if abs(i_onshore_shifted) > i_max_cable:
    #         required_compensation = 'Offshore'
    #         if abs(i_offshore_shifted) < i_max_cable:
    #             while abs(i_onshore_shifted) > i_max_cable:
    #                 i_offshore_shifted += 0.01j
    #                 i_onshore_shifted = i_offshore_shifted * cmath.cosh(
    #                     yp * length_cable) - v_offshore_shifted * cmath.sinh(yp * length_cable) / zc
    #                 if abs(i_offshore_shifted) > i_max_cable:
    #                     required_compensation = 'Mid point'
    #                     break
    #
    #     # ------------------------------#
    #     # range of q_onshore2 with q_wf #
    #     # ------------------------------#
    #
    #     # if abs(i_offshore_shifted.imag) < min(iq_max_wf, iq_max_offshore):
    #     #     i_offshore_shifted2 = ip - iq_max_wf * 1j
    #     #     i_onshore_shifted2 = i_offshore_shifted2 * cmath.cosh(yp * length_cable) - v_offshore_shifted * cmath.sinh(
    #     #         yp * length_cable) / zc
    #     #     v_onshore_shifted2 = v_offshore_shifted * cmath.cosh(yp * length_cable) - i_offshore_shifted2 * cmath.sinh(
    #     #         yp * length_cable) * zc
    #     #
    #     #     s_onshore2 = 3 * v_onshore_shifted2 * np.conj(i_onshore_shifted2)
    #     #     q_onshore2 = -s_onshore2.imag
    #
    #     ######################################
    #     # q_svc_mid (Mid point compensation) #
    #     ######################################
    #
    #     if required_compensation == 'Mid point':
    #
    #         return 'Error: Midpoint compensation is required'
    #
    #         i_offshore_shifted = complex(ip_offshore, iq_offshore_shifted)
    #
    #         v_mid_shifted = v_offshore_shifted * cmath.cosh(
    #             yp * length_cable / 2) - i_offshore_shifted * zc * cmath.sinh(yp * length_cable / 2)
    #         i_offshore_mid = i_offshore_shifted * cmath.cosh(
    #             yp * length_cable / 2) - v_offshore_shifted * cmath.sinh(yp * length_cable / 2) / zc
    #
    #         if abs(i_offshore_mid) > i_max_cable:
    #             return 'Error: the length limit is exceeded (offshore)'
    #
    #         v_mid_shifted_phase = cmath.phase(v_mid_shifted)
    #         v_mid_shifted2 = abs(v_mid_shifted)
    #         s_offshore_mid = 3 * v_mid_shifted * np.conj(i_offshore_mid)
    #         p_offshore_mid = s_offshore_mid.real
    #         ip_offshore_mid = p_offshore_mid / (3 * v_mid_shifted2)
    #
    #         ip_mid_onshore = ip_offshore_mid
    #         iq_mid_onshore = 0
    #
    #         i_mid_onshore = ip_mid_onshore + iq_mid_onshore * 1j
    #         i_onshore_shifted2 = i_mid_onshore * cmath.cosh(yp * length_cable / 2) - v_mid_shifted2 * cmath.sinh(
    #             yp * length_cable / 2) / zc
    #
    #         while abs(i_onshore_shifted2) > i_max_cable:
    #             i_mid_onshore += 0.01j
    #             i_onshore_shifted2 = i_mid_onshore * cmath.cosh(
    #                 yp * length_cable / 2) - v_mid_shifted2 * cmath.sinh(yp * length_cable / 2) / zc
    #             if abs(i_mid_onshore) > i_max_cable:
    #                 return 'Error: the length limit is exceeded (onshore)'
    #
    #         s_offshore_mid = 3 * v_mid_shifted * np.conj(i_offshore_mid)
    #         q_svc_offshore_mid = s_offshore_mid.imag
    #
    #         s_mid_onshore = 3 * v_mid_shifted * np.conj(i_mid_onshore)
    #         q_svc_mid_onshore = s_mid_onshore.imag
    #
    #         q_svc_mid = q_svc_offshore_mid + q_svc_mid_onshore
    #
    #         y_svc_mid = ((3 * v_mid_shifted ** 2) / q_svc_mid) * 1j
    #
    #         i_onshore_shifted = i_onshore_shifted2 * cmath.exp(
    #             1j * v_mid_shifted_phase)
    #         v_onshore_shifted2 = v_mid_shifted * cmath.cosh(
    #             yp * length_cable / 2) - i_mid_onshore * zc * cmath.sinh(
    #             yp * length_cable / 2)
    #         v_onshore_shifted = v_onshore_shifted2 * cmath.exp(
    #             1j * v_mid_shifted_phase)
    #
    #     #############################################################
    #     # No load voltage and current + Maximum voltage and current #
    #     #############################################################
    #
    #     if required_compensation != 'Mid point':
    #         v_onshore_noload = abs(v_offshore_shifted / cmath.cosh(yp * length_cable))
    #         i_offshore_noload = v_onshore_noload * cmath.sinh(yp * length_cable) / zc
    #         i_max = abs(i_onshore_shifted)
    #         v_max = abs(v_onshore_shifted)
    #         v_onshore_noload_limit = abs(v_offshore_shifted * cmath.cosh(yp * length_cable))
    #     else:
    #         matrix_a = np.array([[cmath.cosh(yp * length_cable / 2), zc * cmath.sinh(yp * length_cable / 2)],
    #                              [cmath.sinh(yp * length_cable / 2) / zc, cmath.cosh(yp * length_cable / 2)]])
    #         matrix_b = np.array([[1, 0], [y_svc_mid, 1]])
    #         matrix_equivalent = np.matmul(matrix_a, matrix_b, matrix_a)
    #         v_onshore_noload = v_offshore_shifted / matrix_equivalent[0][0]
    #         i_offshore_noload = v_onshore_noload * matrix_equivalent[0][1]
    #         [v_max, i_max] = [0, 0]  # self.max_vi()
    #
    #     if v_cable_insulation < v_max:
    #         return 'Error: the cable is not fitting (insulation voltage)'
    #
    #     # if v_cable_insulation < v_onshore_noload:
    #     #     return 'The cable is not fitting (insulation voltage - no load)'
    #
    #     # ##########################################
    #     # # Offshore and Onshore Power calculation #
    #     # ##########################################
    #     # s_offshore = 3 * v_offshore_shifted * np.conj(i_offshore_shifted)
    #     # q_svc_offshore = -s_offshore.imag - self.q_wf
    #     #
    #     # s_onshore = 3 * v_onshore_shifted * np.conj(i_onshore_shifted)
    #     # q_onshore = -s_onshore.imag
    #     #
    #     # s_loss = s_offshore - s_onshore
    #     # p_loss = s_loss.real
    #
    #     ###########################################
    #     # Shifting back the voltages and currents #
    #     ###########################################
    #
    #     v_offshore = v_offshore_shifted * cmath.exp(v_offshore_phase * 1j)
    #     v_onshore = v_onshore_shifted * cmath.exp(v_offshore_phase * 1j)
    #     i_offshore = i_offshore_shifted * cmath.exp(v_offshore_phase * 1j)
    #     i_onshore = i_onshore_shifted * cmath.exp(v_offshore_phase * 1j)
    #
    #     v_onshore_min = v_onshore_shifted_min * cmath.exp(v_offshore_phase * 1j)
    #     i_onshore_min = i_onshore_shifted_min * cmath.exp(v_offshore_phase * 1j)
    #     v_onshore_max = v_onshore_shifted_max * cmath.exp(v_offshore_phase * 1j)
    #     i_onshore_max = i_onshore_shifted_max * cmath.exp(v_offshore_phase * 1j)
    #
    #     v_offshore = math.sqrt(3) * v_offshore
    #     v_onshore = math.sqrt(3) * v_onshore
    #     i_offshore = n_cable * (i_offshore)
    #     i_onshore = n_cable * (i_onshore)
    #
    #     v_onshore_min = math.sqrt(3) * v_onshore_min
    #     i_onshore_min = n_cable * np.conj(i_onshore_min)
    #     v_onshore_max = math.sqrt(3) * v_onshore_max
    #     i_onshore_max = n_cable * np.conj(i_onshore_max)
    #
    #     v_onshore_noload = math.sqrt(3) * v_onshore_noload
    #     v_onshore_noload_limit = math.sqrt(3) * v_onshore_noload_limit
    #     v_max = 0  # math.sqrt(3) * v_max
    #
    #     # outputs0 = {'p_loss': round(p_loss, 2),
    #     #             'q_svc_offshore': round(q_svc_offshore, 2),
    #     #             'q_svc_mid': round(q_svc_mid, 2),
    #     #             'q_onshore': round(q_onshore, 2),
    #     #             'q_onshore2': round(q_onshore2, 2),
    #     #             'v_offshore': complex(round(v_offshore.real, 2), round(v_offshore.imag, 2)),
    #     #             'v_onshore': complex(round(v_onshore.real, 2), round(v_onshore.imag, 2)),
    #     #             'i_offshore': complex(round(i_offshore.real, 2), -round(i_offshore.imag, 2)),
    #     #             'i_onshore': complex(round(i_onshore.real, 2), -round(i_onshore.imag, 2)),
    #     #             'i_max_cable': i_max_cable,
    #     #             'v_onshore_noload': round(abs(v_onshore_noload), 2),
    #     #             'i_offshore_noload': round(abs(i_offshore_noload), 2),
    #     #             'v_max': round(v_max, 2),
    #     #             'i_max': round(i_max, 2),
    #     #             'comments': cable_comments}
    #
    #     outputs = {'v_offshore': v_offshore,
    #                'v_onshore': v_onshore,
    #                'i_offshore': i_offshore,
    #                'i_onshore': i_onshore,
    #                'v_onshore_min': v_onshore_min,
    #                'v_onshore_max': v_onshore_max,
    #                'i_onshore_min': abs(i_onshore_min),
    #                'i_onshore_max': abs(i_onshore_max),
    #                'v_onshore_noload': abs(v_onshore_noload),
    #                'i_offshore_noload': abs(i_offshore_noload),
    #                'v_onshore_noload_limit': v_onshore_noload_limit,
    #                'v_max': v_max,
    #                'i_max': i_max,
    #                'q_svc_mid': q_svc_mid,
    #                'y_svc_mid': y_svc_mid,
    #                'comments': cable_comments}
    #
    #     return outputs
    #
    #     #######################
    #     # Power calcaultation #
    #     #######################
    #
    #     # s_offshore = 3 * v_offshore_abs * np.conj(i_offshore_shifted)
    #     # s_onshore = 3 * v_onshore_shifted * np.conj(i_onshore_shifted)
    #     # q_offshore = s_offshore.imag
    #     # q_onshore = s_onshore.imag
    #     # q_onshore2 = 0
    #     # q_svc_offshore = 0
    #     # q_svc_midpoint = 0
    #     #
    #     #
    #     #
    #     # if abs(i_onshore_shifted2) > i_max_cable:
    #     #
    #     #     # q_svc_onshore (Onshore compensation)
    #     #     i_onshore_phase2 = -math.pi / 2
    #     #     i_offshore_shifted2 = 0
    #     #
    #     #     while abs(i_offshore_shifted2.real - ip) > 0.01:
    #     #         i_onshore_shifted2 = i_max_cable * cmath.exp(1j * i_onshore_phase2)
    #     #         i_offshore_shifted2 = (zc * i_onshore_shifted2 + v_offshore_shifted * cmath.sinh(yp * length_cable)) / (
    #     #                 zc * cmath.cosh(yp * length_cable))
    #     #         i_onshore_phase2 += 0.01
    #     #         if i_onshore_phase2 > math.pi / 2:
    #     #             return 'i_offshore is not found'
    #     #
    #     #     i_onshore_shifted2 = i_offshore_shifted2 * cmath.cosh(yp * length_cable) - v_offshore_shifted * cmath.sinh(
    #     #         yp * length_cable) / zc
    #     #     v_onshore_shifted2 = v_offshore_shifted * cmath.cosh(yp * length_cable) - i_offshore_shifted2 * cmath.sinh(
    #     #         yp * length_cable) * zc
    #     #
    #     #     s_offshore2 = 3 * v_offshore_shifted * np.conj(i_offshore_shifted2)
    #     #     s_onshore2 = 3 * v_onshore_shifted2 * np.conj(i_onshore_shifted2)
    #     #     q_onshore2 = s_onshore2.imag
    #     #
    #     #     if abs(i_onshore_shifted) > i_max_cable:
    #     #
    #     #         i_onshore_phase = -math.pi / 2
    #     #         i_offshore_shifted = 0
    #     #
    #     #         while abs(i_offshore_shifted.real - ip) > 0.01:
    #     #             i_onshore_shifted = i_max_cable * cmath.exp(1j * i_onshore_phase)
    #     #             i_offshore_shifted = (zc * i_onshore_shifted + v_offshore_abs * cmath.sinh(yp * length_cable)) / (
    #     #                     zc * cmath.cosh(yp * length_cable))
    #     #             i_onshore_phase += 0.01
    #     #             if i_onshore_phase > math.pi / 2:
    #     #                 return 'i_offshore is not found'
    #     #
    #     #         while abs(i_offshore_shifted.real - ip) > 0.001:
    #     #             i_onshore_shifted = i_max_cable * cmath.exp(1j * i_onshore_phase)
    #     #             i_offshore_shifted = (zc * i_onshore_shifted + v_offshore_abs * cmath.sinh(yp * length_cable)) / (
    #     #                     zc * cmath.cosh(yp * length_cable))
    #     #             i_onshore_phase += 0.001
    #     #
    #     #         i_offshore_shifted = ip + i_offshore_shifted.imag * 1j
    #     #         i_onshore_shifted = i_offshore_shifted * cmath.cosh(yp * length_cable) - v_offshore_abs * cmath.sinh(
    #     #             yp * length_cable) / zc
    #     #         v_onshore_shifted = v_offshore_abs * cmath.cosh(
    #     #             yp * length_cable) - i_offshore_shifted * zc * cmath.sinh(
    #     #             yp * length_cable)
    #     #         abs_i_onshore = abs(i_onshore_shifted)
    #     #         abs_i_offshore = abs(i_offshore_shifted)
    #     #
    #     #         s_offshore = 3 * v_offshore_abs * np.conj(i_offshore_shifted)
    #     #         s_onshore = 3 * v_onshore_shifted * np.conj(i_onshore_shifted)
    #     #         q_offshore = -s_offshore.imag
    #     #         q_onshore = s_onshore.imag
    #     #
    #     #         v_onshore_noload = v_offshore_abs / cmath.cosh(yp * length_cable)
    #     #         i_offshore_noload = v_onshore_noload * cmath.sinh(yp * length_cable) / zc
    #     #         i_max = abs(i_onshore_shifted)
    #     #         v_max = abs(v_onshore_shifted)
    #     #         iq_offshore = i_offshore_shifted.imag
    #     #
    #     #         q_svc_offshore = 0
    #     #
    #     #         if abs(i_offshore_shifted) < i_max_cable:
    #     #             if q_offshore > 0:
    #     #                 q_svc_offshore = q_offshore - q_max_wf
    #     #             else:
    #     #                 q_svc_offshore = q_offshore + q_max_wf
    #     #             q_svc_midpoint = 0
    #     #
    #     #         elif abs_i_offshore > i_max_cable:
    #     #             i_offshore_shifted = complex(ip, iq_max_offshore)
    #     #             s_offshore = 3 * v_offshore_shifted * np.conj(i_offshore_shifted)
    #     #             q_svc_offshore = -s_offshore.imag - self.q_wf
    #     #
    #     #             v_mid_shifted = v_offshore_shifted * cmath.cosh(
    #     #                 yp * length_cable / 2) - i_offshore_shifted * zc * cmath.sinh(yp * length_cable / 2)
    #     #             i_mid_offshore_shifted = i_offshore_shifted * cmath.cosh(
    #     #                 yp * length_cable / 2) - v_offshore_shifted * cmath.sinh(yp * length_cable / 2) / zc
    #     #
    #     #             v_mid_phase = cmath.phase(v_mid_shifted)
    #     #             v_mid_shifted = abs(v_mid_shifted)
    #     #             s_mid_offshore = 3 * v_mid_shifted * np.conj(i_mid_offshore_shifted)
    #     #             p_mid_offshore = s_mid_offshore.real
    #     #             ip_mid_offshore = p_mid_offshore / (3 * v_offshore_shifted)
    #     #
    #     #             i_onshore_phase = -math.pi / 2
    #     #             i_mid_onshore_shifted = 0
    #     #
    #     #             while abs(i_mid_onshore_shifted.real - ip_mid_offshore) > 0.01:
    #     #                 i_onshore_shifted = i_max_cable * cmath.exp(1j * i_onshore_phase)
    #     #                 i_mid_onshore_shifted = (zc * i_onshore_shifted + v_offshore_abs * cmath.sinh(
    #     #                     yp * length_cable / 2)) / (
    #     #                                                 zc * cmath.cosh(yp * length_cable / 2))
    #     #                 i_onshore_phase += 0.01
    #     #                 if i_onshore_phase > math.pi / 2:
    #     #                     return 'i_offshore is not found in the second loop'
    #     #
    #     #             while abs(i_mid_onshore_shifted.real - ip_mid_offshore) > 0.001:
    #     #                 i_onshore_shifted = i_max_cable * cmath.exp(1j * i_onshore_phase)
    #     #                 i_mid_onshore_shifted = (zc * i_onshore_shifted + v_offshore_abs * cmath.sinh(
    #     #                     yp * length_cable / 2)) / (
    #     #                                                 zc * cmath.cosh(yp * length_cable / 2))
    #     #                 i_onshore_phase += 0.001
    #     #
    #     #             i_mid_onshore = i_mid_onshore_shifted * cmath.exp(
    #     #                 1j * v_mid_phase)  # (zc * i_onshore_shifted + v_mid * cmath.sinh(yp * length_cable / 2)) / (zc * cmath.cosh(yp * length_cable / 2))
    #     #             v_onshore_shifted = v_mid_shifted * cmath.cosh(
    #     #                 yp * length_cable / 2) - i_mid_onshore * zc * cmath.sinh(
    #     #                 yp * length_cable / 2)
    #     #
    #     #             if abs(i_mid_offshore_shifted) > i_max_cable or abs(i_mid_onshore) > i_max_cable:
    #     #                 return 'The length limit is exceeded'
    #     #
    #     #             s_mid_onshore = 3 * v_mid_shifted * np.conj(i_mid_onshore)
    #     #             q_svc_mid_onshore = s_mid_onshore.imag
    #     #
    #     #             q_svc_mid_offshore = s_mid_offshore.imag
    #     #             q_svc_midpoint = q_svc_mid_onshore + q_svc_mid_offshore
    #     #
    #     #             s_onshore = 3 * v_onshore_shifted * np.conj(i_onshore_shifted)
    #     #             q_onshore = s_onshore.imag
    #     #
    #     #             y_midsvc = q_svc_midpoint / (math.sqrt(3) * v_mid_shifted)
    #     #             matrix_a = np.array([[cmath.cosh(yp * length_cable / 2), zc * cmath.sinh(yp * length_cable / 2)],
    #     #                                  [cmath.sinh(yp * length_cable / 2) / zc, cmath.cosh(yp * length_cable / 2)]])
    #     #             matrix_b = np.array([[1, 0],
    #     #                                  [y_midsvc, 1]])
    #     #
    #     #             matrix_equivalent = np.matmul(matrix_a, matrix_b, matrix_a)
    #     #
    #     #             v_onshore_noload = v_offshore_abs / matrix_equivalent[0][0]
    #     #             i_offshore_noload = v_onshore_noload * matrix_equivalent[0][1]
    #     #             [v_max, i_max] = self.max_vi()
    #     #
    #     #         q_onshore2 = q_onshore
    #     #
    #     #         if v_cable_insulation < v_max:
    #     #             return 'The cable is not fitting (insulation voltage)'
    #     #
    #     #         if v_cable_insulation < v_onshore_noload:
    #     #             return 'The cable is not fitting (insulation voltage - no load)'
    #     #
    #     # i_offshore = i_offshore_shifted * cmath.exp(v_offshore_phase * 1j)
    #     # v_offshore = v_offshore_shifted * cmath.exp(v_offshore_phase * 1j)
    #     # i_onshore = i_onshore_shifted * cmath.exp(v_offshore_phase * 1j)
    #     # v_onshore = v_onshore_shifted * cmath.exp(v_offshore_phase * 1j)
    #     #
    #     # s_loss = 3 * v_offshore * np.conj(i_offshore) - 3 * v_onshore * np.conj(i_onshore)
    #     # p_loss = s_loss.real
    #     #
    #     # v_offshore = math.sqrt(3) * v_offshore
    #     # v_onshore = math.sqrt(3) * v_onshore
    #     # v_onshore_noload = math.sqrt(3) * v_onshore_noload
    #     # v_max = math.sqrt(3) * v_max
    #     #
    #     # outputs = {'p_loss': round(p_loss, 2),
    #     #            'q_svc_offshore': round(q_svc_offshore, 2),
    #     #            'q_svc_midpoint': round(q_svc_midpoint, 2),
    #     #            'q_onshore': round(q_onshore, 2),
    #     #            'q_onshore2': round(q_onshore2, 2),
    #     #            'v_offshore': complex(round(v_offshore.real, 2), round(v_offshore.imag, 2)),
    #     #            'v_onshore': complex(round(v_onshore.real, 2), round(v_onshore.imag, 2)),
    #     #            'i_offshore': complex(round(i_offshore.real, 2), round(i_offshore.imag, 2)),
    #     #            'i_onshore': complex(round(i_onshore.real, 2), round(i_onshore.imag, 2)),
    #     #            'i_max_cable': i_max_cable,
    #     #            'v_onshore_noload': round(abs(v_onshore_noload), 2),
    #     #            'i_offshore_noload': round(abs(i_offshore_noload), 2),
    #     #            'v_max': round(v_max, 2),
    #     #            'i_max': round(i_max, 2)}
    #     #
    #     # return outputs


class AcTransformer:
    def __init__(self, project_data, df_tr, n_tr):
        self.s = project_data['S (MVA)']
        self.p = int(project_data['P (MW)'])
        self.f = project_data['F (Hz)']
        self.df_tr = df_tr
        capacity_mva = df_tr['Capacity (MVA)']
        v_high_nl = df_tr['High voltage offload (kV)']
        v_low_nl = df_tr['Low voltage offload (kV)']
        self.ratio_nl = (v_high_nl / v_low_nl)
        x_tr = df_tr['x (%)']
        r_tr = df_tr['r (%)']
        z_base_high = (v_high_nl ** 2) / capacity_mva
        z_base_low = (v_low_nl ** 2) / capacity_mva
        self.z_tr_high = (r_tr / 2 + x_tr / 2 * 1j) / n_tr * 0.01 * z_base_high
        self.z_tr_low = (r_tr / 2 + x_tr / 2 * 1j) / n_tr * 0.01 * z_base_low
        p_nl = df_tr['No load losses (kW)'] * 1e-3
        i_nl_percent = df_tr['No load crrent (%)']
        i_rated_low = capacity_mva / (math.sqrt(3) * v_low_nl)
        i_nl_low = i_nl_percent * i_rated_low / 100
        z_nl_abs = v_low_nl / (math.sqrt(3) * (i_nl_low + 1e-19))
        pf_nl_low = p_nl / (math.sqrt(3) * v_low_nl * (i_nl_low + 1e-19))
        z_nl = z_nl_abs * cmath.exp(math.acos(pf_nl_low) * 1j)
        self.z_nl = z_nl / n_tr

    def vi_output(self, v_in, i_in, input_side, tap_changer):
        ratio_nl = self.ratio_nl * (1 + tap_changer / 100)
        z_tr_high = self.z_tr_high
        z_tr_low = self.z_tr_low
        if math.isnan(self.z_nl.real) or math.isnan(self.z_nl.imag):
             z_nl = float('inf')
        else:
             z_nl = self.z_nl

        v_in_lg = v_in / math.sqrt(3) #* cmath.exp(math.pi * 1j / 6)

        if input_side == 'high voltage':
            i_nl = v_in_lg / (z_nl * ratio_nl ** 2)
            i_out = (i_in - i_nl) * ratio_nl
            v_out_lg = (v_in_lg - 1 * z_tr_high * i_in) / ratio_nl - z_tr_low * i_out

        elif input_side == 'low voltage':
            i_nl = v_in_lg / z_nl
            i_out = (i_in - i_nl) / ratio_nl
            v_out_lg = (v_in_lg - 1 * z_tr_low * i_in) * ratio_nl - z_tr_high * i_out

        v_out = v_out_lg * math.sqrt(3) #* cmath.exp(-math.pi * 1j / 6)

        outputs = {'v_output': v_out,
                   'i_output': i_out,
                   'comments': ''}

        return outputs

    def tap_changer(self, v_in, i_in, input_side, v_out):
        ratio_nl = self.ratio_nl
        z_tr_high = self.z_tr_high
        z_tr_low = self.z_tr_low
        v_in_lg = v_in / math.sqrt(3)
        v_out_lg = v_out / math.sqrt(3)

        if input_side == 'high voltage':
            ratio_tap = abs(v_in_lg - 2 * z_tr_high * i_in) / abs(v_out_lg)
        elif input_side == 'low voltage':
            ratio_tap = abs(v_out_lg) / abs(v_in_lg - 2 * z_tr_low * i_in)
        tap_changer_percent = (ratio_tap - ratio_nl) * 100
        v_output = v_out + 10

        while abs(abs(v_out) - abs(v_output)) > 1:
            if input_side == 'high voltage':
                tap_changer_percent -= np.sign(v_out - v_output) * 0.1
            elif input_side == 'low voltage':
                tap_changer_percent += np.sign(v_out - v_output) * 0.1
            vi_output = self.vi_output(v_in, i_in, input_side, tap_changer_percent)
            v_output = vi_output['v_output']

        return tap_changer_percent.real

    def p_loss(self, tr_list, n_tr, i_tr):
        r_tr = tr_list['Resistance 90°C (Ohms/km)'] * self.d
        p = r_tr * i_tr ** 2 * n_tr
        return p

    def q_loss(self, tr_list, v_tr, n_tr, i_tr):
        c_tr = tr_list['Capacitance (µF/km)'] * 1e-6 * self.d
        l_tr = tr_list['Inductance (mH/km)'] * 1e-3 * 2 * math.pi * self.f * self.d
        q = (- c_tr * v_tr ** 2 + l_tr * i_tr ** 2) * n_tr
        return q

    def voltage_drop(self, tr_list, n_cable, i_tr):
        r_tr = tr_list['Resistance 90°C (Ohms/km)'] * self.d
        l_tr = tr_list['Inductance (mH/km)'] * 1e-3 * 2 * math.pi * self.f * self.d
        z_tr = math.sqrt(r_tr ** 2 + l_tr ** 2)
        vd = z_tr * i_tr
        return vd

    def sizing0(self):
        df_tr = self.df_tr
        v_high = df_tr['High voltage offload (kV)']
        v_low = df_tr['Low voltage offload (kV)']
        v_insultaion_high = df_tr['Insulation high voltage (kV)']
        v_insultaion_low = df_tr['Insulation low voltage (kV)']
        tr_v = df_tr['High voltage offload (kV)']
        tr_v = df_tr['High voltage offload (kV)']
        i_tr = self.s / (self.n_tr * self.v_tr)
        s_tr_min = self.s / self.n_tr
        df_tr_v = df_tr['High voltage offload (kV)']
        df_tr_s = df_tr['Capacity (MVA)']
        amp_min = i_tr / self.n_tr
        tr_valid_list = df_tr#[(df_tr_v > self.v_tr) & (df_tr_s > self.s)]
        tr_valid_min = tr_valid_list
        tr_opt = tr_valid_list[tr_valid_list['Capacity (MVA)'] == tr_valid_min]
        return tr_valid_list, tr_opt, self.n_tr, i_tr


class AcGrid:
    def __init__(self, project_data, optimization_data):
        self.f = project_data['F (Hz)']
        #self.pf_min = grid_requirements['minimum power factor']
        #self.pf_max = grid_requirements['maximum power factor']
        self.v_grid = optimization_data['V_grid (kV)']
        self.delta_v_pos = optimization_data['DVmaxOnshore positive (%)']
        self.delta_v_neg = optimization_data['DVmaxOnshore negative (%)']
        self.scr = optimization_data['Grid SCR (MVA)']
        self.x2r = optimization_data['Grid X/R']
        abs_z_grid = (self.v_grid ** 2) / self.scr
        r_grid = math.sqrt(abs_z_grid ** 2 / (1 + self.x2r ** 2))
        x_grid = math.sqrt(abs_z_grid ** 2 - r_grid ** 2)
        self.z_grid = complex(r_grid, x_grid)
        self.v_poc_max = self.v_grid * (1 + self.delta_v_pos/100)
        self.v_poc_min = self.v_grid * (1 - self.delta_v_neg / 100)
        self.q_poc_max = optimization_data['Q_grid max (MVAR)']
        self.q_poc_min = optimization_data['Q_grid min (MVAR)']
        self.q_nl_max_limit = optimization_data['Q_grid no load max (MVAR)']
        self.q_nl_min_limit = optimization_data['Q_grid no load min (MVAR)']

    def get_values(self):
        outputs = {'v_grid': self.v_grid,
                   'delta_v_pos': self.delta_v_pos,
                   'delta_v_neg': self.delta_v_neg,
                   'scr': self.scr,
                   'x2r': self.x2r,
                   'z_grid': self.z_grid,
                   'v_poc_max': self.v_poc_max,
                   'v_poc_min': self.v_poc_min,
                   'q_poc_max': self.q_poc_max,
                   'q_poc_min': self.q_poc_min,
                   'q_nl_max_limit': self.q_nl_max_limit,
                   'q_nl_min_limit': self.q_nl_min_limit,
                   }
        return outputs

    def v_minmax(self):
        outputs = {'v_poc_max': self.v_poc_max,
                   'v_poc_min': self.v_poc_min}
        return outputs

    def check_v_poc_input_versus_limits(self, v_poc_input):
        if self.v_poc_min <= v_poc_input <= self.v_poc_max:
            return True
        else:
            return False

    def v_minmax(self):
        outputs = {'v_poc_max': self.v_poc_max,
                   'v_poc_min': self.v_poc_min}
        return outputs

    def q_poc(self, q_poc):
        q_poc_min = self.q_poc_min
        q_poc_max = self.q_poc_max
        if q_poc < q_poc_min:
            q_poc = q_poc_min
        elif q_poc > q_poc_max:
            q_poc = q_poc_max
        return q_poc

    def v_poc(self, i_poc):
        v_grid = self.v_grid
        z_grid = self.z_grid

        v_poc = v_grid + z_grid * i_poc

        outputs = v_poc

        return outputs

    def vi_poc(self, p_poc, q_poc):
        v_grid = self.v_grid
        z_grid = self.z_grid

        v_poc = v_grid
        v_poc_previous = 0

        while abs(v_poc - v_poc_previous) > 0.01:
            v_poc_previous = v_poc
            v_poc = z_grid * ((p_poc - 1j * q_poc) / np.conj(v_poc_previous)) + v_grid

        i_poc = np.conj((p_poc + 1j * q_poc) / (math.sqrt(3) * v_poc))

        outputs = {'v_poc': v_poc,
                   'i_poc': i_poc}
        return outputs

    def vi_poc2(self, p_poc, q_poc, v_poc_pu):

        z_grid = self.z_grid

        v_poc = self.v_grid * v_poc_pu
        v_grid = v_poc
        v_grid_previous = 0

        while abs(v_grid - v_grid_previous) > 0.01:
            v_grid_previous = v_poc
            v_grid = v_poc - z_grid * ((p_poc - 1j * q_poc) / np.conj(v_grid_previous))

        i_grid = np.conj((p_poc + 1j * q_poc) / (math.sqrt(3) * v_poc))
        i_poc = i_grid

        outputs = {'v_grid': v_grid,
                   'i_grid': i_grid,
                   'v_poc': v_poc,
                   'i_poc': i_poc
                   }
        return outputs


    def qi_limits(self, v_poc):
        q_poc_max = self.q_poc_max
        q_poc_min = self.q_poc_min

        iq_grid_min = q_poc_min / (math.sqrt(3) * abs(v_poc))
        iq_grid_max = q_poc_max / (math.sqrt(3) * abs(v_poc))
        # s_grid = math.sqrt(3) * v_grid * np.conj(i_grid)
        # q_grid_min = np.sign(pf_min) * s_grid.real * math.sqrt(1 / (pf_min ** 2) - 1)
        # q_grid_max = np.sign(pf_max) * s_grid.real * math.sqrt(1 / (pf_max ** 2) - 1)

        outputs = {'iq_grid_min': iq_grid_min,
                   'iq_grid_max': iq_grid_max}

        return outputs

    def voltage_drop(self, tr_list, n_cable, i_tr):
        r_tr = tr_list['Resistance 90°C (Ohms/km)'] * self.d
        l_tr = tr_list['Inductance (mH/km)'] * 1e-3 * 2 * math.pi * self.f * self.d
        z_tr = math.sqrt(r_tr ** 2 + l_tr ** 2)
        vd = z_tr * i_tr
        return vd

class Compensator:
    def __init__(self, project_data):
        self.s = project_data[0]
        self.v = project_data[1]
        self.f = project_data[2]
        self.temp_ambient = project_data[3]
        self.r_g = project_data[4]

    def p_loss(self):
        print('Compensator p_loss ')

    def q_loss(self):
        print('Compensator q_loss ')

    def voltage_drop(self):
        print('Compensator voltage_drop ')


class StatCom:
    def __init__(self, project_data):
        self.s = project_data[0]
        self.v = project_data[1]
        self.f = project_data[2]
        self.temp_ambient = project_data[3]
        self.r_g = project_data[4]

    def p_loss(self):
        print('StatCom p_loss ')

    def q_loss(self):
        print('StatCom q_loss ')

    def voltage_drop(self):
        print('StatCom voltage_drop ')


class DcCable:
    def __init__(self, project_data):
        self.s = project_data[0]
        self.v = project_data[1]
        self.f = project_data[2]
        self.temp_ambient = project_data[3]
        self.r_g = project_data[4]

    def p_loss(self):
        print('DcCable p_loss ')

    def q_loss(self):
        print('DcCable q_loss ')

    def voltage_drop(self):
        print('DcCable voltage_drop ')


class HvdcConv:
    def __init__(self, project_data):
        self.s = project_data[0]
        self.v = project_data[1]
        self.f = project_data[2]
        self.temp_ambient = project_data[3]
        self.r_g = project_data[4]

    def p_loss(self):
        print('HvdcConv p_loss ')

    def q_loss(self):
        print('HvdcConv q_loss ')

    def voltage_drop(self):
        print('HvdcConv voltage_drop ')


class Filter:
    def __init__(self, project_data, n_tr):
        self.n_tr = n_tr
        self.s = project_data[0] / self.n_tr
        self.v = project_data[1]
        self.f = project_data[2]
        self.temp_ambient = project_data[3]
        self.r_g = project_data[4]

    def p_loss(self):
        print('Filter p_loss ')

    def q_loss(self):
        print('Filter q_loss ')

    def voltage_drop(self):
        print('Filter voltage_drop ')
