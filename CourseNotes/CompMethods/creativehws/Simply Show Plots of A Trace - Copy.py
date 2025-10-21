import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np

"""
Change the auto_fit t otake input of list i.e. auto_fit = [2, 5]
This will be a manual check to get those values then fit the rest at those voltages
"""

### Early Jan. 2025 ###
def langmuir_analysis_graphical(filepath, probe_area, auto_fit = True, isat_auto = True, two_fits = False, display_plots = "None", subplot_title = False, sanity_check = False):
    """ DISCLAIMERS
    - auto fit etemp only fits 2 volts worth of data (20 data points cuz typically ΔV = .1)
    - Isat fit is only 10 V worth of data (at the negative end)
    """

    """ Kwargs

    filepath: 
          str of full filepath to .txt format
    probe_area: 
          Exposed area of probe tip in [m^-2]
    isat_auto:
          True(default): Auto fit for isat
          False: manual terminal input for isat
    auto_fit: 
          Whether or not the T_e portion is auto fit (constant ΔV from 0 on the semilog)
          True or False -> False gives terminal input asking for fit inputs
    two_fits: 
          whether or not you do or 1 or 2 fites in the e temp region
          False -> only 1 fit
          True -> 2 fits, disables auto_fit
    display_plots: 
          Options (str)- None, subplots, sequential
          None -> no plots      subplots -> 4 subplots and 1 screen      sequential -> a handful of plots, 1 per screen , sequentially
    sanity_check: 
          True/False -> terminal output of outputs"""

    """ Output
    [floating_pot, ion_saturation_current, etemp_ev, ne, plasma_potential]
    """


    #-------------#
    # Import data #
    #-------------#

    # This pulls from txt and removes header by just looking at all rows below header
    rawdatatxtmess = pd.read_csv(filepath).iloc[20:].reset_index() #it's header until row 20
                                # if this is  r"C:etcetc", the direction of the slashes don't matter
    # This brick below: makes the one column ^^ into a list, 
    #                   then takes components of each into voltage or curr column and makes them floats

    splitted_column = rawdatatxtmess['LabVIEW Measurement\t'].apply(lambda x: x.split('\t')[- (len(x.split('\t'))- 1):])    # Going from a txt file
    goodrawdata = pd.DataFrame({})                                                                                          # 
    goodrawdata['Applied_Voltage'] = splitted_column.apply(lambda x: x[0]).astype(float)                                    # To two columns in a dataframe
    goodrawdata['Probe_Current'] = splitted_column.apply(lambda x: x[1]).astype(float)                                      # 

    raw_voltage = np.array(goodrawdata.Applied_Voltage) # NumPy array for voltage
    raw_current = np.array(goodrawdata.Probe_Current)   # NumPy array for the current


    #-------------------#
    # Misc/HouseKeeping #
    #-------------------#

    # Font Sizes for Plots #
    axis_label_fontsize = 20
    legend_fontsize = 20
    title_fontsize = 24
    tick_fontsize = 20
    linewidth = 3
    spot_size = 5

    if two_fits == True:
        auto_fit = False

    display_plots = display_plots.lower()

    # Handy #
    voltage_stepsize = abs(round(raw_voltage[0] - raw_voltage[1], 2)) # ΔV between each datapoint
    steps_per_volt = int(round(1 / voltage_stepsize, 1))              # X amount of values/steps per 1 Volt


    #--------------------#
    # Floating Potential # 
    #--------------------#

    if raw_current[-1] > 0:     # Important cause magentic fields can give an all negative current langmuir trace
        last_neg = np.where(raw_current < 0)[0][-1]     # This pull out index of last negative current value
        barely_neg_volt = raw_voltage[last_neg]         # last NEGATIVE  current value's corresponding voltage
        barely_pos_volt = raw_voltage[last_neg + 1]     # first POSITIVE current value's corresponding voltage

        floating_pot = round((barely_pos_volt + barely_neg_volt) / 2, 6)    # floating potential to 6th decimal (keithley report voltages to the 6th decimal)


    #------------------------#
    # Ion Saturation Current #
    #------------------------#

    if isat_auto == True:
        intro_skip = 2 * steps_per_volt     # skip the first 2 V of the IV trace, but in num of steps, 2V to skip just feels right. No hard science here
        end_10V_later = intro_skip + 20 * steps_per_volt    # avoid the exponential increase

        ion_portion_volt = raw_voltage[intro_skip:end_10V_later]    # Section of the IV trace that skips the most negative volts (startup)
        ion_portion_curr = raw_current[intro_skip:end_10V_later]    # and ends before the exponential bit

        isat_fit_par = np.polyfit(ion_portion_volt, ion_portion_curr, 1)    # Linear fit for Isat
        isat_regress = isat_fit_par[0] * raw_voltage + isat_fit_par[1]      # Regression for a big straight line

        ion_saturation_current = round(isat_regress[0], 11)     # What will be used as the isat is the most negative point of the regression

    if isat_auto == False:

        ### Raw IV ###  
        plt.scatter(raw_voltage, raw_current * 10**6, s = spot_size, marker = ',')  # 

        plt.title("Raw Data: IV Trace", fontsize = title_fontsize)                  # Gotta see the plot if manual fitting
        plt.xlabel("V" + r"$_{bias}$" + " [V]", fontsize = axis_label_fontsize)     # 
        plt.ylabel("I" + r"$_{probe}$" + " [μA]", fontsize = axis_label_fontsize)   # 
        plt.gca().tick_params(axis='both', labelsize = tick_fontsize)               # 
        plt.grid(True)

        plt.show()

        begin_voltage = int(input('Pick starting voltage: '))   # User input beginning voltage
        end_voltage = int(input('Pick ending voltage: '))       # Use input ending voltage
        print()

        # does not play well with enagtive values of voltage (become negative index)

        begin_index = np.where(raw_voltage < begin_voltage)
        end_index = np.where(raw_voltage < end_voltage)

        start = begin_index[0][-1]
        end = end_index[0][-1]

        isat_fit_volt = raw_voltage[start:end]  # Chopping up the voltage and current
        isat_fit_curr = raw_current[start:end]  # to make the fit prettier

        isat_fit_par = np.polyfit(isat_fit_volt, isat_fit_curr, 1)    # Linear fit for Isat
        isat_regress = isat_fit_par[0] * raw_voltage + isat_fit_par[1]      # Regression for a big straight line

        ion_saturation_current = round(isat_regress[0], 11)     # What will be used as the isat is the most negative point of the regression

        ### Isat Fit ###
        plt.scatter(raw_voltage, raw_current * 10**6, s = spot_size)
        plt.scatter(raw_voltage, isat_regress * 10**6, s = spot_size)

        plt.legend(["Raw IV", "I" + r"$_{sat}$" + " Fit"], fontsize = legend_fontsize)
        plt.title("IV Trace & I" + r"$_{sat}$" + " Fit", fontsize = title_fontsize)
        plt.xlabel("V" + r"$_{bias}$" + " [V]", fontsize = axis_label_fontsize)
        plt.ylabel("I" + r"$_{probe}$" + " [μA]", fontsize = axis_label_fontsize)
        plt.gca().tick_params(axis='both', labelsize = tick_fontsize)
        plt.grid(True)

        plt.show()



    #----------------------#
    # Semilog the IV trace #
    #----------------------#

    delta_curr = raw_current - ion_saturation_current   # Shift the current data up by the isat
    valid_semilog = np.where(delta_curr > 0)[0][5]      # Excluding the first couple points to I don't ln(neg number)

    delta_voltage = (raw_voltage - floating_pot)[valid_semilog:]  # This is because the equation to be fitted is actually w.r.t ΔV
    semilog_curr = np.log(delta_curr[valid_semilog:])             # Just log of the currents (Only care about >0 volts, but the whole for the sake of a plot)


    #--------------------------#
    # Electron Temperature Fit #
    #--------------------------#

    if auto_fit == True:
        center_0Volts = np.where(delta_voltage >= 0)[0][0]              # finding where ΔV = 0 (index)
        end_index = center_0Volts + int(round(2 * steps_per_volt, 1))   # end index for etemp fit. 2V is just a guess of good region

        etemp_volt = delta_voltage[center_0Volts:end_index] # Voltage section to be fit for electron temp
        etemp_curr = semilog_curr[center_0Volts:end_index]  # Current section to be fit for electron temp

        etemp_fit_param = np.polyfit(etemp_volt, etemp_curr, 1)                 # linear fit for my etemp region ^^
        etemp_regress = etemp_fit_param[0] * delta_voltage + etemp_fit_param[1] # Regression for the whol semilog plot

        etemp_ev = round(1 / etemp_fit_param[0], 3) # The electron temp [eV] from slope of etemp line


    if auto_fit == False:

        plt.scatter(delta_voltage, semilog_curr, s = spot_size)         # Gotta see the plot
        plt.title("Semilog IV", fontsize = title_fontsize)              # if we're doing a manual fit
        plt.xlabel("ΔV [V]", fontsize = axis_label_fontsize)            # 
        plt.ylabel("ln(ΔI) [Arb]", fontsize = axis_label_fontsize)      # 
        plt.gca().tick_params(axis='both', labelsize = tick_fontsize)   # 
        plt.grid(True)                                                  # 

        plt.show()

        begin_voltage = int(input('Pick starting voltage: '))   # User input beginning voltage
        end_voltage = int(input('Pick ending voltage: '))       # Use input ending voltage

        print()

        start = np.where(delta_voltage < begin_voltage)[0][-1]  # start index
        end = np.where(delta_voltage < end_voltage)[0][-1]      # end index

        etemp_fit_param = np.polyfit(delta_voltage[start:end], semilog_curr[start:end], deg = 1)    # getting slope and intercept

        etemp_regress = delta_voltage * etemp_fit_param[0] + etemp_fit_param[1]     # y = m x + b for a regression

        etemp_ev = round(1 / etemp_fit_param[0], 3) # The electron temp [eV] from slope of etemp line


        plot_list = [semilog_curr, etemp_regress]   # Future proofing plot production if two fits is true
        legend_list = ["Semilog", "T_e Fit"]        # 

        if two_fits == True:    # Same as ^^ but again
            begin_voltage = int(input('Pick starting voltage: '))   # User input beginning voltage
            end_voltage = int(input('Pick ending voltage: '))       # Use input ending voltage

            print()

            start = np.where(delta_voltage < begin_voltage)[0][-1]  # start index
            end = np.where(delta_voltage < end_voltage)[0][-1]      # end index

            etemp_fit_param2 = np.polyfit(delta_voltage[start:end], semilog_curr[start:end], deg = 1)   # getting slope and intercept
            etemp_regress2 = delta_voltage * etemp_fit_param2[0] + etemp_fit_param2[1]  # y = m x + b for a regression

            etemp_ev2 = round(1 / etemp_fit_param2[0], 3) # The electron temp [eV] from slope of etemp line

            plot_list.append(etemp_regress2)        # Future proofing plotting
            legend_list.append("Second T_e Fit")    # 

        for series in plot_list:
            plt.scatter(delta_voltage, series, s = spot_size)
            
        plt.title("Semilog Temp Fit", fontsize = title_fontsize)              # It is nice to visualise the fit 
        plt.legend(legend_list, fontsize = legend_fontsize)             # you just manually entered
        plt.xlabel("ΔV [V]", fontsize = axis_label_fontsize)            # 
        plt.ylabel("ln(ΔI) [Arb]", fontsize = axis_label_fontsize)      # 
        plt.gca().tick_params(axis='both', labelsize = tick_fontsize)   # 
        plt.grid(True)

        plt.show()

    #------#
    # Math #
    #------#

    q_e = 1.6021766E-19 # [C]
    k_b = 1.3806482E-23 # [J/K]
    kb_ev = 8.6173324E-5# [eV/K]
    m_argon = 6.62E-26  # [kg]

    etemp_kelvin = etemp_ev / kb_ev   # Getting T_e [K] to make equation below kind

    firstbit = abs(ion_saturation_current) / (q_e * probe_area * np.exp(- 1/ 2) )
    secondbit = np.sqrt(m_argon / (k_b * etemp_kelvin))  # 
    ne = "{:e}".format(firstbit * secondbit)       # this is  m^-3. the e stuff is to make it scientific notation

    plasma_potential = round(floating_pot + 5.4 * etemp_ev, 6)   # I'm assuming this is in Volts. Eq'n courtesy of Saikat



    #-------#
    # Plots #
    #-------#

    ### This is to plot raw IV trace ###
    if display_plots == "sequential": # Sequential

        ### Raw IV ###
        plt.scatter(raw_voltage, raw_current * 10**6, s = spot_size, marker = ',')

        plt.title("Raw Data: IV Trace", fontsize = title_fontsize)
        plt.xlabel("V" + r"$_{bias}$" + " [V]", fontsize = axis_label_fontsize)
        plt.ylabel("I" + r"$_{probe}$" + " [μA]", fontsize = axis_label_fontsize)
        plt.gca().tick_params(axis='both', labelsize = tick_fontsize)
        plt.grid(True)

        plt.show()

        ### Isat Fit ###
        plt.scatter(raw_voltage, raw_current * 10**6, s = spot_size)
        plt.scatter(raw_voltage, isat_regress * 10**6, s = spot_size)

        plt.legend(["Raw IV", "I" + r"$_{sat}$" + " Fit"], fontsize = legend_fontsize)
        plt.title("IV Trace & I" + r"$_{sat}$" + " Fit", fontsize = title_fontsize)
        plt.xlabel("V" + r"$_{bias}$" + " [V]", fontsize = axis_label_fontsize)
        plt.ylabel("I" + r"$_{probe}$" + " [μA]", fontsize = axis_label_fontsize)
        plt.gca().tick_params(axis='both', labelsize = tick_fontsize)
        plt.grid(True)

        plt.show()

        ### Semilog ###
        plt.scatter(delta_voltage, semilog_curr, s = spot_size)

        plt.title("Semilog", fontsize = title_fontsize)
        plt.xlabel("ΔV [V]", fontsize = axis_label_fontsize)
        plt.ylabel("ln(ΔI) [Arb]", fontsize = axis_label_fontsize)
        plt.gca().tick_params(axis='both', labelsize = tick_fontsize)
        plt.grid(True)
        
        plt.show()

        ### Etemp fits ###
        plot_list = [semilog_curr, etemp_regress]
        legend_list = ["Semilog IV", "T" + r"$_{e}$" + " Fit"]

        if two_fits == True:
            plot_list.append(etemp_regress2)
            legend_list.append("T" + r"$_{e}$" + "2 Fit")

        for series in plot_list:
            plt.scatter(delta_voltage, series, s = spot_size)
            
        plt.title("Semilog Temp Fit", fontsize = title_fontsize)              # It is nice to visualise the fit 
        plt.legend(legend_list, fontsize = legend_fontsize * .7)             # you just manually entered
        plt.xlabel("ΔV [V]", fontsize = axis_label_fontsize)            # 
        plt.ylabel("ln(ΔI) [Arb]", fontsize = axis_label_fontsize)      # 
        plt.gca().tick_params(axis='both', labelsize = tick_fontsize)   # 
        plt.grid(True)

        plt.show()

    if display_plots == "subplots":

        if subplot_title != False:
            plt.suptitle(subplot_title, fontsize = 1.5 * title_fontsize)

        ### Raw IV ###
        plt.subplot(221)

        plt.scatter(raw_voltage, raw_current * 10**6, s = spot_size, marker = ',')
        
        plt.title("Raw Data: IV Trace", fontsize = title_fontsize)
        plt.xlabel("V" + r"$_{bias}$" + " [V]", fontsize = axis_label_fontsize)
        plt.ylabel("I" + r"$_{probe}$" + " [μA]", fontsize = axis_label_fontsize)
        plt.gca().tick_params(axis='both', labelsize = tick_fontsize)
        plt.grid(True)

        ### Raw IV & Isat Regression ###
        plt.subplot(223)

        plt.scatter(raw_voltage, raw_current * 10**6, s = spot_size)
        plt.scatter(raw_voltage, isat_regress * 10**6, s = spot_size)

        plt.legend(["Raw IV", "I" + r"$_{sat}$" + " Fit"], fontsize = legend_fontsize)
        plt.title("IV Trace & I" + r"$_{sat}$" + " Fit", fontsize = title_fontsize)
        plt.xlabel("V" + r"$_{bias}$" + " [V]", fontsize = axis_label_fontsize)
        plt.ylabel("I" + r"$_{probe}$" + " [μA]", fontsize = axis_label_fontsize)
        plt.gca().tick_params(axis='both', labelsize = tick_fontsize)
        plt.grid(True)

        ### Semilog ###
        plt.subplot(222)

        plt.scatter(delta_voltage, semilog_curr, s = spot_size)

        plt.title("Semilog", fontsize = title_fontsize)
        plt.xlabel("ΔV [V]", fontsize = axis_label_fontsize)
        plt.ylabel("ln(ΔI) [Arb]", fontsize = axis_label_fontsize)
        plt.gca().tick_params(axis='both', labelsize = tick_fontsize)
        plt.grid(True)
        
        ### Semilog and Etemp fit ###
        plt.subplot(224)

        plot_list = [semilog_curr, etemp_regress]
        legend_list = ["Semilog IV", "T" + r"$_{e}$" + " Fit"]

        if two_fits == True:
            plot_list.append(etemp_regress2)
            legend_list.append("T" + r"$_{e}$" + "2 Fit")
            
        for series in plot_list:
            plt.scatter(delta_voltage, series, s = spot_size)
            
        plt.title("Semilog IV", fontsize = title_fontsize)              # It is nice to visualise the fit 
        plt.legend(legend_list, fontsize = legend_fontsize)             # you just manually entered
        plt.xlabel("ΔV [V]", fontsize = axis_label_fontsize)            # 
        plt.ylabel("ln(ΔI) [Arb]", fontsize = axis_label_fontsize)      # 
        plt.gca().tick_params(axis='both', labelsize = tick_fontsize)   # 
        plt.grid(True)

        plt.tight_layout(pad = .01)
        plt.show()


    #--------------#
    # Sanity Check #
    #--------------#

    if sanity_check == True:
        print()
        print(filepath)
        print(f"  Floating Potential [V]: {floating_pot}", 
                "\n", f" Isat [A]: {ion_saturation_current}", 
                "\n", f" Electron Temp [eV]: {etemp_ev}", 
                "\n", f" Electron Density [m^-3]: {ne}", 
                "\n", f" Plasma Pot [V]: {plasma_potential}")

        if two_fits == True:
            print(f"  Electron Temp 2[eV]: {etemp_ev2}",)
        print()


    #------------------#
    # Values to Return #
    #------------------#

    relevant_vals = [floating_pot, ion_saturation_current, etemp_ev, ne, plasma_potential]

    if two_fits == True:
        relevant_vals.append(etemp_ev2)
    return relevant_vals

# Dirty Single Tip Feb 7th #
full_filepath = r"CompMethods/creativehws/LP_r_0_p_85_B_0_Num_1.txt"

# Dimensions and plot function # 
probe_rad = .00318 / 2  # m
probe_len = .0045 # m

# probe_rad = .002 / 2
# probe_len = .002

exposed_area = np.pi * probe_rad **2 + 2 * np.pi * probe_rad * probe_len

langmuir_analysis_graphical(full_filepath, probe_area=exposed_area, display_plots= 'subplots', sanity_check= True, auto_fit= True)