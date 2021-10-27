# Jack C. Cook
# Tuesday, October 26, 2021

import ghedt
import ghedt.PLAT as PLAT
import gFunctionDatabase as gfdb
import ghedt.PLAT.pygfunction as gt
import pandas as pd


def main():
    # Borehole dimensions
    # -------------------
    H = 96.  # Borehole length (m)
    D = 2.  # Borehole buried depth (m)
    r_b = 0.075  # Borehole radius]
    B = 5.  # Borehole spacing (m)

    # Pipe dimensions
    # ---------------
    r_out = 26.67 / 1000. / 2.  # Pipe outer radius (m)
    r_in = 21.6 / 1000. / 2.  # Pipe inner radius (m)
    s = 32.3 / 1000.  # Inner-tube to inner-tube Shank spacing (m)
    epsilon = 1.0e-6  # Pipe roughness (m)

    # Pipe positions
    # --------------
    # Single U-tube [(x_in, y_in), (x_out, y_out)]
    pos = PLAT.media.Pipe.place_pipes(s, r_out, 1)
    # Single U-tube BHE object
    bhe_object = PLAT.borehole_heat_exchangers.SingleUTube

    # Thermal conductivities
    # ----------------------
    k_p = 0.4  # Pipe thermal conductivity (W/m.K)
    k_s = 2.0  # Ground thermal conductivity (W/m.K)
    k_g = 1.0  # Grout thermal conductivity (W/m.K)

    # Volumetric heat capacities
    # --------------------------
    rhoCp_p = 1542. * 1000.  # Pipe volumetric heat capacity (J/K.m3)
    rhoCp_s = 2343.493 * 1000.  # Soil volumetric heat capacity (J/K.m3)
    rhoCp_g = 3901. * 1000.  # Grout volumetric heat capacity (J/K.m3)

    # Thermal properties
    # ------------------
    # Pipe
    pipe = PLAT.media.Pipe(pos, r_in, r_out, s, epsilon, k_p, rhoCp_p)
    # Soil
    ugt = 18.3  # Undisturbed ground temperature (degrees Celsius)
    soil = PLAT.media.Soil(k_s, rhoCp_s, ugt)
    # Grout
    grout = PLAT.media.ThermalProperty(k_g, rhoCp_g)

    # Eskilson's original ln(t/ts) values
    log_time = ghedt.utilities.Eskilson_log_times()

    # Inputs related to fluid
    # -----------------------
    # Fluid properties
    mixer = 'MEG'  # Ethylene glycol mixed with water
    percent = 0.  # Percentage of ethylene glycol added in
    fluid = gt.media.Fluid(mixer=mixer, percent=percent)

    # Coordinates
    Nx = 12
    Ny = 13
    coordinates = gfdb.coordinates.rectangle(Nx, Ny, B, B)

    # Fluid properties
    V_flow_borehole = 0.2  # System volumetric flow rate (L/s)
    V_flow_system = V_flow_borehole * float(Nx * Ny)
    # Total fluid mass flow rate per borehole (kg/s)
    m_flow_borehole = V_flow_borehole / 1000. * fluid.rho

    # Define a borehole
    borehole = gt.boreholes.Borehole(H, D, r_b, x=0., y=0.)

    # Simulation start month and end month
    # --------------------------------
    # Simulation start month and end month
    start_month = 1
    n_years = 20
    end_month = n_years * 12
    # Maximum and minimum allowable fluid temperatures
    max_EFT_allowable = 35  # degrees Celsius
    min_EFT_allowable = 5  # degrees Celsius
    # Maximum and minimum allowable heights
    max_Height = 150  # in meters
    min_Height = 60  # in meters
    sim_params = PLAT.media.SimulationParameters(
        start_month, end_month, max_EFT_allowable, min_EFT_allowable,
        max_Height, min_Height)

    # Process loads from file
    # -----------------------
    # read in the csv file and convert the loads to a list of length 8760
    hourly_extraction: dict = \
        pd.read_csv('../GHE/Atlanta_Office_Building_Loads.csv').to_dict('list')
    # Take only the first column in the dictionary
    hourly_extraction_ground_loads: list = \
        hourly_extraction[list(hourly_extraction.keys())[0]]

    # Calculate a g-function for uniform inlet fluid temperature with
    # 8 unequal segments using the equivalent solver
    g_function = ghedt.gfunction.compute_live_g_function(
        B, [H], [r_b], [D], m_flow_borehole, bhe_object, log_time,
        coordinates, fluid, pipe, grout, soil)

    # --------------------------------------------------------------------------

    # Initialize the GHE object
    ghe = ghedt.ground_heat_exchangers.GHE(
        V_flow_system, B, bhe_object, fluid, borehole, pipe, grout, soil,
        g_function, sim_params, hourly_extraction_ground_loads)

    # Simulate after computing just one g-function
    max_HP_EFT, min_HP_EFT = ghe.simulate()

    print('Min EFT: {}\nMax EFT: {}'.format(min_HP_EFT, max_HP_EFT))

    # Compute a range of g-functions for interpolation
    H_values = [24., 48., 96., 192., 384.]
    r_b_values = [r_b] * len(H_values)
    D_values = [2.] * len(H_values)

    g_function = ghedt.gfunction.compute_live_g_function(
        B, H_values, r_b_values, D_values, m_flow_borehole, bhe_object,
        log_time, coordinates, fluid, pipe, grout, soil)

    # Re-Initialize the GHE object
    ghe = ghedt.ground_heat_exchangers.GHE(
        V_flow_system, B, bhe_object, fluid, borehole, pipe, grout, soil,
        g_function, sim_params, hourly_extraction_ground_loads)

    ghe.size(method='hybrid')

    print('Height of boreholes: {}'.format(ghe.bhe.b.H))


if __name__ == '__main__':
    main()