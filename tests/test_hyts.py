import unittest
import ghedt.peak_load_analysis_tool as plat
import os

import numpy as np
import pandas as pd
import pygfunction as gt

"""
Unit test the Hybrid time step object.
"""

TESTDATA_FILENAME = os.path.join(os.path.dirname(__file__),
                                 'Atlanta_Office_Building_Loads.csv')


class TestHYTS(unittest.TestCase):

    def setUp(self) -> None:
        hourly_extraction: dict = \
            pd.read_csv(TESTDATA_FILENAME).to_dict('list')
        # Take only the first column in the dictionary
        self.hourly_extraction_ground_loads: list = \
            hourly_extraction[list(hourly_extraction.keys())[0]]

    def test_hybrid_time_step(self):
        hourly_rejection_loads, hourly_extraction_loads = \
            plat.ground_loads.HybridLoad.split_heat_and_cool(
                self.hourly_extraction_ground_loads)

        # Borehole dimensions
        H = 100.  # Borehole length (m)
        D = 2.  # Borehole buried depth (m)
        r_b = 150. / 1000. / 2.  # Borehole radius

        # Pipe dimensions
        r_out = 26.67 / 1000. / 2.  # Pipe outer radius (m)
        r_in = 21.6 / 1000. / 2.  # Pipe inner radius (m)
        s = 32.3 / 1000.  # Inner-tube to inner-tube Shank spacing (m)
        epsilon = 1.0e-6  # Pipe roughness (m)

        # Pipe positions
        # Single U-tube [(x_in, y_in), (x_out, y_out)]
        pos = plat.media.Pipe.place_pipes(s, r_out, 1)

        # Thermal conductivities
        k_p = 0.4  # Pipe thermal conductivity (W/m.K)
        k_s = 2.0  # Ground thermal conductivity (W/m.K)
        k_g = 1.0  # Grout thermal conductivity (W/m.K)

        # Volumetric heat capacities
        rhoCp_p = 1542. * 1000.  # Pipe volumetric heat capacity (J/K.m3)
        rhoCp_s = 2343.493 * 1000.  # Soil volumetric heat capacity (J/K.m3)
        rhoCp_g = 3901. * 1000.  # Grout volumetric heat capacity (J/K.m3)

        # Thermal properties
        # Pipe
        pipe = plat.media.Pipe(pos, r_in, r_out, s, epsilon, k_p, rhoCp_p)
        # Soil
        ugt = 18.3  # Undisturbed ground temperature (degrees Celsius)
        soil = plat.media.Soil(k_s, rhoCp_s, ugt)
        # Grout
        grout = plat.media.Grout(k_g, rhoCp_g)

        # Fluid properties
        fluid = gt.media.Fluid('MEG', 0.)
        V_flow_borehole = 0.2  # Volumetric flow rate per borehole (L/s)
        # Total fluid mass flow rate per borehole (kg/s)
        m_flow_borehole = V_flow_borehole / 1000. * fluid.rho

        # Define a borehole
        borehole = gt.boreholes.Borehole(H, D, r_b, x=0., y=0.)

        single_u_tube = plat.borehole_heat_exchangers.SingleUTube(
            m_flow_borehole, fluid, borehole, pipe, grout, soil)

        radial_numerical = plat.radial_numerical_borehole.RadialNumericalBH(
            single_u_tube)
        radial_numerical.calc_sts_g_functions(single_u_tube)

        # Simulation start month and end month
        # --------------------------------
        # Simulation start month and end month
        start_month = 1
        n_years = 2
        end_month = n_years * 12
        # Maximum and minimum allowable fluid temperatures
        max_EFT_allowable = 35  # degrees Celsius
        min_EFT_allowable = 5  # degrees Celsius
        # Maximum and minimum allowable heights
        max_Height = 150  # in meters
        min_Height = 60  # in meters
        sim_params = plat.media.SimulationParameters(
            start_month, end_month, max_EFT_allowable, min_EFT_allowable,
            max_Height, min_Height)

        hybrid_load = plat.ground_loads.HybridLoad(
            hourly_rejection_loads, hourly_extraction_loads, single_u_tube,
            radial_numerical, sim_params)

        # Rejection
        monthly_cl = [0.0, 2985.279414548951, 7348.185264939054,
                      24242.4729437719, 38355.355848082654, 66819.12789289287,
                      92585.00987606372, 101198.14499566889, 113087.48169935631,
                      72806.64725332448, 36683.97513419995, 14823.049038232912,
                      6905.23604019025]
        monthly_avg_cl = [0.0, 4.012472331382999, 10.934799501397402,
                          32.583969010446104, 53.271327566781466,
                          89.81065577001729, 128.59029149453295,
                          136.0190120909528, 151.99930335934988,
                          101.12034340739511, 49.30641819112897,
                          20.58756810865682, 9.281231236814852]
        monthly_peak_cl = [0.0, 101.65464560000001, 168.583588315556,
                           249.85718963555598, 254.22101848, 337.1167684,
                           357.063188222222, 372.182726844444, 379.231087066667,
                           369.814712266667, 319.93987057777804, 187.06670408,
                           164.354865253333]
        monthly_peak_cl_day = [0.0, 31, 28, 30, 17, 30, 27, 6, 7, 7, 10, 2, 12]
        monthly_peak_cl_duration = [0, 3.1631360357235083, 4.536722758635145,
                                    5.275275024136286, 5.884801974995289,
                                    6.566516432677471, 6.2120005864026115,
                                    9.064031315164565, 7.022117117872036,
                                    5.294648619721669, 6.014175916357069,
                                    4.304879113701285, 3.887673847693884]

        assert np.allclose(monthly_avg_cl,
                           hybrid_load.rejection.monthly_average)
        assert np.allclose(monthly_cl, hybrid_load.rejection.monthly_total)
        assert np.allclose(monthly_peak_cl, hybrid_load.rejection.monthly_peak)
        assert np.allclose(monthly_peak_cl_day,
                           hybrid_load.rejection.monthly_peak_day)
        assert np.allclose(
            monthly_peak_cl_duration,
            hybrid_load.rejection.monthly_peak_duration)

        # Extraction
        monthly_avg_hl = [0, 15.459590754237224, 9.060390078577253,
                          0.6989300720989248, 0.2025507049590962,
                          0.0005207361557825568, 0.0, 0.0, 0.0, 0.0,
                          0.26110093017241937, 1.923909858140038,
                          6.383868764643993]
        monthly_hl = [0.0, 11501.935521152494, 6088.582132803914,
                      520.0039736416, 145.83650757054926, 0.3874276999022223,
                      0.0, 0.0, 0.0, 0.0, 194.25909204828, 1385.2150978608274,
                      4749.598360895131]
        monthly_peak_hl = [0, 175.98041008888902, 180.00339723111102,
                           44.3234887022222, 45.4052141733333,
                           0.300107748488889, 0.0, 0.0, 0.0, 0.0,
                           44.4398379333333, 111.84531427555599,
                           149.047468048889]
        monthly_peak_hl_day = [0.0, 9, 13, 23, 6, 10, 0, 0, 0, 0, 26, 22, 18]
        monthly_peak_hl_duration = [0, 1.8092592525113897, 1.7442685745218776,
                                    1.0020826443712465, 1.0, 1.0, 1e-06, 1e-06,
                                    1e-06, 1e-06, 1.0, 1.2259173461572768,
                                    1.380120484557252]

        assert np.allclose(monthly_avg_hl, hybrid_load.extraction.monthly_average)
        assert np.allclose(monthly_hl, hybrid_load.extraction.monthly_total)
        assert np.allclose(monthly_peak_hl, hybrid_load.extraction.monthly_peak)
        assert np.allclose(monthly_peak_hl_day, hybrid_load.extraction.monthly_peak_day)
        assert np.allclose(
            monthly_peak_hl_duration, hybrid_load.extraction.monthly_peak_duration)

        _load = list(hybrid_load.load)
        _hour = list(hybrid_load.hour)
        _sfload = list(hybrid_load.sfload)

        load = [0.0, 0.0, -11.528404811628684, -175.98041008888902, -11.528404811628684, 101.65464560000001, -11.528404811628684, 1.214867517465956, -180.00339723111102, 1.214867517465956, 168.583588315556, 1.214867517465956, 30.42989070508177, -44.3234887022222, 30.42989070508177, 249.85718963555598, 30.42989070508177, 51.54690897747506, -45.4052141733333, 51.54690897747506, 254.22101848, 51.54690897747506, 87.72734974048161, -0.300107748488889, 87.72734974048161, 337.1167684, 87.72734974048161, 126.60192289907663, 357.063188222222, 126.60192289907663, 133.10638396532786, 372.182726844444, 133.10638396532786, 149.83417997461356, 379.231087066667, 149.83417997461356, 99.1298137893778, 369.814712266667, 99.1298137893778, 46.961529088849524, 319.93987057777804, 46.961529088849524, -44.4398379333333, 46.961529088849524, 17.872915222508293, 187.06670408, 17.872915222508293, -111.84531427555599, 17.872915222508293, 2.3315391132531786, 164.354865253333, 2.3315391132531786, -149.047468048889, 2.3315391132531786, -11.528404811628684, -175.98041008888902, -11.528404811628684, 101.65464560000001, -11.528404811628684, 1.214867517465956, -180.00339723111102, 1.214867517465956, 168.583588315556, 1.214867517465956, 30.42989070508177, -44.3234887022222, 30.42989070508177, 249.85718963555598, 30.42989070508177, 51.54690897747506, -45.4052141733333, 51.54690897747506, 254.22101848, 51.54690897747506, 87.72734974048161, -0.300107748488889, 87.72734974048161, 337.1167684, 87.72734974048161, 126.60192289907663, 357.063188222222, 126.60192289907663, 133.10638396532786, 372.182726844444, 133.10638396532786, 149.83417997461356, 379.231087066667, 149.83417997461356, 99.1298137893778, 369.814712266667, 99.1298137893778, 46.961529088849524, 319.93987057777804, 46.961529088849524, -44.4398379333333, 46.961529088849524, 17.872915222508293, 187.06670408, 17.872915222508293, -111.84531427555599]
        hour = [0.0, 0.0, 204.09537037374432, 205.9046296262557, 731.4184319821383, 734.5815680178617, 744.0, 1044.127865712739, 1045.872134287261, 1402.7316386206824, 1407.2683613793174, 1416.0, 1956.4989586778145, 1957.5010413221858, 2122.362362487932, 2127.637637512068, 2160.0, 2292.5, 2293.5, 2554.0575990125026, 2559.942400987498, 2880.0, 3108.5, 3109.5, 3585.7167417836613, 3592.2832582163387, 3624.0, 4257.893999706798, 4264.106000293201, 4344.0, 4472.467984342417, 4481.532015657582, 5088.0, 5241.488941441064, 5248.5110585589355, 5832.0, 5986.352675690139, 5991.647324309861, 6552.0, 6777.992912041822, 6784.007087958179, 7164.5, 7165.5, 7296.0, 7330.84756044315, 7335.152439556851, 7812.3870413269215, 7813.6129586730785, 8016.0, 8291.056163076153, 8294.943836923847, 8436.30993975772, 8437.690060242277, 8760.0, 8964.095370373745, 8965.904629626257, 9491.418431982138, 9494.581568017862, 9504.0, 9804.127865712739, 9805.872134287261, 10162.731638620682, 10167.268361379318, 10176.0, 10716.498958677814, 10717.501041322184, 10882.362362487931, 10887.637637512067, 10920.0, 11052.5, 11053.5, 11314.057599012502, 11319.942400987497, 11640.0, 11868.5, 11869.5, 12345.71674178366, 12352.28325821634, 12384.0, 13017.8939997068, 13024.106000293203, 13104.0, 13232.467984342418, 13241.532015657584, 13848.0, 14001.488941441065, 14008.511058558937, 14592.0, 14746.35267569014, 14751.64732430986, 15312.0, 15537.992912041822, 15544.007087958178, 15924.5, 15925.5, 16056.0, 16090.84756044315, 16095.15243955685, 16572.38704132692, 16573.61295867308]
        sfload = [0.0, 0.0, -11.528404811628684, -164.45200527726033, 164.45200527726033, 113.1830504116287, -113.1830504116287, 12.74327232909464, -181.21826474857698, 181.21826474857698, 167.36872079809004, -167.36872079809004, 29.215023187615817, -74.75337940730397, 74.75337940730397, 219.42729893047422, -219.42729893047422, 21.11701827239329, -96.95212315080836, 96.95212315080836, 202.67410950252494, -202.67410950252494, 36.18044076300655, -88.0274574889705, 88.0274574889705, 249.3894186595184, -249.3894186595184, 38.87457315859501, 230.46126532314537, -230.46126532314537, 6.504461066251238, 239.07634287911614, -239.07634287911614, 16.727796009285697, 229.39690709205343, -229.39690709205343, -50.704366185235756, 270.6848984772892, -270.6848984772892, -52.16828470052828, 272.9783414889285, -272.9783414889285, -91.40136702218282, 91.40136702218282, -29.08861386634123, 169.1937888574917, -169.1937888574917, -129.7182294980643, 129.7182294980643, -15.541376109255115, 162.02332614007983, -162.02332614007983, -151.37900716214216, 151.37900716214216, -13.859943924881863, -164.45200527726033, 164.45200527726033, 113.1830504116287, -113.1830504116287, 12.74327232909464, -181.21826474857698, 181.21826474857698, 167.36872079809004, -167.36872079809004, 29.215023187615817, -74.75337940730397, 74.75337940730397, 219.42729893047422, -219.42729893047422, 21.11701827239329, -96.95212315080836, 96.95212315080836, 202.67410950252494, -202.67410950252494, 36.18044076300655, -88.0274574889705, 88.0274574889705, 249.3894186595184, -249.3894186595184, 38.87457315859501, 230.46126532314537, -230.46126532314537, 6.504461066251238, 239.07634287911614, -239.07634287911614, 16.727796009285697, 229.39690709205343, -229.39690709205343, -50.704366185235756, 270.6848984772892, -270.6848984772892, -52.16828470052828, 272.9783414889285, -272.9783414889285, -91.40136702218282, 91.40136702218282, -29.08861386634123, 169.1937888574917, -169.1937888574917, -129.7182294980643]
        
        assert np.allclose(load, list(hybrid_load.load)[0:len(load)])
        assert np.allclose(hour, list(hybrid_load.hour)[0:len(hour)])
        assert np.allclose(sfload, list(hybrid_load.sfload)[0:len(sfload)])
