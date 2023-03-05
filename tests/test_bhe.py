import unittest
import ghedt.peak_load_analysis_tool as plat
import pygfunction as gt

"""
Test the single U-tube, double U-tube and coaxial borehole heat exchangers.
"""


class TestBHE(unittest.TestCase):

    def setUp(self) -> None:
        # Borehole dimensions
        self.H = 100.  # Borehole length (m)
        self.D = 2.  # Borehole buried depth (m)
        self.r_b = 150. / 1000. / 2.  # Borehole radius

        # Thermal conductivities
        self.k_p = 0.4  # Pipe thermal conductivity (W/m.K)
        k_s = 2.0  # Ground thermal conductivity (W/m.K)
        k_g = 1.0  # Grout thermal conductivity (W/m.K)

        # Volumetric heat capacities
        self.rhoCp_p = 1542. * 1000.  # Pipe volumetric heat capacity (J/K.m3)
        rhoCp_s = 2343.493 * 1000.  # Soil volumetric heat capacity (J/K.m3)
        rhoCp_g = 3901. * 1000.  # Grout volumetric heat capacity (J/K.m3)

        # Thermal properties
        # Soil
        ugt = 18.3  # Undisturbed ground temperature (degrees Celsius)
        self.soil = plat.media.Soil(k_s, rhoCp_s, ugt)
        # Grout
        self.grout = plat.media.Grout(k_g, rhoCp_g)

        # Fluid properties
        self.fluid = gt.media.Fluid('MEG', 0.)
        V_flow_borehole = 0.2  # Volumetric flow rate per borehole (L/s)
        # Total fluid mass flow rate per borehole (kg/s)
        self.m_flow_borehole = V_flow_borehole / 1000. * self.fluid.rho

    def test_single_u_tube(self):
        # Pipe dimensions
        r_out = 26.67 / 1000. / 2.  # Pipe outer radius (m)
        r_in = 21.6 / 1000. / 2.  # Pipe inner radius (m)
        s = 32.3 / 1000.  # Inner-tube to inner-tube Shank spacing (m)
        epsilon = 1.0e-6  # Pipe roughness (m)

        # Pipe positions
        # Single U-tube [(x_in, y_in), (x_out, y_out)]
        pos = plat.media.Pipe.place_pipes(s, r_out, 1)

        # Thermal properties
        # Pipe
        pipe = plat.media.Pipe(
            pos, r_in, r_out, s, epsilon, self.k_p, self.rhoCp_p)

        # Define a borehole
        borehole = gt.boreholes.Borehole(self.H, self.D, self.r_b, x=0., y=0.)

        single_u_tube = plat.borehole_heat_exchangers.SingleUTube(
            self.m_flow_borehole, self.fluid, borehole, pipe, self.grout, self.soil)

        Re = plat.borehole_heat_exchangers.compute_Reynolds(
            single_u_tube.m_flow_pipe, r_in, epsilon, self.fluid)

        # Reynolds number
        self.assertAlmostEqual(Re, 11667.361867903752, places=8)
        # Pipe resistance (K/(W/m))
        self.assertAlmostEqual(single_u_tube.R_p, 0.08389296717975074, places=8)
        # Convection coefficient (W/(m2.K))
        self.assertAlmostEqual(single_u_tube.h_f, 2522.2943235234375, places=8)
        # Convective resistance (K/(W/m))
        self.assertAlmostEqual(
            single_u_tube.R_fp, 0.08973549264054603, places=8)
        # Borehole thermal resistance ((m.K)/W)
        self.assertAlmostEqual(
            single_u_tube.compute_effective_borehole_resistance(),
            0.20730076327484714, places=8)

    def test_double_u_tube(self):
        # Pipe dimensions
        r_out = 26.67 / 1000. / 2.  # Pipe outer radius (m)
        r_in = 21.6 / 1000. / 2.  # Pipe inner radius (m)
        s = 32.3 / 1000.  # Inner-tube to inner-tube Shank spacing (m)
        epsilon = 1.0e-6  # Pipe roughness (m)

        # Pipe positions
        # Double U-tube [(x_in, y_in), (x_out, y_out), (x_in, y_in), (x_out, y_out)]
        pos = plat.media.Pipe.place_pipes(s, r_out, 2)

        # Thermal properties
        # Pipe
        pipe = plat.media.Pipe(pos, r_in, r_out, s, epsilon, self.k_p, self.rhoCp_p)

        # Define a borehole
        borehole = gt.boreholes.Borehole(self.H, self.D, self.r_b, x=0., y=0.)

        double_u_tube_parallel = plat.borehole_heat_exchangers.MultipleUTube(
            self.m_flow_borehole, self.fluid, borehole, pipe, self.grout,
            self.soil, config='parallel')

        # Intermediate variables
        Re = plat.borehole_heat_exchangers.compute_Reynolds(
            double_u_tube_parallel.m_flow_pipe, r_in, epsilon, self.fluid)

        # Reynolds number
        self.assertAlmostEqual(Re, 5833.680933951876, places=8)
        # Pipe resistance (K/(W/m))
        self.assertAlmostEqual(
            double_u_tube_parallel.R_p, 0.08389296717975074, places=8)
        # Convection coefficient (W/(m2.K))
        self.assertAlmostEqual(
            double_u_tube_parallel.h_f, 1292.3277737336607, places=8)
        # Convective resistance (K/(W/m))
        self.assertAlmostEqual(
            double_u_tube_parallel.R_fp, 0.09529608727383472, places=8)
        # Borehole thermal resistance ((m.K)/W)
        self.assertAlmostEqual(
            double_u_tube_parallel.compute_effective_borehole_resistance(),
            0.15970364062604958, places=8)

        double_u_tube_series = plat.borehole_heat_exchangers.MultipleUTube(
            self.m_flow_borehole, self.fluid, borehole, pipe, self.grout,
            self.soil, config='series')

        # Intermediate variables
        Re = plat.borehole_heat_exchangers.compute_Reynolds(
            double_u_tube_series.m_flow_pipe, r_in, epsilon, self.fluid)

        # Reynolds number
        self.assertAlmostEqual(Re, 11667.361867903752, places=8)
        # Pipe resistance (K/(W/m))
        self.assertAlmostEqual(double_u_tube_series.R_p, 0.08389296717975074, places=8)
        # Convection coefficient (W/(m2.K))
        self.assertAlmostEqual(double_u_tube_series.h_f, 2522.2943235234375, places=8)
        # Convective resistance (K/(W/m))
        self.assertAlmostEqual(
            double_u_tube_series.R_fp, 0.08973549264054603, places=8)
        # Borehole thermal resistance ((m.K)/W)
        self.assertAlmostEqual(
            double_u_tube_series.compute_effective_borehole_resistance(),
            0.16239288378526318, places=8)

    def test_coaxial(self):
        # Pipe dimensions
        # Inner pipe radii
        r_in_in = 44.2 / 1000. / 2.
        r_in_out = 50. / 1000. / 2.
        # Outer pipe radii
        r_out_in = 97.4 / 1000. / 2.
        r_out_out = 110. / 1000. / 2.
        # Pipe radii
        # Note: This convention is different from pygfunction
        r_inner = [r_in_in,
                   r_in_out]  # The radii of the inner pipe from in to out
        r_outer = [r_out_in,
                   r_out_out]  # The radii of the outer pipe from in to out
        epsilon = 1.0e-6  # Pipe roughness (m)

        # Pipe positioning
        pos = (0, 0)
        s = 0

        # Thermal properties
        k_p = [0.4, 0.4]  # Inner and outer pipe thermal conductivity (W/m.K)

        # Volumetric heat capacities
        rhoCp_p = 1542. * 1000.  # Pipe volumetric heat capacity (J/K.m3)

        # Thermal properties
        # Pipe
        pipe = plat.media.Pipe(pos, r_inner, r_outer, s, epsilon, k_p, rhoCp_p)

        # Define a borehole
        borehole = gt.boreholes.Borehole(self.H, self.D, self.r_b, x=0., y=0.)

        coaxial = plat.borehole_heat_exchangers.CoaxialPipe(
            self.m_flow_borehole, self.fluid, borehole, pipe, self.grout,
            self.soil)

        Re = plat.borehole_heat_exchangers.compute_Reynolds(
            coaxial.m_flow_pipe, coaxial.pipe.r_out[1], epsilon, self.fluid)

        # Reynolds number
        self.assertAlmostEqual(Re, 2291.0456031520093, places=8)
        # Convection coefficient (W/(m2.K))
        self.assertAlmostEqual(coaxial.h_fluid_a_out, 57.129718138096315, places=8)
        # Convective resistance (K/(W/m))
        self.assertAlmostEqual(
            coaxial.R_fp, 0.10560900480844831, places=8)
        # Borehole thermal resistance ((m.K)/W)
        self.assertAlmostEqual(
            coaxial.compute_effective_borehole_resistance(),
            0.19286205271006363, places=8)

    def test_coaxial_to_single_u_tube(self):
        # Pipe dimensions
        # Inner pipe radii
        r_in_in = 44.2 / 1000. / 2.
        r_in_out = 50. / 1000. / 2.
        # Outer pipe radii
        r_out_in = 97.4 / 1000. / 2.
        r_out_out = 110. / 1000. / 2.
        # Pipe radii
        # Note: This convention is different from pygfunction
        r_inner = [r_in_in,
                   r_in_out]  # The radii of the inner pipe from in to out
        r_outer = [r_out_in,
                   r_out_out]  # The radii of the outer pipe from in to out
        epsilon = 1.0e-6  # Pipe roughness (m)

        # Pipe positioning
        pos = (0, 0)
        s = 0

        # Thermal properties
        k_p = [0.4, 0.4]  # Inner and outer pipe thermal conductivity (W/m.K)

        # Thermal properties
        # Pipe
        pipe = plat.media.Pipe(
            pos, r_inner, r_outer, s, epsilon, k_p, self.rhoCp_p)

        # Define a borehole
        borehole = gt.boreholes.Borehole(self.H, self.D, self.r_b, x=0., y=0.)

        coaxial = plat.borehole_heat_exchangers.CoaxialPipe(
            self.m_flow_borehole, self.fluid, borehole, pipe, self.grout,
            self.soil)

        # Test intermediate variables
        V_fluid, V_pipe, R_conv, R_pipe = \
            plat.equivalance.concentric_tube_volumes(coaxial)
        self.assertAlmostEqual(V_fluid, 0.007021773740038546, places=8)
        self.assertAlmostEqual(V_pipe, 0.0024815440370705762, places=8)
        self.assertAlmostEqual(R_conv, 0.04531412438364096, places=8)
        self.assertAlmostEqual(R_pipe, 0.048404650347060686, places=8)

        single_u_tube = plat.equivalance.compute_equivalent(coaxial)

        # Test the single U-tube equivalent parameters
        V_flow = single_u_tube.m_flow_pipe * 1000. / single_u_tube.fluid.rho
        self.assertAlmostEqual(V_flow, 0.2)
        self.assertAlmostEqual(single_u_tube.r_in, 0.03342977714553299)
        self.assertAlmostEqual(single_u_tube.r_out, 0.03889087296526011)
        self.assertAlmostEqual(single_u_tube.pipe.s, 0.005370899457402695)
        self.assertAlmostEqual(single_u_tube.h_f, 247.4249791690743)
        self.assertAlmostEqual(single_u_tube.pipe.k, 0.3233497585764886)
        self.assertAlmostEqual(single_u_tube.grout.k, 0.5354034529364631)
        self.assertAlmostEqual(
            single_u_tube.compute_effective_borehole_resistance(),
            0.19286205242585608)

    def test_double_to_single_u_tube(self):
        # Pipe dimensions
        r_out = 26.67 / 1000. / 2.  # Pipe outer radius (m)
        r_in = 21.6 / 1000. / 2.  # Pipe inner radius (m)
        s = 32.3 / 1000.  # Inner-tube to inner-tube Shank spacing (m)
        epsilon = 1.0e-6  # Pipe roughness (m)

        # Pipe positions
        # Double U-tube [(x_in, y_in), (x_out, y_out), (x_in, y_in), (x_out, y_out)]
        pos = plat.media.Pipe.place_pipes(s, r_out, 2)

        # Thermal properties
        # Pipe
        pipe = plat.media.Pipe(pos, r_in, r_out, s, epsilon, self.k_p, self.rhoCp_p)

        # Define a borehole
        borehole = gt.boreholes.Borehole(self.H, self.D, self.r_b, x=0., y=0.)

        # Double U-tube defaults to parallel
        double_u_tube = plat.borehole_heat_exchangers.MultipleUTube(
            self.m_flow_borehole, self.fluid, borehole, pipe, self.grout,
            self.soil)

        # Test intermediate variables
        V_fluid, V_pipe, R_conv, R_pipe = \
            plat.equivalance.u_tube_volumes(double_u_tube)
        self.assertAlmostEqual(V_fluid, 0.001465741468458854)
        self.assertAlmostEqual(V_pipe, 0.0007688385143611112)
        self.assertAlmostEqual(R_conv, 0.13198055664449054)
        self.assertAlmostEqual(R_pipe, 0.020973241794937685)

        single_u_tube = plat.equivalance.compute_equivalent(double_u_tube)

        # Test the single U-tube equivalent parameters
        V_flow = single_u_tube.m_flow_pipe * 1000. / single_u_tube.fluid.rho
        self.assertAlmostEqual(V_flow, 0.2)
        self.assertAlmostEqual(single_u_tube.r_in, 0.015273506473629427)
        self.assertAlmostEqual(single_u_tube.r_out, 0.01885853785424522)
        self.assertAlmostEqual(single_u_tube.pipe.s, 0.02485528286100637)
        self.assertAlmostEqual(single_u_tube.h_f, 1287.1059002600034)
        self.assertAlmostEqual(single_u_tube.pipe.k, 0.2316554011453535)
        self.assertAlmostEqual(single_u_tube.grout.k, 1.6966493055905383)
        self.assertAlmostEqual(
            single_u_tube.compute_effective_borehole_resistance(),
            0.15970357274051686)
