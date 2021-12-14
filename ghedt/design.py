# Jack C. Cook
# Friday, December 10, 2021

import ghedt as dt
import ghedt.PLAT as PLAT
import ghedt.PLAT.pygfunction as gt


# Common design interface
class Design:
    def __init__(self, V_flow: float, borehole: gt.boreholes.Borehole,
                 bhe_object: PLAT.borehole_heat_exchangers,
                 fluid: gt.media.Fluid, pipe: PLAT.media.Pipe,
                 grout: PLAT.media.ThermalProperty, soil: PLAT.media.Soil,
                 sim_params: PLAT.media.SimulationParameters,
                 geometric_constraints: dt.media.GeometricConstraints,
                 coordinates_domain: list, hourly_extraction_ground_loads: list,
                 routine: str = 'near-square', flow: str = 'borehole'):
        self.V_flow = V_flow  # volumetric flow rate, m3/s
        self.borehole = borehole
        self.bhe_object = bhe_object  # a borehole heat exchanger object
        self.fluid = fluid  # a fluid object
        self.pipe = pipe
        self.grout = grout
        self.soil = soil
        self.sim_params = sim_params
        self.geometric_constraints = geometric_constraints
        self.coordinates_domain = coordinates_domain
        self.hourly_extraction_ground_loads = hourly_extraction_ground_loads

        # Check the routine parameter
        self.routine = routine
        available_routines = ['near-square']
        if routine in available_routines:
            # If a near-square design routine is requested, then we go from a
            # 1x1 to 32x32 at the B-spacing
            if routine == 'near-square':
                self.coordinates_domain = \
                    dt.domains.square_and_near_square(
                        1, 32, self.geometric_constraints.B_max_x)
        else:
            raise ValueError('The requested routine is not available. '
                             'The currently available routines are: '
                             '`near-square`.')

        # Check the flow rate parameter
        if flow == 'borehole':
            self.V_flow_borehole = V_flow  # borehole volumetric flow rate, m3/s
        elif flow == 'system':
            self.V_flow_system = V_flow  # system volumetric flow rate, m3/s
        else:
            raise ValueError('The flow rate should be on a `borehole` or '
                             '`system` basis.')

    def find_design(self, disp=False):
        # Find near-square
        if self.routine == 'near-square':
            bisection_search = dt.search_routines.Bisection1D(
                self.coordinates_domain, self.V_flow_borehole, self.borehole,
                self.bhe_object, self.fluid, self.pipe, self.grout,
                self.soil, self.sim_params, self.hourly_extraction_ground_loads,
                disp=disp)
        else:
            raise ValueError('The requested routine is not available. '
                             'The currently available routines are: '
                             '`near-square`.')

        return bisection_search
