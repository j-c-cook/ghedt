import ghedt as dt
import ghedt.peak_load_analysis_tool as plat
import pygfunction as gt
import numpy as np
import textwrap


# Common design interface
class Design:
    def __init__(self, V_flow: float, borehole: gt.boreholes.Borehole,
                 bhe_object: plat.borehole_heat_exchangers,
                 fluid: gt.media.Fluid, pipe: plat.media.Pipe,
                 grout: plat.media.Grout, soil: plat.media.Soil,
                 sim_params: plat.media.SimulationParameters,
                 geometric_constraints: dt.media.GeometricConstraints,
                 hourly_extraction_ground_loads: list, method: str = 'hybrid',
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
        self.hourly_extraction_ground_loads = hourly_extraction_ground_loads
        self.method = method
        if self.method == 'hourly':
            msg = 'Note: It is not recommended to perform a field selection ' \
                  'with the hourly simulation due to computation time. If ' \
                  'the goal is to validate the selected field with the ' \
                  'hourly simulation, the better solution is to utilize the ' \
                  'hybrid simulation to automatically select the field. Then ' \
                  'perform a sizing routine on the selected GHE with the ' \
                  'hourly simulation.'
            # Wrap the text to a 50 char line width and print it
            wrapper = textwrap.TextWrapper(width=72)
            word_list = wrapper.wrap(text=msg)
            for element in word_list:
                print(element)
            print('\n')

        # Check the routine parameter
        self.routine = routine
        available_routines = ['near-square', 'rectangle', 'bi-rectangle',
                              'bi-zoned']
        self.geometric_constraints.check_inputs(self.routine)
        gc = self.geometric_constraints
        if routine in available_routines:
            # If a near-square design routine is requested, then we go from a
            # 1x1 to 32x32 at the B-spacing
            # The lower end of the near-square routine is always 1 borehole.
            # There would never be a time that a user would __need__ to give a
            # different lower range. The upper number of boreholes range is
            # calculated based on the spacing and length provided.
            if routine == 'near-square':
                number_of_boreholes = \
                    dt.utilities.number_of_boreholes(
                        gc.length, gc.B, func=np.floor)
                self.coordinates_domain = \
                    dt.domains.square_and_near_square(
                        1, number_of_boreholes, self.geometric_constraints.B)
            elif routine == 'rectangle':
                self.coordinates_domain = dt.domains.rectangular(
                    gc.length, gc.width, gc.B_min, gc.B_max_x, disp=False
                )
            elif routine == 'bi-rectangle':
                self.coordinates_domain_nested = dt.domains.bi_rectangle_nested(
                    gc.length, gc.width, gc.B_min, gc.B_max_x, gc.B_max_y,
                    disp=False)
            elif routine == 'bi-zoned':
                self.coordinates_domain_nested = \
                    dt.domains.bi_rectangle_zoned_nested(
                        gc.length, gc.width, gc.B_min, gc.B_max_x, gc.B_max_y)
        else:
            raise ValueError('The requested routine is not available. '
                             'The currently available routines are: '
                             '`near-square`.')
        self.flow = flow

    def find_design(self, disp=False):
        if disp:
            title = 'Find {}...'.format(self.routine)
            print(title + '\n' + len(title) * '=')
        # Find near-square
        if self.routine == 'near-square':
            bisection_search = dt.search_routines.Bisection1D(
                self.coordinates_domain, self.V_flow, self.borehole,
                self.bhe_object, self.fluid, self.pipe, self.grout,
                self.soil, self.sim_params, self.hourly_extraction_ground_loads,
                method=self.method, flow=self.flow, disp=disp)
        # Find a rectangle
        elif self.routine == 'rectangle':
            bisection_search = dt.search_routines.Bisection1D(
                self.coordinates_domain, self.V_flow, self.borehole,
                self.bhe_object, self.fluid, self.pipe, self.grout, self.soil,
                self.sim_params, self.hourly_extraction_ground_loads,
                method=self.method, flow=self.flow, disp=disp)
        # Find a bi-rectangle
        elif self.routine == 'bi-rectangle':
            bisection_search = dt.search_routines.Bisection2D(
                self.coordinates_domain_nested, self.V_flow,
                self.borehole, self.bhe_object, self.fluid, self.pipe,
                self.grout, self.soil, self.sim_params,
                self.hourly_extraction_ground_loads, method=self.method,
                flow=self.flow, disp=disp)
        # Find bi-zoned rectangle
        elif self.routine == 'bi-zoned':
            bisection_search = dt.search_routines.BisectionZD(
                self.coordinates_domain_nested, self.V_flow, self.borehole,
                self.bhe_object, self.fluid, self.pipe, self.grout, self.soil,
                self.sim_params, self.hourly_extraction_ground_loads,
                method=self.method, flow=self.flow, disp=disp)
        else:
            raise ValueError('The requested routine is not available. '
                             'The currently available routines are: '
                             '`near-square`.')

        return bisection_search
