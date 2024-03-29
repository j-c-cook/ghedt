import copy
from typing import List

from calendar import monthrange

import scipy.optimize.optimize

import ghedt.peak_load_analysis_tool as plat
from scipy import interpolate
import math
import pandas as pd
import numpy as np
import pygfunction as gt


def synthetic_load_function(
        hourly_loads, avg_hourly_load, peak_load, peak_load_duration, peak_day,
        hours_in_previous_months, hours_in_month) -> None:
    # This function can be used to help generate synthetic loads. Note that the
    # lists in python (as well as dictionaries and user defined objects) are
    # passed by reference, and do not need to be returned. Therefore, the only
    # variable here that needs returned is the integer we need to update,
    # because integers, floats and strings standing alone do not fall into the
    # same category for passing by reference. They can be updated locally, and
    # the change will have no effect upon return

    hours_in_day: int = 24

    # The hourly loads for the current month are set to avg_hourly_load
    hourly_loads[hours_in_previous_months:
                 hours_in_previous_months + hours_in_month] = \
        [avg_hourly_load] * hours_in_month

    # The peak load will split the noon hour
    peak_day_start_hour = peak_day * hours_in_day + hours_in_previous_months
    peak_day_noon = peak_day_start_hour + 12
    # If the peak load is odd, then the peak load is shifted one to the left
    # set the peak_load for a duration of peak_load_duration
    peak_load_start_hour = \
        math.floor(peak_day_noon - (float(peak_load_duration) / 2.))
    hourly_loads[peak_load_start_hour:
                 peak_load_start_hour + peak_load_duration] = \
        [peak_load] * peak_load_duration

    return


def create_synthetic_doubling_load_profile(units='W', year=2019) -> tuple:
    # Create a synthetic load doubling profile: The first month has an average
    # load of 1 (Watt or kiloWatt based on units), with a peak of 2 units. Month
    # 2 has 2 units with a 2 hour peak of 4 units. These are applied to the
    # extraction loads. The rejection loads get the opposite treatment (nearly).
    # Simply reversing the order of te list of extraction loads gives different
    # peak load days, so the following is done with rejection: the first month
    # has an average load of 12 units and a peak duration of 24 units.

    # Check the units
    if units == 'W':
        scale = 1000.
    elif units == 'kW':
        scale = 1.
    else:
        raise ValueError('Units provided are not an option.')

    # Get the number of days in each month for a given year (make 0 NULL)
    days_in_month = \
        [0] + [monthrange(year, i)[1] for i in range(1, 13)]
    days_in_year = sum(days_in_month)
    hours_in_day: int = 24

    # Heat rejection in the ground occurs when buildings are in cooling
    # mode, these loads appear negative on Ground extraction loads plots
    # hourly_rejection_loads: list = [0.] * days_in_year * hours_in_day
    # Heat extraction in the ground occurs when buildings are in heating
    # mode, these loads appear positive on Ground extraction load plots
    hourly_extraction_loads: list = [0.] * days_in_year * hours_in_day
    hourly_rejection_loads: list = [0.] * days_in_year * hours_in_day

    # the peak extraction will occur on day 14
    peak_day_ext: int = 7
    # the peak rejection will occur on day 17
    peak_day_rej: int = 19

    # Store the index of the last months hours
    hours_in_previous_months: int = 0
    for i in range(1, len(days_in_month)):
        # The total number of hours in this month
        hours_in_month: int = hours_in_day * days_in_month[i]
        # The extraction monthly average load
        avg_hourly_load_ext: float = float(i) * scale
        # The rejection monthly average load
        avg_hourly_load_rej: float = float(len(days_in_month) - i)

        # The peak load occurs in day peak_day and for time peak_load_duration
        peak_load_ext: float = 2. * float(i) * scale  # The monthly peak load
        peak_load_duration_ext: int = i  # The peak load duration

        # The peak load and peak load duration for rejection
        peak_load_rej: float = 2. * float(len(days_in_month) - i) * scale
        peak_load_duration_rej: int = len(days_in_month) - i

        synthetic_load_function(
            hourly_extraction_loads, avg_hourly_load_ext, peak_load_ext,
            peak_load_duration_ext, peak_day_ext, hours_in_previous_months,
            hours_in_month)
        synthetic_load_function(
            hourly_rejection_loads, avg_hourly_load_rej, peak_load_rej,
            peak_load_duration_rej, peak_day_rej, hours_in_previous_months,
            hours_in_month)

        # Count up on extraction days from 7 to 18
        peak_day_ext += 1
        # Count down on rejection days from 18 to 7
        peak_day_rej -= 1

        # track the hours that have passed
        hours_in_previous_months += hours_in_month

    return hourly_rejection_loads, hourly_extraction_loads


class Loads:
    def __init__(self, hourly_loads: List[float], g_sts, ts, two_pi_k, Rb, year=2019):
        self.hourly = hourly_loads

        self.g_sts = g_sts
        self.ts = ts
        self.two_pi_k = two_pi_k
        self.Rb = Rb

        # Get the number of days in each month for a given year (make 0 NULL)
        self.days_in_month = \
            [0] + [monthrange(year, i)[1] for i in range(1, 13)]
        assert len(hourly_loads) == sum(self.days_in_month) * 24., \
            "The total number of hours in the year are not equal. Is this a" \
            " leap year?"

        # This block of data holds the compact monthly representation of the
        # loads. The intention is that these loads will usually repeat. It's
        # possible that for validation or design purposes, users may wish to
        # specify loads that differ from year to year. For these arrays,
        # January is the second item (1) and December the last (12)
        # We'll reserve the first item (0) for an annual total or peak
        self.monthly_total = [0.] * 13
        self.monthly_peak = [0.] * 13
        self.monthly_average = [0.] * 13
        self.monthly_peak_day = [0.] * 13
        self.split_load_by_month()

        # 48 hour loads are going to be necessary for the hourly simulation for
        # finding the peak load duration
        # These will be a 2D list, a list of 48 hour loads in each index
        # Make 0 position NULL
        self.monthly_two_day_hourly_peak = [[0]]
        self.process_two_day_loads()

        # Now we need to perform 48-hour simulations to determine the
        # monthly peak load hours
        self.monthly_two_day_fluid_temps_nm = [[0]]
        # Stores two day (48 hour) fluid temperatures with peak load
        self.monthly_two_day_fluid_temps_pk = [[0]]

        # duration of monthly peak clg load in hours
        self.monthly_peak_duration = [0] * 13
        self.find_peak_durations()

    def split_load_by_month(self):
        # Split the loads into peak, total and average loads for each month

        hours_in_day = 24
        # Store the index of the last months hours
        hours_in_previous_months = 0
        for i, days_in_month in enumerate(self.days_in_month[1:], start=1):
            hours_in_month = hours_in_day * days_in_month
            # Slice the hours in this current month
            month_slice = slice(hours_in_previous_months,
                                hours_in_previous_months + hours_in_month)
            month_loads = self.hourly[month_slice]
            assert len(month_loads) == hours_in_month

            # Total loads for the month in kW
            self.monthly_total[i] = sum(month_loads)

            # Peak load in kW
            # monthly peak cooling
            self.monthly_peak[i] = max(month_loads)

            # Average monthly load in kW
            self.monthly_average[i] = self.monthly_total[i] / len(month_loads)

            # Day of month the peak load occurs (e.g. 1-31)
            self.monthly_peak_day[i] = math.ceil(month_loads.index(
                self.monthly_peak[i]) / hours_in_day)

            hours_in_previous_months += hours_in_month

        return

    def process_two_day_loads(self) -> None:
        # The two day (48 hour) two day loads are selected by locating the day
        # the peak load of the month occurs on, and pulling a 48-hour load
        # profile -- the day before and the day of

        hours_in_day = 24
        hours_in_year = len(self.hourly)

        # Add the last day of the year to the beginning of the loads to account
        # for the possibility that a peak load occurs on the first day of the
        # year

        hourly_loads = \
            self.hourly[hours_in_year-hours_in_day:hours_in_year] + \
            self.hourly + self.hourly[0:hours_in_day]

        # Keep track of how many hours are in
        # start at 24 since we added the last day of the year to the beginning
        hours_in_previous_months = hours_in_day
        # loop over all 12 months
        for i, days_in_month in enumerate(self.days_in_month[1:], start=1):
            hours_in_month = hours_in_day * days_in_month

            monthly_peak_day = self.monthly_peak_day[i]

            # Get the starting hour of the day before the peak cooling load day
            monthly_peak_hour_start = hours_in_previous_months + (monthly_peak_day-1) * hours_in_day

            # monthly cooling loads (or heat rejection) in kWh
            two_day_hourly_peak_load = \
                hourly_loads[monthly_peak_hour_start:monthly_peak_hour_start+2*hours_in_day]

            assert len(two_day_hourly_peak_load) == 2*hours_in_day

            # Double check ourselves
            monthly_peak_day_start = \
                int((monthly_peak_hour_start-hours_in_day) / hours_in_day)
            monthly_peak_cl_hour_month = \
                int(monthly_peak_day_start - sum(self.days_in_month[0:i]))
            assert monthly_peak_cl_hour_month == monthly_peak_day - 1

            # monthly loads in kWh
            self.monthly_two_day_hourly_peak.append(two_day_hourly_peak_load)

            hours_in_previous_months += hours_in_month

        return

    def find_peak_durations(self) -> None:
        # Find the peak durations using hourly simulations for 2 days

        for i in range(1, len(self.days_in_month)):
            # Scale all the loads by the peak load
            # Perform an hourly simulation with the scaled loads
            # Perform an hourly simulation with a load of 1, or the peak loads
            # divided by the peak

            # two day cooling loads (or heat rejection) in kWh
            current_two_day_load = [0.] + self.monthly_two_day_hourly_peak[i]

            # This tolerance applies to the difference between the current
            # months peak load and the maximum of the two-day load. If the
            # absolute value of the difference between the current months
            # peak load and the current two-day peak load is within this
            # tolerance, then the maximum of the two-day load is equal to the
            # maximum of the current month. If the absolute difference is
            # greater than the tolerance, then the two-day peak load contains
            # a load greater than the current months peak load. The tolerance
            # could ONLY be exceeded when the first 24 hours is located in the
            # previous month.
            tol = 0.1

            # Ensure the peak load for the two-day load profile is the same or
            # greater than the monthly peak load. This check is done in case
            # the previous month contains a higher load than the current month.
            load_diff = self.monthly_peak[i] - max(current_two_day_load)
            # monthly peak cooling load (or heat rejection) in kW
            if abs(load_diff) < tol:
                current_month_peak_cl = self.monthly_peak[i]
            else:
                current_month_peak_cl = max(current_two_day_load)

            # monthly average cooling load (or heat rejection) in kW
            current_month_avg_cl = self.monthly_average[i]

            if current_month_peak_cl != 0.0:
                peak_duration, q_peak, q_nominal = \
                    self.perform_current_month_simulation(
                        current_two_day_load, current_month_peak_cl,
                        current_month_avg_cl, self.monthly_two_day_fluid_temps_pk,
                        self.monthly_two_day_fluid_temps_nm)
            else:
                peak_duration = 1.0e-6

            # Set the monthly cooling load duration
            self.monthly_peak_duration[i] = peak_duration

        return

    @staticmethod
    def simulate_hourly(hour_time, q, g_sts, Rb, two_pi_k, ts):
        # An hourly simulation for the fluid temperature
        # Chapter 2 of Advances in Ground Source Heat Pumps

        q_dt = np.hstack((q[1:] - q[:-1]))

        dT_fluid = [0]
        for n in range(1, len(hour_time)):
            # Take the last i elements of the reversed time array
            _time = hour_time[n] - hour_time[0:n]
            # _time = time_values_reversed[n - i:n]
            g_values = g_sts(np.log((_time * 3600.) / ts))
            # Tb = Tg + (q_dt * g)  (Equation 2.12)
            delta_Tb_i = (q_dt[0:n] / two_pi_k).dot(g_values)
            # Delta mean HPEFT fluid temperature
            Tf_mean = delta_Tb_i + q[n] * Rb
            dT_fluid.append(Tf_mean)

        return dT_fluid

    def perform_current_month_simulation(
            self, two_day_hourly_peak_load, peak_load, avg_load,
            two_day_fluid_temps_pk, two_day_fluid_temps_nm):
        g_sts = self.g_sts
        ts = self.ts
        two_pi_k = self.two_pi_k
        Rb = self.Rb

        hours_in_day = 24
        hour_time = np.array(list(range(0, 2 * hours_in_day + 1)))
        # Two day peak cooling load scaled down by average (q_max - q_avg)
        q_peak = np.array([0.] + [peak_load - avg_load] * (2 * hours_in_day))
        # Two day nominal cooling load (q_i - q_avg) / q_max * q_i
        q_nominal = np.array(
            [0.] + [(two_day_hourly_peak_load[i] - avg_load) / peak_load *
                    two_day_hourly_peak_load[i] for i in range(1, len(q_peak))])
        # Get peak fluid temperatures using peak load
        dT_fluid_pk = self.simulate_hourly(
            hour_time, q_peak, g_sts, Rb, two_pi_k, ts)
        two_day_fluid_temps_pk.append(dT_fluid_pk)
        # Get nominal fluid temperatures using nominal load
        dT_fluid_nm = self.simulate_hourly(
            hour_time, q_nominal, g_sts, Rb, two_pi_k, ts)
        two_day_fluid_temps_nm.append(dT_fluid_nm)

        dT_fluid_nm_max = max(dT_fluid_nm)

        if dT_fluid_nm_max > 0.0:
            f = scipy.interpolate.interp1d(dT_fluid_pk, hour_time)
            peak_duration = f(dT_fluid_nm_max).tolist()
        else:
            peak_duration = 1.0e-6

        return peak_duration, q_peak, q_nominal

    def peak_energy(self, m_idx):
        # gives htg load pk energy in kWh
        return self.monthly_peak[m_idx] * self.monthly_peak_duration[m_idx]

    def peak_hour_range(self, m_idx, total_months):
        # Place the peaks roughly midway through the day they occur on.
        # (In JDS's opinion, this should be amply accurate for the
        # hybrid time step.)
        # Catch the first and last peak hours to make sure they aren't 0
        # Could only be 0 when the first month has no load.
        _hours_to_mid_day = (self.monthly_peak_day[m_idx] - 1) * 24 + 12
        _half_peak_duration = (self.monthly_peak_duration[m_idx] / 2)
        first_hour_peak = \
            firstmonthhour(total_months) + _hours_to_mid_day - _half_peak_duration
        if first_hour_peak < 0.:
            first_hour_peak = 1.0e-6
        last_hour_peak = first_hour_peak + self.monthly_peak_duration[m_idx]
        return first_hour_peak, last_hour_peak


class HybridLoad:
    def __init__(self,
                 hourly_rejection_loads: list, hourly_extraction_loads: list,
                 bhe: plat.borehole_heat_exchangers.SingleUTube,
                 radial_numerical: plat.radial_numerical_borehole.RadialNumericalBH,
                 sim_params: plat.media.SimulationParameters,
                 COP_rejection=None, COP_extraction=None, year=2019):
        # Store the borehole heat exchanger
        self.bhe = bhe
        # Store the radial numerical g-function value
        # Note: this is intended to be a scipy.interp1d object
        self.radial_numerical = radial_numerical

        _g_sts = self.radial_numerical.g_sts
        _ts = self.radial_numerical.t_s
        _two_pi_k = 2. * np.pi * self.bhe.soil.k
        _Rb = self.bhe.compute_effective_borehole_resistance()

        # monthly cooling loads (or heat rejection)
        self.rejection = Loads(
            hourly_rejection_loads, _g_sts, _ts, _two_pi_k, _Rb, year=year)
        # monthly heating loads (or heat extraction)
        self.extraction = Loads(
            hourly_extraction_loads, _g_sts, _ts, _two_pi_k, _Rb, year=year)

        if COP_extraction is None:
            self.COP_extraction = 2.5  # When the building is heating mode
        else:
            self.COP_extraction = COP_extraction
        if COP_rejection is None:
            self.COP_rejection = 4.    # When the building is in cooling mode
        else:
            self.COP_rejection = COP_rejection

        # Simulation start and end month
        self.startmonth = sim_params.start_month
        self.endmonth = sim_params.end_month

        self.peakretainstart = 12  # use peak laods for first 12 months
        self.peakretainend = 12  # use peak loads for last 12 months
        # This block of data holds the sequence of loads. This is an
        # intermediate form, where the load values hold the actual loads,
        # not the the devoluted loads
        self.load = np.array(0)  # holds the load during the period
        self.hour = np.array(0)  # holds the last hour of a period
        self.sfload = np.array(0)  # holds the load in terms of step functions

        self.process_monthly_loads()

    def __repr__(self):
        output = str(self.__class__) + '\n'

        output += self.create_dataframe_of_peak_analysis().to_string()

        return output

    @staticmethod
    def split_heat_and_cool(hourly_heat_extraction, units='W'):
        """
         JCC 02.16.2020
         Split the provided loads into heating and cooling. Heating is positive,
         cooling is negative.
         :return: Loads split into heating and cooling
         """
        # Expects hourly_heat_extraction to be in Watts

        # Heat rejection in the ground occurs when buildings are in cooling
        # mode, these loads appear negative on Ground extraction loads plots
        hourly_rejection_loads: list = [0.] * len(hourly_heat_extraction)
        # Heat extraction in the ground occurs when buildings are in heating
        # mode, these loads appear positive on Ground extraction load plots
        hourly_extraction_loads: list = [0.] * len(hourly_heat_extraction)

        if units == 'W':
            scale = 1000.
        elif units == 'kW':
            scale = 1.
        else:
            raise ValueError('Units provided are not an option.')

        for i, l_hour in enumerate(hourly_heat_extraction):
            if l_hour >= 0.0:
                # Heat is extracted from ground when > 0
                hourly_extraction_loads[i] = l_hour / scale
            else:
                # Heat is rejected to ground when < 0
                hourly_rejection_loads[i] = l_hour / -scale

        return hourly_rejection_loads, hourly_extraction_loads

    def create_dataframe_of_peak_analysis(self) -> pd.DataFrame:
        # cooling loads (or heat rejection)
        cl = self.rejection
        # heating loads (or heat extraction)
        hl = self.extraction

        # The fields are: sum, peak, avg, peak day, peak duration
        hybrid_time_step_fields = {'Total': {}, 'Peak': {},
                  'Average': {}, 'Peak Day': {},
                  'Peak Duration': {}}

        d: dict = {}
        # For all the months, create dictionary of fields
        for i in range(1, 13):
            month_name = number_to_month(i)
            d[month_name] = copy.deepcopy(hybrid_time_step_fields)

            # set total
            d[month_name]['Total']['rejection'] = cl.monthly_total[i]
            d[month_name]['Total']['extraction'] = hl.monthly_total[i]
            # set peak
            d[month_name]['Peak']['rejection'] = cl.monthly_peak[i]
            d[month_name]['Peak']['extraction'] = hl.monthly_peak[i]
            # set average
            d[month_name]['Average']['rejection'] = cl.monthly_average[i]
            d[month_name]['Average']['extraction'] = hl.monthly_average[i]
            # set peak day
            d[month_name]['Peak Day']['rejection'] = cl.monthly_peak_day[i]
            d[month_name]['Peak Day']['extraction'] = hl.monthly_peak_day[i]
            # set peak duration
            d[month_name]['Peak Duration']['rejection'] = \
                cl.monthly_peak_duration[i]
            d[month_name]['Peak Duration']['extraction'] = \
                hl.monthly_peak_duration[i]

        # Convert the dictionary into a multi-indexed pandas dataframe
        arrays = [[], []]
        for field in hybrid_time_step_fields:
            arrays[0].append(field)
            arrays[0].append(field)
            arrays[1].append('rejection')
            arrays[1].append('extraction')
        tuples = list(zip(*arrays))
        index = pd.MultiIndex.from_tuples(tuples, names=['Fields', 'Load Type'])

        res = []
        for month in d:
            tmp = []
            for field in d[month]:
                for load_type in d[month][field]:
                    tmp.append(d[month][field][load_type])
            res.append(tmp)
        res = np.array(res)

        df = pd.DataFrame(res, index=list(d.keys()), columns=index)

        return df

    def process_monthly_loads(self):
        # Converts monthly load format to sequence of loads needed for
        # simulation
        # This routine is taking loads applied to the ground NOT to a heat pump.

        # cooling loads (or heat rejection)
        cl = self.rejection
        # heating loads (or heat extraction)
        hl = self.extraction

        # First, begin array with zero load before simulation starts.
        self.load = np.append(self.load, 0)
        #        self.sfload = np.append(self.sfload,0)
        lastzerohour = firstmonthhour(self.startmonth) - 1
        self.hour = np.append(self.hour, lastzerohour)

        for i in range(self.startmonth, self.endmonth + 1):
            mi = i % 12  # month index
            if mi == 0:
                mi = 12

            # set flag that determines if peak load will be included
            if i < self.startmonth + self.peakretainstart:
                include_peak = True
            elif i > self.endmonth - self.peakretainend:
                include_peak = True
            else:
                include_peak = False

            def _place_loads(order: List[Loads], sign: List[float], rate):
                for k, _ in enumerate(order):
                    if (order[k].monthly_peak[mi] > 0 and include_peak):
                        _first_hour_peak, _last_hour_peak = \
                            order[k].peak_hour_range(mi, i)
                        self.load = np.append(self.load, rate)
                        self.hour = np.append(self.hour, _first_hour_peak)
                        # peak
                        self.load = np.append(
                            self.load, sign[k] * order[k].monthly_peak[mi])
                        self.hour = np.append(self.hour, _last_hour_peak)
                # rest of month
                lastavghour = lastmonthhour(i)
                self.load = np.append(self.load, rate)
                self.hour = np.append(self.hour, lastavghour)

            # There may be a more sophisticated way to do this, but I will loop
            # through the lists mduration is the number of hours over which to
            # calculate the average value for the month (edited)
            if include_peak:
                mduration = \
                    monthdays(mi) * 24 - cl.monthly_peak_duration[mi] - \
                    hl.monthly_peak_duration[mi]
                mpeak_cl = cl.peak_energy(mi)  # peak clg energy
                mpeak_hl = hl.peak_energy(mi)  # peak htg energy
                mload = cl.monthly_total[mi] - hl.monthly_total[mi] - mpeak_cl + mpeak_hl
                mrate = mload / mduration
                peak_day_diff = cl.monthly_peak_day[mi] - hl.monthly_peak_day[mi]
            else:
                # peak load not used this month
                mduration = monthdays(mi) * 24
                mload = cl.monthly_total[mi] - hl.monthly_total[mi]
                mrate = mload / mduration
                peak_day_diff = 0

            if peak_day_diff < 0:
                order = [cl, hl]
                sign = [1.0, -1.0]
                _place_loads(order, sign, mrate)
            elif peak_day_diff > 0:
                order = [hl, cl]
                sign = [-1.0, 1.0]
                _place_loads(order, sign, mrate)
            else:
                # monthly peak heating day and cooling day are the same
                # in this case, we are ignoring the peaks
                # This is also used for the case where ipf[i] is False
                # A more sophisticated default could be use, like placing the
                # peaks on the 10th and 20th
                lastavghour = lastmonthhour(i)
                self.load = np.append(self.load, mrate)
                self.hour = np.append(self.hour, lastavghour)

        # Now fill array containing step function loads
        # Note they are paired with the ending hour, so the ith load will start
        # with the (i-1)th time
        n = self.hour.size
        # Note at this point the load and hour np arrays contain zeroes in
        # indices zero and one, then continue from there.
        for i in range(1, n):
            step_load = self.load[i] - self.load[i - 1]
            self.sfload = np.append(self.sfload, step_load)

    def hourly_load_representation(self, year=2019):

        # cooling loads (or heat rejection)
        cl = self.rejection
        # heating loads (or heat extraction)
        hl = self.extraction

        # Get the number of days in each month for a given year (make 0 NULL)
        days_in_month = \
            [0] + [monthrange(year, i)[1] for i in range(1, 13)]
        days_in_year = sum(days_in_month)
        hours_in_day: int = 24

        # Heat rejection in the ground occurs when buildings are in cooling
        # mode, these loads appear negative on Ground extraction loads plots
        # hourly_rejection_loads: list = [0.] * days_in_year * hours_in_day
        # Heat extraction in the ground occurs when buildings are in heating
        # mode, these loads appear positive on Ground extraction load plots
        hourly_extraction_loads: list = [0.] * days_in_year * hours_in_day
        hourly_rejection_loads: list = [0.] * days_in_year * hours_in_day

        # Store the index of the last months hours
        hours_in_previous_months: int = 0
        for i in range(1, len(days_in_month)):
            # The total number of hours in this month
            hours_in_month: int = hours_in_day * days_in_month[i]
            # monthly average cooling load (or heat rejection) in kW
            avg_hourly_load_rej: float = cl.monthly_average[i]
            # monthly average heating load (or heat extraction) in kW
            avg_hourly_load_ext: float = hl.monthly_average[i]

            # monthly peak heating load (or heat extraction) in kW
            peak_load_ext: float = hl.monthly_peak[i]
            # duration of monthly peak htg load in hours
            peak_load_duration_ext: int \
                = int(round(hl.monthly_peak_duration[i]))
            # day of the month on which peak htg load occurs (e.g. 1-31)
            peak_day_ext: int = hl.monthly_peak_day[i]

            peak_load_rej: float = cl.monthly_peak[i]
            peak_load_duration_rej: int \
                = int(round(cl.monthly_peak_duration[i]))
            peak_day_rej: int = cl.monthly_peak_day[i]

            synthetic_load_function(
                hourly_extraction_loads, avg_hourly_load_ext, peak_load_ext,
                peak_load_duration_ext, peak_day_ext, hours_in_previous_months,
                hours_in_month)
            synthetic_load_function(
                hourly_rejection_loads, avg_hourly_load_rej, peak_load_rej,
                peak_load_duration_rej, peak_day_rej, hours_in_previous_months,
                hours_in_month)

            # track the hours that have passed
            hours_in_previous_months += hours_in_month

        return hourly_rejection_loads, hourly_extraction_loads

    def visualize_hourly_heat_extraction(self):

        # The rejection and extraction lists will only contain loads in kW
        # that are not 0.0
        rejection: list = []
        extraction: list = []

        # cooling loads (or heat rejection)
        cl = self.rejection
        # heating loads (or heat extraction)
        hl = self.extraction

        for i in range(len(hl.hourly)):
            if self.hourly_rejection_loads[i] != 0.0:
                rejection.append((i, -hl.hourly[i]))
            if cl.hourly[i] != 0.0:
                # Flip the sign back negative because split_heat_and_cool made
                # all loads positive
                extraction.append((i, cl.hourly[i]))

        fig = gt.gfunction._initialize_figure()
        ax = fig.add_subplot(111)

        # Catch the error in cases where there are no loads
        if len(extraction) > 0:
            h_hour, h_load = zip(*extraction)
            ax.scatter(
                h_hour, h_load, c='red', s=7,
                label='Building in heating mode (ground heat extraction)')
        if len(rejection) > 0:
            c_hour, c_load = zip(*rejection)
            ax.scatter(c_hour, c_load, c='blue', s=7,
                       label='Building in cooling mode (ground heat rejection)')

        ax.set_ylabel('Ground extraction loads (kW)')
        ax.set_xlabel('Hours in one year')

        ax.grid()
        ax.set_axisbelow(True)

        fig.legend()

        fig.tight_layout(rect=(0, 0, 1, .9))

        return fig

    def visualize_fluid_simulation_results(self, month_num, extraction=True,
                                           rejection=True):
        # Plots the simulation result by month
        # rejection on ground is cooling load in building in kWh
        # extraction on ground is heating load in building in kWh

        fig = gt.gfunction._initialize_figure()
        ax = fig.add_subplot(111)

        hours_in_day = 24
        two_day_hour_time = list(range(0, 2 * hours_in_day+1))

        def draw_horizontal(two_day_fluid_temps_nm, peak_duration):
            # Find the place where the peak occurs in two day fluid temp nominal
            two_day_fluid_nm_max = max(two_day_fluid_temps_nm)
            two_day_fluid_nm_max_idx = two_day_fluid_temps_nm.index(
                two_day_fluid_nm_max)

            points = [(peak_duration, two_day_fluid_nm_max),
                      (two_day_fluid_nm_max_idx, two_day_fluid_nm_max)]
            x, y = list(zip(*points))
            ax.plot(x, y, color='orange', zorder=0)
            # Plot x marks the spot
            ax.scatter(peak_duration, two_day_fluid_nm_max, marker='X', s=50,
                       color='purple', zorder=2)
            return

        new_ticks = []

        # cooling loads (or heat rejection)
        cl = self.rejection
        # heating loads (or heat extraction)
        hl = self.extraction

        if extraction:
            # plot two day (48 hour) fluid temperatures for cooling with peak
            # load
            ax.plot(two_day_hour_time,
                    cl.monthly_two_day_fluid_temps_pk[month_num],
                    'k--', label='Peak Extraction', zorder=0)
            # plot two day (48 hour) fluid temperatures for cooling with nominal
            # load
            ax.plot(two_day_hour_time,
                    cl.monthly_two_day_fluid_temps_nm[month_num],
                    'b-', label='Extraction', zorder=1)
            # Draw horizontal line from max time to fluid sim
            draw_horizontal(cl.monthly_two_day_fluid_temps_nm[month_num],
                            cl.monthly_peak_duration[month_num])
            pk_cl_duration = cl.monthly_peak_duration[month_num]
            pk_cl_duration_rounded = round(pk_cl_duration, 1)
            new_ticks.append(pk_cl_duration_rounded)

        if rejection:
            # plot two day (48 hour) fluid temperatures for heating with peak
            # load
            ax.plot(two_day_hour_time,
                    hl.monthly_two_day_fluid_temps_pk[month_num],
                    'k.', label='Peak Rejection', zorder=0)
            # plot two day (48 hour) fluid temperatures for heating with nominal
            # load
            ax.plot(two_day_hour_time,
                    hl.monthly_two_day_fluid_temps_nm[month_num],
                    'r-.', label='Rejection', zorder=1)
            draw_horizontal(hl.monthly_two_day_fluid_temps_nm[month_num],
                            hl.monthly_peak_duration[month_num])
            # duration of monthly peak htg load in hours
            pk_hl_duration = hl.monthly_peak_duration[month_num]
            pk_hl_duration_rounded = round(pk_hl_duration, 1)
            new_ticks.append(pk_hl_duration_rounded)

        ax.set_xticks(list(ax.get_xticks()) + new_ticks)

        ax.set_xlabel('Time (hours)')
        ax.set_ylabel('Change in fluid temperature, $\Delta$T$_f$ ($\degree$C)')

        ax.set_xlim([-2, 50])

        ax.grid()
        ax.set_axisbelow(True)

        fig.tight_layout()

        fig.legend()

        return fig


def number_to_month(x):
    # Convert a numeric 1-12 to a month name
    if int(x) <= 12 and int(x) > 0:

        list_of_months = {'1': 'January', '2': 'February', '3': 'March',
                          '4': 'April', '5': 'May', '6': 'June', '7': 'July',
                          '8': 'August', '9': 'September', '10': 'October',
                          '11': 'November', '12': 'December'}

        return list_of_months[str(x)]

    else:
        print('num_to_month function error: "num=' + str(x) + '"')


def monthindex(mname):
    months = {'January': 1, 'February': 2, 'March': 3, 'April': 4, 'May': 5,
              'June': 6, 'July': 7, 'August': 8, 'September': 9, 'October': 10,
              'November': 11, 'December': 12}
    mi = months.get(mname)
    return mi


def monthdays(month):
    if month > 12:
        md = month % 12
    else:
        md = month
    ndays = [31, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    monthdays = ndays[md]
    return monthdays


def firstmonthhour(month):
    fmh = 1
    if month > 1:
        for i in range(1,month):
            mi = i % 12
            fmh = fmh + 24 * monthdays(mi)
    return fmh


def lastmonthhour(month):
    lmh = 0
    for i in range(1, month + 1):
        lmh = lmh + monthdays(i) * 24
    if month == 1:
        lmh = 31 * 24
    return lmh
