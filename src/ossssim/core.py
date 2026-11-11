"""
The OSSOS Survey Simulator module.  This is front-end for the OSSSim.  The primary OSSSim is writen in f95.
The OSSOS Survey Simulator module.  This is front-end for the OSSSim.  The primary OSSSim is writen in f95.
"""
import logging
from astropy import units as u
from astropy.units import Quantity
from astropy.table import Row
from astropy.time import Time

from .color import PhotSpec
import ossssimlib
import rebound
import os
import numpy
from . import definitions
T_ORB_M_UNITS = definitions.T_ORB_M_UNITS
T_ORB_M_FIELDS = ('a', 'e', 'inc', 'node', 'peri', 'm')
# Characters per detection key returned by detos_batch; at least key_len in F95/parameters.f95.
_SURVEY_KEY_CODES = 128

_UNIT_SCALES = {}


def _to_value(value, unit):
    """
    Value of a Quantity in unit (plain numbers are returned unchanged).

    Same result as value.to(unit).value, but the unit scale factor is looked up
    once per pair of units instead of on every call.
    """
    if not isinstance(value, Quantity):
        return value
    from_unit = value.unit
    if from_unit is unit:
        return value.value
    scale = _UNIT_SCALES.get((from_unit, unit))
    if scale is None:
        scale = _UNIT_SCALES[(from_unit, unit)] = from_unit.to(unit)
    return value.value * scale


_LIGHT_CURVE_UNITS = {'gb': u.mag, 'phase': u.radian, 'period': u.day, 'amplitude': u.mag}
_LIGHT_CURVE_DEFAULTS = {key: definitions.LIGHT_CURVE_PARAMS[key].to(unit).value
                         for key, unit in _LIGHT_CURVE_UNITS.items()}


def _light_curve_value(row, key):
    """Light-curve parameter from row, or the default when absent or zero."""
    if key in row:
        value = _to_value(row[key], _LIGHT_CURVE_UNITS[key])
        if value:
            return value
    return _LIGHT_CURVE_DEFAULTS[key]


class Cartesian(object):
    """
    A class to provide the cartesian state vector of a set of objects.
    """
    def __init__(self, epoch=definitions.Neptune['Epoch']):
        self.epoch = epoch
        self.rebound_archive = f"Rebound_Archive{epoch.jd}.bin"
        self._base_sim = None
        self.initialize_sim()

    def _make_rebound_archive_file(self):
        _sim = rebound.Simulation()
        _sim.add("Sun",  date=self.epoch.to_datetime())
        _sim.add("Jupiter", date=self.epoch.to_datetime())
        _sim.add("Saturn", date=self.epoch.to_datetime())
        _sim.add("Uranus", date=self.epoch.to_datetime())
        _sim.add("Neptune", date=self.epoch.to_datetime())
        _sim.move_to_com()
        _sim.save_to_file(self.rebound_archive)
        del _sim

    def initialize_sim(self):
        """
        Start/open a rebound Simulation object used to compute the cartesian locations of the particles in this model.

        The archive is read once; each call returns a copy of that simulation.
        """
        if self._base_sim is None:
            if not os.access(self.rebound_archive, os.F_OK):
                self._make_rebound_archive_file()
            self._base_sim = rebound.Simulation(self.rebound_archive)
        return self._base_sim.copy()

    def __call__(self, **kwargs) -> dict:
        """
        Provide the state vector of the orbits at Mt or M if Mt is not provided.


        Returns:
            dictionary of x/y/z/vx/vy/vz
        """
        if "rows" not in kwargs:
            return self.single_row(**kwargs)

        sim = self.initialize_sim()
        rows = kwargs['rows']
        M = 'Mt' in rows and 'Mt' or 'M'
        for idx in range(len(rows['a'])):
            sim.add(a=rows['a'][idx].to('au').value,
                    e=rows['e'][idx],
                    inc=rows['inc'][idx].to('rad').value,
                    Omega=rows['node'][idx].to('rad').value,
                    omega=rows['peri'][idx].to('rad').value,
                    M=rows[M][idx].to('rad').value,
                    )

        for s in ['x', 'y', 'z', 'vx', 'vy', 'vz']:
            rows[s] = [p.__getattribute__(s) * definitions.column_unit[s] for p in sim.particles[5:]]
        return rows

    def single_row(self, **kwargs):
        """
        Provide the state vector of the orbits at Mt or M if Mt is not provided.

        Returns:
            dictionary of x/y/z/vx/vy/vz
        """
        sim = self.initialize_sim()
        sim.add(a=kwargs['a'].to('au').value,
                e=kwargs['e'],
                inc=kwargs['inc'].to('rad').value,
                Omega=kwargs['node'].to('rad').value,
                omega=kwargs['peri'].to('rad').value,
                M=kwargs.get('Mt', kwargs.get('M')).to('rad').value
                )
        cartesian = {}
        for s in ['x', 'y', 'z']:
            cartesian[s] = sim.particles[-1].__getattribute__(s) * u.au
        for s in ['vx', 'vy', 'vz']:
            cartesian[s] = sim.particles[-1].__getattribute__(s) * (2 * numpy.pi) * u.au / u.year
        return cartesian


class OSSSSim:
    """
    Outer Solar System Survey Simulator.

    This class simulates the process of observing a model of the solar system using a set of characterized observations.
    """

    def __init__(self, characterization_directory, seed):
        """
        Args:
            characterization_directory (str): path to a survey *root*.
                If ``root/pointings.list`` exists it is a single-epoch survey;
                otherwise each child directory with ``pointings.list`` is an
                epoch. Optional ``survey.conf`` sets ``detections_required``
                (default 1 = any epoch; set to n_epochs for AND). 
                Format docs are under docs/; small example characterizations
                live in data/Surveys/.
        """

        self.characterization_directory = characterization_directory
        self.cartesian = Cartesian(epoch=definitions.Neptune['Epoch'])
        self.seed = seed 
        ossssimlib.surveysub.reset_simulator()
        # One t_orb_m reused for every call; simulate() sets all of its fields each time.
        self._o_m = ossssimlib.datadec.t_orb_m()
        self._epoch_jd = (None, None)


    def simulate(self, row: dict, colors: PhotSpec, model_band, seed=None, epoch=None, debug=False) -> dict:
        """
        Pass the target elements to detos1 and determine if this target would be detected.

        Pass the survey root once; Detos1 loads all epochs and applies
        ``detections_required`` from ``survey.conf``. Changing
        ``self.characterization_directory`` between calls reloads that survey
        root without resetting the ran3 stream. Do not construct a new
        OSSSSim per epoch: ``__init__`` calls ``reset_simulator()`` and would
        reseed the RNG.

        Args:
            row (Row or dict): elements of the target to simulate
            colors (PhotSpec): colors of the target index by 'g-x' etc, see below.
            model_band (str): a single letter designator of the bandpass that the H values are provided in.
            seed (int): a seed to pass to the fortran code to allow reproducible simulations
            epoch (Quantity or float or Time): JD of elements, can also be provided for each row sent to simulate
            debug (bool): if True, have Fortran detos routine create a log of detection process.

        Returns:
            Row or dict: the target elements/values at time of simulated detection.

        row should have values for a, e, inc, node, peri, M, H, [epoch]
        can, optionally, also define: gb, phase, period, amplitude

        the result row has, in addition to above, the following items:

        flag: 0 - not detected, 1 - detected, 2 - tracked, 3 - characterized and lost, 4 - tacked and characterised

        When comparing to a detected sample you likely want to only use those with flag == 4

        following are values at time of detection, None if not detected

        RA: RA of target

        DEC: DEC of target

        d_ra: rate of RA sky motion

        d_dec: rate of DEC sky motion

        delta: distance from Earth

        r: distance from Sun

        m_int: intrinsic magnitude in filter of target H (ie. the model filter)

        m_rand: magnitude at detection, includes scatter due to flux measurement uncertainty (ie. the survey filter)

        h_rand: inferred absolute magnitude based on m_rand and detection circumstances

        eff: the detection field's efficiency of detection for a source of m_rand

        M: the mean anomaly at detection

        Survey: detection key "survey/block" for the characterizing block

        The colors list declares the color of the KBO in multiple filters. The value of model_band_pass is
        used to transform the color dictionary into colors in the bandpass of the model.

        A default light_curve_params dictionary is used if target provided to simulate doesn't have any, set in
        ossssim.definitions.LIGHT_CURVE_PARAMS
        """

        if seed is not None:
            logging.warning('seed is now set in the constructor. ignored here')

        # pack the orbit into a t_orb_m object to pass to fortran module.
        o_m = self._o_m
        row = dict(row)
        elements = {}
        for colname, value in row.items():
            field = colname.lower()
            if field in T_ORB_M_FIELDS:
                elements[field] = _to_value(value, T_ORB_M_UNITS.get(colname))
        for field in T_ORB_M_FIELDS:
            setattr(o_m, field, elements.get(field, 0.0))

        # attempt to detect this object
        color_offset_array = colors.colors_array(row['comp'], model_band)
        if 'epoch' in row.keys():
            epoch_jd = row['epoch']
        else:
            epoch_jd = epoch
        if isinstance(epoch_jd, Quantity):
            epoch_jd = _to_value(epoch_jd, u.day)
        if isinstance(epoch_jd, Time):
            if self._epoch_jd[0] is not epoch_jd:
                self._epoch_jd = (epoch_jd, epoch_jd.jd)
            epoch_jd = self._epoch_jd[1]
        gb = _light_curve_value(row, 'gb')
        phase = _light_curve_value(row, 'phase')
        period = _light_curve_value(row, 'period')
        amplitude = _light_curve_value(row, 'amplitude')
        H = None
        if 'H' in row:
            H = _to_value(row['H'], u.mag) or None
        if H is None and 'h' in row:
            H = _to_value(row['h'], u.mag) or None
        if H is None:
            raise ValueError("H or h must be provided to simulate")

        row['flag'], row['RA'], row['DEC'], row['d_ra'], row['d_dec'], row['r'], row['delta'], \
            row['m_int'], row['m_rand'], row['eff'], isur, row['Mt'], jdayp, ic, row['Survey'], \
            row['h_rand'], ierr = \
            ossssimlib.surveysub.detos1(o_m,
                                           epoch_jd,
                                           H,
                                           color_offset_array,
                                           gb,
                                           phase,
                                           period,
                                           amplitude,
                                           self.characterization_directory,
                                           self.seed,
                                           debug)
        if ierr != 0:
            raise IOError(f"SSim failed with error code: {ierr}")
        row['j'] = row.get('j', 0)
        row['k'] = row.get('k', 0)
        row['delta'] = Quantity(row['delta'], u.au)
        row['r'] = Quantity(row['r'], u.au)
        row['m_int'] = Quantity(row['m_int'], u.mag)
        row['m_rand'] = Quantity(row['m_rand'], u.mag)
        row['h_rand'] = Quantity(row['h_rand'], u.mag)
        # ic is the 1-based Fortran filter_to_index of the discovery bandpass.
        # Convert to band letter and look up the matching 0-based colors_list slot.
        if ic > 0:
            band = chr(ic + ord('A') - 1)
            row['band'] = band
            row['color'] = Quantity(color_offset_array[ord(band) - ord('A')], u.mag)
        else:
            row['band'] = ''
            row['color'] = Quantity(0, u.mag)
        row['q'] = row['a'] * (1 - row['e'])
        row['Mt'] = Quantity(row['Mt'], u.rad)
        row['RA'] = Quantity(row['RA'], u.rad)
        row['DEC'] = Quantity(row['DEC'], u.rad)
        row['Survey'] = row['Survey'].decode('utf-8')
        if row['flag'] > 0:
            row.update(self.cartesian(**row))

        # m_int and h are in "x" band (filter of object creation)
        # m_rand and h_rand are in discovery filter
        return row

    def simulate_batch(self, columns, colors: PhotSpec, model_band, epoch=None, debug=False) -> dict:
        """
        Run detos1 on many targets with one call into the Fortran library.

        Args:
            columns (dict or QTable): equal-length columns with the names simulate() takes for a row:
                a, e, inc, node, peri, M, H (or h) and comp, and optionally epoch, gb, phase, period,
                amplitude, j and k. Quantity columns are converted to the units the Fortran code uses;
                plain arrays are taken to be in those units already (au, rad, mag, day).
            colors (PhotSpec): colors of the targets, see simulate().
            model_band (str): bandpass of the H values.
            epoch (Quantity or float or Time): JD of the elements when there is no epoch column.
            debug (bool): if True, have the Fortran detos routine log the detection process.

        Returns:
            dict: the input columns plus the columns simulate() adds, as arrays. x, y, z, vx, vy, vz
            are computed for detected targets (flag > 0) and are NaN otherwise.

        The values are the same as calling simulate() on each row in order.
        """
        result = {name: columns[name] for name in columns.keys()}
        n = len(result['a'])
        if n == 0:
            raise ValueError("simulate_batch needs at least one target")

        names = {}
        for name in result:
            if name.lower() in T_ORB_M_FIELDS:
                names[name.lower()] = name
        elements = []
        for field in T_ORB_M_FIELDS:
            if field in names:
                value = _to_value(result[names[field]], T_ORB_M_UNITS.get(names[field]))
                elements.append(numpy.ascontiguousarray(numpy.broadcast_to(value, (n,)), dtype=numpy.float64))
            else:
                elements.append(numpy.zeros(n))

        epoch_jd = result['epoch'] if 'epoch' in result else epoch
        if isinstance(epoch_jd, Quantity):
            epoch_jd = _to_value(epoch_jd, u.day)
        if isinstance(epoch_jd, Time):
            epoch_jd = epoch_jd.jd
        jday = numpy.ascontiguousarray(numpy.broadcast_to(epoch_jd, (n,)), dtype=numpy.float64)

        light_curve = []
        for key in ('gb', 'phase', 'period', 'amplitude'):
            default = _LIGHT_CURVE_DEFAULTS[key]
            if key in result:
                value = numpy.asarray(_to_value(result[key], _LIGHT_CURVE_UNITS[key]), dtype=numpy.float64)
                value = numpy.where(value != 0, value, default)
            else:
                value = default
            light_curve.append(numpy.ascontiguousarray(numpy.broadcast_to(value, (n,)), dtype=numpy.float64))

        H = numpy.zeros(n)
        if 'H' in result:
            H = numpy.array(_to_value(result['H'], u.mag), dtype=numpy.float64)
        if 'h' in result:
            h = numpy.asarray(_to_value(result['h'], u.mag), dtype=numpy.float64)
            H = numpy.where(H != 0, H, h)
        if numpy.any(H == 0):
            raise ValueError("H or h must be provided to simulate")

        comps, comp_index = numpy.unique(numpy.asarray(result['comp']).astype(str), return_inverse=True)
        comp_colors = numpy.stack([colors.colors_array(str(comp), model_band) for comp in comps], axis=1)
        color = numpy.asfortranarray(comp_colors[:, comp_index.reshape(-1)])

        flag, isur, ic = (numpy.zeros(n, dtype=numpy.int32) for _ in range(3))
        ra, dec, d_ra, d_dec, r, delta, m_int, m_rand, eff, mt, jdayp, h_rand = \
            (numpy.zeros(n) for _ in range(12))
        survey_codes = numpy.zeros((_SURVEY_KEY_CODES, n), dtype=numpy.int32, order='F')
        ierr = ossssimlib.surveysub.detos_batch(n, _SURVEY_KEY_CODES, *elements, jday, H, color,
                                                *light_curve, self.characterization_directory,
                                                self.seed, debug,
                                                flag, ra, dec, d_ra, d_dec, r, delta, m_int, m_rand,
                                                eff, isur, mt, jdayp, ic, survey_codes, h_rand)
        if ierr != 0:
            raise IOError(f"SSim failed with error code: {ierr}")

        result['flag'] = flag
        result['RA'] = Quantity(ra, u.rad)
        result['DEC'] = Quantity(dec, u.rad)
        result['d_ra'] = d_ra
        result['d_dec'] = d_dec
        result['r'] = Quantity(r, u.au)
        result['delta'] = Quantity(delta, u.au)
        result['m_int'] = Quantity(m_int, u.mag)
        result['m_rand'] = Quantity(m_rand, u.mag)
        result['eff'] = eff
        result['Mt'] = Quantity(mt, u.rad)
        result['h_rand'] = Quantity(h_rand, u.mag)
        result['Survey'] = numpy.array([bytes(code).decode('utf-8').rstrip(' \x00')
                                        for code in survey_codes.T.astype(numpy.uint8)])
        result['j'] = result.get('j', numpy.zeros(n, dtype=int))
        result['k'] = result.get('k', numpy.zeros(n, dtype=int))
        detected_band = ic > 0
        result['band'] = numpy.where(detected_band, [chr(c + ord('A') - 1) if c > 0 else '' for c in ic], '')
        result['color'] = Quantity(numpy.where(detected_band, color[numpy.maximum(ic - 1, 0), numpy.arange(n)], 0.0),
                                   u.mag)
        result['q'] = result['a'] * (1 - result['e'])

        cartesian = {s: numpy.full(n, numpy.nan) for s in ('x', 'y', 'z', 'vx', 'vy', 'vz')}
        for i in numpy.flatnonzero(flag > 0):
            state = self.cartesian(a=result['a'][i], e=result['e'][i], inc=result['inc'][i],
                                   node=result['node'][i], peri=result['peri'][i], Mt=result['Mt'][i])
            for s in cartesian:
                cartesian[s][i] = state[s].value
        for s in ('x', 'y', 'z'):
            result[s] = Quantity(cartesian[s], u.au)
        for s in ('vx', 'vy', 'vz'):
            result[s] = Quantity(cartesian[s], u.au / u.year)
        return result
