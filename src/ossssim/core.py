"""
The OSSOS Survey Simulator module.  This is front-end for the OSSSim.  The primary OSSSim is writen in f95.
The OSSOS Survey Simulator module.  This is front-end for the OSSSim.  The primary OSSSim is writen in f95.
"""
import logging
from astropy import units as u
from astropy.units import Quantity
from astropy.table import Column, Row, Table
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


def _as_quantity(value):
    """A Table Column that carries a unit as a Quantity; anything else unchanged."""
    if isinstance(value, Column) and value.unit is not None:
        return value.quantity
    return value


def _to_value(value, unit):
    """
    Value of a Quantity (or a Column with a unit) in unit; plain numbers are returned unchanged.

    Same result as value.to(unit).value, but the unit scale factor is looked up
    once per pair of units instead of on every call.
    """
    value = _as_quantity(value)
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


_OUTPUT_UNITS = {'RA': u.rad, 'DEC': u.rad, 'r': u.au, 'delta': u.au, 'm_int': u.mag,
                 'm_rand': u.mag, 'Mt': u.rad, 'h_rand': u.mag}
_CARTESIAN_UNITS = {'x': u.au, 'y': u.au, 'z': u.au,
                    'vx': u.au / u.year, 'vy': u.au / u.year, 'vz': u.au / u.year}


def _row_units(row) -> dict:
    """Units of the unit-bearing plain Columns behind an astropy Row ({} for anything else)."""
    if not isinstance(row, Row):
        return {}
    return {name: column.unit for name, column in row.table.columns.items()
            if isinstance(column, Column) and column.unit is not None}


def _survey_key(codes) -> str:
    """Detection key from the character codes Detos_batch returns for one object."""
    return codes.astype(numpy.uint8).tobytes().decode('utf-8').rstrip(' \x00')


def _is_chunk(targets) -> bool:
    """True for a table or a mapping of arrays, False for a single row."""
    if isinstance(targets, Table):
        return True
    if isinstance(targets, Row):
        return False
    try:
        return numpy.ndim(targets['a']) > 0
    except (KeyError, TypeError, IndexError):
        return False


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
        a = rows['a'].to_value(u.au)
        inc = rows['inc'].to_value(u.rad)
        node = rows['node'].to_value(u.rad)
        peri = rows['peri'].to_value(u.rad)
        mean_anomaly = rows[M].to_value(u.rad)
        for idx in range(len(a)):
            sim.add(a=a[idx], e=rows['e'][idx], inc=inc[idx], Omega=node[idx], omega=peri[idx],
                    M=mean_anomaly[idx])

        particles = sim.particles[5:]
        for s in ['x', 'y', 'z', 'vx', 'vy', 'vz']:
            rows[s] = Quantity([p.__getattribute__(s) for p in particles], definitions.column_unit[s])
        return rows

    def single_row(self, **kwargs):
        """
        Provide the state vector of the orbits at Mt or M if Mt is not provided.

        Returns:
            dictionary of x/y/z/vx/vy/vz
        """
        state = self.state(kwargs['a'].to('au').value,
                           kwargs['e'],
                           kwargs['inc'].to('rad').value,
                           kwargs['node'].to('rad').value,
                           kwargs['peri'].to('rad').value,
                           kwargs.get('Mt', kwargs.get('M')).to('rad').value)
        cartesian = {}
        for s in ['x', 'y', 'z']:
            cartesian[s] = Quantity(state[s], u.au)
        for s in ['vx', 'vy', 'vz']:
            cartesian[s] = Quantity(state[s], u.au / u.year)
        return cartesian

    def state(self, a, e, inc, node, peri, M) -> dict:
        """
        Cartesian state of one orbit as plain floats: x/y/z [au] and vx/vy/vz [au/yr].

        Elements are in au and radians.
        """
        sim = self.initialize_sim()
        sim.add(a=a, e=e, inc=inc, Omega=node, omega=peri, M=M)
        particle = sim.particles[-1]
        state = {s: particle.__getattribute__(s) for s in ['x', 'y', 'z']}
        for s in ['vx', 'vy', 'vz']:
            state[s] = particle.__getattribute__(s) * (2 * numpy.pi)
        return state


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
        self._epoch_jd = (None, None)
        self._row_buffers = None


    def simulate(self, row: dict, colors: PhotSpec, model_band, seed=None, epoch=None, debug=False) -> dict:
        """
        Pass the target elements to detos1 and determine if this target would be detected.

        Pass the survey root once; Detos1 loads all epochs and applies
        ``detections_required`` from ``survey.conf``. Changing
        ``self.characterization_directory`` between calls reloads that survey
        root without resetting the ran3 stream. Do not construct a new
        OSSSSim per epoch: ``__init__`` calls ``reset_simulator()`` and would
        reseed the RNG.

        row can also be a chunk of targets: a Table/QTable (e.g. ``table[0:1000]`` or a chunk from
        ``model.iter_chunks()``) or a dict of equal-length arrays. A chunk is simulated with a single
        call into the Fortran code, which is much faster than one call per row, and gives the same
        values as simulating its rows one at a time in order. The result is the same kind of object:
        a QTable chunk returns a QTable, a dict of arrays returns a dict of arrays. In a chunk,
        x/y/z/vx/vy/vz of undetected targets keep their input values, or are NaN if the input had none.

        Args:
            row (Row or dict or Table): elements of the target (or targets) to simulate
            colors (PhotSpec): colors of the target index by 'g-x' etc, see below.
            model_band (str): a single letter designator of the bandpass that the H values are provided in.
            seed (int): a seed to pass to the fortran code to allow reproducible simulations
            epoch (Quantity or float or Time): JD of elements, can also be provided for each row sent to simulate
            debug (bool): if True, have Fortran detos routine create a log of detection process.

        Returns:
            dict (or chunk): the target elements/values at time of simulated detection.

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

        if _is_chunk(row):
            return self._simulate_chunk(row, colors, model_band, epoch, debug)
        return self._simulate_row(row, colors, model_band, epoch, debug)

    def _simulate_row(self, row, colors: PhotSpec, model_band, epoch, debug) -> dict:
        """simulate() for a single row: the row is passed to Fortran as a chunk of one."""
        units = _row_units(row)
        row = dict(row)
        # A Row of a plain Table gives bare numbers; put its columns' units back for the conversions
        # below, while the returned row keeps the values it was given.
        values = row if not units else {name: Quantity(value, units[name]) if name in units else value
                                        for name, value in row.items()}
        elements = {}
        for colname, value in values.items():
            field = colname.lower()
            if field in T_ORB_M_FIELDS:
                elements[field] = _to_value(value, T_ORB_M_UNITS.get(colname))
        H = None
        if 'H' in values:
            H = _to_value(values['H'], u.mag) or None
        if H is None and 'h' in values:
            H = _to_value(values['h'], u.mag) or None
        if H is None:
            raise ValueError("H or h must be provided to simulate")

        # One-element input and output arrays, reused for every single-row call.
        if self._row_buffers is None:
            self._row_buffers = ([numpy.zeros(1) for _ in T_ORB_M_FIELDS], numpy.zeros(1), numpy.zeros(1),
                                 [numpy.zeros(1) for _ in _LIGHT_CURVE_UNITS], self._output_buffers(1))
        element_buffers, jday, H_buffer, light_curve, out = self._row_buffers
        for buffer, field in zip(element_buffers, T_ORB_M_FIELDS):
            buffer[0] = elements.get(field, 0.0)
        jday[0] = self._jd(values['epoch'] if 'epoch' in values else epoch)
        H_buffer[0] = H
        for buffer, key in zip(light_curve, _LIGHT_CURVE_UNITS):
            buffer[0] = _light_curve_value(values, key)
        color = colors.colors_array(row['comp'], model_band).reshape(-1, 1)

        self._detos(element_buffers, jday, H_buffer, color, light_curve, debug, out)
        flag = int(out['flag'][0])
        ic = int(out['ic'][0])
        row['flag'] = flag
        for name in ('RA', 'DEC', 'd_ra', 'd_dec', 'r', 'delta', 'm_int', 'm_rand', 'eff', 'Mt'):
            unit = _OUTPUT_UNITS.get(name)
            row[name] = float(out[name][0]) if unit is None else Quantity(out[name][0], unit)
        row['Survey'] = _survey_key(out['survey_codes'][:, 0])
        row['h_rand'] = Quantity(out['h_rand'][0], u.mag)
        row['j'] = row.get('j', 0)
        row['k'] = row.get('k', 0)
        # ic is the 1-based Fortran filter_to_index of the discovery bandpass: convert it to the
        # band letter and look up the matching 0-based colors_list slot.
        row['band'] = chr(ic + ord('A') - 1) if ic > 0 else ''
        row['color'] = Quantity(color[ic - 1, 0] if ic > 0 else 0.0, u.mag)
        row['q'] = row['a'] * (1 - row['e'])
        if flag > 0:
            state = self.cartesian.state(*(buffer[0] for buffer in element_buffers[:5]), out['Mt'][0])
            for s, unit in _CARTESIAN_UNITS.items():
                row[s] = Quantity(state[s], unit)

        # m_int and h are in "x" band (filter of object creation)
        # m_rand and h_rand are in discovery filter
        return row

    def _simulate_chunk(self, chunk, colors: PhotSpec, model_band, epoch, debug):
        """simulate() for a chunk of rows; returns the same kind of container it was given."""
        result = {name: chunk[name] for name in chunk.keys()}
        n = len(result['a'])
        if n == 0:
            raise ValueError("simulate needs at least one target")

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
        jday = numpy.ascontiguousarray(numpy.broadcast_to(self._jd(epoch_jd), (n,)), dtype=numpy.float64)

        light_curve = []
        for key, unit in _LIGHT_CURVE_UNITS.items():
            default = _LIGHT_CURVE_DEFAULTS[key]
            if key in result:
                value = numpy.asarray(_to_value(result[key], unit), dtype=numpy.float64)
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

        out = self._output_buffers(n)
        self._detos(elements, jday, H, color, light_curve, debug, out)
        flag, ic = out['flag'], out['ic']
        detected = flag > 0
        result['flag'] = flag
        for name in ('RA', 'DEC', 'd_ra', 'd_dec', 'r', 'delta', 'm_int', 'm_rand', 'eff', 'Mt'):
            unit = _OUTPUT_UNITS.get(name)
            result[name] = out[name] if unit is None else Quantity(out[name], unit)
        result['Survey'] = numpy.array([_survey_key(codes) for codes in out['survey_codes'].T])
        result['h_rand'] = Quantity(out['h_rand'], u.mag)
        result['j'] = result.get('j', numpy.zeros(n, dtype=int))
        result['k'] = result.get('k', numpy.zeros(n, dtype=int))
        result['band'] = numpy.array([chr(c + ord('A') - 1) if c > 0 else '' for c in ic])
        result['color'] = Quantity(numpy.where(ic > 0, color[numpy.maximum(ic - 1, 0), numpy.arange(n)], 0.0),
                                   u.mag)
        result['q'] = result['a'] * (1 - result['e'])

        state = {s: numpy.full(n, numpy.nan) for s in _CARTESIAN_UNITS}
        for i in numpy.flatnonzero(detected):
            for s, value in self.cartesian.state(*(element[i] for element in elements[:5]),
                                                 out['Mt'][i]).items():
                state[s][i] = value
        for s, unit in _CARTESIAN_UNITS.items():
            values = state[s]
            if s in result:
                # Keep the model's own state vector for objects that were not detected.
                values = numpy.where(detected, values, numpy.asarray(_to_value(result[s], unit), dtype=numpy.float64))
            result[s] = Quantity(values, unit)

        if isinstance(chunk, Table):
            return type(chunk)(result, meta=chunk.meta)
        return result

    def _jd(self, epoch):
        """JD of an epoch given as Time, Quantity or a plain number (or arrays of these)."""
        epoch = _as_quantity(epoch)
        if isinstance(epoch, Quantity):
            return _to_value(epoch, u.day)
        if isinstance(epoch, Time):
            if self._epoch_jd[0] is not epoch:
                self._epoch_jd = (epoch, epoch.jd)
            return self._epoch_jd[1]
        return epoch

    @staticmethod
    def _output_buffers(n) -> dict:
        """Arrays for the Detos_batch outputs of n objects."""
        out = {name: numpy.zeros(n, dtype=numpy.int32) for name in ('flag', 'isur', 'ic')}
        for name in ('RA', 'DEC', 'd_ra', 'd_dec', 'r', 'delta', 'm_int', 'm_rand', 'eff', 'Mt', 'jdayp', 'h_rand'):
            out[name] = numpy.zeros(n)
        out['survey_codes'] = numpy.zeros((_SURVEY_KEY_CODES, n), dtype=numpy.int32, order='F')
        return out

    def _detos(self, elements, jday, H, color, light_curve, debug, out):
        """Run Detos_batch on plain arrays (au, radians, mag, days), filling the arrays in out."""
        ierr = ossssimlib.surveysub.detos_batch(len(jday), _SURVEY_KEY_CODES, *elements, jday, H, color,
                                                *light_curve, self.characterization_directory,
                                                self.seed, debug,
                                                out['flag'], out['RA'], out['DEC'], out['d_ra'], out['d_dec'],
                                                out['r'], out['delta'], out['m_int'], out['m_rand'], out['eff'],
                                                out['isur'], out['Mt'], out['jdayp'], out['ic'],
                                                out['survey_codes'], out['h_rand'])
        if ierr != 0:
            raise IOError(f"SSim failed with error code: {ierr}")
