import os
import pathlib
import unittest
from tempfile import NamedTemporaryFile

import numpy
from astropy import units
from astropy.table import QTable
from astropy.units import Quantity

import ossssim
from ossssim.models import Resonant

script_directory = pathlib.Path(__file__).parent.resolve()
CFEPS_DIR = script_directory / 'data' / 'Surveys' / 'CFEPS'
SEED = 123456789
COMPARED = ['flag', 'RA', 'DEC', 'd_ra', 'd_dec', 'r', 'delta', 'm_int', 'm_rand', 'eff', 'Mt',
            'h_rand', 'color', 'q', 'x', 'y', 'z', 'vx', 'vy', 'vz']


def _values(value):
    return value.value if isinstance(value, Quantity) else numpy.asarray(value)


class ChunkedSimulateTest(unittest.TestCase):
    """simulate() on a chunk must give the same results as simulate() on its rows in order."""

    def setUp(self):
        # The two test_model.dat orbits at many mean anomalies, so some of them are detected.
        header = [line for line in open(script_directory / 'data' / 'test_model.dat') if line.startswith('#')]
        rows = [line.split() for line in open(script_directory / 'data' / 'test_model.dat')
                if not line.startswith('#')]
        tmp = NamedTemporaryFile('w', suffix='.dat', delete=False)
        tmp.writelines(header)
        for i in range(300):
            for row in rows:
                row = list(row)
                row[5] = f"{(float(row[5]) + 1.2 * i) % 360:.3f}"
                tmp.write(" ".join(row) + "\n")
        tmp.close()
        self.model_filename = tmp.name

    def tearDown(self):
        os.remove(self.model_filename)

    def _per_row(self):
        model = ossssim.ModelFile(self.model_filename)
        sim = ossssim.OSSSSim(characterization_directory=str(CFEPS_DIR), seed=SEED)
        results = [sim.simulate(row, colors=model.colors, model_band=model.model_band, epoch=model.epoch)
                   for row in model]
        model.close()
        return results

    def _chunked(self, as_dict=False):
        model = ossssim.ModelFile(self.model_filename)
        sim = ossssim.OSSSSim(characterization_directory=str(CFEPS_DIR), seed=SEED)
        results = []
        for chunk in model.iter_chunks(250):
            self.assertIsInstance(chunk, QTable)
            if as_dict:
                chunk = {name: chunk[name] for name in chunk.colnames}
            result = sim.simulate(chunk, colors=model.colors, model_band=model.model_band, epoch=model.epoch)
            self.assertIsInstance(result, dict if as_dict else QTable)
            results.append(result)
        model.close()
        return results

    def _assert_same(self, per_row, chunks):
        self.assertEqual(len(chunks), 3)
        flags = numpy.concatenate([_values(chunk['flag']) for chunk in chunks])
        self.assertEqual(len(flags), len(per_row))
        self.assertGreater((flags > 0).sum(), 0)
        for name in COMPARED:
            expected = numpy.array([float(_values(row.get(name, numpy.nan))) for row in per_row])
            got = numpy.concatenate([_values(chunk[name]) for chunk in chunks]).astype(float)
            numpy.testing.assert_array_equal(got, expected, err_msg=name)
        for name in ('Survey', 'band', 'comp'):
            expected = [row[name] for row in per_row]
            got = [str(value) for chunk in chunks for value in chunk[name]]
            self.assertEqual(got, expected, msg=name)

    def test_qtable_chunks_match_rows(self):
        self._assert_same(self._per_row(), self._chunked())

    def test_dict_chunks_match_rows(self):
        self._assert_same(self._per_row(), self._chunked(as_dict=True))

    def test_single_row_types(self):
        types = {'flag': int, 'd_ra': float, 'd_dec': float, 'eff': float, 'Survey': str, 'band': str,
                 'RA': Quantity, 'DEC': Quantity, 'r': Quantity, 'delta': Quantity, 'm_int': Quantity,
                 'm_rand': Quantity, 'h_rand': Quantity, 'Mt': Quantity, 'color': Quantity, 'q': Quantity}
        per_row = self._per_row()
        for row in per_row:
            for name, kind in types.items():
                self.assertIs(type(row[name]), kind, msg=name)
            self.assertEqual('x' in row, row['flag'] > 0)
        detected = next(row for row in per_row if row['flag'] > 0)
        self.assertEqual(detected['x'].unit, units.au)
        self.assertEqual(detected['vx'].unit, units.au / units.year)

    def test_chunks_match_rows(self):
        model = ossssim.ModelFile(self.model_filename)
        rows = list(model)
        model.close()
        model = ossssim.ModelFile(self.model_filename)
        chunks = list(model.iter_chunks(250))
        model.close()
        for name in rows[0]:
            expected = [_values(row[name]).item() for row in rows]
            got = [value.item() for chunk in chunks for value in numpy.atleast_1d(_values(chunk[name]))]
            self.assertEqual(got, expected, msg=name)


class ParametricChunksTest(unittest.TestCase):

    def test_chunks_keep_model_state_vectors(self):
        model = Resonant(j=3, k=2, seed=42, size=500)
        sim = ossssim.OSSSSim(characterization_directory=str(CFEPS_DIR), seed=SEED)
        chunks = model.iter_chunks(200)
        sizes = [len(next(chunks)) for _ in range(4)]
        self.assertEqual(sizes, [200, 200, 100, 200])

        chunk = next(Resonant(j=3, k=2, seed=42, size=500).iter_chunks(200))
        result = sim.simulate(chunk, colors=model.colors, model_band=model.model_band, epoch=model.epoch)
        undetected = result['flag'] == 0
        numpy.testing.assert_array_equal(result['x'][undetected].to_value('au'),
                                         chunk['x'][undetected].to_value('au'))


class StreamedWriteTest(unittest.TestCase):
    """write_row streams rows after the first; the file must read back with masked non-detections."""

    def _results(self):
        model_file = script_directory / 'data' / 'test_model.dat'
        model = ossssim.ModelFile(str(model_file))
        sim = ossssim.OSSSSim(characterization_directory=str(CFEPS_DIR), seed=SEED)
        written = []
        for i in range(20):
            for row in ossssim.ModelFile(str(model_file)):
                result = sim.simulate(row, colors=model.colors, model_band=model.model_band,
                                      epoch=model.epoch)
                for name in ('x', 'y', 'z'):
                    result.setdefault(name, Quantity(0.0, 'au'))
                written.append(result)
        return model, written

    def _read_back(self, model, write):
        with NamedTemporaryFile(suffix='.ecsv') as fobj:
            filename = fobj.name
        writer = ossssim.DetectFile(filename, seed=SEED, epoch=model.epoch,
                                    longitude_neptune=model.longitude_neptune,
                                    colors=model.colors, model_band=model.model_band)
        write(writer)
        writer.write_footer(n_iter=0, n_hits=0, n_track=0)
        table = ossssim.ModelFile(filename).table
        lines = open(filename).readlines()
        os.remove(filename)
        return table, [line for line in lines if 'Creation_time' not in line]

    def test_detect_file_round_trip(self):
        model, written = self._results()
        table, _ = self._read_back(model, lambda writer: [writer.write_row(row) for row in written])

        self.assertEqual(len(table), len(written))
        self.assertGreater(sum(result['flag'] > 0 for result in written), 0)
        a = table['a'].unmasked.to_value('au')
        ra = table['RA'].unmasked.to_value('rad')
        for i, result in enumerate(written):
            self.assertEqual(table['flag'][i], result['flag'])
            self.assertAlmostEqual(a[i], result['a'].to_value('au'), places=4)
            if result['flag'] > 0:
                self.assertEqual(table['Survey'][i], result['Survey'])
                self.assertAlmostEqual(ra[i], result['RA'].to_value('rad'), places=4)
            else:
                self.assertTrue(table['RA'].mask[i])

    def test_write_table_matches_rows(self):
        model, written = self._results()
        names = list(written[0])
        chunk = QTable({name: [row[name] for row in written] for name in names})
        _, by_row = self._read_back(model, lambda writer: [writer.write_row(row) for row in written])
        _, by_table = self._read_back(model, lambda writer: writer.write_row(chunk))
        self.assertEqual(by_table, by_row)


if __name__ == '__main__':
    unittest.main()
