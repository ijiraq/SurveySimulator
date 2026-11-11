import os
import pathlib
import unittest
from tempfile import NamedTemporaryFile

import numpy
from astropy.units import Quantity

import ossssim

script_directory = pathlib.Path(__file__).parent.resolve()
CFEPS_DIR = script_directory / 'data' / 'Surveys' / 'CFEPS'
SEED = 123456789
COMPARED = ['flag', 'RA', 'DEC', 'd_ra', 'd_dec', 'r', 'delta', 'm_int', 'm_rand', 'eff', 'Mt',
            'h_rand', 'color', 'q', 'x', 'y', 'z', 'vx', 'vy', 'vz']


def _values(value):
    return value.value if isinstance(value, Quantity) else numpy.asarray(value)


class SimulateBatchTest(unittest.TestCase):
    """simulate_batch over iter_chunks must give the same results as simulate() row by row."""

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

    def test_batch_matches_simulate(self):
        model = ossssim.ModelFile(self.model_filename)
        sim = ossssim.OSSSSim(characterization_directory=str(CFEPS_DIR), seed=SEED)
        per_row = [sim.simulate(row, colors=model.colors, model_band=model.model_band, epoch=model.epoch)
                   for row in model]
        model.close()

        model = ossssim.ModelFile(self.model_filename)
        sim = ossssim.OSSSSim(characterization_directory=str(CFEPS_DIR), seed=SEED)
        chunks = [sim.simulate_batch(columns, colors=model.colors, model_band=model.model_band,
                                     epoch=model.epoch)
                  for columns in model.iter_chunks(250)]
        model.close()

        self.assertEqual(len(chunks), 3)
        flags = numpy.concatenate([chunk['flag'] for chunk in chunks])
        self.assertEqual(len(flags), len(per_row))
        self.assertGreater((flags > 0).sum(), 0)
        for name in COMPARED:
            expected = numpy.array([float(_values(row.get(name, numpy.nan))) for row in per_row])
            batch = numpy.concatenate([_values(chunk[name]) for chunk in chunks]).astype(float)
            numpy.testing.assert_array_equal(batch, expected, err_msg=name)
        for name in ('Survey', 'band', 'comp'):
            expected = [row[name] for row in per_row]
            batch = list(numpy.concatenate([chunk[name] for chunk in chunks]))
            self.assertEqual(batch, expected, msg=name)

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


class StreamedWriteTest(unittest.TestCase):
    """write_row streams rows after the first; the file must read back with masked non-detections."""

    def test_detect_file_round_trip(self):
        model_file = script_directory / 'data' / 'test_model.dat'
        model = ossssim.ModelFile(str(model_file))
        sim = ossssim.OSSSSim(characterization_directory=str(CFEPS_DIR), seed=SEED)
        with NamedTemporaryFile(suffix='.ecsv') as fobj:
            filename = fobj.name
        writer = ossssim.DetectFile(filename, seed=SEED, epoch=model.epoch,
                                    longitude_neptune=model.longitude_neptune,
                                    colors=model.colors, model_band=model.model_band)
        written = []
        for i in range(20):
            for row in ossssim.ModelFile(str(model_file)):
                result = sim.simulate(row, colors=model.colors, model_band=model.model_band,
                                      epoch=model.epoch)
                for name in ('x', 'y', 'z'):
                    result.setdefault(name, Quantity(0.0, 'au'))
                writer.write_row(result)
                written.append(result)
        writer.write_footer(n_iter=len(written), n_hits=0, n_track=0)
        table = ossssim.ModelFile(filename).table
        os.remove(filename)

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


if __name__ == '__main__':
    unittest.main()
