import unittest
from astropy.units import Quantity
import ossssim
from tempfile import NamedTemporaryFile


class OSSSIMTest(unittest.TestCase):
    # def setUp(self):

    def test_simulate(self):
        self.model = ossssim.ModelFile('data/test_model.dat')
        result = ossssim.ModelFile('data/test_detect.dat')
        self.result_row = next(iter(result))
        result.close()
        self.seed = int(result.header['Seed'][0])
        self.osssim = ossssim.OSSSSim('data/Surveys/CFEPS')

        with NamedTemporaryFile() as fobj:
            self.detect_filename = fobj.name
            temp_results = ossssim.DetectFile(self.detect_filename,
                                              seed=self.seed,
                                              epoch=self.model.epoch,
                                              longitude_neptune=self.model.longitude_neptune,
                                              colors=self.model.colors)

        # loop over the model file until we have a detection and then compare the detected row values to the test row
        for row in self.model:
            result_row = self.osssim.simulate(row, seed=self.seed, epoch=self.model.epoch,
                                              colors=self.model.colors, model_band=self.model.model_band_pass)
            if result_row['flag'] > 0:
                temp_results.write_row(result_row)
                for key in row:
                    test_value = row[key]
                    result_value = result_row[key]
                    if isinstance(test_value, Quantity):
                        test_value = test_value.to(ossssim.definitions.column_unit[key]).value
                        result_value = result_value.to(ossssim.definitions.column_unit[key]).value
                    self.assertAlmostEqual(test_value, result_value, 4)
                break
        self.model.close()
        temp_results.write_footer(n_iter=1, n_hits=1, n_track=0)
        results2 = ossssim.ModelFile(self.detect_filename)
        for row in results2:
            pass
        self.assertEqual(len(results2.table), 1)
        self.assertAlmostEqual(results2.colors.colors['default']['r-g'],
                               result.colors.colors['default']['r-g'])
        self.assertEqual(results2.seed, result.seed)
        results2.close()

if __name__ == '__main__':
    unittest.main()
