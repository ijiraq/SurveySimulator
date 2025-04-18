import unittest
from ossssim.models import ModelFile
from astropy import units
from astropy.time import Time


class SSimModelFileTest(unittest.TestCase):

    def setUp(self):
        self.epoch = Time(val=2453157.50000, format='jd')
        self.longitude_neptune = 5.489 * units.radian
        self.colors = [float(x.replace("d","e")) for
                       x in "0.0d0 -0.70d0 -1.2d0 -1.7d0 0.8d0 0.5d0  0.1d0 -0.8d0 -1.2d0 0.0d0".split()] * units.mag
        self.colnames = "  a      e     inc     node    peri    M       H      dist     comp      j  k".split()
        self.row = {'a':     88.253*units.au,
                    'e':      0.564,
                    'inc':    8.807*units.degree,
                    'node': 316.126*units.degree,
                    'peri': 286.639*units.degree,
                    'M':     60.537*units.degree,
                    'H':      8.05*units.mag,
                    'dist':  90.70*units.au,
                    'comp': 'resonant',
                    'j':     5,
                    'k':     1}
        self.model = ModelFile('data/test_model.dat')

    def test_epoch(self):
        self.assertAlmostEqual(self.epoch, self.model.epoch)

    def test_lambda_neptune(self):
        self.assertAlmostEqual(self.longitude_neptune, self.model.longitude_neptune)

    def test_colors(self):
        for idx, key in enumerate(list(self.model.colors.colors['default'])):
            if idx > len(self.colors)-1:
                break
            self.assertEqual(self.colors[idx], self.model.colors.colors['default'][key])

    def test_colnames(self):
        self.assertEqual(self.colnames, self.model.colnames)

    def test_row(self):
        for row in self.model:
            for key in self.row:
                self.assertEqual(self.row[key], row[key])
            break

    def test_targets(self):
        self.assertEqual(len(self.model.targets), 2)
        for row in self.model:
            for key in self.row:
                self.assertEqual(self.row[key], row[key])
            break

    def tearDown(self):
        self.model.close()


if __name__ == '__main__':
    unittest.main()
