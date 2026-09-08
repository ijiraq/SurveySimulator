"""
Example: basic RosePlot of a model and detections.
"""
from pathlib import Path

import ossssim
from ossssim import plotter
from astropy.time import Time

REPO = Path(__file__).resolve().parents[1]


def run():
    model_path = REPO / 'F95' / 'tests' / 'Models' / 'L7model-3.0-9.0'
    if not model_path.is_file():
        # Fall back to sibling SurveySimulator-Data if present
        model_path = REPO.parent / 'SurveySimulator-Data' / 'Models' / 'L7model-3.0-9.0'
    model = ossssim.ModelFile(str(model_path))

    plot_driver = plotter.RosePlot(epoch=Time(model.epoch))
    plot_driver.add_scale_rings()
    plot_driver.add_model(model, sample_size=10**4)
    plot_driver.show()


if __name__ == '__main__':
    run()
