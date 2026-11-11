"""
Expose a model file to a survey characterization, simulating the model in chunks.

simulate() accepts a single row or a chunk of rows (a QTable such as those from
ModelFile.iter_chunks, a table slice, or a dict of arrays). A chunk is simulated with
one call into the Fortran code, which is about ten times faster than calling
simulate() on each row, and gives the same results.

Run from the top of the repository:

    python examples/read_model_in_chunks.py data/Models/L7model-3.0-9.0 data/Surveys/CFEPS
"""
import argparse

import ossssim

TRACKED = 4


def run(model_filename: str, characterization_dir: str, output: str, seed: int = 123456789,
        chunk_size: int = 10000) -> None:
    model = ossssim.ModelFile(model_filename)
    ssim = ossssim.OSSSSim(characterization_directory=characterization_dir, seed=seed)
    detect_file = ossssim.DetectFile(output, seed=seed, epoch=model.epoch,
                                     longitude_neptune=model.longitude_neptune,
                                     colors=model.colors, model_band=model.model_band)

    n_iter = n_hits = n_track = 0
    for chunk in model.iter_chunks(chunk_size):
        result = ssim.simulate(chunk, colors=model.colors, model_band=model.model_band, epoch=model.epoch)
        detected = result[result['flag'] > 0]
        detect_file.write_row(detected)
        n_iter += len(result)
        n_hits += len(detected)
        n_track += int((detected['flag'] == TRACKED).sum())
    detect_file.write_footer(n_iter=n_iter, n_hits=n_hits, n_track=n_track)
    model.close()
    print(f"{n_iter} objects, {n_hits} detected, {n_track} tracked: {output}")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('model', help="SSim model file with a commented column-name header line")
    parser.add_argument('characterization', help="Survey characterization directory")
    parser.add_argument('--output', default='detect.ecsv', help="File for the detected objects")
    parser.add_argument('--seed', type=int, default=123456789)
    parser.add_argument('--chunk-size', type=int, default=10000)
    args = parser.parse_args()
    run(args.model, args.characterization, args.output, seed=args.seed, chunk_size=args.chunk_size)


if __name__ == '__main__':
    main()
