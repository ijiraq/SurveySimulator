# setup.py
from distutils.dir_util import remove_tree
from glob import glob
import os
import pathlib
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import shutil
import subprocess
import sys

# Remove the build directory at the start

if pathlib.PosixPath("build").is_dir():
    remove_tree('build', verbose=True)


ROOT = pathlib.Path(__file__).parent.resolve()
FORTRAN_DIR = ROOT / "F95"
make = shutil.which('make')
ext_name = "ossssimlib"


class BuildExtWithMake(build_ext):
    def run(self):
        # Drop stale objects/mods and f2py/meson tree before every extension build.
        # Stale debug.mod previously left f90wrap wrappers with ambiguous names.
        subprocess.check_call([make, "clean"], cwd=FORTRAN_DIR)

        # 1. Run make to build the Fortran-based extension
        module_name = ext_name.split(".")[-1]   # "simsubs"
        try:
            subprocess.check_call([make, f"MODULE={module_name}"], cwd=FORTRAN_DIR)
        except subprocess.CalledProcessError:
            log = FORTRAN_DIR / "f2py_f90wrap.log"
            if log.is_file():
                lines = log.read_text(errors="replace").splitlines()
                print(f"==== {log} (last 100 lines) ====", file=sys.stderr)
                print("\n".join(lines[-100:]), file=sys.stderr)
            raise

        # 2. Continue normal extension build process
        super().run()


    def build_extension(self, ext):
        # We assume the Makefile produced a shared library for this extension.
        assert ext_name == ext.name # "ossssimlib"
        target_path = pathlib.Path(self.get_ext_fullpath(ext.name)).parent
        built_so = pathlib.Path(self.get_ext_fullpath(ext.name)).name
        built_so = FORTRAN_DIR / f"_{built_so}"
        built_py = FORTRAN_DIR / f"{ext_name}"
        built_so.rename(target_path/built_so.name)
        built_py.rename(target_path / built_py.name)


setup(
    cmdclass={"build_ext": BuildExtWithMake},
)

