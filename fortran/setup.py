# setup.py
from glob import glob
import os
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import pathlib
import subprocess
import shutil
import sys


ROOT = pathlib.Path(__file__).parent.resolve()
FORTRAN_DIR = ROOT / "src/F95"
make = shutil.which('make')
LOG_FILE = open(f"{ROOT.joinpath("log.txt")}", "w")

ext_name = "ossssimlib"


class BuildExtWithMake(build_ext):
    def run(self):
        # 1. Run make to build the Fortran-based extension
        module_name = ext_name.split(".")[-1]   # "simsubs"
        subprocess.check_call([make, f"MODULE={module_name}"], cwd=FORTRAN_DIR)
        
        # 2. Continue normal extension build process
        super().run()


    def build_extension(self, ext):
        # We assume the Makefile produced a shared library for this extension.
        # Name it so that we can locate it deterministically.
        # Example: for ext.name == "surveysim._simsub", Makefile builds:
        #   fortran/_simsub.so
        ext_name = ext.name.split(".")[-1]   # "_simsubs"
        target_path = pathlib.Path(self.get_ext_fullpath(ext.name)).parent
        built_so = pathlib.Path(self.get_ext_fullpath(ext.name)).name
        built_so = FORTRAN_DIR / f"_{built_so}"
        built_py = FORTRAN_DIR / f"{ext_name}"
        LOG_FILE.write(f"{built_so} -> {target_path}\n")
        built_so.rename(target_path/built_so.name)
        LOG_FILE.write(f"{built_py} -> {target_path}\n")
        built_py.rename(target_path / built_py.name)


ext_modules = [
    # No sources: our build_ext copies in the prebuilt .so
    Extension(ext_name, sources=[]),
]

setup(
    ext_modules=ext_modules,
    cmdclass={"build_ext": BuildExtWithMake},
    package_dir={"": "src"},
    packages=["ossssim"],
)
