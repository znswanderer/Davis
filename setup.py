from distutils.core import setup, Extension

# https://stackoverflow.com/questions/58797673/how-to-compile-init-py-file-using-cython-on-windows
from distutils.command.build_ext import build_ext
"""
def get_export_symbols_fixed(self, ext):
      print("ext.name: {}".format(ext.name))
      names = ext.name.split('.')
      if names[-1] != "__init__":
            initfunc_name = "PyInit_" + names[-1]
      else:
            # take name of the package if it is an __init__-file
            initfunc_name = "PyInit_" + names[-2]
      print("initfunc_name: {}".format(initfunc_name))
      if initfunc_name not in ext.export_symbols:
            ext.export_symbols.append(initfunc_name)
      return ext.export_symbols
"""

# This is a BAD quick fix for problems with the windows 
# compilation
def get_export_symbols_fixed(self, ext):
    return [
          "Cells_free", 
          "Cells_new",
          "Cells_clear",
          "dvs_advance",
          "dvs_correct",
          "dvs_calc_force",
          "dvs_populate_cells",
          "dvs_calc_forces",
          "dvs_calc_brute_forces",
          "dvs_visualise_positions",
          "dvs_copy_particles",
          "dvs_collect_forces"
          ] 

# replace wrong version with the fixed:
build_ext.get_export_symbols = get_export_symbols_fixed

module1 = Extension("davis_clib", 
                    sources = ["davis.c"],
                    extra_compile_args = ["-O2"])

setup(name = "davis", 
      version = "1.0",
      ext_modules = [module1])
