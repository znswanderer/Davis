from distutils.core import setup, Extension

module1 = Extension("davis_clib", 
                    sources = ["davis.c"],
                    extra_compile_args = ["-O2", "-std=c99"])

setup(name = "davis", 
      version = "1.0",
      ext_modules = [module1])
