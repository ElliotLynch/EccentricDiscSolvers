#!/usr/bin/env python

from distutils.core import setup

setup(name="eccentric_disc",
      version='1.0',
      description='Eccentric Disc Solver',
      author='Elliot Lynch',
      #packages=['distutils', 'distutils.command'],
      packages=["eccentric_disc", "eccentric_disc.VerticalStructure", "eccentric_disc.ModeSolvers","eccentric_disc.plottools","eccentric_disc.geometry","eccentric_disc.eos","eccentric_disc.datatools"],
      #packages=["eccentric_disc"],
     )



