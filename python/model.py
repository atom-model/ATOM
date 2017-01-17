from pyatom import Atmosphere, Hydrosphere


class Model(object):
    """
    ATOM Model Object
    """
    def __init__(self):
        self.atm = Atmosphere()
        self.hyd = Hydrosphere()

    def run(self, t, output_dir=None):
        print 'TODO RUN %s %s' % (t, output_dir)
