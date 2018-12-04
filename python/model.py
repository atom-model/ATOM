#!/usr/bin/env python

from pyatom import Atmosphere, Hydrosphere


class Model ( object ):
    """
    ATOM Model Object
    """
    def __init__( self ):
        self.atm = Atmosphere()
        self.hyd = Hydrosphere()
        self.output_path = ""
        self.time_slice = 0
        self.config_xml = ""


    def load_config_atm ( self, cfg_xml ):
        self.config_xml = cfg_xml
        self.atm.load_config ( self.config_xml )
        self.print_config_atm()



    def load_config_hyd ( self, cfg_xml ):
        self.config_xml = cfg_xml
        self.hyd.load_config ( self.config_xml )
        self.print_config_hyd()



    def print_config_atm ( self ):
        print ( "\n   atmosphere configuration file name = %s\n" % ( self.config_xml ) )
        print ( "   output path is           %s" % ( self.atm.output_path ) )
        print ( "   bathymetry path is at    %s" % ( self.atm.bathymetry_path ) )
        print ( "\n" )



    def print_config_hyd ( self ):
        print ( "\n   hydrosphere configuration file name = %s\n" % ( self.config_xml ) )
        print ( "   output path is           %s" % ( self.hyd.output_path ) )
        print ( "   bathymetry path is at    %s" % ( self.hyd.bathymetry_path ) )
        print ( "\n" )



    def run_Model_atm ( self, t_s ):
        self.time_slice = t_s
        print ( "\n   run_Model for the Atmosphere code prepared for time-slice    Ma = %s\n" % ( self.time_slice ) )
        self.atm.run_time_slice( self.time_slice )
        print ( "\n    successfully terminated Atmosphere code for time-slice    Ma = %s\n" % ( self.time_slice ) )



    def run_Model_hyd ( self, t_s ):
        self.time_slice = t_s
        print ( "\n   run_Model for the Hydrosphere code prepared for time-slice    Ma = %s\n" % ( self.time_slice ) )
        self.hyd.run_time_slice( self.time_slice )
        print ( "\n    successfully terminated Hydrosphere code for time-slice    Ma = %s\n" % ( self.time_slice ) )



model = Model()

model.load_config_atm ( "config_atm.xml" )
model.run_Model_atm ( 5 )

#model.load_config_hyd ( "config_hyd.xml" )
#model.run_Model_hyd ( 5 )
