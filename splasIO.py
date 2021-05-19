import json
from splaspy import wrapper_py as splaspyw
import numpy as np
import os.path
import logging

class io_json:

    """
    This function reads the json input to set all parameter values
    """

    def __init__(self, json_input):
        """ Constructor: reads json string and calls set functions """
        if os.path.isfile(json_input):
            self.init_configure = json.loads(open(json_input).read())
            logging.debug(" INPUT = json file ")
        else :
            self.init_configure = json.loads(json_input)
            logging.debug(" INPUT = json string ")



    def read(self):
        "Set input parameter from json file"
        self._read_name_filepath()
        self._read_cell()
        self._read_scalar()
        self._read_sim_par()

    def _read_name_filepath(self):
        """Sets simulation name and file output directory"""
        self.name = self.init_configure["splas_simulation"]
        outpath = self.init_configure["filepath"]
        self.outdir = outpath["dir_out"]

    def _read_cell(self):
        """ Sets parameters for each cell """
        cellpar = self.init_configure["cell"]
        self.nut_supply = cellpar["input_supply_nut"]
        self.time_flag = cellpar["time-varying"]
        self.nuts = cellpar["nutrients"]
        #nuts = cellpar["nutrients"]
        self.nuts_name = self.nuts["name"]
        self.nuts_value = self.nuts["value"]

    def _read_scalar(self):
        """ Sets scalar parameters of the model's equations """
        scal = self.init_configure["scalar"]
        self.nut_number = scal["number_nutrients"]
        self.phyto_num_class = scal["phyto_class_size"]
        self.zoo_num_class = scal["zoo_class_size"]
        self.phyto_min_size = scal["phyto_x_min"]
        self.phyto_max_size = scal["phyto_x_max"]
        self.zoo_min_size = scal["zoo_x_min"]
        self.zoo_max_size = scal["zoo_x_max"]
        self.zoo_growth = scal["zoo_growth_eff"]
        self.zoo_egestion = scal["zoo_eg_frac"]
        self.phyto_mortality = scal["phyto_mort"]
        self.phyto_sink = scal["phyto_sinking"]
        self.zoo_mortality = scal["zoo_mort"]
        self.zoo_hs = scal["zoo_half_sat"]
        self.zoo_dx = scal["zoo_prey_size_tolerance"]
        self.init_p = scal["initial_p"]
        self.init_z = scal["initial_z"]
        self.init_pquota = scal["initial_p_quota"]
        self.init_diff = scal["initial_diff"]

    def _read_sim_par(self):
        """ Sets simulation parameters (length and steps) """
        sim_par = self.init_configure["simulation_par"]
        self.simulation_length = sim_par["sim_length"]
        self.time_step = sim_par["delta_t"]
        self.print_step = sim_par["print_step"]
        self.dz = sim_par["delta_z"]
        self.m = sim_par["n_layers"]
        self.i = sim_par["izero"]
        self.sd = sim_par["sd"]
        self.delta = sim_par["delta"]
