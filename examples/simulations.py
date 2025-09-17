import numpy as np
import os
import h5py
from scipy import interpolate

base_path = os.getcwd()


class fit:
    def __init__(self, halo_id, data_path="/data/EAGLE-50-data/", model="SIDM1b"):
        self.data_path = base_path + data_path
        self.halo_id = halo_id
        self.model = model
        self._model_list = ["SIDM1b", "CDMb", "vdSIDMb"]
        self.hid = int(halo_id) - 1  # halo index in the data files (0-indexed)

        # Check is data path exists
        if not os.path.exists(self.data_path):
            raise Exception(f"Data path {self.data_path} does not exist.")

        # Check if model is valid
        if self.model not in self._model_list:
            raise Exception(
                f"Model {self.model} not recognized. Choose from {self._model_list}."
            )

        sph_filename = (
            self.data_path + f"{self.model}_sphericallyAveraged_density_profiles.hdf5"
        )
        cyl_filename = (
            self.data_path + f"{self.model}_cylindrical_density_and_potential.hdf5"
        )
        axi_filename = self.data_path + f"{self.model}_axisymmetric_shape_profiles.hdf5"

        # Check if files exist
        for filename in [sph_filename, cyl_filename, axi_filename]:
            if not os.path.isfile(filename):
                print(
                    f"File {filename} not found. If you do not need this file, ignore this message..."
                )

        # Load spherical data
        self.sph_data = {}
        try:
            with h5py.File(sph_filename, "r") as f:
                for key in f.keys():
                    self.sph_data[key] = f[key][self.hid]
        except:
            pass

        self.cyl_data = {}
        try:
            with h5py.File(cyl_filename, "r") as f:
                for key in f.keys():
                    self.cyl_data[key] = f[key][self.hid]
        except:
            pass

        self.shape_data = {}
        try:
            with h5py.File(axi_filename, "r") as f:
                for key in f.keys():
                    self.shape_data[key] = f[key][self.hid]
        except:
            pass

        # Load halo properties as attributed of the class
        self.M200 = self.sph_data.get("M200", None)
        self.R200 = self.sph_data.get("R200", None) * 1e3  # convert to kpc

        # Class Methods -- easy access to common profiles
        def sph_avg_dm_density(self):
            r_list = self.sph_data["rs"] * 1e3  # convert to kpc
            rho = self.sph_data["dm_rho"] * 1e-9  # convert to Msol/kpc^3
            return {"r": r_list, "rho": rho}

        def sph_avg_star_density(self):
            r_list = self.sph_data["rs"] * 1e3  # convert to kpc
            rho = self.sph_data["star_rho"] * 1e-9  # convert to Msol/kpc^3
            return {"r": r_list, "rho": rho}

        def sph_avg_gas_density(self):
            r_list = self.sph_data["rs"] * 1e3  # convert to kpc
            rho = self.sph_data["gas_rho"] * 1e-9  # convert to Msol/kpc^3
            return {"r": r_list, "rho": rho}

        def sph_avg_bh_density(self):
            r_list = self.sph_data["rs"] * 1e3  # convert to kpc
            rho = self.sph_data["bh_rho"] * 1e-9  # convert to Msol/kpc^3
            return {"r": r_list, "rho": rho}

        def sph_avg_baryon_density(self):
            r_list = self.sph_data["rs"] * 1e3  # convert to kpc
            rho_star = self.sph_data["star_rho"] * 1e-9  # convert to Msol/kpc^3
            rho_gas = self.sph_data["gas_rho"] * 1e-9  # convert to Msol/kpc^3
            rho_bh = self.sph_data["bh_rho"] * 1e-9  # convert to Msol/kpc^3
            rho = rho_star + rho_gas + rho_bh
            return {"r": r_list, "rho": rho}
