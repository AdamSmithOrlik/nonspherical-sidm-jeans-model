import sys
import os
import numpy as np
import jeans
import datetime as dt
import time as t

sys.path.append("../configs")
from run_dict import run_dictionary as rd


def main(rd, filename=None):
    start = t.time()
    now_str = dt.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")

    # check is save directory exists, if not create it
    if rd["save_profile"]:
        os.makedirs(rd["save_dir"], exist_ok=True)

    ac_input_dict = {}
    if rd["AC_prescription"] is not None:
        ac_input_dict["AC_prescription"] = rd["AC_prescription"]
        if rd["AC_prescription"] == "Gnedin":
            ac_input_dict["Gnedin_params"] = rd["Gnedin_params"]
        else:
            print(
                "Warning: AC_prescription not recognized. No adiabatic contraction will be applied."
            )

    if rd["model"] == "spherical":

        print("Running spherical model...")
        profile = jeans.spherical(
            rd["r1"],
            rd["M200"],
            rd["c"],
            Phi_b=rd["Phi_b"],
            verbose=rd["verbose"],
            **ac_input_dict,
        )
        if profile:
            print(f"Profile computed in {t.time() - start:.2f} seconds.")
        else:
            print("Profile computation failed.")
            return None
    elif rd["model"] == "cdm":
        rd["r1"] = 0.0  # override r1 so that CDM halo is assumed

        print("Running CDM model...")
        profile = jeans.cdm(
            rd["r1"],
            rd["M200"],
            rd["c"],
            q0=rd["q0"],
            Phi_b=rd["Phi_b"],
            verbose=rd["verbose"],
            **ac_input_dict,
        )
        if profile:
            print(f"Profile computed in {t.time() - start:.2f} seconds.")
        else:
            print("Profile computation failed.")
            return None

    elif rd["model"] == "squashed":
        print("Running squashed model...")
        profile = jeans.squashed(
            rd["r1"],
            rd["M200"],
            rd["c"],
            q0=rd["q0"],
            alpha=rd["alpha"],
            Phi_b=rd["Phi_b"],
            q_mode=rd["q_mode"],
            verbose=rd["verbose"],
            **ac_input_dict,
        )
        if profile:
            print(f"Profile computed in {t.time() - start:.2f} seconds.")
        else:
            print("Profile computation failed.")
            return None
    elif rd["model"] == "isothermal":
        print("Running isothermal model...")
        profile = jeans.isothermal(
            rd["r1"],
            rd["M200"],
            rd["c"],
            q0=rd["q0"],
            alpha=rd["alpha"],
            Phi_b=rd["Phi_b"],
            L_list=rd["L_list"],
            M_list=rd["M_list"],
            verbose=rd["verbose"],
            **ac_input_dict,
        )
        if profile:
            print(f"Profile computed in {t.time() - start:.2f} seconds.")
        else:
            print("Profile computation failed.")
            return None
    else:
        print(
            "Error: model not recognized. Choose from 'spherical', 'cdm', 'squashed' or 'isothermal'."
        )
        return None

    if rd["save_profile"]:
        if filename is None:
            filename = f"profile_{rd['model']}_r1_{rd['r1']}_logM200_{np.log10(rd['M200']):.1f}_c_{rd['c']}_{now_str}.npz"
            filepath = os.path.join(rd["save_dir"], filename)
            profile.save(filepath)
            print(f"Profile saved to {filepath}")
        else:
            filepath = os.path.join(rd["save_dir"], filename)
            profile.save(filepath)
            print(f"Profile saved to {filepath}")

    return None


def run_jeans_multi(rd, key="r1", values=[0, 5, 10]):
    for val in values:
        rd[key] = val
        print(f"Running model with {key}={val}")
        main(rd)


if __name__ == "__main__":
    filename = None  # specify for a custom filename

    # main(rd)
    run_jeans_multi(rd, key="r1", values=[0, 5, 10])
