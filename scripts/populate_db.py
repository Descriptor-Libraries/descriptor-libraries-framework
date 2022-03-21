import os
import sys
from argparse import ArgumentParser
from glob import glob

import yaml
from tqdm import tqdm

sys.path.append("../backend/app")
from app.db.session import SessionLocal, models

parser = ArgumentParser()
parser.add_argument("--xtb", action="store_true", default=False)
parser.add_argument("--xtb-ni", action="store_true", default=False)
parser.add_argument("--ml", action="store_true", default=False)
parser.add_argument("--dft", action="store_true", default=False)

parser.add_argument("--data-folder", type=str)


def extract_conformer_data(filename):
    data = yaml.load(open(filename, "r"), Loader=yaml.CLoader)
    data_to_store = []
    for key, conformer_data in data.items():
        if key == "smiles":
            continue
        conformer_data_to_store = {}
        # naming is different in xtb file.
        if "coords" in conformer_data.keys():
            conformer_data_to_store["coords"] = conformer_data["coords"]
        else:
            conformer_data_to_store["coords"] = conformer_data["coords_extended"]
        if "elements" in conformer_data.keys():
            conformer_data_to_store["elements"] = conformer_data["elements"]
        else:
            conformer_data_to_store["elements"] = conformer_data["sterimol_parameters"][
                "elements_extended"
            ]
        data_to_store.append(conformer_data_to_store)
    return data_to_store, data["smiles"]


def extract_molecule_data(filename, field):
    data = yaml.load(open(filename, "r"), Loader=yaml.CLoader)
    data_to_store = {field: {}}
    keys = [
        "boltzmann_averaged_data",
        "max_data",
        "min_data",
        "delta_data",
        "vburminconf_data",
    ]

    for key in keys:
        if key in data.keys():
            data_to_store[field][key] = data[key]

    if "smi" in data.keys():
        data_to_store["smiles"] = data["smiv"]
    elif "smiles" in data.keys():
        data_to_store["smiles"] = data["smiles"]
    return data_to_store


def populate_ml_data(path, session):
    data = open(path, "r").read().split("\n")[:-1]
    headers = data[0].split(";")[1:]
    for line in tqdm(data[1:]):
        max_data = {}
        min_data = {}
        delta_data = {}
        boltz_data = {}
        vburminconf_data = {}
        out = line.split(";")
        smiles = out[0]
        for name, value in zip(headers, out[1:]):
            name = "".join(name)
            if name.endswith("_boltz"):
                name = name[: -len("_boltz")]
                boltz_data[name] = value
            elif name.endswith("_delta"):
                name = name[: -len("_delta")]
                delta_data[name] = value
            elif name.endswith("_min"):
                name = name[: -len("_min")]
                min_data[name] = value
            elif name.endswith("_max"):
                name = name[: -len("_max")]
                max_data[name] = value
            elif name.endswith("_vburminconf"):
                name = name[: -len("_vburminconf")]
                vburminconf_data[name] = value
            else:
                raise ValueError("Something is wrong!")
        ml_data = {}
        ml_data["smiles"] = smiles
        ml_data["ml_data"] = {}
        ml_data["ml_data"]["max_data"] = max_data
        ml_data["ml_data"]["min_data"] = min_data
        ml_data["ml_data"]["vburminconf_data"] = vburminconf_data
        ml_data["ml_data"]["delta_data"] = delta_data
        ml_data["ml_data"]["boltzmann_averaged_data"] = boltz_data

        molecule = (
            session.query(models.molecule)
            .filter(models.molecule.smiles == smiles)
            .one_or_none()
        )
        if molecule is None:
            molecule = models.molecule(**ml_data)
            session.add(molecule)
            session.commit()
            continue
        setattr(molecule, "ml_data", ml_data["ml_data"])
        session.commit()


def main(args):
    field = ""
    if args.dft:
        field = "dft_data"
    elif args.xtb:
        field = "xtb_data"
    elif args.ml:
        field = "ml_data"
    elif args.xtb_ni:
        field = "xtb_ni_data"
    else:
        raise ValueError("Datatype not set! Please set one of --dft/--ml/--xtb!")

    session = SessionLocal()
    if field == "ml_data":
        populate_ml_data(args.data_folder, session)
        return

    data_files = glob(os.path.join(args.data_folder, "*.yml"))
    # DFT naming convention is different
    data_files.sort(reverse=True if args.dft else False)

    for i, filename in enumerate(tqdm(data_files)):
        if filename.endswith("_combined.yml"):
            continue

        elif filename.endswith("_confdata.yml") or filename.endswith("_confs.yml"):
            if args.xtb or args.xtb_ni:
                continue
            data, smiles = extract_conformer_data(filename)
            molecule = (
                session.query(models.molecule)
                .filter(models.molecule.smiles == smiles)
                .one_or_none()
            )
            if molecule is None:
                raise ValueError("Could not add conformer: molecule not found")
            for confo in data:
                confo["molecule_id"] = molecule.molecule_id
                session.add(models.conformer(**confo))
            session.commit()

        else:
            data = extract_molecule_data(filename, field)
            molecule = (
                session.query(models.molecule)
                .filter(models.molecule.smiles == data["smiles"])
                .one_or_none()
            )
            if molecule is None:
                molecule = models.molecule(**data)
                session.add(molecule)
                session.commit()
                continue
            setattr(molecule, field, data[field])
            session.add(molecule)
            session.commit()


if __name__ == "__main__":
    args = parser.parse_args()
    if (args.dft and args.ml) or (args.dft and args.xtb) or (args.xtb and args.ml) or (args.xtb_ni and args.ml):
        raise ValueError("only one of --dft/--xtb/--ml can be set at the same time!")
    main(args)
