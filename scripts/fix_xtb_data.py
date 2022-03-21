from glob import glob
import os
import yaml
from tqdm import tqdm


def main(data_folder):
    data_files = glob(os.path.join(data_folder, '*.yml'))
    data_files.sort()

    smiles = None
    import ipdb; ipdb.set_trace()
    for filename in tqdm(data_files):
        if not '_' in os.path.basename(filename):
            new_data = {}
            data = yaml.load(open(filename, 'r'), Loader=yaml.CLoader)
            for section_name, section in data.items():
                if isinstance(section, str) or isinstance(section, int) or isinstance(section, float):
                    new_data[section_name] = section
                    continue
                new_data[section_name] = {}
                if not isinstance(section, dict):
                    continue
                for name, value in section.items():
                    if isinstance(value, list):
                        continue
                    new_data[section_name][name] = value
            smiles = data['smiles']
            with open(filename, 'w') as f:
                yaml.dump(new_data, f, Dumper=yaml.CDumper)
            continue
        data = yaml.load(open(filename, 'r'), Loader=yaml.CLoader)
        data['smiles'] = smiles
        with open(filename, 'w') as f:
            yaml.dump(data, f, Dumper=yaml.CDumper)


if __name__ == '__main__':
    main('./xtb/results_all_noNi/')
