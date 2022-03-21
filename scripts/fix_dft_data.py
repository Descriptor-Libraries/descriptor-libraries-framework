from glob import glob
import os
import yaml
from tqdm import tqdm


def main(smiles_filename, data_folder):
    data_files = glob(os.path.join(data_folder, '*.yml'))
    data_files.sort()

    raw_smiles_data = open(smiles_filename, 'r').read().split('\n')[:-1]
    smiles_data = {}
    for line in raw_smiles_data:
        num, smiles = line.split('    ')
        smiles_data[num] = smiles

    for filename in tqdm(data_files):
        if not '_' in os.path.basename(filename):
            continue
        num = os.path.basename(filename).split('_')[0]
        data = yaml.load(open(filename, 'r'), Loader=yaml.CLoader)
        smiles = smiles_data[num]
        data['smiles'] = smiles
        with open(filename, 'w') as f:
            yaml.dump(data, f, Dumper=yaml.CDumper)


if __name__ == '__main__':
    main('/Users/tga/Downloads/dft_data_210125.smi',
         './xtb/results_all_noNi/')
