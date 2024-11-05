mass = {
    'MET': 131.199997,
    'GLY': 57.049999,
    'LYS': 128.199997,
    'THR': 101.099998,
    'ARG': 156.199997,
    'ALA': 71.080002,
    'ASP': 115.099998,
    'GLU': 129.100006,
    'TYR': 163.199997,
    'VAL': 99.07,
    'LEU': 113.199997,
    'GLN': 128.100006,
    'TRP': 186.199997,
    'PHE': 147.199997,
    'SER': 87.080002,
    'HIS': 137.100006,
    'ASN': 114.099998,
    'PRO': 97.120003,
    'CYS': 103.099998,
    'ILE': 113.199997
}

amino_acids = {
    'ALA': 6,
    'ARG': 5,
    'ASN': 17,
    'ASP': 7,
    'CYS': 19,
    'GLU': 8,
    'GLY': 2,
    'HIS': 16,
    'ILE': 20,
    'LEU': 11,
    'LYS': 3,
    'MET': 1,
    'PHE': 14,
    'PRO': 18,
    'GLN': 12,
    'SER': 15,
    'THR': 4,
    'TRP': 13,
    'TYR': 9,
    'VAL': 10
}

single_letter = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}

reverse = {v: k for k, v in amino_acids.items()}
reverse_single = {v: k for k, v in single_letter.items()}