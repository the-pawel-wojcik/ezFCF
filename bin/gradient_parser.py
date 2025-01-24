#!/usr/bin/env python

import sys
import argparse


def get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "gradient_file",
        help="Output of a Q-Chem or CFOUR gradient calculation."
    )
    args = parser.parse_args()
    return args


def find_program_name(lines: list[str]) -> str | None:
    for line in lines[:10]:
        if line.strip() == "Welcome to Q-Chem":
            return "qchem"

    for line in lines[:10]:
        if line.strip() == "--invoking executable--":
            return "cfour"

    return None


def get_analytical_gradient_range(lines):
    START_STR = " Full Analytical Gradient: (in a.u., or Hartree/Bohr)"
    END_STR = " Gradient time:  CPU "

    start_index = -1
    end_index = -1

    # Find the section of the output file that contains the gradient
    for n, line in enumerate(lines):
        if start_index == -1 and line.startswith(START_STR):
            start_index = n
        if end_index == -1 and line.startswith(END_STR):
            end_index = n + 1

    if start_index == -1:
        print("Warning! No analytical gradient section was found.",
              file=sys.stderr)

    if end_index == -1:
        print("Warning! No end of the analytical gradient section was found.",
              file=sys.stderr)

    return (start_index, end_index)


def parse_analytical_gradient_QChem(
        qchemoutput: list[str],
) -> dict:
    grad_start_index, grad_end_index = get_analytical_gradient_range(
        qchemoutput
    )

    if grad_start_index == -1 or grad_end_index == -1:
        print("Info: Unable to parse analytical gradient", file=sys.stderr)
        return None

    # Discard the header and time
    grad_input = qchemoutput[grad_start_index+1:grad_end_index-1]

    N_LINES_PER_BLOCK = 4
    n_grad_lines = len(grad_input)
    if n_grad_lines % N_LINES_PER_BLOCK != 0:
        print("Error: The number of gradient lines is NOT divisible by the\n"
              f"\texpected number of lines per block ({N_LINES_PER_BLOCK}).",
              file=sys.stderr)

        print(f"N grad lines: {n_grad_lines}", file=sys.stderr)
        for n, line in enumerate(grad_input):
            print(f"{n}: {line[:-1]}", file=sys.stderr)

        raise RuntimeError()

    gradient = []
    n_grad_blocks = n_grad_lines // N_LINES_PER_BLOCK
    for block_id in range(n_grad_blocks):
        block_header = grad_input[block_id * N_LINES_PER_BLOCK]
        atom_idxs = [int(idx) for idx in block_header.split()]

        block_data = grad_input[block_id * N_LINES_PER_BLOCK + 1:
                                (block_id + 1) * N_LINES_PER_BLOCK]

        gradient_qchem_form = []
        # import re
        # n_atoms_in_block = len(atom_idxs)
        # x_data = re.match(r'\s*1\s*(-?\d+.\d+)'+n_atoms_in_block, block_data)
        # y_data = re.match(r'\s*2\s*(-?\d+.\d+)'+n_atoms_in_block, block_data)
        # z_data = re.match(r'\s*3\s*(-?\d+.\d+)'+n_atoms_in_block, block_data)

        for xyz in range(3):
            data_only = block_data[xyz].split()[1:]
            grad_components = [float(g_xyz) for g_xyz in data_only]
            gradient_qchem_form += [grad_components]

        for idx, atom_idx in enumerate(atom_idxs):
            atom = {
                "atom_idx": atom_idx,
                "x": gradient_qchem_form[0][idx],
                "y": gradient_qchem_form[1][idx],
                "z": gradient_qchem_form[2][idx]
            }

            gradient += [atom]

    return {"unit": "a.u.", "data": gradient}


def skip_to_cfour_gradient(it):
    line = next(it)
    while True:
        if line.strip() == "gradient from JOBARC":
            break
        line = next(it)


def parse_gradient_CFOUR(cfourout: list[str]) -> dict:
    it = iter(cfourout)
    try:
        skip_to_cfour_gradient(it)
    except StopIteration:
        raise RuntimeError("No gradient in the CFOUR output.")

    gradient = list()
    atom_idx = 0
    while True:
        try:
            line = next(it)
        except StopIteration:
            raise RuntimeError("Malformed gradient in the CFOUR output.")

        if line.strip() == '':
            break

        xyz = [float(xyz) for xyz in line.strip().split()]
        gradient += [{
            'atom_idx': atom_idx,
            'x': xyz[0],
            'y': xyz[1],
            'z': xyz[2],
        }]
        atom_idx += 1

    return {"unit": "a.u.", "data": gradient}


def parse_findiff_gradient(qchemoutput):
    grad_start_index, grad_end_index = get_findiff_gradient_range(qchemoutput)

    if grad_start_index == -1 or grad_end_index == -1:
        print("Info: Unable to finite difference gradient", file=sys.stderr)
        return None

    grad_lines = qchemoutput[grad_start_index:grad_end_index]

    gradient = []
    for line in grad_lines:
        sline = line.split()
        atom = {
            "atom_idx": int(sline[0]),
            "x": float(sline[1]),
            "y": float(sline[2]),
            "z": float(sline[3]),
        }

        gradient += [atom]

    return {"unit": "a.u.", "data": gradient}


def get_findiff_gradient_range(lines):
    """ Some methods do not have analytical gradient implement. They can still
    get the gradient by finite difference. The output look like this:
FINAL TENSOR RESULT:
Order 1, Length 78
    Atom      X          Y          Z
     1    3.212317155650272e-05  -1.157126743522474e-06  -2.589693857926722e-02
     2   -1.202336143814783e-04  -3.419040609746764e-06   1.208464090684659e-04
     3    3.079124251153599e-04  -1.214038561104478e-06   2.380809836314850e-02
    ...
    24    1.719858282392718e-04  -2.543103967047002e-06  -2.620156864839180e-02
    25    5.022606270724348e-05  -2.294280206455916e-06   1.727509613030894e-03
    26   -1.412507996837069e-05  -1.066043771197788e-06   2.342512899664083e-02



 ###########################################
 # Finishing finite difference calculation #
 ###########################################
"""
    START_STR = "FINAL TENSOR RESULT:"
    END_STR = " # Finishing finite difference calculation #"

    start_index = -1
    end_index = -1

    # Find the section of the output file that contains the gradient
    for n, line in enumerate(lines):
        if start_index == -1 and line.startswith(START_STR):
            start_index = n
        if end_index == -1 and line.startswith(END_STR):
            end_index = n

    if start_index == -1:
        print("Warning! No finite difference gradient section was found.",
              file=sys.stderr)

    if end_index == -1:
        print("Warning!"
              " No end of the findinte difference gradient section was found.",
              file=sys.stderr)

    return (start_index + 3, end_index - 4)


def pretty_print_gradient(gradient):
    print(f"Parsed gradient ({gradient['unit']}):")
    print(" " * 4 + " " + "    x   " + "    y   " + " " + "    z   ")
    for component in gradient['data']:
        print(f"{component['atom_idx']:3d}:", end='')
        for xyz in ["x", "y", "z"]:
            print(f" {component[xyz]:-8.3f}", end='')
        print("")


GRADIENT_ezFCF = '''
<target_state>
  <vertical_excitation_energy units="eV"> 1.00 </vertical_excitation_energy>
  <gradient
'''
GRADIENT_ezFCF_end = '''
  >
  </gradient>
</target_state>
'''


def print_ezFCF_gradient(gradient, filename: str):
    print(f'<!-- THIS TARGET STATE IS FROM "{filename}" FILE -->')
    print(GRADIENT_ezFCF, end='')
    print('    units = "' + gradient['unit'] + '"')
    print('    text = "', end='')
    for atom in gradient['data']:
        print("\n" + " " * 5, end='')
        for xyz in ['x', 'y', 'z']:
            print(f" {atom[xyz]:-12.5f}", end='')
    print('"', end='')
    print(GRADIENT_ezFCF_end)


def main():
    args = get_args()
    filename = args.gradient_file

    with open(filename, 'r') as f:
        lines = f.readlines()

    program_name = find_program_name(lines)
    parsers = {
        'qchem': parse_analytical_gradient_QChem,
        'cfour': parse_gradient_CFOUR,
    }

    if program_name is None or program_name not in parsers:
        raise RuntimeError(
            "The input file is of an unknown filetype. "
            f"The supported programs are {' '.join(parsers)}."
        )

    print(
        f"Info: the input file format recognized as {program_name}",
        file=sys.stderr,
    )
    parser = parsers[program_name]
    gradient = parser(lines)

    if gradient is None and program_name == "qchem":
        print("Info: Trying finite difference gradient.", file=sys.stderr)
        gradient = parse_findiff_gradient(lines)

    if gradient is None:
        print("Error in parsing gradient.", file=sys.stderr)
        return -1

    # pretty_print_gradient(gradient)
    print_ezFCF_gradient(gradient, filename)

    return 0


if __name__ == "__main__":
    main()
