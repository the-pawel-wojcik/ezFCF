#!/usr/bin/env python3

"""
This script creates an input XML file for the 'ezSpectrum 3.0' programm

Script supports the following ab initio outputs:

Q-Chem 
ACES II (versions 0.3 [old] and 2.5.0 [new])
Molpro
GAMESS (linear molecules are not supported)
Gaussian (we can not guarantee it works, since we are not allowed to touch it)

Script requires a set of ab initio outputs with frequency jobs in
the optimized geometries. First file in the argument's list is
the initial state, others are the target states.
From each file this script copies a geometry, frequencies and
normal modes into the XML file.
Default job parameters and energy difference between the electronic
states (IP, EA, etc) should be adjusted as necessary in the created XML file.
"""

import sys
import numpy as np

def parse_qchem(StateF, data: dict):
    Line = StateF.readline()
    while Line:
        if (Line.find('Standard Nuclear Orientation') >= 0) and (data['ifGeometryIsLoaded'] is False):
            StateF.readline()
            StateF.readline()
            Line = StateF.readline()
            while Line.find('----') == -1:
                data['NAtoms'] += 1
                data['Geometry'] += Line[5:]
                data['atoms_list'] += Line[8:14]
                Line = StateF.readline()
            data['ifGeometryIsLoaded'] = True

        if Line.find('Raman Active: ') >= 0:
            StateF.readline()
            for i in range(data['NAtoms']):
                Line = StateF.readline()
                data['NormalModes'] += Line[2:]
            data['NormalModes'] += '\n'

        # TODO: Pawel: I do not understand the comment below. It comes from the previous version of the script
        # remove end of the line symbols!!!
        if Line.find('Frequency: ') >= 0:
            data['Frequencies'] += Line.replace('Frequency:','')
            data['NLinesWithFrequencies'] += 1

        Line = StateF.readline()

    if data['NLinesWithFrequencies'] == data['NAtoms']-1:
        data['ifLinear'] = True

    data['if_normal_modes_weighted'] = True
    data['geometry_units'] = "angstr"

    # END Q-Chem
    # ================================================================================

def parse_aces_old(StateF, data: dict):
    Line = StateF.readline()
    while Line:
        # TODO: Pawel: I do not understand the comment below. It comes from the previous version of the script
        # in FREQ job only one such line. Check it if it will be from OPT job....
        if (Line.find('Coordinates (in bohr)') >= 0) and (data['ifGeometryIsLoaded'] is False):
            StateF.readline()
            StateF.readline()
            Line = StateF.readline()
            while Line.find('----') == -1:
                data['NAtoms'] += 1
                data['Geometry'] += "      " + Line[5:9] + Line[20:]
                data['atoms_list'] += Line[5:9]
                Line = StateF.readline()
                data['ifGeometryIsLoaded'] = True
    
        if Line.find('Normal Coordinates') >= 0:
            Line = StateF.readline()
            while (Line == "\n"):
                StateF.readline()
                Line = StateF.readline()
                data['NLinesWithFrequencies'] += 1
                data['Frequencies'] += Line
                StateF.readline()
                StateF.readline()
                for i in range(data['NAtoms']):
                    Line = StateF.readline()
                    data['NormalModes'] += Line[4:]
                data['NormalModes'] += '\n'
                Line = StateF.readline()
        Line = StateF.readline()
    if data['NLinesWithFrequencies'] == data['NAtoms']-1:
        data['ifLinear'] = True
    
    data['if_normal_modes_weighted'] = False
    data['geometry_units'] = "au"
    
    # END ACES_OLD
    # ================================================================================


def parse_aces_new(StateF, data: dict):
    Line = StateF.readline()
    while Line:
        # TODO: Pawel: I do not understand the comment below. It comes from the previous version of the script
        # in FREQ job only one such line. Check it if it will be from OPT job....
        if (Line.find('C o o r d i n a t e s') >= 0) and (data['ifGeometryIsLoaded'] is False):
            StateF.readline()
            StateF.readline()
            Line = StateF.readline()
            while Line.find('----') == -1:
                data['NAtoms'] += 1
                data['Geometry'] += "      " + Line[5:9] + Line[20:]
                data['atoms_list'] += Line[5:9]
                Line = StateF.readline()
            data['ifGeometryIsLoaded'] = True

        if Line.find('Normal Coordinates') >= 0:
            Line = StateF.readline()
            Line = StateF.readline()
            while (Line == "\n"):
                StateF.readline()
                Line = StateF.readline()
                data['NLinesWithFrequencies'] += 1
                data['Frequencies'] += Line
                StateF.readline()
                StateF.readline()
                for i in range(data['NAtoms']):
                    Line = StateF.readline()
                    data['NormalModes'] += Line[4:]
                data['NormalModes'] += '\n'
                Line = StateF.readline()
        Line = StateF.readline()
    if data['NLinesWithFrequencies'] == data['NAtoms']-1:
        data['ifLinear'] = True

    data['if_normal_modes_weighted'] = False
    data['geometry_units'] = "au"

# END ACES_NEW
# ================================================================================


def parse_molpro(StateF, data: dict):
    normal_coordinates = []
    n_normal_modes = 0
    
    Line = StateF.readline()
    while Line:
        # geometry
        if (Line.find('ATOMIC COORDINATES') >= 0) and (data['ifGeometryIsLoaded'] is False):
            StateF.readline()
            StateF.readline()
            StateF.readline()
            Line = StateF.readline()
            while Line.find('Bond lengths in Bohr (Angstrom)') == -1:
                data['NAtoms'] += 1
                data['Geometry'] += "      " + Line[5:8] + Line[19:]
                data['atoms_list'] += Line[5:8]
                Line = StateF.readline()
            data['ifGeometryIsLoaded'] = True
            data['NAtoms'] -= 1
            # create an empty list of coordinates for every atom in form "X_nm1 Y_nm1 Z_nm1 X_nm2 Y_nm2 Z_nm2 X_nm3..."
            for i in range(data['NAtoms']):
                normal_coordinates.append([])
     
        if (Line.find('Intensities [relative]') >= 0) and (data['ifNormalModesLoaded'] is False):
            for i in range(data['NAtoms']):
                LineX = StateF.readline()
                LineY = StateF.readline()
                LineZ = StateF.readline()
                
                # analyse the length of the line 20+12*j:
                nm_per_line = (len(LineX)-23)//12
                
                n_normal_modes += nm_per_line
                
                for j in range(nm_per_line):
                    normal_coordinates[i].append( LineX[ 23+12*j : 23+12*(j+1) ]+LineY[ 23+12*j : 23+12*(j+1) ]+LineZ[ 23+12*j : 23+12*(j+1)] )
    
            n_normal_modes //= data['NAtoms'] # repeated "NAtoms"-times
            
        if (Line.find('Wavenumbers [cm-1]') >= 0) and (data['ifFrequenciesLoaded'] is False):
            data['Frequencies'] += Line.replace('Wavenumbers [cm-1]','')
            data['NLinesWithFrequencies'] += 1
    
        if Line.find('Normal Modes of low/zero frequencies') >= 0:
            data['ifFrequenciesLoaded'] = True
            data['ifNormalModesLoaded'] = True
    
        Line = StateF.readline()
    
    # now create normal modes in q-chem format:
    printed_normal_modes = 0
    for j in range( len(normal_coordinates[0])//3 ):
        for i in range(data['NAtoms']):
            data['NormalModes'] += normal_coordinates[i][j*3] + "   " + normal_coordinates[i][j*3+1] + "   " + normal_coordinates[i][j*3+2] + '\n'
        printed_normal_modes += 3
        data['NormalModes'] += '\n'
    
    # add the reminder of 3
    if ( (len(normal_coordinates[0]) % 3) > 0 ):
        for i in range(data['NAtoms']):
            for j in range( len(normal_coordinates[0]) % 3 ):
                data['NormalModes'] += normal_coordinates[i][printed_normal_modes+j]+"   "
            data['NormalModes'] += '\n'
        data['NormalModes'] += '\n'
    
    if (n_normal_modes == (3*data['NAtoms']-5)):
        data['ifLinear'] = True
    
    data['if_normal_modes_weighted'] = True
    data['geometry_units'] = "au"
    
    # END MOLPRO
    # ================================================================================


def parse_gamess(StateF, data: dict, run_type):
    if run_type != "web":
        print("\nWarning! All geometries from GAMESS' outputs treated as non-linear.")
    
    data['ifLinear'] = False
    
    normal_coordinates = []
    n_normal_modes = 0
    
    Line = StateF.readline()
    while Line:
        # geometry
        if (Line.find('COORDINATES (BOHR)') >= 0) and (data['ifGeometryIsLoaded'] is False):
            StateF.readline()
            Line = StateF.readline()
            while Line.find('INTERNUCLEAR DISTANCES (ANGS.)') == -1:
                data['NAtoms'] += 1
                data['Geometry'] += "     " + Line[0:5] + Line[19:]
                data['atoms_list'] += Line[0:5]
                Line = StateF.readline()
            data['ifGeometryIsLoaded'] = True
            data['NAtoms'] -= 1
            # create an empty list of coordinates for every atom in form "X_nm1 Y_nm1 Z_nm1 X_nm2 Y_nm2 Z_nm2 X_nm3..."
            for i in range(data['NAtoms']):
                normal_coordinates.append([])
     
        if (Line.find(' IR INTENSITY:') >= 0) and (data['ifNormalModesLoaded'] is False):
            Line = StateF.readline()
           
            for i in range(data['NAtoms']):
                LineX = StateF.readline()
                LineY = StateF.readline()
                LineZ = StateF.readline()
                
                # analyse the length of the line 20+12*j:
                nm_per_line = (len(LineX)-20)//12
                
                for j in range(nm_per_line):
                    normal_coordinates[i].append( LineX[ 20+12*j : 20+12*(j+1) ]+LineY[ 20+12*j : 20+12*(j+1) ]+LineZ[ 20+12*j : 20+12*(j+1)] )
    
            n_normal_modes //= data['NAtoms'] # repeated "NAtoms"-times
    
        if (Line.find('   FREQUENCY:') >= 0) and (data['ifFrequenciesLoaded'] is False):
            data['Frequencies'] += Line.replace('   FREQUENCY:','')
            data['NLinesWithFrequencies'] += 1
    
        if Line.find('REFERENCE ON SAYVETZ CONDITIONS') >= 0:
            data['ifFrequenciesLoaded'] = True
            data['ifNormalModesLoaded'] = True
    
        Line = StateF.readline()
    
    # now create normal modes in q-chem format:
    printed_normal_modes = 0
    for j in range( len(normal_coordinates[0])//3 ):
        for i in range(data['NAtoms']):
            data['NormalModes'] += normal_coordinates[i][j*3]+"   "+normal_coordinates[i][j*3+1]+"   "+normal_coordinates[i][j*3+2]+'\n'
        printed_normal_modes += 3
        data['NormalModes'] += '\n'
    
    # add the reminder of 3
    if ( (len(normal_coordinates[0]) % 3) > 0 ):
        for i in range(data['NAtoms']):
            for j in range( len(normal_coordinates[0]) % 3 ):
                data['NormalModes'] += normal_coordinates[i][printed_normal_modes+j]+"   "
            data['NormalModes'] += '\n'
        data['NormalModes'] += '\n'
    
    # remove first 6 normal modes (two lines of three normal modes)
    # i.e. rotation+translation (works only for non-linear molecules!):
    for i in range((data['NAtoms']+1)*2):
        data['NormalModes'] = data['NormalModes'][data['NormalModes'].find('\n')+1:]
    # remove first 6 frequencies:
    Frequencies_set = data['Frequencies'].split()
    Frequencies_set = Frequencies_set[6:]
    data['Frequencies'] = " ".join([f"{f:s}" for f in Frequencies_set])
    data['Frequencies'] = "       " + data['Frequencies'] + "\n"
    
    data['if_normal_modes_weighted'] = True
    data['geometry_units'] = "au"
    
    # END GAMESS
    # ================================================================================


def parse_other(StateF, data: dict):
    ifAtomsLoaded = False
    
    Line = StateF.readline()
    while Line:
        if (Line.find('Distance matrix') >= 0) and (ifAtomsLoaded is False):
            StateF.readline()
            while Line.find('Stoichiometry') == -1:
                data['atoms_list'] += Line[8:12]
                Line = StateF.readline()
            atoms_set=data['atoms_list'].split()
            ifAtomsLoaded = True
    
        if (Line.find('Standard orientation') >= 0) and (data['ifGeometryIsLoaded'] is False):
            for L in range(4):
                StateF.readline()
            Line = StateF.readline()
            while Line.find('----') == -1:
                data['Geometry'] += "      " + atoms_set[data['NAtoms']] + Line[32:]
                Line = StateF.readline()
                data['NAtoms'] += 1
            data['ifGeometryIsLoaded'] = True
    
        if Line.find('Atom AN') >= 0:
            for i in range(data['NAtoms']):
                Line = StateF.readline()
                data['NormalModes'] += Line[10:]
            data['NormalModes'] += '\n'
    
        # TODO: I do not understand the comment below. It comes from the previous version of the script. Pawel
        # remove end of the line symbols!!!
        if Line.find('Frequencies --') >= 0:
            data['Frequencies'] += Line.replace('Frequencies --', '')
            data['NLinesWithFrequencies'] += 1
    
        Line = StateF.readline()
    
    if data['NLinesWithFrequencies'] == data['NAtoms']-1:
        data['ifLinear'] = True
    
    data['if_normal_modes_weighted'] = True
    data['geometry_units'] = "angstr"
    
    # END OTHER
    # ================================================================================


def parse_orca(StateF, data: dict):
    data['if_normal_modes_weighted'] = True                 # this is what the ORCA output is saying
    data['geometry_units'] = "angstr"                       # the same with previous

    Lines = [ line.rstrip() for line in StateF.readlines()] # load file with newline symbols ('\n') stripped
    IndNormModes = Lines.index('$normal_modes')             # index of the line where normal modes start
    IndVibFreq = Lines.index('$vibrational_frequencies')    # index of the line where vibrational frequencies start
    IndXYZ = Lines.index('$atoms')                          # index of the XYZ coordinates

    NAt = int(Lines[IndXYZ+1].split()[0])                   # number of atoms in the molecule (next line after $atoms)
    for a in Lines[(IndXYZ+2):(IndXYZ+2+NAt)]:              # parse the XYZ block
        words = a.split()
        # Geometry is the XYZ geometry in the format <Atom Label>  <X> <Y> <Z>, atomic masses (2nd column) are ignored
        data['Geometry'] += f"{words[0]:6s}  "
        for word in words[2:5]:
            coordinate = 0.529177210903 * float(word)
            data['Geometry'] += f" {coordinate:12.6f}"
        data['Geometry'] += "\n" 
        # atom list is the list of atomic labels
        data['atoms_list'] += "   " + words[0] + " "
    
    NInRow = len(Lines[IndNormModes+2].split())
    NVib = 3*NAt
    NormModes = np.zeros( (NVib,NAt,3), dtype=float )
    data['NAtoms'] = NAt
    
    BlockNum = 0
    while BlockNum*NInRow < 3*NAt:
        NModes = np.array(Lines[IndNormModes + 2 + BlockNum*(NVib + 1)].split(), dtype=int)
        for i in range(0, NVib):
            words = Lines[IndNormModes + 2 + BlockNum*(NVib + 1) + 1 + i].split()
            tmp = np.array(words[1:], dtype=float)
            for imode, nmode in enumerate(NModes):
                NormModes[nmode][i / 3][i % 3] = tmp[imode]
        BlockNum += 1

    for nmodeblock in range(6, NVib, 3):
        for nat in range(0, NAt):
            for nmib in range(0, 3):
                data['NormalModes'] += "  "
                for mode in NormModes[nmodeblock+nmib][nat]:
                    data['NormalModes'] += f" {mode:7.3f}"
            data['NormalModes'] += "\n"
        data['NormalModes'] += "\n"

    count = 0
    for i in range(6, NVib):
        freq = float(Lines[IndVibFreq + 2 + i].split()[1])
        data['Frequencies'] += "    " + f"{freq:7.2f}"
        count += 1
        if count%3 == 0 and count>=3:
            data['Frequencies'] += "\n"

    # END ORCA
    # ================================================================================


def ReadAndWriteState(FileName, which_state, run_type, xmlF):
    """Reads one state from an output file. Writes to the XML file"""

    # check which ab initio program created the input file
    input_type = ""
    with open(FileName, 'r') as StateF:
        Line = StateF.readline()
        while Line and (input_type == ""):
            if (Line.find('Welcome to Q-Chem') >= 0) and (input_type == ""):
                input_type="q-chem"
            if (Line.find('* ACES2:') >= 0) and (input_type == ""):
                input_type="aces_old"
            if (Line.find('* ACES :') >= 0) and (input_type == ""):
                input_type="aces_new"
            if (Line.find('PROGRAM SYSTEM MOLPRO') >= 0) and (input_type == ""):
                input_type="molpro"
            if (Line.find('GAMESS VERSION = ') >= 0) and (input_type == ""):
                input_type="gamess"
            if (Line.find('$orca_hessian_file') >= 0) and (input_type == ""):
                input_type="orca"
            if (Line.find('RESTRICTED RIGHTS') >= 0) and (input_type == ""):
                input_type="other"

            Line = StateF.readline()

    if input_type == "":
        error_message = f'Error: File "{FileName}" has an unknown format.'
        error_message += ' Q-Chem, Molpro, ACES II, GAMESS, or Gaussian frequency jobs are supported'
        if run_type != "web":
            print(error_message)
            sys.exit(2)
        else:
            return error_message

    # a container for data
    data = {'NAtoms': 0,
            'NLinesWithFrequencies': 0,
            'ifLinear': False,
            'Geometry': '',
            'atoms_list': '',
            'NormalModes': '',
            'Frequencies': '',
            'if_normal_modes_weighted': None,
            'geometry_units': None,
            'ifGeometryIsLoaded': False,
            'ifFrequenciesLoaded': False,
            'ifNormalModesLoaded': False} 

    with open(FileName, 'r') as StateF:
        if (input_type == "q-chem"):
            parse_qchem(StateF, data)
               
        if (input_type == "aces_old"):
            parse_aces_old(StateF, data)
           
        if (input_type == "aces_new"):
            parse_aces_new(StateF, data)

        if (input_type == "molpro"):
            parse_molpro(StateF, data)

        if (input_type == "gamess"):
            parse_gamess(StateF, data, run_type)

        if (input_type == "other"):
            parse_other(StateF, data)

        if (input_type == "orca"):
            parse_orca(StateF, data)

    # Write the state to the xml file
    
    xmlF.write('  <geometry\n')
    xmlF.write(f'    number_of_atoms = "{str(data["NAtoms"])}"\n')
    # linear?
    if data['ifLinear'] is True:
        xmlF.write('    linear = "true"\n')
    else:
        xmlF.write('    linear = "false"\n')

    xmlF.write(f'    units   = "{data["geometry_units"]}"\n')
        
    xmlF.write('    text   = "\n')
    xmlF.write(data['Geometry'])
    xmlF.write('             ">\n')
    xmlF.write('  </geometry>\n\n')

    atoms_order = " ".join(["{str(nm)}" for nm in range(data['NAtoms'])])

    xmlF.write('  <OPT_manual_atoms_reordering\n')
    xmlF.write(f'     new_order="{atoms_order}">\n')
    xmlF.write('  </OPT_manual_atoms_reordering>\n\n')

    xmlF.write('  <normal_modes\n')
    xmlF.write(f'    if_mass_weighted="{data["if_normal_modes_weighted"]}"\n')
    xmlF.write('    text = "\n')
    xmlF.write(data['NormalModes'])
    xmlF.write('           "\n')
    xmlF.write('   atoms = "')
    xmlF.write(data['atoms_list'])
    xmlF.write('           ">\n')
    xmlF.write('  </normal_modes>\n\n')

    if data['ifLinear'] is True:
        normal_modes = " ".join([f"{str(nm)}" for nm in range(3*data['NAtoms']-5)])
    else:
        normal_modes = " ".join([f"{str(nm)}" for nm in range(3*data['NAtoms']-6)])
        
    if which_state == "target":
        if run_type != "web":
            xmlF.write('  <OPT_manual_normal_modes_reordering\n')
            xmlF.write(f'     new_order="{normal_modes}">\n')
            xmlF.write('  </OPT_manual_normal_modes_reordering>\n\n')
        else: # web version:
            xmlF.write('  <manual_normal_modes_reordering\n')
            xmlF.write(f'     new_order="{normal_modes}">\n')
            xmlF.write('  </manual_normal_modes_reordering>\n\n')

    xmlF.write('  <frequencies\n')
    xmlF.write('    text = "\n')
    xmlF.write(data['Frequencies'])
    xmlF.write('             ">\n')
    xmlF.write('  </frequencies>\n\n')

    return ""


def main(xml_filename, ai_filenames, run_type):
    if len(ai_filenames) == 0:
        if run_type != "web":
            print('\nTo create an XML file from ab initio outputs type:\n'
                  '"make_xml.py  <filename.xml> <initial_state.out> <target_state_1.out> <target_state_2.out> etc..."\n\n')
            sys.exit(2)
        else:
            return "Error. make_xml.py: no ab-initio file names found"

    xmlF = open(xml_filename, 'w')

    # Write default job parameters
    xmlF.write("""<input
  job = "harmonic_pes">

<job_parameters 
        temperature                              = "300"
        spectrum_intensity_threshold             = "0.001" >
</job_parameters>

<!-- 
  ______________________________________________________________________ 

    Tags which start with "OPT_" will be ignored.
    To include these optional keywords please "uncomment" by removing 
    "OPT_" from the start and the corresponding end tag (if present)
  ______________________________________________________________________ 

 -->

<parallel_approximation
        max_vibr_excitations_in_initial_el_state = "1"
        max_vibr_excitations_in_target_el_state  = "4"
        combination_bands                        = "true"
        use_normal_coordinates_of_target_states  = "true"
 >

  <OPT_do_not_excite_subspace size = "0" normal_modes = " " >
  </OPT_do_not_excite_subspace>

  <OPT_energy_thresholds  units="eV, K, cm-1">
    <initial_state   units="K">      1000    </initial_state>
    <target_state   units="eV">      0.25    </target_state>
  </OPT_energy_thresholds> 

  <OPT_print_franck_condon_matrices flag="true">
  </OPT_print_franck_condon_matrices>

</parallel_approximation>

<!-- 
  ______________________________________________________________________ 

 -->

<OPT_dushinsky_rotations target_state="1"
      max_vibr_excitations_in_initial_el_state = "1"
      max_vibr_excitations_in_target_el_state  = "4"
      > 
  <OPT_max_vibr_to_store  target_el_state="4">
  </OPT_max_vibr_to_store>

  <OPT_do_not_excite_subspace size = "2" normal_modes = "0 1">
  </OPT_do_not_excite_subspace>

  <OPT_energy_thresholds  units="eV, K, cm-1">
    <initial_state   units="K">      1000    </initial_state>
    <target_state   units="eV">      0.25    </target_state>
  </OPT_energy_thresholds> 

  <OPT_single_excitation 
       ini="0" 
       targ="1v1">
  </OPT_single_excitation>

</OPT_dushinsky_rotations>

<!-- 
  ______________________________________________________________________ 

 -->\n\n""")
  
    # web version
    if run_type == "web": 
        xmlF.write('<if_web_version flag="true">\n</if_web_version>\n\n')
 
    xmlF.write('<initial_state>\n')
    xmlF.write(f'  <!-- THIS INITIAL STATE IS FROM "{ai_filenames[0]}" FILE -->\n\n')

    return_string = ReadAndWriteState(ai_filenames[0], "initial", run_type, xmlF)
    # TODO: This is quite unique error handling mechanism, in later versions substitute with python error mechanism 
    if not (return_string == ""):
        return return_string
    
    xmlF.write('</initial_state>\n\n')
    xmlF.write("""<!-- 
  ______________________________________________________________________ 

 -->\n\n""")

    state_n = 0
    for targetStateFileName in ai_filenames[1:]:
        state_n += 1
        xmlF.write('<target_state>\n\n')
        xmlF.write('  <ip units="eV"> {str(state_n)} </ip>\n\n')
        xmlF.write('  <!-- THIS TARGET STATE IS FROM "{targetStateFileName}" FILE -->\n')
        return_string = ReadAndWriteState(targetStateFileName, "target", run_type, xmlF)
        if not(return_string == ""):
            return return_string
        xmlF.write('</target_state>\n\n')
        xmlF.write("""<!-- 
  ______________________________________________________________________ 

 -->\n\n""")
    xmlF.write('</input>\n')

    return ""


if __name__ == "__main__":
    if len(sys.argv) < 2:
        error_message = '\nTo create an XML file from ab initio outputs type:\n\
        "make_xml.py  <filename.xml> <initial_state.out> <target_state_1.out> <target_state_2.out> etc..."\n\n'
        print(error_message)
        sys.exit(2)
    main(sys.argv[1], sys.argv[2:], "command_line")
