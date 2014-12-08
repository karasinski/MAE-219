import subprocess
import os
from PrettyPlots import *


def subprocess_cmd(command):
    process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
    proc_stdout = process.communicate()[0].strip()
    # print proc_stdout


def generate_folders(ARs, Res):
    for AR, Re in zip(ARs, Res):
        run = "Run" + str(AR) + '-' + str(Re)
        if not os.path.exists(run):
            command = "cp -rf base/ " + run + "/; "
            subprocess_cmd(command)

    print ('Folders generated.')


def create_mesh_file(AR, Re):
    mesh = 200.

    config = '''
    /*--------------------------------*- C++ -*----------------------------------*\
    | =========                 |                                                 |
    | \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
    |  \\    /   O peration     | Version:  2.3.0                                 |
    |   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
    |    \\/     M anipulation  |                                                 |
    \*---------------------------------------------------------------------------*/
    FoamFile
    {
        version     2.0;
        format      ascii;
        class       dictionary;
        object      blockMeshDict;
    }
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    convertToMeters 0.1;

    vertices
    (
        (0 0 0)
        (1 0 0)
        (1 ''' + str(AR) + ''' 0)
        (0 ''' + str(AR) + ''' 0)
        (0 0 0.1)
        (1 0 0.1)
        (1 ''' + str(AR) + ''' 0.1)
        (0 ''' + str(AR) + ''' 0.1)
    );

    blocks
    (
        hex (0 1 2 3 4 5 6 7) (''' + str(int(mesh)) + ' ' + str(int(mesh * AR)) + ''' 1) simpleGrading (1 1 1)
    );

    edges
    (
    );

    boundary
    (
        movingWall
        {
            type wall;
            faces
            (
                (3 7 6 2)
            );
        }
        fixedWalls
        {
            type wall;
            faces
            (
                (0 4 7 3)
                (2 6 5 1)
                (1 5 4 0)
            );
        }
        frontAndBack
        {
            type empty;
            faces
            (
                (0 3 2 1)
                (4 5 6 7)
            );
        }
    );

    mergePatchPairs
    (
    );

    // ************************************************************************* //
    '''

    return config


def create_properties_file(AR, Re):
    nu = 0.1 / Re

    config = '''
    /*--------------------------------*- C++ -*----------------------------------*\
    | =========                 |                                                 |
    | \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
    |  \\    /   O peration     | Version:  2.3.0                                 |
    |   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
    |    \\/     M anipulation  |                                                 |
    \*---------------------------------------------------------------------------*/
    FoamFile
    {
        version     2.0;
        format      ascii;
        class       dictionary;
        location    "constant";
        object      transportProperties;
    }
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    nu              nu [ 0 2 -1 0 0 0 0 ] ''' + str(nu) + ''';


    // ************************************************************************* //
    '''

    return config


def update_dimensions(ARs, Res):
    for AR, Re in zip(ARs, Res):
        run = "Run" + str(AR) + '-' + str(Re)
        path = run + '/constant/polyMesh/blockMeshDict'
        with open(path, 'w') as config_file:
            config_file.write(create_mesh_file(AR, Re))

        path = run + '/constant/transportProperties'
        with open(path, 'w') as config_file:
            config_file.write(create_properties_file(AR, Re))

    print ('Config generated.')


def run_simulations(ARs, Res):
    for AR, Re in zip(ARs, Res):
        run = "Run" + str(AR) + '-' + str(Re)
        if not os.path.exists(run + '/log'):
            print(run + ' running now.')
            command = "hdiutil attach -quiet -mountpoint $HOME/OpenFOAM OpenFOAM.sparsebundle; "
            command += "sleep 1; "
            command += "source $HOME/OpenFOAM/OpenFOAM-2.3.0/etc/bashrc; "
            command += "cd " + run + "; "
            command += "blockMesh; "
            command += "icoFoam > log; "
            command += "streamFunction"
            subprocess_cmd(command)
        print(run + ' complete.')

    print('Simulations complete.')


def main(ARs, Res):
    print('Running ARs ' + str(ARs) + ' with Res ' + str(Res) + '.')
    generate_folders(ARs, Res)
    update_dimensions(ARs, Res)
    run_simulations(ARs, Res)
    # generate_plots(ARs, Res)
    print('Done!')

if __name__ == "__main__":
    # Base case
    ARs = [  0.5]
    Res = [100.0]
    #          o                      Broken=x Working=o
    main(ARs, Res)

    # Additional cases
    # ARs = [0.5,    0.5,   2.0,   5.0]
    # Res = [1.0, 2000.0, 100.0, 100.0]
    #        x       o      o      o  Broken=x Working=o

    main(ARs, Res)


