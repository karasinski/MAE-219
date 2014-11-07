import subprocess


def subprocess_cmd(command):
    process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
    proc_stdout = process.communicate()[0].strip()
    print proc_stdout


def generate_folders(widths):
    command = "rm -rf Run*/; "
    subprocess_cmd(command)

    for width in widths:
        run = "Run" + str(width)
        command = "cp -rf base/ " + run + "/; "
        subprocess_cmd(command)


def create_config_file(width):
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

    convertToMeters 1;


    vertices
    (
        (0.5 0 0)
        (1 0 0)
        (''' + str(width) + ''' 0 0)
        (''' + str(width) + ''' 0.707107 0)
        (0.707107 0.707107 0)
        (0.353553 0.353553 0)
        (''' + str(width) + ''' 2 0)
        (0.707107 2 0)
        (0 2 0)
        (0 1 0)
        (0 0.5 0)
        (0.5 0 0.5)
        (1 0 0.5)
        (''' + str(width) + ''' 0 0.5)
        (''' + str(width) + ''' 0.707107 0.5)
        (0.707107 0.707107 0.5)
        (0.353553 0.353553 0.5)
        (''' + str(width) + ''' 2 0.5)
        (0.707107 2 0.5)
        (0 2 0.5)
        (0 1 0.5)
        (0 0.5 0.5)
    );

    blocks
    (
        hex (5 4 9 10 16 15 20 21) (''' + str(width * 5) + ' ' + str(width * 5) + ''' 1) simpleGrading (1 1 1)
        hex (0 1 4 5 11 12 15 16) (''' + str(width * 5) + ' ' + str(width * 5) + ''' 1) simpleGrading (1 1 1)
        hex (1 2 3 4 12 13 14 15) (''' + str(width * 10) + ' ' + str(width * 5) + ''' 1) simpleGrading (1 1 1)
        hex (4 3 6 7 15 14 17 18) (''' + str(width * 10) + ' ' + str(width * 10) + ''' 1) simpleGrading (1 1 1)
        hex (9 4 7 8 20 15 18 19) (''' + str(width * 5) + ' ' + str(width * 10) + ''' 1) simpleGrading (1 1 1)
    );

    edges
    (
        arc 0 5 (0.469846 0.17101 0)
        arc 5 10 (0.17101 0.469846 0)
        arc 1 4 (0.939693 0.34202 0)
        arc 4 9 (0.34202 0.939693 0)
        arc 11 16 (0.469846 0.17101 0.5)
        arc 16 21 (0.17101 0.469846 0.5)
        arc 12 15 (0.939693 0.34202 0.5)
        arc 15 20 (0.34202 0.939693 0.5)
    );

    boundary
    (
        left
        {
            type symmetryPlane;
            faces
            (
                (8 9 20 19)
                (9 10 21 20)
            );
        }
        right
        {
            type patch;
            faces
            (
                (2 3 14 13)
                (3 6 17 14)
            );
        }
        down
        {
            type symmetryPlane;
            faces
            (
                (0 1 12 11)
                (1 2 13 12)
            );
        }
        up
        {
            type patch;
            faces
            (
                (7 8 19 18)
                (6 7 18 17)
            );
        }
        hole
        {
            type patch;
            faces
            (
                (10 5 16 21)
                (5 0 11 16)
            );
        }
        frontAndBack
        {
            type empty;
            faces
            (
                (10 9 4 5)
                (5 4 1 0)
                (1 4 3 2)
                (4 7 6 3)
                (4 9 8 7)
                (21 16 15 20)
                (16 11 12 15)
                (12 13 14 15)
                (15 14 17 18)
                (15 18 19 20)
            );
        }
    );

    mergePatchPairs
    (
    );

    // ************************************************************************* //

    '''

    return config


def update_dimensions(widths):
    for width in widths:
        run = "Run" + str(width)
        path = run + '/constant/polyMesh/blockMeshDict'
        with open(path, 'w') as config_file:
            config_file.write(create_config_file(width))


def run_simulations(widths):
    for width in widths:
        run = "Run" + str(width)
        command = "hdiutil attach -quiet -mountpoint $HOME/OpenFOAM OpenFOAM.sparsebundle; "
        command += "sleep 1; "
        command += "source $HOME/OpenFOAM/OpenFOAM-2.3.0/etc/bashrc; "
        command += "cd " + run + "; "
        command += "blockMesh; "
        command += "solidDisplacementFoam > log; "
        command += "foamCalc components sigma; "
        command += "sample"
        subprocess_cmd(command)


def main():
    widths = [2, 4, 6, 8, 10, 20]

    generate_folders(widths)
    update_dimensions(widths)
    run_simulations(widths)

if __name__ == "__main__":
    main()