import subprocess
import numpy as np
import matplotlib.pyplot as plt
import os
import time
import datetime


# Configure figures for production
WIDTH = 495.0  # the number latex spits out
FACTOR = 1.0   # the fraction of the width the figure should occupy
fig_width_pt  = WIDTH * FACTOR

inches_per_pt = 1.0 / 72.27
golden_ratio  = (np.sqrt(5) - 1.0) / 2.0      # because it looks good
fig_width_in  = fig_width_pt * inches_per_pt  # figure width in inches
fig_height_in = fig_width_in * golden_ratio   # figure height in inches
fig_dims      = [fig_width_in, fig_height_in] # fig dims as a list


def subprocess_cmd(command):
    process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
    proc_stdout = process.communicate()[0].strip()
    # print proc_stdout


def generate_folders(widths):
    command = "rm -rf Run*/; "
    subprocess_cmd(command)

    for width in widths:
        run = "Run" + str(width)
        command = "cp -rf base/ " + run + "/; "
        subprocess_cmd(command)

    print ('Folders generated.')


def create_config_file(width, mesh):
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
        hex (5 4 9 10 16 15 20 21) (''' + str(mesh) + ' ' + str(mesh) + ''' 1) simpleGrading (1 1 1)
        hex (0 1 4 5 11 12 15 16) (''' + str(mesh) + ' ' + str(mesh) + ''' 1) simpleGrading (1 1 1)
        hex (1 2 3 4 12 13 14 15) (''' + str(mesh * 2) + ' ' + str(mesh) + ''' 1) simpleGrading (1 1 1)
        hex (4 3 6 7 15 14 17 18) (''' + str(mesh * 2) + ' ' + str(mesh * 2) + ''' 1) simpleGrading (1 1 1)
        hex (9 4 7 8 20 15 18 19) (''' + str(mesh) + ' ' + str(mesh * 2) + ''' 1) simpleGrading (1 1 1)
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


def update_dimensions(widths, meshes):
    for width, mesh in zip(widths, meshes):
        run = "Run" + str(width)
        path = run + '/constant/polyMesh/blockMeshDict'
        with open(path, 'w') as config_file:
            config_file.write(create_config_file(width, mesh))

    print ('Config generated.')


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
        print(run + ' complete.')

    print('Simulations complete.')


def generate_plots(widths):
    # Format our plots
    plt.figure(figsize=fig_dims)
    plt.xlabel('Distance, y (m)')
    plt.ylabel('Stress ($\sigma_{xx}$)$_{x=0}$(kPa)')
    plt.title('Normal stress along the vertical symmetry')

    x = np.linspace(0.5, 2)
    sigmaxx = 1E4*(1+(0.125/(x**2))+(0.09375/(x**4)))

    plt.plot(x, sigmaxx, '-k', label='Analytic Solution')
    plt.xlim(0.5, 2)

    for width in widths:
        path = 'Run' + str(width) + '/postProcessing/sets/100/leftPatch_sigmaxx.xy'
        data = np.loadtxt(path)
        label = 'Explicit Solution (y=' + str(2*width) + ')'
        plt.plot(data[:, 0], data[:, 1], '--', markersize=5, label=label)

    plt.legend(loc='best')

    # Save plots
    ts = time.time()
    st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d-%H-%M-%S')
    save_name = 'test-' + st + '.pdf'
    try:
        os.mkdir('figures')
    except Exception:
        pass

    plt.savefig('figures/' + save_name, bbox_inches='tight')
    plt.clf()

    print('Plots generated.')


def main():
    widths = [2, 3, 4, 5]
    meshes = [10, 10, 10, 10]

    print('Running widths ' + str(widths) + ' with meshes ' + str(meshes) + '.')
    generate_folders(widths)
    update_dimensions(widths, meshes)
    run_simulations(widths)
    generate_plots(widths)
    print('Done!')

if __name__ == "__main__":
    main()
