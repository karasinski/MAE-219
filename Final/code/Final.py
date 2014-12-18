import subprocess
import os


def inplace_change(filename, old_string, new_string):
    with open(filename) as f:
        s = f.read()

    if old_string in s:
        # print 'Changing "{old_string}" to "{new_string}"'.format(**locals())
        s = s.replace(old_string, new_string)
        with open(filename, 'w') as f:
            f.write(s)
    # else:
        # print 'No occurances of "{old_string}" found.'.format(**locals())


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


def create_mesh_file(path, AR, Re):
    M = 100.

    YMESH       = str(int(M * AR))
    INV_GRADING = str(0.1)
    GRADING     = str(10.)
    MESH        = str(int(M))
    AR          = str(-AR)

    inplace_change(path, 'AR',                   AR)
    inplace_change(path, 'XMESH',              MESH)
    inplace_change(path, 'YMESH',             YMESH)
    inplace_change(path, 'INV_GRADING', INV_GRADING)
    inplace_change(path, 'GRADING',         GRADING)
    inplace_change(path, 'MESH',               MESH)


def create_properties_files(path, AR, Re):
    d = 1.           # characteristic size of domain
    NU = str(d / Re)

    inplace_change(path, 'NU_VAR', NU)


def update_dimensions(ARs, Res):
    for AR, Re in zip(ARs, Res):
        run = "Run" + str(AR) + '-' + str(Re)
        path = run + '/constant/polyMesh/blockMeshDict'
        create_mesh_file(path, AR, Re)

        path = run + '/constant/transportProperties'
        create_properties_files(path, AR, Re)

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
            command += "decomposePar; "
            command += "mpirun -np 4 icoFoam -parallel > log; "
            command += "reconstructPar; "
            command += "foamCalc components U; "
            command += "sample; "

            subprocess_cmd(command)
        print(run + ' complete.')

    print('Simulations complete.')


def main(ARs, Res):
    print('Running ARs ' + str(ARs) + ' with Res ' + str(Res) + '.')
    generate_folders(ARs, Res)
    update_dimensions(ARs, Res)
    run_simulations(ARs, Res)
    print('Done!')

if __name__ == "__main__":
    # Base case
    ARs = [  0.5]
    Res = [100.0]
    #          o                      Broken=x Working=o
    main(ARs, Res)

    # Additional cases
    ARs = [0.5,    0.5,   2.0,   5.0]
    Res = [1.0, 2000.0, 100.0, 100.0]
    #        x       o      o      o  Broken=x Working=o
    main(ARs, Res)
