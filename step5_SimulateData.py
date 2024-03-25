import os
import itertools
from functions.create_config import generate_config
import tqdm
import joblib
import argparse
from functions.utilities import str2bool


def main(tetrises, output_path, GPU_ID, dimensions, phase_plates, phase_shift, defoci, total_doses, tilt_steps,
         start_angles, ice_densities):
    print("Generating configurations from tetrises")
    if not os.path.exists(output_path):
        os.mkdir(output_path)

    # repeat the GPU_IDs so that each simulation is done on a new GPU
    GPU_list = GPU_ID * 1000

    tetris_id = 0
    for combination in itertools.product(total_doses, tilt_steps, start_angles, ice_densities, defoci):
        tetris_file = '{}/{}/atomic_angposfile.txt'.format(tetrises, tetris_id)
        total_dose, tilt_step, start_angle, ice_density, defocus = combination
        eletrons_per_angstrom = total_dose / len(range(start_angle, -start_angle + tilt_step, tilt_step))
        config_dir = '{}/{}'.format(output_path, tetris_id)
        os.mkdir(config_dir)
        config_path = '{}/config.yaml'.format(config_dir)
        generate_config(eletrons_per_angstrom, dimensions, phase_plates, phase_shift, defocus, tetris_file, start_angle, tilt_step,
                        ice_density, config_path, GPU_list[tetris_id])
        tetris_id += 1

    simulation_list = ["{}/{}".format(output_path, i) for i in range(tetris_id)]

    parent_dir = os.path.abspath(os.curdir)
    print(simulation_list)

    # Running the simulation
    def simulate(simulation, parent_dir):
        os.chdir(simulation)
        os.system('parakeet.run -c config.yaml --steps all')
        os.chdir(parent_dir)


    joblib.Parallel(n_jobs=GPU_ID.__len__())(joblib.delayed(simulate)(simulation, parent_dir)
                                             for simulation in tqdm.tqdm(simulation_list))




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run simulations with varying parameters.",
                                     epilog="Example of use: python %(prog)s --tetrises tetrises --output parakeet"
                                            " --GPU_ID 0 1 --dimensions 2048 2048 1024 --phase_plates --phase_shift 90"
                                            " --defoci 0 -0.5 1 --total_doses 75 150 --tilt_steps 2 4"
                                            " --start_angles -60 -40 --ice_densities 0.9 1.1")

    parser.add_argument('--tetrises', type=str, default='tetrises', help='Default: %(default)s. Directory where all the tetrises exist.')
    parser.add_argument('--output', type=str, default='parakeet', help='Default: %(default)s. Directory where the output simulations will be stored.')

    parser.add_argument('--GPU_ID', nargs='+', default=['0', '1'],
                        help="Default: %(default)s. GPU IDs to use for simulation. Example: --GPU_ID 0 1")
    parser.add_argument('--dimensions', type=int, nargs=3, default=[2048, 2048, 1024],
                        help="Default: %(default)s. Dimensions in Angstroms.")
    parser.add_argument('--phase_plates', type=str2bool, default=True, help="Default: %(default)s."
                                                                             " Enable phase plates if set.")
    parser.add_argument('--phase_shift', type=int, default=90, help="Default: %(default)s. Phase shift "
                                                                     "(will not be effective if phase plate are not used).")
    parser.add_argument('--defoci', nargs='+', type=float, default=[0, -0.5, -1],
                        help=" In microns. If you are not using VPP, use values like [-2.5, -3.25, -4]")
    parser.add_argument('--total_doses', nargs='+', type=int, default=[75, 150],
                        help="Total doses to simulate. Example: --total_doses 75 150")
    parser.add_argument('--tilt_steps', nargs='+', type=int, default=[2, 4],
                        help="Tilt steps to simulate. Example: --tilt_steps 2 4")
    parser.add_argument('--start_angles', nargs='+', type=int, default=[-60, -40],
                        help="Start angles to simulate. Example: --start_angles -60 -40")
    parser.add_argument('--ice_densities', nargs='+', type=float, default=[0.9, 1.1],
                        help="Ice densities to simulate. Example: --ice_densities 0.9 1.1")

    args = parser.parse_args()

    main(args.tetrises, args.output, args.GPU_ID, tuple(args.dimensions), args.phase_plates, args.phase_shift, args.defoci,
         args.total_doses, args.tilt_steps, args.start_angles, args.ice_densities)
