
#!/usr/bin/env python3

import numpy as np
import elfi
import argparse
import pickle
import multiprocessing
import logging
import scipy
import math
import matplotlib.pyplot as plt
import glob
from time import gmtime, strftime

logging.basicConfig(level=logging.DEBUG)

logger = logging.getLogger(__name__)

##################################
##handles command line arguements#
##################################

def process_inputs():
    #defines parser
    parser = argparse.ArgumentParser(description='Fits NFDS model of bacterial evolution using BOLFI.')

    #adds individual commmand line arguements
    parser.add_argument('-t',  type= str, required=True,
                        help= 'immigration type [0 - by strain, 1 - by SC]')
    parser.add_argument('-n',  type= str, required=True,
                        help= 'population carrying capacity  [integer]')
    parser.add_argument('-g',  type= str, required=True,
                        help= 'number of generations [integer]')
    parser.add_argument('-u',  type= str, required=True,
                        help= 'upper gene frequency limit [double between 0 and 1]')
    parser.add_argument('-l',  type= str, required=True,
                        help= 'lower gene frequency limit [double between 0 and 1]')
    parser.add_argument('-q',  type= int, required=False,
                        help= 'generation in which vaccine formulation is changed')
    parser.add_argument('-c',  type= str, required=True,
                        help= 'vaccine target COG name [string]')
    parser.add_argument('-f',  type= str, required=True,
                        help= 'input file')
    parser.add_argument('-p',  type= str, required=True,
                        help= 'prefix for output files')

    #converts command line arguements to dictionary
    args = vars(parser.parse_args())

    return(args)
##############################################################################

########################################################################################
#creates elfi simulatr using shell command which is constructed from command line inputs
########################################################################################

#inputs are shell command and the number of simulations to run, returns numpy array of parameter estimates
def create_model(shell_command,m, cog_catergories_dict):

        nfds_sim = elfi.tools.external_operation(shell_command,stdout=True, prepare_inputs = prepare_inputs) #wraps command in a python function
        nfds_sim_vector = elfi.tools.vectorize(nfds_sim) # vectorizes command


        #Defines 3 priors for the model
        elfi.Prior('uniform', -3.5, 3.4, model = m, name = 'v')
        elfi.Prior('uniform', -6, 5.9, model = m, name = 's')
        elfi.Prior('uniform', -3.5, 2.5, model = m, name = 'i')

        #
        for i in cog_catergories_dict:
            elfi.Prior('uniform', 0.01, 0.98, model = m, name = i)

        # create central node that runs the model, I think
        nfds_simulator = elfi.Simulator(nfds_sim_vector, m['v'], m['s'], m['i'],  m['Bacteriocin'], m['BCA'], m['CPS_BCA'], m['CPS'],
                                        m['GeneralMetabolism'], m['GeneralMetabolism_NutrientTransporters'], m['NutrientTransporters'],
                                        m['Prophage_ICE_PRCI_ProphageRemnant_OtherMGE'], m['Resistance'],
                                        name = 'nfds', observed = 0.0)

        return(nfds_simulator)
##############################################################################

########################################################################################
# creates shell command to call C exe for model
########################################################################################

#input is args ductionary from command line arguments
def create_shell_command(args):
    #creates list off all command line arguments
    all_arguments = []
    for i in args:
        if args[i] != None and i != 'p': #assert that commanmd line argument was added
            command_argument = '-' + i + ' ' + args[i]
            all_arguments.append(command_argument)

    all_arguments_string = ' '.join(all_arguments) #converts list to space seperated string

    prefix = args['p']

    shell_command = 'freqDepSelect.exe ' + all_arguments_string #creates full shell command to be passed to elfi # removed by Nick for Mac
    shell_command = shell_command + ' -p f -o ' + prefix +  '_temp_output -v {0} -s {1} -i {2} -w {3} | perl -lane "print @F[2];"'

    return(shell_command)
##############################################################################

##################################################
#prepres inputs so logging can be done on log scale
##################################################
def prepare_inputs(*inputs, **kwinputs):

    input_list = list(inputs)
    for i in range(3):
        input_list[i] = 10**(input_list[i])

    with open('cog_categories_dict.pkl', 'rb') as a:
        cogs_dict = pickle.load(a)
    all_categories_list = list(cogs_dict.keys())

    with open('weighting_file.tsv', 'w') as a:
        for i in cogs_dict:
    	    for j in all_categories_list:
                    if i == j:
                        for cog in cogs_dict[i]:
                            position_in_inputs = all_categories_list.index(j) + 3
                            a.write(cog + '\t' + str(input_list[position_in_inputs]) + '\n')

    input_list = input_list[:3]
    input_list.append('weighting_file.tsv')

    inputs = tuple(input_list)


    # inputs_check = 'Inputs: ' + str(inputs)
    # logger.debug('{}'.format(inputs_check))

    return inputs, kwinputs
##############################################################################
##############################
#returns current time
############################
def get_time():
    return strftime("%Y-%m-%d %H:%M:%S", gmtime())
##############################################################################

# def make_weighting_dict():
#     cog_catergories_dict = {}
#     for i in glob.glob("Q:\\Project2\\nfds_model\\mass_input_files\\mass\\functional_catergories\\*.cogOrdering"):
#         group_name = i.split('Fake_')[1]
#         group_name = group_name.split('_mass.cogOrdering')[0]
#
#         full_path = "Q:\\Project2\\nfds_model\\mass_input_files\\mass\\functional_catergories\\Fake_" + group_name + "_mass.cogOrdering"
#         cogs_list = []
#         with open(full_path, 'r') as a:
#             for line in a:
#                 line = line.split("\n")[0]
#                 cogs_list.append(line)
#         cog_catergories_dict[group_name] = cogs_list
#
#     return cog_catergories_dict

# def make_weighting_file(cog_catergories_dict):
#
#     closest_catergory = 'NutrientTransporters'
#     weighting_file = 'cog_weighting.csv'
#     with open(weighting_file, 'w') as a:
#
#         for i in cog_catergories_dict:
#             if i == closest_catergory:
#                 for cog in cog_catergories_dict[closest_catergory]:
#                     line = cog + ',' + '1' + '\n'
#                         a.write(line)
#
#             else:

if __name__ == '__main__':
    fitting_steps = 200
    samples = 100000

    current_time = get_time()
    logger.debug('Started at {}'.format(current_time))

    args = process_inputs()

    prefix = args['p']

    with open('cog_categories_dict.pkl', 'rb') as a:
        cog_catergories_dict = pickle.load(a)

    elfi.set_client('multiprocessing')

    shell_command = create_shell_command(args) #builds shell command

    m = elfi.ElfiModel(name = 'nfds') #creates model
    nfds_simulator = create_model(shell_command,m, cog_catergories_dict) #creates model


    # Set an arbitrary global seed to keep the randomly generated quantities the same
    seed = 1
    np.random.seed(seed)
    ######################################

    # rejection
    initial_evidence = 5
    rej = elfi.Rejection(nfds_simulator, batch_size = 1)
    n_rejection_samples = 2*initial_evidence
    rej_result = rej.sample(n_samples=n_rejection_samples,quantile = 1)

    update_interval = 1
    # BOLFI
    noise = 0.01
    n_evidence = 10
    parameter_bounds = {'v':(-3.5, -0.1), 's':(-6, -0.1), 'i':(-3.5, -1)}
    for i in cog_catergories_dict:
        parameter_bounds[i] = (0.01, 0.99)

    try:
        with open(prefix + '_run_details.txt', 'w') as a:
            a.write('parameter_bounds: ' + str(parameter_bounds) + '\n')
            a.write('fitting_steps: ' + str(fitting_steps) + '\n')
            a.write('samples: ' + str(samples) + '\n')
            a.close()
    except:
        pass
    # gp = elfi.GPyRegression(m.parameter_names, parameter_bounds)
    # prior = elfi.methods.utils.ModelPrior(m)
    # acq = elfi.methods.bo.acquisition.RandMaxVar(model=gp, prior=prior)
    # bolfi = elfi.BOLFI(nfds_simulator, batch_size=1, initial_evidence=rej_result.outputs, update_interval=update_interval,
    #                     bounds=parameter_bounds, seed=scipy.random.seed(), async = True, target_model=gp, acquisition_method=acq,
    #                      acq_noise_var=[0.1, 0.01, 0.1], batches_per_acquisition=1)

    bolfi = elfi.BOLFI(nfds_simulator, batch_size=1, initial_evidence=rej_result.outputs, update_interval=update_interval,
                        bounds=parameter_bounds, seed=scipy.random.seed(), async = True, batches_per_acquisition=1)


    post = bolfi.fit(n_evidence=fitting_steps) # fits hypermodel to data by minimizing the strain divergence metric - generates very wide probabiolity distributions

    try:
        with open(prefix + '_bolfi_with_posterior.pkl', 'wb') as a:
            pickle.dump(bolfi, a, pickle.HIGHEST_PROTOCOL)
    except:
        pass

    elfi.set_client('native')

    result_BOLFI = bolfi.sample(samples, algorithm = 'metropolis') # samples from probability distributions for each param to get narrower credible intervals

    try:
        with open(prefix + '_bolfi_result.pkl', 'wb') as output:
            pickle.dump(result_BOLFI, output, pickle.HIGHEST_PROTOCOL)
    except:
        pass
    try:
        with open(prefix + '_bolfi_sample_means.txt', 'w') as a:
            a.write(str(result_BOLFI.sample_means))
            a.close()
    except:
        pass

    try:
        # plots
        result_BOLFI.plot_traces();
        plt.savefig(prefix + '_bolfi_sample_traces.png')
        result_BOLFI.plot_marginals();
        plt.savefig(prefix + '_bolfi_sample_marginals.png')
        plt.gcf().clear()
    except:
        pass

    current_time = get_time()
    logger.debug('Finished at {}'.format(current_time))
