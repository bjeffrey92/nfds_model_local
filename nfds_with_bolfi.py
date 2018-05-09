#!/usr/bin/env python3

import numpy as np
import elfi
import argparse
import pickle

results_file = '3test_bolfi_results.pkl'
posterior_file = '3test_posterior_results.pkl'

##################################
##handles command line arguements#
##################################

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
parser.add_argument('-r',  type= str, required=True,
                    help= 'cog reordering file')

#converts command line arguements to dictionary
args = vars(parser.parse_args())


########################################################################################
#creates elfi simulatr using shell command which is constructed from command line inputs
########################################################################################

#inputs are shell command and the number of simulations to run, returns numpy array of parameter estimates
def create_model(shell_command):

        m = elfi.ElfiModel(name = 'nfds') #creates model

        # shell_command = './updated_nfds_model/freqDepSelect -f input_files/mass/mass.input -c rmlB -v {0} -u 0.95 -n 50000 -s {1} -g 73 -t 1 -l 0.05 -i {2} -p f -o temp_output | perl -lane "print $F[2];"'
        nfds_sim = elfi.tools.external_operation(shell_command,stdout=True) #wraps command in a python function
        nfds_sim_vector = elfi.tools.vectorize(nfds_sim) # vectorizes command

        # m = elfi.ElfiModel(name = 'nfds') #creates model

        #Defines 3 priors for the model
        elfi.Prior('uniform', 0, 1, model = m, name = 't1')
        elfi.Prior('uniform', 0, 0.05, model = m, name = 't2')
        elfi.Prior('uniform', 0, 1, model = m, name = 't3')

        # create central node that runs the model, I think
        nfds_simulator = elfi.Simulator(nfds_sim_vector, m['t1'], m['t2'], m['t3'],  name = 'nfds', observed = 0.0)

        return(nfds_simulator)

########################################################################################
# creates shell command to call C exe for model
########################################################################################

#input is args ductionary from command line arguments
def create_shell_command(args):

    #creates list off all command line arguments
    all_arguments = []
    for i in args:
        if args[i] != None: #assert that commanmd line argument was added
            command_argument = '-' + i + ' ' + args[i]
            all_arguments.append(command_argument)

    all_arguments_string = ' '.join(all_arguments) #converts list to space seperated string

    # shell_command = './linux_exe/freqDepSelect ' + all_arguments_string #creates full shell command to be passed to elfi
    shell_command = 'freqDepSelect.exe ' + all_arguments_string #creates full shell command to be passed to elfi
    shell_command = shell_command + ' -p f -o 2temp_output -v {0} -s {1} -i {2} | perl -lane "print @F[2];"'
    # shell_command = shell_command + ' -p f -o temp_output -v {0} -s {1} -i {2}'
    # shell_command = 'for /f "tokens=2" %a in ("' + shell_command + '-v' + '{0}' + '-s' + '{1}' + '-i' + '{2}' + '") do echo %a'

    return(shell_command)


shell_command = create_shell_command(args) #builds shell command

nfds_simulator = create_model(shell_command) #creates model


# Set an arbitrary global seed to keep the randomly generated quantities the same
seed = 1
np.random.seed(seed)

bolfi = elfi.BOLFI(nfds_simulator, batch_size=1, initial_evidence=20, update_interval=10,
            bounds={'t1':(0, 1), 't2':(0, 0.05), 't3':(0, 1)}, acq_noise_var=[0.1, 0.01, 0.1], seed=seed) #creates bolfi object which will be fited to data

print('fitting')
post = bolfi.fit(n_evidence=2000) # fits hypermodel to data by minimizing the strain divergence metric, generates very wide probabiolity distributions
print('fitted')
with open(posterior_file, 'wb') as output:
    pickle.dump(post, output, pickle.HIGHEST_PROTOCOL)

result_BOLFI = bolfi.sample(4000) # samples from probability distributions for each param to get narrower credible intervals
with open(results_file, 'wb') as output:
    pickle.dump(result_BOLFI, output, pickle.HIGHEST_PROTOCOL)
