# SCRIPT FROM Hamish McGregor to run the simulation for power analysis, received 07/08/2024
# I adapted file paths to point to files / folders on my computer (local)

# 2024-08-06
# simulation of uk Biobank without ageing or competition - single mutants only?
# this version without double mutants (for now) and unlimited cell pop size
# produces only biobank-realistic file to save space

# adapted from 2022-11-07 clonal_competition_simple background

# with a significant amount of cleanup - removed broken plots at the end and various redundancies such as
# caroline's long lists of colours and plotting options

# how to run the simulation
# cd ~/Desktop/msc_thesis/task1_predict_binding_to_HLA/scripts02
# python3 9a_simulation_code.py --name 'DNMT3A_R882H_015' --mu 1e-08 --s 0.15 --N 100000 --depth 70 --lifespan 74

import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import scipy.integrate as it
import scipy
import copy
# import datetime
import csv
import os
from argparse import ArgumentParser


# arguments
parser = ArgumentParser()
parser.add_argument("--name", dest='sim_name', help='simulation_name', required=True)
parser.add_argument("--mu", dest="mu", help="variant mutation rate", type=float, required=True)
parser.add_argument("--s", dest="s", help="variant fitness", type=float, required=True)
#parser.add_argument("--tau-params", nargs="+", dest="tau_params", help='tau_start tau_end tau_age1 tau_age2 space-separated', required=True)
parser.add_argument("--N", dest="N", help="stem cell initial pop size", type=int, required=True)
#parser.add_argument("--sample-ages", nargs="+", dest="sample_ages", help="ages at blood draw space-separated", required=True)
parser.add_argument("--depth", dest="depth", help="sequencing depth", type=int, required=True)
parser.add_argument("--lifespan", dest="lifespan", help="simulation lifetime (years)", type=int, required=True)
#parser.add_argument("--quantiles", dest="number_of_quantiles", help="number of quantiles for plotting", type=int, required=True)
#parser.add_argument("--plot-save-name", dest="plot_save_name", help="save name for final plot", required=True)
o = parser.parse_args()

seed = 3
np.random.seed(seed)

# read biobank data
# import curated phenotype data (only includes individuals present in the exome data)
print('Importing biobank phenotypes...')
biobank_healthy_pheno = pd.read_csv('/Users/barbarawalkowiak/Desktop/msc_thesis/data/2024-04-24_healthy_pheno_more_cancer_cols.tsv', sep='\t', low_memory=False)

print('Total individuals in phenotype data:', len(biobank_healthy_pheno))
biobank_ages = biobank_healthy_pheno[['ID_v0', 'Age.when.attended.assessment.centre_v0']].copy().reset_index(drop=True)

# age bins
biobank_ages_bins = biobank_ages['Age.when.attended.assessment.centre_v0'].value_counts().reset_index().sort_values('Age.when.attended.assessment.centre_v0').rename(
    columns={'count': 'age_count'}).reset_index(drop=True)

# define function 'divide'
def divide(clone_info, dt, B0, D0):
    # dt = a small interval of time,
    # D0 = symmetric rate to differentiated
    # B0 = symmetric rate to self renewal
    # B  = modified birth rate including fitness advantage

    mutations = clone_info[0]  # [(ID1, s1)]
    clone_size = clone_info[1]  # n
    fitness_effects = [i[1] for i in mutations]
    F = sum(fitness_effects)
    B = B0 * (
            1.0 + F)  # B = birth rate, which is 1 + the overall fitness effect of the mutations within a particular clone
    number_of_births = np.random.poisson(
        clone_size * B * dt)  # pulls a random number from the Poisson distribution with/
    # mean of clone_size x birth rate x interval of time
    number_of_deaths = np.random.poisson(clone_size * D0 * dt)
    new_clone_size = clone_size + number_of_births - number_of_deaths
    return [mutations, new_clone_size]

# define function 'mutation_fitness' of the form:
# sb=0.05 #beneficial fitness effect
# sd=0.00 #deleterious fitness effect
# Pb=2/3 #probability of beneficial fitness effect
# Pn=1/3 #probability of neutral fitness effect
# Pd=0.00 #probability of deleterious fitness effect
# DFE = [[sd, Pd], [0.0, Pn], [sb, Pb]] #distribution of fitness effects doesn't need to be normalized
# (because normalized in function below)
def mutation_fitness(DFE):
    total_prob = sum([i[1] for i in DFE])
    normalized_DFE = [[i[0], i[1] / total_prob] for i in DFE]

    random_number = np.random.random()
    cumulative = 0.0
    for [s, P] in normalized_DFE:
        cumulative = cumulative + P
        if random_number < cumulative:
            fitness_effect = s
            break
    return fitness_effect


# define function 'mutate'
def mutate(clone_info, dt, u, last_mutation_ID, DFE): #u=mutation rate
    mutations=clone_info[0] #(ID, s)
    clone_size=clone_info[1]

    # prevent a second mutation in an already mutated clone
    if len(mutations) > 1:
        return([], last_mutation_ID)

    number_of_mutations = np.random.poisson(clone_size*u*dt)
    list_of_new_clones=[]
    for i in range(number_of_mutations):
        last_mutation_ID = last_mutation_ID+1
        new_fitness_effect = mutation_fitness(DFE)
        new_clone = [mutations+[(last_mutation_ID, new_fitness_effect)],1]
        list_of_new_clones.append(new_clone)
    return (list_of_new_clones, last_mutation_ID)


# define function 'pop_size'
def pop_size(pop):
    total=0
    for clone in pop:
        sub_pop_size=clone[-1]
        total=total+sub_pop_size
    return total

# define function 'sequence' (to include R882 variants)
# change so now there are two depths: one deep and one shallow, to compare
def sequence(pop, mean_depth1, depth_SD1, mean_depth2, depth_SD2):
    total_pop_size = pop_size(pop)  # this is confusing: used to be called N (which is a constant =1e5 in main code). Now called total_pop_size
    variant_dict = {}
    list_of_variants = []

    # first determine total freq of variant across all clones and record clones in which it arises
    for (clone1, size) in pop:
        if len(clone1) > 1:  # ensures non-WT clones only
            clone = copy.deepcopy(clone1)
            clone.pop(0)

            multi_mutant = len(clone) # number of mutations in clone (not including WT)
            for mutation in clone:

                # update mutations frequency
                mID = mutation[0]

                if mID in variant_dict:

                    # need to check what this actually does. Don't think it's required with single mutants only
                    variant_dict[mID]['true_freq'] += 0.5 * (size / total_pop_size)
                    # variant_dict[mID]['true_VAF'] += 0.5*(size/total_pop_size)

                else:
                    variant_dict[mID] = {
                        'true_freq': 0.5 * (size / total_pop_size),
                        'true_clone_size': size,
                        'total_pop_size': total_pop_size}  # , 'true_VAF':0.5*(size/total_pop_size)}

                    # mutation type - modified (see below)
                    # if 0.00<mutation[1]<R882:
                    #    mutation_type="n" # non-synonymous (with a range of s)
                    # elif mutation[1]==0.0:
                    #    mutation_type="s" #synonymous
                    # elif mutation[1]==R882:
                    #    mutation_type="R882" #R882
                    # variant_dict[mID]['mutation_type']=mutation_type

                    if mutation[1] == s_R882H:
                        mutation_type = 'R882H'
                    elif mutation[1] > 0:
                        mutation_type = 'driver'
                    variant_dict[mID]['mutation_type'] = mutation_type

                    variant_dict[mID]['fitness'] = mutation[1]
                    variant_dict[mID]['multi_mutant'] = multi_mutant

    for m, e in variant_dict.items():

        mutation_type = e['mutation_type']
        true_freq = e['true_freq']
        fitness = e['fitness']
        multi_mutant = e['multi_mutant']
        true_clone_size = e['true_clone_size']
        popsize_temp = e['total_pop_size']
        # true_VAF=e['true_VAF']

        # depth1 = int(mean_depth1*np.exp(np.random.normal(0.0, depth_SD1)))
        # if depth1 <= 0:
        #    depth1 = 1
        # expected_reads1=true_freq*depth1
        ##reads1=np.random.poisson(expected_reads1)
        # reads1=np.random.binomial(depth1, true_freq)
        #
        # depth2 = int(mean_depth2*np.exp(np.random.normal(0.0, depth_SD2)))
        # if depth2 <= 0:
        #    depth2 = 1
        # expected_reads2=true_freq*depth2
        ##reads2=np.random.poisson(expected_reads2)
        # reads2=np.random.binomial(depth2, true_freq)

        # if reads1<4:
        #    reads1=0
        # if reads2<4:
        #    reads2=0

        # new method for making normal distribution of depths
        depth1 = int(np.random.normal(mean_depth1, depth_SD1))
        if depth1 <= 1:
            depth1 = 2
        reads1 = np.random.binomial(depth1, true_freq)

        depth2 = int(np.random.normal(mean_depth2, depth_SD2))
        if depth2 <= 1:
            depth2 = 2
        reads2 = np.random.binomial(depth2, true_freq)

        freq1 = reads1 / depth1
        freq2 = reads2 / depth2

        list_of_variants.append({"mutation_ID": m,
                                 "mutation_type": mutation_type,
                                 "fitness": fitness,
                                 "true_clone_size": true_clone_size,
                                 "popsize_temp": popsize_temp,
                                 "true_freq": true_freq,
                                 "multi_mutant":multi_mutant,
                                 # "true_VAF": true_VAF,
                                 "freq1": freq1,
                                 "reads1": reads1,
                                 "depth1": depth1,
                                 "freq2": freq2,
                                 "reads2": reads2,
                                 "depth2": depth2
                                 })
    return list_of_variants


# region PARAMETERS
# distribution of fitness: delta function at s_R882H and second delta function at specified 'background' mu/s
# try a wide variety of different 'background' mutations

# mutation rates
u_R882H =  o.mu             #15e-9  #18.82e-9  #  all per year
s_R882H = o.s

u_background = 0            # fixed at zero for variable tau model

lam = 5
DFE = ((s_R882H, u_R882H), (0.15, 0)) # fake s_background and u_background with zero mutatiion rate - otherwise the dfe doesn't work
u = u_R882H/lam #+ u_background/lam

# other parameters
eq_population_size = o.N  # number of stem cells (HSCs)
lifespan = o.lifespan*lam  # measured in cell divisions (e.g. 100 = 100 cell divisions) - try not to extend too much otherwise sim stops often (overloaded with variants)
dt = 1  # measured in units of the overall cell division rate

last_mutation_ID = 0

depth1 = 500
depth2 = o.depth
depthsd1 = 0
depthsd2 = 0

# endregion

# region NAMES

extended_dirname = f'/Users/barbarawalkowiak/Desktop/msc_thesis/results/simulation/{o.sim_name}'
os.mkdir(extended_dirname)

biobank_variants_name = extended_dirname+'/biobank_variants.csv'
# endregion


# region SIMULATION

population_dict = {}
clone_histories = {}

with open(biobank_variants_name, 'w') as biobank_file:

    writer_biobank = csv.writer(biobank_file, delimiter=',')
    writer_biobank.writerow(["person_ID",
                             "age",
                             'mutation_ID',
                             "mutation_type",
                             "fitness",
                             "multi_mutant",
                             "true_clone_size",
                             "popsize_temp",
                             "true_freq",
                             # "true_VAF",
                             "freq1",
                             "var_reads1",
                             "depth1",
                             "freq2",
                             "var_reads2",
                             "depth2"])
    number_alive_at_end = []
    for person_index in range(len(biobank_ages)):
        if person_index % 1000 == 0:
            print(person_index)

        person_ID = biobank_ages.ID_v0[person_index]
        ukb_blood_draw = biobank_ages['Age.when.attended.assessment.centre_v0'][person_index]

        n0 = eq_population_size  # starting number of cells in first clone
        founding_pop = [[[(0, 0.0), ], n0]]
        development_phase = 1
        t = 0.0

        population = founding_pop
        population_size = n0

        while t < lifespan:
            if population_size < eq_population_size and development_phase == 1:
                t = t + dt
                B0 = 1
                D0 = 0

                population_size = pop_size(population)

                new_population = []
                for clone_info in population:
                    # first put the ID and clone size in our 'clone histories' dictionary
                    clone_ID = tuple(clone_info[0])
                    clone_size = clone_info[1]
                    if clone_ID not in clone_histories:
                        clone_histories[clone_ID] = {'sizes': [clone_size], 'times': [t]}
                    else:
                        clone_histories[clone_ID]['sizes'].append(clone_size)
                        clone_histories[clone_ID]['times'].append(t)
                    # then update clone sizes/ population list
                    updated_clone = divide(clone_info, dt, B0, D0)
                    (new_clones, new_latest_ID) = mutate(clone_info, dt, u,
                                                         last_mutation_ID, DFE)
                    last_mutation_ID = new_latest_ID
                    new_population = new_population + new_clones
                    if updated_clone[1] > 0:
                        new_population.append(updated_clone)
                population = new_population

            else:
                development_phase = 0  # stops the dev phase first time popsize > eq_pop_size

                t = t + dt
                tau = 1 #tau_of_t(t / 5, *tau_params)
                B0 = 0.2 / tau
                D0 = 0.2 / tau

                population_size = pop_size(population)
                # print(population_size)

                new_population = []
                for clone_info in population:

                    clone_ID = tuple(clone_info[0])
                    clone_size = clone_info[1]

                    if person_ID < 500:
                        # first put the ID and clone size in our 'clone histories' dictionary
                        if t % 1 < dt:
                            if clone_ID not in clone_histories:
                                clone_histories[clone_ID] = {'sizes': [clone_size], 'times': [t],
                                                             'VAF': [0.5 * clone_size / population_size]}
                            else:
                                clone_histories[clone_ID]['sizes'].append(clone_size)
                                clone_histories[clone_ID]['times'].append(t)
                                clone_histories[clone_ID]['VAF'].append(0.5 * clone_size / population_size)

                    # then update clone sizes/ population list
                    updated_clone = divide(clone_info, dt, B0, D0)
                    (new_clones, new_latest_ID) = mutate(clone_info, dt, u, last_mutation_ID, DFE)
                    last_mutation_ID = new_latest_ID
                    new_population = new_population + new_clones
                    if updated_clone[1] > 0:
                        new_population.append(updated_clone)
                population = new_population


                ukb_blood_draw_time = int(ukb_blood_draw * lam / dt)
                #time_5 = lifespan / dt


                # sequence for biobank file
                if int(t / dt) == ukb_blood_draw_time:
                    biobank_variants = sequence(population, depth1, depthsd1, depth2, depthsd2)
                    for var in biobank_variants:
                        writer_biobank.writerow([person_ID,
                                                 int(t/lam),
                                                 var['mutation_ID'],
                                                 var['mutation_type'],
                                                 np.round(var['fitness'], 3),
                                                 var['multi_mutant'],
                                                 var['true_clone_size'],
                                                 var['popsize_temp'],
                                                 var['true_freq'],
                                                 var['freq1'],
                                                 var['reads1'],
                                                 var['depth1'],
                                                 var['freq2'],
                                                 var['reads2'],
                                                 var['depth2']])

        number_alive_at_end.append(len(population))
# endregion
