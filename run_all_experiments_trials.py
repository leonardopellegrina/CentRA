import os
from multiprocessing import Pool
import random

num_tests = 10
ms = [500000]
ks = [50]
ts = [1 , 10 , 50 , 100 , 250 , 500]
num_cores = 10

def run_command(cmd):
    os.system("sleep "+str(random.randint(0, 5)))
    os.system(cmd)

pool_ = Pool(processes=num_cores)

graph_paths = ["com-youtube.ungraph_pre.txt",
    "soc-LiveJournal1_pre.txt",
    "wiki-Vote_pre.txt",
    "email-Enron_pre.txt"]
directed_graphs = ["soc-LiveJournal1_pre.txt","wiki-Vote_pre.txt"]

graphs_rel_dir = ""

list_commands = []
for i in range(num_tests):
    for m in ms:
        for k in ks:
            for graph_path in graph_paths:
                for t in ts:
                    cmd = "python run_centra.py -db "+graphs_rel_dir+str(graph_path)+" -k "+str(k)+" -p 0 -m "+str(m)+" -c "+str(t)
                    cmd = cmd + " -r results_trials.csv"
                    if graph_path in directed_graphs:
                        cmd = cmd + " -t 1"
                    list_commands.append(cmd)

res_list = []
for i in range(len(list_commands)):
    cmd_i = list_commands[i]
    out_path_i = "output_trials_exp_"+str(i)+".txt"
    cmd_i = cmd_i + " -o "+out_path_i
    res = pool_.apply_async(run_command, (cmd_i,))
    res_list.append(res)

for res in res_list:
    res.get()
