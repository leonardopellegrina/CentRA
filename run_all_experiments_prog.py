import os
from multiprocessing import Pool
import random

num_tests = 10
eps = [0.2 , 0.1 , 0.05]
ks = [10 , 25 , 50 , 75 , 100]
num_cores = 20

def run_command(cmd):
    os.system("sleep "+str(random.randint(0, 9)))
    os.system(cmd)

pool_ = Pool(processes=num_cores)

graph_paths = ["actor-collaboration_pre.txt",
    "com-amazon.ungraph_pre.txt",
    "com-dblp.ungraph_pre.txt",
    "com-youtube.ungraph_pre.txt",
    "email-Enron_pre.txt",
    "soc-LiveJournal1_pre.txt",
    "soc-pokec-relationships_pre.txt",
    "wiki-Talk_pre.txt",
    "wiki-Vote_pre.txt",
    "wiki-topcats_pre.txt"]
directed_graphs = [
    "soc-LiveJournal1_pre.txt",
    "soc-pokec-relationships_pre.txt",
    "wiki-Talk_pre.txt",
    "wiki-topcats_pre.txt",
    "wiki-Vote_pre.txt"
    ]



list_commands = []
for i in range(num_tests):
    for u in [1]:#[0,1]:
        for k in ks:
            for graph_path in graph_paths:
                for epsilon in eps:
                    cmd = "python run_centra.py -db "+str(graph_path)+" -k "+str(k)+" -p 1 -e "+str(epsilon)
                    if graph_path in directed_graphs:
                        cmd = cmd + " -t 1"
                    cmd = cmd + " -u "+str(u)
                    list_commands.append(cmd)

res_list = []
for i in range(len(list_commands)):
    cmd_i = list_commands[i]
    out_path_i = "output_prog_exp_"+str(i)+".txt"
    cmd_i = cmd_i + " -o "+out_path_i
    res = pool_.apply_async(run_command, (cmd_i,))
    res_list.append(res)

for res in res_list:
    res.get()
