import os
from multiprocessing import Pool
from tqdm import tqdm

def run_command(cmd):
    os.system(cmd)

num_cores = 32
pool_ = Pool(processes=num_cores)

num_tests = 10
ms = [50000 , 100000 , 200000 , 500000 , 1000000]
ks = [10 , 25, 50, 75, 100]
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
    for graph_path in graph_paths:
        for m in ms:
            for k in ks:
                cmd = "python run_centra.py -db "+str(graph_path)+" -k "+str(k)+" -p 0 -m "+str(m)
                if graph_path in directed_graphs:
                    cmd = cmd + " -t 1"
                #os.system(cmd)
                list_commands.append(cmd)

res_list = []
for i in range(len(list_commands)):
    cmd_i = list_commands[i]
    out_path_i = "output_fixed_exp_"+str(i)+".txt"
    cmd_i = cmd_i + " -o "+out_path_i
    res = pool_.apply_async(run_command, (cmd_i,))
    res_list.append(res)

for res in tqdm(res_list):
    res.get()
