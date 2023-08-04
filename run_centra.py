import os
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-db", help="path to graph input file")
parser.add_argument("-e","--epsilon", type=float ,help="approximation accuracy parameter (in (0,1))",default=0.1)
parser.add_argument("-d","--delta", type=float ,help="approximation confidence parameter (in (0,1), def. 0.1)",default=0.05)
parser.add_argument("-k","--k", type=int ,help="number of nodes to select in the greedy algorithm",default=10)
parser.add_argument("-c","--mctrials", type=int ,help="number of mc trials to do to estimate the mcera",default=100)
parser.add_argument("-p","--progressive", type=int ,help="1 means progressive sampling, 0 fixed number of samples (def. progressive)",default=1)
parser.add_argument("-m","--numsamples", type=int ,help="number of samples to take when progressuve = 0",default=10000)
parser.add_argument("-t","--type", type=int ,help="type of graph. 1 means directed, 0 undirected (def. undirected)",default=0)
parser.add_argument("-u","--unionbound", type=int ,help="1 use the union bound to check the stopping condition",default=0)
parser.add_argument("-o","--outputpath" ,help="output file path to use",default="")
parser.add_argument("-r","--resultspath" ,help="output file path to use for results",default="")
args = parser.parse_args()

delta = args.delta
file_output_path = "output_fixed.txt"
path_executable = "./centra"
output_path = "results_fixed.csv"
directed_flag = ""
if args.type == 1:
    directed_flag = "-d "

if args.progressive == 1:
    file_output_path = "output_prog_centra.txt"
    output_path = "results_prog_centra.csv"
    cmd = path_executable+" "+str(args.k)+" "+str(args.delta)+" "+str(args.db)+" "+directed_flag+"-p -e "+str(args.epsilon)+" -t "+str(args.mctrials)
    if args.unionbound == 1:
        cmd = cmd+" -u "
        file_output_path = "output_prog_unionbound.txt"
        output_path = "results_prog_unionbound.csv"
else:
    cmd = path_executable+" "+str(args.k)+" "+str(args.delta)+" "+str(args.db)+" "+directed_flag+"-m "+str(args.numsamples)+" -t "+str(args.mctrials)
if len(args.outputpath) > 0:
    file_output_path = args.outputpath
cmd = cmd+" > "+file_output_path

if len(args.resultspath) > 0:
    output_path = args.resultspath

print(cmd)
os.system(cmd)


if os.path.isfile(output_path) == False:
    # write the header for the results
    out = "graph_name;epsilon;delta;k;type;numtrials;numsamples;time;diameter;iterations;est_opt;sd_bound_centra;sd_bound_unionbound;mcera_est\n"
    out_file = open(output_path,"w")
    out_file.write(out)
    print(out)
    out_file.close()

def get_result(pattern , path ,  verbose=1):
    fin = open(path,'r')
    to_return = ""
    line_to_print = ""
    for line in fin:
        if pattern in line:
            line = line.replace('\n','')
            line_to_print = line
            to_return = line[len(pattern):]
    fin.close()
    if verbose == 1 and len(line_to_print) > 0:
        print(line_to_print)
    return to_return


# results to gather from output file
results_strings = ["time for second pass ","estimated diameter of the graph: ","iteration_index ", "covered_samples_fraction ", "sd_bound ","sd_bound_ub ","final nmcera: "]

out_file = open(output_path,"a")
line_out = str(args.db)+";"+str(args.epsilon)+";"+str(args.delta)+";"+str(args.k)+";"+str(args.type)+";"+str(args.mctrials)+";"

if args.progressive == 0:
    line_out = line_out + str(args.numsamples)+";"
else:
    value_result = get_result("FINISHED " , file_output_path ,  1)
    line_out = line_out + str(value_result)+";"

for idx, result_string in enumerate(results_strings):
    value_result = get_result(result_string , file_output_path ,  1)
    line_out = line_out+value_result
    if idx < len(results_strings)-1:
        line_out = line_out+";"
line_out = line_out+"\n"
print(line_out)
out_file.write(line_out)
