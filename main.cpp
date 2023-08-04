#include <iostream>
#include <dirent.h>
#include <limits.h>
#include <math.h>
#include <ctime>
#include <cstddef>
#include <unistd.h>
#include <ctype.h>

#include "Rand_gen.h"
#include "Graph.h"
#include "utilities.h"
#include "Probabilistic.h"

extern char *optarg;
static const std::string ERROR_HEADER = "ERROR: ";
using namespace std;

bool directed = false;
bool progressive = false;
bool progressive_union_bound = false;
double verb = 60;
double delta;
double err;
char *graph_file;
std::string output_file;
int64_t k = 0;
int64_t fixed_num_samples = 0;
double sampling_rate = 2.3;
bool alpha_given = false;
bool m_given = false;
bool eps_given = false;
double empirical_peeling_param = 2.0;

// mcrade
int num_mc = 10;

/**
 * Print usage on stderr.
 */
void usage(const char *binary_name) {
    std::cerr << binary_name
        << ": compute betweenness centrality approximations for all nodes"
        << std::endl;
    std::cerr << "USAGE: " << binary_name << " [-dh] [-v verbosity] [-o output] [-a a_emp_peeling] [-s alpha] [-e epsilon] k delta graph"
        << std::endl;
    std::cerr << "\t-d: consider the graph as directed" << std::endl;
    std::cerr << "\t-k: compute the top-k betweenness centralities (if 0, compute all of them with absolute error) " << std::endl;
    std::cerr << "\t-h: print this help message" << std::endl;
    std::cerr << "\t-v: print additional messages (verbosity is the time in second between subsequent outputs)" << std::endl;
    std::cerr << "\t-o: path for the output file (if empty, do not write the output)" << std::endl;
    std::cerr << "\t-a: parameter a for empirical peeling (def. = 2)" << std::endl;
    std::cerr << "\t-s: parameter alpha for sampling shortest paths (def. = 2.3)" << std::endl;
    std::cerr << "\terr: accuracy (0 < epsilon < 1), relative accuracy if k > 0" << std::endl;
    std::cerr << "\tdelta: confidence (0 < delta < 1)" << std::endl;
    std::cerr << "\tgraph: graph edge list file" << std::endl;
}

/**
 * Parse command line options.
 * Return 0 if everything went well, 1 if there were errors, 2 if -h was specified.
 */
int parse_command_line(int& argc, char *argv[]) {
    int opt;
    while ((opt = getopt(argc, argv, "dpuho:s:a:v:m:t:e:")) != -1) {
        switch (opt) {
        case 'd':
            directed = true;
            break;
        case 'p':
            progressive = true;
            break;
        case 'u':
            progressive_union_bound = true;
            break;
        case 'h':
            return 2;
            break;
        case 'o':
            std::cerr << "Writing output to " << optarg << std::endl;
            output_file = optarg;
            break;
        case 's':
            sampling_rate = std::strtod(optarg, NULL);
            alpha_given = true;
            break;
        case 'a':
            empirical_peeling_param = std::strtod(optarg, NULL);
            if (errno == ERANGE || empirical_peeling_param <= 1) {
                std::cerr << ERROR_HEADER
                    << "The value a should be >= 1."
                    << std::endl;
                return 1;
            }
            alpha_given = true;
            break;
        case 'e':
            err = std::strtod(optarg, NULL);
            if (errno == ERANGE || err >= 1.-exp(-1.) || err <= 0.) {
                std::cerr << ERROR_HEADER
                    << "The value of epsilon should be >= 0 and <= 1 - 1/e"
                    << std::endl;
                return 1;
            }
            eps_given = true;
            break;
        case 't':
            num_mc = std::strtod(optarg, NULL);
            if (errno == ERANGE || num_mc <= 0 || num_mc > UINT_MAX) {
                std::cerr << ERROR_HEADER
                    << "The value of MC trials t should be between 1 and 2^32-1."
                    << std::endl;
                return 1;
            }
            break;
        case 'm':
            fixed_num_samples = std::strtod(optarg, NULL);
            m_given = true;
            if (errno == ERANGE || fixed_num_samples <= 0 || fixed_num_samples > UINT_MAX) {
                std::cerr << ERROR_HEADER
                    << "The value m should be between 1 and 2^32-1."
                    << std::endl;
                return 1;
            }
            break;
        case 'v':
            verb = std::strtod(optarg, NULL);
            if (errno == ERANGE || verb < 0) {
                std::cerr << ERROR_HEADER
                    << "The verbosity should be a positive number, or 0 to produce no output."
                    << std::endl;
                return 1;
            }
            break;
        }
    }

    if (optind != argc - 3) {
        std::cerr << ERROR_HEADER << "Wrong number of arguments" << std::endl;
        return 1;
    } else {
        k = std::strtod(argv[argc - 3], NULL);
        if (errno == ERANGE || k <= 0 ) {
            std::cerr << ERROR_HEADER <<
                "The parameter k should be greater than 0"
                << std::endl;
            return 1;
        }
        delta = std::strtod(argv[argc - 2], NULL);
        if (errno == ERANGE || delta >= 1.0 || delta <= 0.0) {
            std::cerr << ERROR_HEADER <<
                "Delta should be greater than 0 and smaller than 1"
                << std::endl;
            return 1;
        }
        graph_file = argv[argc - 1];
    }
    if(progressive == false && m_given == false){
      std::cerr << "Fixed size sample mode selected, but no value of sample size m is given " << std::endl;
      return 1;
    }
    if(progressive == true && eps_given == false){
      std::cerr << "Progressive sampling mode selected, but no value of epsilon is given " << std::endl;
      return 1;
    }

    return 0;
}

int main(int argc, char *argv[]){
    int correct_parse = parse_command_line(argc, argv);

    if (correct_parse != 0) {
        usage(argv[0]);
        return correct_parse!=2;
    }

    Probabilistic G( graph_file, directed, verb , num_mc , sampling_rate , alpha_given , empirical_peeling_param , output_file);
    G.run((uint32_t) k, delta, err, fixed_num_samples, progressive , progressive_union_bound);
    std::cout << "run finished" << std::endl;
    return 0;
}
