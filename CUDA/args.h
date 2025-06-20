#ifndef ARGS_H
#define ARGS_H

#ifdef __cplusplus
extern "C" {
#endif

extern int verbose;
extern int input_problem;
extern int no_output;
extern int output_freq;
extern int enable_checkpoints;

void parse_args(int argc, char *argv[]);
void print_opts();

#ifdef __cplusplus
}
#endif

#endif
