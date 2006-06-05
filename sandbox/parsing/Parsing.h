#ifndef __PARSING_H__
#define __PARSING_H__

/* Parsing.c */
void setup_parser(char *filename, int verbose_output);
void shutdown_parser(void);
int add_to_parameter_list(register char *name, register char *value);
int compute_parameter_hash_table(register char *s);
int input_string(char *name, char *value, char *Default);
int input_boolean(char *name, int *value, char *interpret);
int input_float(char *name, float *value, char *interpret);
int input_int(char *name, int *value, char *interpret);
int input_double(char *name, double *value, char *interpret);
int input_char_vector(char *name, int number, char *value);
int input_int_vector(char *name, int number, int *value);
int input_float_vector(char *name, int number, float *value);
int input_double_vector(char *name, int number, double *value);
int interpret_control_string(char *interpret, int *essential, double *Default, double *minvalue, double *maxvalue);
/* Generated with: /usr/bin/cproto -q -p -f 3 Parsing.c 2>/dev/null */

#endif
