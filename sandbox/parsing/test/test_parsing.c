#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "../Parsing.h"

void test_string(int verbose)
{
    char mystring[2][100];
    
    input_string("mystring0", mystring[0], NULL);
    if(verbose)
        printf("mystring0 is %s\n", mystring[0]);
    assert(strcmp(mystring[0], "thisisatest") == 0);

    input_string("mystring1", mystring[1], "notininputfile");
    if(verbose)
        printf("mystring1 is %s\n", mystring[1]);
    assert(strcmp(mystring[1], "notininputfile") == 0);
}

void test_boolean(int verbose)
{
    int i;
    int mybool[4];

    input_boolean("mybool0", &mybool[0], "essential");
    input_boolean("mybool1", &mybool[1], "essential");
    input_boolean("mybool2", &mybool[2], "off");
    input_boolean("mybool3", &mybool[3], "on");
    
    if(verbose)
        for(i = 0; i < 4; i++)
            printf("mybool%d is %d\n", i, mybool[i]);
    
    assert(!mybool[0]);  /* off in file */
    assert(mybool[1]);   /* on in file */
    assert(!mybool[2]);  /* not in file, but defaults to off */
    assert(mybool[3]);   /* not in file, but defaults to on */
}

void test_float(int verbose)
{
    int i, len;
    char variable[24][300];
    float myfloat[24];
    float expect[] = {   123.45,
                         123.45,
                        -123.45,
                         123.45,
                         123.45,
                         0.012345,
                         123.45,
                         123.45,
                         0.012345,
                        -123.45,
                        -123.45,
                        -0.012345,
                         0.12345,
                         0.12345,
                        -0.12345,
                         123.45,
                         12.345,
                         0.0012345,
                         12.345,
                         12.345,
                         0.0012345,
                        -12.345,
                        -12.345,
                        -0.0012345 };

    for(i = 0; i < 24; i++)
    {
        len = sprintf(variable[i], "myfloat%02d", i);
        variable[i][len] = '\0';
    }

    for(i = 0; i < 24; i++)
        input_float(variable[i], &myfloat[i], "1");

    if(verbose)
        for(i = 0; i < 24; i++)
            printf("myfloat%02d is %g\n", i, myfloat[i]);

    for(i = 0; i < 24; i++)
        assert(myfloat[i] == expect[i]);
}

void test_double(int verbose)
{
    int i, len;
    char variable[24][300];
    double mydouble[24];
    double expect[] = {  123.4567,
                         123.4567,
                        -123.4567,
                         123.4567,
                         123.4567,
                         0.01234567,
                         123.4567,
                         123.4567,
                         0.01234567,
                        -123.4567,
                        -123.4567,
                        -0.01234567,
                         0.1234567,
                         0.1234567,
                        -0.1234567,
                         123.4567,
                         12.34567,
                         0.001234567,
                         12.34567,
                         12.34567,
                         0.001234567,
                        -12.34567,
                        -12.34567,
                        -0.001234567 };
    
    
    for(i = 0; i < 24; i++)
    {
        len = sprintf(variable[i], "mydouble%02d", i);
        variable[i][len] = '\0';
    }

    for(i = 0; i < 24; i++)
        input_double(variable[i], &mydouble[i], "1");

    if(verbose)
        for(i = 0; i < 24; i++)
            printf("mydouble%02d is %g\n", i, mydouble[i]);

    for(i = 0; i < 24; i++)
        assert(mydouble[i] == expect[i]);
}

void test_vectors(int verbose)
{
    int i;

    char   my_char_vector[3];
    int    my_int_vector[3];
    float  my_float_vector[3];
    double my_double_vector[4];

    char   expect_char[]   = {'a', 'b', 'c'};
    int    expect_int[]    = {1, 2, 3};
    float  expect_float[]  = {0.0, 0.5, 1.0};
    double expect_double[] = {0.0, 0.33333, 0.66666, 1.0};

    input_char_vector("my_char_vector", 3, my_char_vector);
    if(verbose)
        for(i = 0; i < 3; i++)
            printf("my_char_vector[%d] is %c\n", i, my_char_vector[i]);
    for(i = 0; i < 3; i++)
        assert(my_char_vector[i] == expect_char[i]);

    input_int_vector("my_int_vector", 3, my_int_vector);
    if(verbose)
        for(i = 0; i < 3; i++)
            printf("my_int_vector[%d] is %d\n", i, my_int_vector[i]);
    for(i = 0; i < 3; i++)
        assert(my_int_vector[i] == expect_int[i]);

    input_float_vector("my_float_vector", 3, my_float_vector);
    if(verbose)
        for(i = 0; i < 3; i++)
            printf("my_float_vector[%d] is %g\n", i, my_float_vector[i]);
    for(i = 0; i < 3; i++)
        assert(my_float_vector[i] == expect_float[i]);

    input_double_vector("my_double_vector", 4, my_double_vector);
    if(verbose)
        for(i = 0; i < 4; i++)
            printf("my_double_vector[%d] is %g\n", i, my_double_vector[i]);
    for(i = 0; i < 4; i++)
        assert(my_double_vector[i] == expect_double[i]);

}


int main(int argc, char **argv)
{
    int mybool[2];
    int myint[2];
    float myfloat[24];
    double mydouble[24];
    char mystring[100];

    int i;
    char variable[300];

    int v = 0;

    /* Initialize the input file parser */
    
    if(argc == 1)
        setup_parser("input1", 0);
    else if(argc == 2)
    {
        setup_parser(argv[1], 0);
    }
    else if(argc == 3)
    {
        v = (strcmp(argv[2],"on") == 0);
        setup_parser(argv[1], 0);
    }
    else if(argc == 4)
    {
        v = (strcmp(argv[2], "on") == 0);
        setup_parser(argv[1], (strcmp(argv[3], "on") == 0));
    }
    
    test_string(v);
    test_boolean(v);
    test_float(v);
    test_double(v);
    test_vectors(v);

    return 0;
}
