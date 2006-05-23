/* .indent.pro -- a configuration file for /usr/bin/indent */
/* Refer to indent(1) for an explanation of each option */

/* Setting the maximum line length for code and comments */
--line-length500
--comment-line-length80

/* Setting indentation options */
--no-tabs
--tab-size4
--indent-level4 
--parameter-indentation4 /* Indent parameter types in old-style func defs */
--continue-at-parentheses /* Line up continued lines at parentheses */

/* Controlling the placement of braces */
--brace-indent0
--braces-after-if-line
--braces-after-struct-decl-line
--dont-cuddle-else
--cuddle-do-while

/* Controlling the placement of spaces */
--space-special-semicolon
--no-space-after-function-call-names
--no-space-after-parentheses
--no-space-after-casts
--no-space-after-while
--no-space-after-for
--no-space-after-if
--dont-break-procedure-type

/* Controlling the placement of blank lines */
--no-blank-lines-before-block-comments
--no-blank-lines-after-declarations
--blank-lines-after-procedures

/* Put the '*' character at the left of comments. */
--start-left-side-of-comments

/* Preserve access and modification times on output files */
--preserve-mtime
