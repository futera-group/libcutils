/******************************************************************************\
 *                                                                            * 
 *  Libcutils - library of C function                                         * 
 *                                                                            *
 *  Version:             3.4                                                  * 
 *  Date:                28/01/2017                                           *
 *                                                                            * 
 *  Author:              Zdenek Futera                                        * 
 *                                                                            * 
 *  Address:             University College London (UCL)                      * 
 *                       Department of Physics and Astronomy                  * 
 *                       Gower Street, London WC1E 6BT                        * 
 *                       United Kingdom                                       * 
 *                                                                            * 
 *  E-Mail:              z.futera@ucl.ac.uk                                   * 
 *                                                                            * 
\******************************************************************************/

#include <stdlib.h>
#include "cmn/arg.h"
#include "cmn/string.h"
#include "cmn/types.h"

/* -------------------------------------------------------------------------- */

/* return pointer of argument structure with given flag name
   
   a      - command-line arguments
   name   - flag name of the argument
   single - single switch / item */
struct arg_dat* arg_read_get_ptr(struct arg *a, char name, short single) {
  struct arg_dat *t = NULL;
  struct ldata *p;
  if (a && a->data && a->data->first) { 
    for (p=a->data->first; p; p=p->l_next) {
      t = (struct arg_dat*)p->l_data;
      /* single argument */
      if (single) {
        if (!t->swtch && t->simple && !t->set)
          return(t);
        }
      /* named option */
      else {
        if (t->name && t->name[0]==name)
          return(t);
        }
      }
    }
  return(NULL);
  }

/* read and analyze command line arguments
   
   a    - list of possible and required argument
   argc - number of given arguments
   argv - list of given arguments */
void arg_read(struct arg *a, unsigned argc, char **argv) {
  char *program_name,*msg;
  size_t arg_length;
  unsigned i,j;
  struct arg_dat *t;
  /* name of the program */
  program_name = str_file_name_new(argv[0]);
  for (i=1; i<argc; i++) {
    /* switches and options */
    if (argv[i][0]=='-') {
      arg_length = str_length(argv[i]);
      /* single switch or option */
      if (arg_length==2) {
        /* print help */
        if (a->help && (argv[i][1]=='h' || argv[i][1]=='H')) {
          arg_print_usage(a,program_name,NULL,0);
          a->help_fce();
          exit(0);
          }
        /* assing name to argument definition */
        t = arg_read_get_ptr(a,argv[i][1],0);
  	if (!t) {
          msg = str_new(18+str_length(argv[i]));
          sprintf(msg,"unknown option \"%s\"",argv[i]);
  	  arg_print_usage(a,program_name,msg,1);
  	  }
        /* switch */
  	else if (t->swtch) {
  	  (*(short*)t->val) = 1;
          if (t->simple)
            return;
          }
        /* repeated option */
  	else if (t->set) {
          msg = str_new(38+str_length(argv[i]));
          sprintf(msg,"more than one occurence of option \"%s\"",argv[i]);
  	  arg_print_usage(a,program_name,msg,1);
  	  }
        else {
          if ((i+1)==argc) {
            msg = str_new(24+str_length(argv[i]));
            sprintf(msg,"option %s without value",argv[i]);
            arg_print_usage(a,program_name,msg,1);
            }
          else if (!type_read(argv[i+1],t->val,t->type)) {
            msg = str_new(38+str_length(argv[i])+str_length(argv[i+1]));
  	    sprintf(msg,"cannot assign value \"%s\" to option \"%s\"",
              argv[i+1],argv[i]);
            arg_print_usage(a,program_name,msg,1);
            }
          else
            i++;
          }
  	t->set = 1;
        }
      /* grouped switches */
      else if (arg_length>2) {
        for (j=1; j<arg_length; j++) {
          t = arg_read_get_ptr(a,argv[i][j],0);
          if (!t) {
            msg = str_new(20);
            sprintf(msg,"unknown option \"-%c\"",argv[i][j]);
            arg_print_usage(a,program_name,msg,1);
            }
          else if (!t->swtch)
            arg_print_usage(a,program_name,NULL,1);
          else {
  	    (*(short*)t->val) = 1;
            if (t->simple)
              return;
            }
          }
        }
      /* unrecognized switches */
      else {
        msg = str_new(18+str_length(argv[i]));
        sprintf(msg,"unknown option \"%s\"",argv[i]);
        arg_print_usage(a,program_name,msg,1);
        }
      }
    /* simple arguments */
    else {
      t = arg_read_get_ptr(a,0,1);
      if (!t) {
        msg = str_new(22+str_length(argv[i]));
        sprintf(msg,"excessive argument \"%s\"",argv[i]);
        arg_print_usage(a,program_name,msg,1);
        }
      else if (!type_read(argv[i],t->val,t->type)) {
        msg = str_new(24+str_length(argv[i]));
        sprintf(msg,"cannot assign value \"%s\"",argv[i]);
        arg_print_usage(a,program_name,msg,1);
        }
      t->set = 1;
      }
    }
  /* check presence and compatibility */
  arg_check(a,program_name);
  str_free(program_name);
  }

/* -------------------------------------------------------------------------- */
