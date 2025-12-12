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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cmn/message.h"
#include "cmn/queue.h"
#include "cmn/string.h"

/* -------------------------------------------------------------------------- */

/* read line form input stream (unspecified lenght)

   file - the input stream */
char* str_read_line_new(FILE *file) {
  static char buf[STR_BUFFER_LENGTH];
  unsigned i,n,length;
  char *t,*s = NULL;
  struct queue *q;
  /* first buffer reading */
  if (fgets(buf,STR_BUFFER_LENGTH,file)) {
    n = str_length(buf);
    /* line fitted to buffer */
    if (n<(STR_BUFFER_LENGTH-1) || buf[STR_BUFFER_LENGTH-1]=='\n') {
      s = str_copy_new(buf);
      return(s);
      }
      /* buffer is too small, save segments to queue */
    length = n;
    q = queue_alloc();
    queue_sadd(q,buf);
    /* read line segments */
    while (fgets(buf,STR_BUFFER_LENGTH,file)) {
      if (feof(file)) break;
      if (ferror(file)) break;
      queue_sadd(q,buf);
      n = str_length(buf);
      length += n;
      if (n<(STR_BUFFER_LENGTH-1) || buf[STR_BUFFER_LENGTH-1]=='\n')
        break;
      }
    /* merge segments */
    if (q->num) {
      s = str_new(length);
      n = 0;
      while (q->num) {
        t = queue_get(q);
        for (i=0; i<str_length(t); i++) 
          s[n++] = t[i];
        }
      s[n] = '\0';
      }
    /* clean memory */
    queue_free(q);
    }
  return(s);
  }

/* -------------------------------------------------------------------------- */

/* read line from standard input and return first word

   word - the word which will be returned */
void str_read_word(char *word) {
  int i,j,n;
  char *line = NULL;
  /* read line */
  line = str_read_line_new(stdin);
  if (!line)
    msg_error("unexpected error while reading from standard input",1);
  /* copy first word */
  i = j = 0;
  n = str_length(line);
  while (i<n && (line[i]==' ' || line[i]=='\t'))
    i++;
  while (i<n && line[i]!='\n' && line[i]!=' ' && line[i]!='\t') {
    word[j] = line[i];
    j++;
    i++;
    }
  word[j] = '\0';
  /* clean memory */
  str_free(line);
  }

/* read line from standard input and return new string with first word */
char* str_read_word_new(void) {
  int i,j,k,n;
  char *line = NULL,*word = NULL;
  /* read line */
  line = str_read_line_new(stdin);
  if (!line)
    msg_error("unexpected error while reading from standard input",1);
  /* locate the first word */
  n = str_length(line);
  for (i=0; i<n && (line[i]==' ' || line[i]=='\t'); i++);
  for (j=i; j<n && line[j]!='\n' && line[j]!=' ' && line[j]!='\t'; j++);
  /* copy the word */
  if (i<n && j>i) {
    word = str_new(j-i);
    for (k=i; k<j; k++)
      word[k-i] = line[k];
    word[k-i] = '\0';
    }
  /* clean memory */
  str_free(line);
  return(word);
  }

/* read one word from the given string

   str  - the string from where the word is read
   word - the word which is read
   del  - deliminer between words
   id   - starting position for reading */
unsigned str_read_word_s(char *str, char *word, char del, unsigned id) {
  unsigned i,j,n;
  n = str_length(str);
  for (i=id; i<n && str[i]==del; i++);
  for (j=0; i<n && str[i]!=del; i++,j++)
    word[j] = str[i];
  word[j] = '\0';
  return(i);
  }

/* -------------------------------------------------------------------------- */

/* read yes or no answer to given question

   s - the question string */
short str_read_yesno(char *s) {
  char *w;
  short v = -1;
  while (v<0) {
    printf("%s [y/n] ",s);
    w = str_read_word_new();
    if (w && str_length(w)==1) 
      v = ((w[0]=='y' || w[1]=='Y') ? 1 : ((w[0]=='n' || w[1]=='N') ? 0 : -1));
    str_free(w);
    }
  return(v);
  }

/* read yes or no answer to given question with pre-set default 

   s - the question string
   d - the default answer */
short str_read_yesno_def(char *s, short d) {
  char *w;
  short v=-1;
  while (v<0) {
    printf("%s [%c/%c] ",s,(d ? 'Y' : 'y'),(!d ? 'N' : 'n'));
    w = str_read_word_new();
    if (!w)
      v = d;
    else if (str_length(w)==1) 
      v = ((w[0]=='y' || w[1]=='Y') ? 1 : ((w[0]=='n' || w[1]=='N') ? 0 : -1));
    str_free(w);
    }
  return(v);
  }

/* -------------------------------------------------------------------------- */

/* ask for and integer number

   s - question that will be printed before reading */
int str_read_inum(char *s) {
  char *t;
  int i,val,valid;
  short ok = 0;
  while (!ok) {
    printf("%s ",s);
    /* read word */
    t = str_read_word_new();
    /* check characters */
    valid = 0;
    if (t) {
      for (i=0,valid=1; i<str_length(t); i++)
        if ((t[i]<'0' || t[i]>'9') && 
           !(t[i]=='-' && i==0) &&
           !(t[i]=='+' && i==0)) {
          printf("Invalid value, integer number expected\n");
          valid = 0;
          break;
          }
      }
    /* convert the value */
    if (valid) {
      if (sscanf(t,"%d",&val)!=1)
        printf("Invalid value, integer number expected\n");
      else
        ok=1;
      }
    str_free(t);
    }
  return(val);
  }

/* ask for and read unsigned integer number

   s - question that will be printed before reading */
unsigned str_read_unum(char *s) {
  char *t;
  int i,valid;
  unsigned val;
  short ok = 0;
  while (!ok) {
    printf("%s ",s);
    /* read word */
    t = str_read_word_new();
    /* check characters */
    valid = 0;
    if (t) {
      for (i=0,valid=1; i<str_length(t); i++)
        if ((t[i]<'0' || t[i]>'9') && 
           !(t[i]=='+' && i==0)) {
          printf("Invalid value, unsigned integer number expected\n");
          valid = 0;
          break;
          }
      }
    /* convert the value */
    if (valid) {
      if (sscanf(t,"%u",&val)!=1)
        printf("Invalid value, unsigned integer number expected\n");
      else
        ok=1;
      }
    str_free(t);
    }
  return(val);
  }

/* ask for and read unsigned integer number, offer default

   s - question that will be printed before reading
   d - default value  */
unsigned str_read_unum_def(char *s, unsigned d) {
  char *t;
  int i,valid;
  unsigned val;
  short ok = 0;
  while (!ok) {
    printf("%s [%d] ",s,d);
    /* read word */
    t = str_read_word_new();
    /* default value */
    if (!t)
      return(d);
    /* check characters */
    for (i=0,valid=1; i<str_length(t); i++)
      if ((t[i]<'0' || t[i]>'9') && 
         !(t[i]=='+' && i==0)) {
        printf("Invalid value, unsigned integer number expected\n");
        valid = 0;
        break;
        }
    /* convert the value */
    if (valid) {
      if (sscanf(t,"%u",&val)!=1)
        printf("Invalid value, unsigned integer number expected\n");
      else
        ok=1;
      }
    str_free(t);
    }
  return(val);
  }

/* ask for and read unsigned integer number from given interval

   s - question that will be printed before reading
   a - lower boundary of the interval (greater or equal)
   b - upper boundary of the interval (less or equal) */
unsigned str_read_unum_int(char *s, unsigned a, unsigned b) {
  char *t;
  int i,valid;
  unsigned val;
  short ok = 0;
  while (!ok) {
    printf("%s ",s);
    /* read word */
    t = str_read_word_new();
    /* check characters */
    valid = 0;
    if (t) {
      for (i=0,valid=1; i<str_length(t); i++)
        if ((t[i]<'0' || t[i]>'9') && 
           !(t[i]=='+' && i==0)) {
          printf("Invalid value, unsigned integer number expected\n");
          valid = 0;
          break;
          }
      }
    /* convert the value */
    if (valid) {
      if (sscanf(t,"%u",&val)!=1)
        printf("Invalid value, unsigned integer number expected\n");
      else if (val<a || val>b)
        printf("Value from interval [%d,%d] expected\n",a,b);
      else
        ok=1;
      }
    str_free(t);
    }
  return(val);
  }

/* ask for and read double precision real number

   s - question that will be printed before reading */
double str_read_fnum(char *s) {
  char *t;
  int i,valid;
  double val;
  short ok = 0;
  while (!ok) {
    printf("%s ",s);
    /* read word */
    t = str_read_word_new();
    /* check characters */
    valid = 0;
    if (t) {
      for (i=0,valid=1; i<str_length(t); i++)
        if (!((t[i]>='0' && t[i]<='9') || t[i]=='e' || t[i]=='E' ||
             (t[i]=='.' && (i==0 || (t[i-1]>='0' && t[i-1]<='9'))) ||
             ((t[i]=='+' || t[i]=='-') && 
              (i==0 || t[i-1]=='e' || t[i-1]=='E')))) {
          printf("Invalid value, real number expected\n");
          valid = 0;
          break;
          }
      }
    /* convert the value */
    if (valid) {
      if (sscanf(t,"%lf",&val)!=1)
        printf("Invalid value, real number expected\n");
      else
        ok=1;
      }
    str_free(t);
    }
  return(val);
  }

/* -------------------------------------------------------------------------- */
