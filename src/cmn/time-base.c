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
#include <time.h>
#include "cmn/string.h"

/* -------------------------------------------------------------------------- */

long int tm_start=0;

/* -------------------------------------------------------------------------- */

/* return actual time */
long int time_actual(void) {
  return(time(NULL));
  }

/* -------------------------------------------------------------------------- */

/* set starting time and return its value */
long int time_start(void) {
  return(time(&tm_start));
  }

/* return time from starting point */
long int time_stop(void) {
  return(time(NULL)-tm_start);
  }

/* -------------------------------------------------------------------------- */

/* convert absolute time to formatted string

   tsec  - absolute time in second from 00:00:00 1.1.1970
   stime - storage for the formatted time string */
void time_str(long int tsec, char *stime) {
  struct tm *data;
  data = localtime(&tsec);
  sprintf(stime,"%02d.%02d.%4d %02d:%02d:%02d",data->tm_mday,data->tm_mon+1,
    data->tm_year+1900,data->tm_hour,data->tm_min,data->tm_sec);
  }

/* convert absolute time to formatted time string

   tsec  - absolute time in second from 00:00:00 1.1.1970
   stime - storage for the formatted time string */
void time_str_time(long int tsec, char *stime) {
  struct tm *data;
  data = localtime(&tsec);
  sprintf(stime,"%02d:%02d:%02d",data->tm_hour,data->tm_min,data->tm_sec);
  }

/* convert absolute time to formatted date string

   tsec  - absolute time in second from 00:00:00 1.1.1970
   sdate - storage for the formatted time string */
void time_str_date(long int tsec, char *sdate) {
  struct tm *data;
  data = localtime(&tsec);
  sprintf(sdate,"%02d.%02d.%4d",data->tm_mday,data->tm_mon+1,
    data->tm_year+1900);
  }

/* convert time interval in seconds to formatted string 

   tsec  - time in seconds
   stime - place for formatted time string */
char* time_str_int(long int tsec) {
  static char stime[80];
  short space = 0;
  unsigned i,data[6] = {0,0,0,0,0,0};
  unsigned tdiv[6] = {365*24*60*60,30*24*60*60,24*60*60,60*60,60,1};
  char desc[5][7] = {"year\0","month\0","day\0","hour\0","minute\0"};
  /* prepare data */
  for (i=0; i<6; i++) {
    data[i] = tsec/tdiv[i];
    tsec %= tdiv[i];
    }
  /* create string */
  stime[0]='\0';
  for (i=0; i<5; i++) {
    if (data[i]) {
      if (space)
        sprintf(stime+str_length(stime)," ");
      sprintf(stime+str_length(stime),"%d %s",data[i],desc[i]);
      if (data[i]>1)
        sprintf(stime+str_length(stime),"s");
      space=1;
      }
    }
  if (space)
    sprintf(stime+str_length(stime)," ");
  sprintf(stime+str_length(stime),"%d second",data[5]);
  if (data[5]!=1)
    sprintf(stime+str_length(stime),"s");
  return(stime);
  }

/* -------------------------------------------------------------------------- */
