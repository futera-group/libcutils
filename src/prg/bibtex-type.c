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
#include <cmn/string.h>
#include "prg/bibtex.h"

/* -------------------------------------------------------------------------- */

/* return internal ID of specified bibtex entry type
 
   type - name of the entry (type) */
short bibtex_type_entry_id(char *type) {
  if (str_compare(type,"comment"))
    return(BIB_TYPE_COMMENT);
  if (str_compare(type,"ARTICLE"))
    return(BIB_TYPE_ARTICLE);
  if (str_compare(type,"BOOK"))
    return(BIB_TYPE_BOOK);
  if (str_compare(type,"BOOKLET"))
    return(BIB_TYPE_BOOKLET);
  if (str_compare(type,"CONFERENCE"))
    return(BIB_TYPE_CONF);
  if (str_compare(type,"INBOOK"))
    return(BIB_TYPE_INBOOK);
  if (str_compare(type,"INCOLLECTION"))
    return(BIB_TYPE_INCOLL);
  if (str_compare(type,"INPROCEEDINGS"))
    return(BIB_TYPE_INPROC);
  if (str_compare(type,"MANUAL"))
    return(BIB_TYPE_MANUAL);
  if (str_compare(type,"MASTERSTHESIS"))
    return(BIB_TYPE_MTHESIS);
  if (str_compare(type,"MISC"))
    return(BIB_TYPE_MISC);
  if (str_compare(type,"PHDTHESIS"))
    return(BIB_TYPE_PHDTHESIS);
  if (str_compare(type,"PROCEEDINGS"))
    return(BIB_TYPE_PROC);
  if (str_compare(type,"TECHREPORT"))
    return(BIB_TYPE_TECHREP);
  if (str_compare(type,"UNPUBLISHED"))
    return(BIB_TYPE_UNPUB);
  return(BIB_TYPE_UNKNOWN);
  }

/* return name of speciefied bibtex entry type
 
   id - internal ID of the entry type  */
char* bibtex_type_entry_name(short id) {
  static char name[80]="\0";
  switch (id) {
    case BIB_TYPE_ARTICLE:   sprintf(name,"ARTICLE");       break;
    case BIB_TYPE_BOOK:      sprintf(name,"BOOK");          break;
    case BIB_TYPE_BOOKLET:   sprintf(name,"BOOKLET");       break;
    case BIB_TYPE_CONF:      sprintf(name,"CONFERENCE");    break;
    case BIB_TYPE_INBOOK:    sprintf(name,"INBOOK");        break;
    case BIB_TYPE_INCOLL:    sprintf(name,"INCOLLECTION");  break;
    case BIB_TYPE_INPROC:    sprintf(name,"INPROCEEDINGS"); break;
    case BIB_TYPE_MANUAL:    sprintf(name,"MANUAL");        break;
    case BIB_TYPE_MTHESIS:   sprintf(name,"MASTERSTHESIS"); break;
    case BIB_TYPE_MISC:      sprintf(name,"MISC");          break;
    case BIB_TYPE_PHDTHESIS: sprintf(name,"PHDTHESIS");     break;
    case BIB_TYPE_PROC:      sprintf(name,"PROC");          break;
    case BIB_TYPE_TECHREP:   sprintf(name,"TECHREPORT");    break;
    case BIB_TYPE_UNPUB:     sprintf(name,"UNPUBLISHED");   break;
    default: sprintf(name,"MISC"); break;
    }
  return(name);
  }

/* -------------------------------------------------------------------------- */

/* return internal ID of specified bibtex item type
 
   type - name of the itey */
short bibtex_type_item_id(char *type) {
  if (str_compare(type,"author"))
    return(BIB_ITEM_AUTHOR);
  if (str_compare(type,"title"))
    return(BIB_ITEM_TITLE);
  if (str_compare(type,"booktitle"))
    return(BIB_ITEM_BOOKTITLE);
  if (str_compare(type,"journal"))
    return(BIB_ITEM_JOURNAL);
  if (str_compare(type,"year"))
    return(BIB_ITEM_YEAR);
  if (str_compare(type,"volume"))
    return(BIB_ITEM_VOLUME);
  if (str_compare(type,"chapter"))
    return(BIB_ITEM_CHAPTER);
  if (str_compare(type,"number"))
    return(BIB_ITEM_NUMBER);
  if (str_compare(type,"pages"))
    return(BIB_ITEM_PAGES);
  if (str_compare(type,"file"))
    return(BIB_ITEM_FILE);
  if (str_compare(type,"editor"))
    return(BIB_ITEM_EDITOR);
  if (str_compare(type,"publisher"))
    return(BIB_ITEM_PUBLISH);
  if (str_compare(type,"organization"))
    return(BIB_ITEM_ORGANIZ);
  if (str_compare(type,"school"))
    return(BIB_ITEM_SCHOOL);
  if (str_compare(type,"address"))
    return(BIB_ITEM_ADDRESS);
  if (str_compare(type,"keywords"))
    return(BIB_ITEM_KEYWORD);
  if (str_compare(type,"owner"))
    return(BIB_ITEM_OWNER);
  if (str_compare(type,"timestamp"))
    return(BIB_ITEM_TIMESTAMP);
  return(BIB_ITEM_UNKNOWN);
  }

/* return name of specified bibtex entry item

   id - internal ID of file type */
char* bibtex_type_item_name(short id) {
  static char name[80]="\0";
  switch (id) {
    case BIB_ITEM_LABEL:     sprintf(name,"label"); break;
    case BIB_ITEM_TYPE:      sprintf(name,"type"); break;
    case BIB_ITEM_AUTHOR:    sprintf(name,"author"); break;
    case BIB_ITEM_TITLE:     sprintf(name,"title"); break;
    case BIB_ITEM_JOURNAL:   sprintf(name,"journal"); break;
    case BIB_ITEM_YEAR:      sprintf(name,"year"); break;
    case BIB_ITEM_PAGES:     sprintf(name,"pages"); break;
    case BIB_ITEM_OWNER:     sprintf(name,"owner"); break;
    case BIB_ITEM_TIMESTAMP: sprintf(name,"timestamp"); break;
    case BIB_ITEM_VOLUME:    sprintf(name,"volume"); break;
    case BIB_ITEM_FILE:      sprintf(name,"file"); break;
    case BIB_ITEM_NUMBER:    sprintf(name,"number"); break;
    case BIB_ITEM_ORGANIZ:   sprintf(name,"organization"); break;
    case BIB_ITEM_SCHOOL:    sprintf(name,"school"); break;
    case BIB_ITEM_ADDRESS:   sprintf(name,"address"); break;
    case BIB_ITEM_CHAPTER:   sprintf(name,"chapter"); break;
    case BIB_ITEM_PUBLISH:   sprintf(name,"publisher"); break;
    case BIB_ITEM_EDITOR:    sprintf(name,"editor"); break;
    case BIB_ITEM_BOOKTITLE: sprintf(name,"book title"); break;
    case BIB_ITEM_KEYWORD:   sprintf(name,"keywords"); break;
    case BIB_ITEM_UNKNOWN:   sprintf(name,"unknown"); break;
    default: sprintf(name,"uknown"); break;
    }
  return(name);
  }

/* -------------------------------------------------------------------------- */

/* return internal ID of specified bibtex file type
 
   type - name of the file type */
short bibtex_type_file_id(char *type) {
  if (str_compare(type,"PDF"))
    return(BIB_FILE_PDF);
  return(BIB_FILE_UNKNOWN);
  }

/* return name of specified bibtex file type

   id - internal ID of file type */
char* bibtex_type_file_name(short id) {
  static char name[80]="\0";
  switch (id) {
    case BIB_FILE_PDF: sprintf(name,"PDF"); break;
    default: sprintf(name,"TXT"); break;
    }
  return(name);
  }

/* -------------------------------------------------------------------------- */
