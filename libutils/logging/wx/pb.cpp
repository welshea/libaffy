#include <wx/wx.h>
#include <wx/progdlg.h>
#include <stdarg.h>
#include "util_log.h"

extern wxApp *wxapp;

#define MAXBUF 1024
#define PB_NUM_TICKS 20

extern "C"
{
 void pb_init(LIBUTILS_PB_STATE *pbs)
 {
   if ( pbs == NULL )
      return;

   pbs->depth=0;
 }
}

/* This will involve C++ calls, but the function should not be mangled */
extern "C"
{
void pb_begin(LIBUTILS_PB_STATE *pbs, unsigned int max, char *title, ...)
{
  char buf[1024];
  va_list ap;
  int i;

  if ( pbs == NULL ) 
     return;
  assert(pbs->depth < LIBUTILS_MAX_PB_DEPTH);
  pbs->cur_ticks[pbs->depth] =0;
  pbs->tick_interval[pbs->depth]=ceil((double)max/PB_NUM_TICKS);
  pbs->max[pbs->depth]=PB_NUM_TICKS;

  /* Handle the optional printf-like arguments */
  if ( title != NULL )
  {
     va_start(ap, title);
     vsprintf(buf, title, ap);
     va_end(ap);
  }
  for (i=strlen(buf); i < 120; i++) {
     strcat(buf," ");
  }

  /* Create the progress bar */
  pbs->ptr[pbs->depth] = (void *) new wxProgressDialog(wxString(buf, wxConvUTF8), wxT("Initializing..."), pbs->max[pbs->depth],NULL,wxPD_AUTO_HIDE|wxPD_APP_MODAL|wxRESIZE_BORDER);

  pbs->depth++;

  wxapp->Yield();
  return;
}
}

/* This will involve C++ calls, but the function should not be mangled */
extern "C"
{
void pb_cleanup(LIBUTILS_PB_STATE *pbs)
{
   return;
}
}
/* This will involve C++ calls, but the function should not be mangled */
extern "C"
{
/*
 * Execute a tick of the progress bar. More than a single tick can be
 * accomplished by specifying a tick_sz > 1. An optional msg can be
 * provided ("" for nothing)
 */
void pb_tick(LIBUTILS_PB_STATE *pbs, unsigned int tick_sz, char *msg, ...)
{
  int i;
  wxProgressDialog *p;
  char buf[MAXBUF]="";
  va_list ap;
  unsigned int oldticks=0;
  unsigned int oldinterval, newinterval;

  if ( pbs == NULL )
     return;
  assert(pbs->depth > 0);
  /* Current progress bar is -1 from depth */
  i=pbs->depth-1;
  
  /* Update current tick count */ 
  oldinterval=pbs->cur_ticks[i]/pbs->tick_interval[i];
  pbs->cur_ticks[i]+=tick_sz;
  newinterval=pbs->cur_ticks[i]/pbs->tick_interval[i];

  /* Only update message if there is one. */
  if ( msg != NULL && msg[0] != '\0' ) {
     va_start(ap, msg);
     vsprintf(buf, msg, ap);
     va_end(ap);
  }
  /* Only run wx update if needed */
  if ( buf[0] != '\0' && oldinterval != newinterval ) {
     /* Get the wx ptr */
     p=(wxProgressDialog *)pbs->ptr[i];
     p->Update(newinterval,wxString(buf, wxConvUTF8));
     wxapp->Yield();
  }

}

void pb_msg(LIBUTILS_PB_STATE *pbs, char *msg, ...)
{
  va_list ap;
  char update_msg[MAXBUF]="";
  wxProgressDialog *p;
  int i;

  assert(msg != NULL);
  if ( pbs == NULL )
     return;
  assert(pbs->depth > 0);

  /* Handle the optional printf-like arguments */
  if ( msg[0] != '\0' ) {
     va_start(ap, msg);
     vsprintf(update_msg, msg, ap);
     va_end(ap);
  }

  i=pbs->depth-1;
  p=(wxProgressDialog *)pbs->ptr[i];
  p->Update(pbs->cur_ticks[i]/pbs->tick_interval[i], wxString(update_msg,wxConvUTF8));
  
  wxapp->Yield();
}


void pb_finish(LIBUTILS_PB_STATE *pbs, char *msg, ...)
{
  int i;
  va_list ap;
  wxProgressDialog *p;
  char buf[MAXBUF]="";

  assert (msg!=NULL);
  if ( pbs==NULL )
     return;
  assert(pbs->depth>0);

  if ( msg[0] != '\0')
  { 
    va_start(ap,msg);
    vsprintf(buf, msg, ap);
    va_end(ap);
  }
  i=pbs->depth-1;
  p=(wxProgressDialog *)(pbs->ptr[i]);

  p->Update(pbs->max[i], wxString(buf,wxConvUTF8));
  wxapp->Yield();

  delete (wxProgressDialog *)(pbs->ptr[i]);
  pbs->ptr[i]=NULL;
  pbs->depth--;

}
}

