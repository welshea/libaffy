#include <wx/wx.h>
#include <wx/log.h>
#include <stdarg.h>
#include "util_log.h"

extern wxApp *wxapp;

/* This will involve C++ calls, but the function should not be mangled */
extern "C"
{
void info(const char *msg, ...)
{

  va_list ap;

  va_start(ap, msg);
  wxVLogMessage(wxString(msg, wxConvUTF8), ap);
  va_end(ap);

  wxapp->Yield();
  return;
}
}

