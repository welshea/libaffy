#include <wx/wx.h>
#include <wx/log.h>
#include <stdarg.h>
#include <stdlib.h>
#include "util_log.h"



/* This will involve C++ calls, but the function should not be mangled */
extern "C"
{
/*
 * This routine actually calls the wxWidget function, it needed
 * to be separate for the variable arguments.
 */
void die(const char *msg, ...)
{
	va_list ap;
	va_start(ap,msg);
	wxVLogFatalError(wxString(msg, wxConvUTF8), ap);
	va_end(ap);
	
        exit(1);
}
}

