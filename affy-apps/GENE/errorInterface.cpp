#include <wx/wx.h>
#include <wx/msgdlg.h>
#include "affy.h"

extern wxApp *wxapp;

void wx_handle_error(AFFY_ERROR *err)
{

   char buf[1024];

   sprintf(buf,"ERROR: %s (%s) [%s:%d]\n",
		err->descr,
		affy_strerror(err->type),
		err->module,
		err->location);
   wxMessageDialog *d=new wxMessageDialog(wxapp->GetTopWindow(), wxString(buf), wxString("Error"),
	wxOK | wxICON_ERROR);
   d->ShowModal();
}

