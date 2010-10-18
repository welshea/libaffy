#include "gene.h"

// ----------------------------------------------------------------------------
// wxLogTextCtrl implementation
// ----------------------------------------------------------------------------

LogTextCtrl::LogTextCtrl(wxTextCtrl *pTextCtrl)
{
    m_pTextCtrl = pTextCtrl;
}

void LogTextCtrl::DoLogString(const wxChar *szString, time_t WXUNUSED(t))
{
    wxString msg;
    TimeStamp(&msg);

#if defined(__WXMAC__)
    // VZ: this is a bug in wxMac, it *must* accept '\n' as new line, the
    //     translation must be done in wxTextCtrl, not here! (FIXME)
    msg << szString << wxT('\r');
#else
    msg << szString;
#endif

    m_pTextCtrl->AppendText(msg);
}
