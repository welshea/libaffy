// -*- C++ -*-

#include "gene.h"

/**
 * The RMA C interface must be defined here as C-code.
 * This section handles getter functions, the setters
 * are in the parent (MainFrame).
 */
extern "C" {
	void ri_init();
	int ri_getBackground();
	int ri_isAFFXProbeNormalization();
	char *ri_getNormalization();
	char *ri_getOutputFile();
	void ri_callRMA();
};

BEGIN_EVENT_TABLE(RMAPanel, wxPanel)
EVT_BUTTON(ID_ChooseRMAOutputFile, RMAPanel::ChooseOutputFile)
END_EVENT_TABLE()

RMAPanel::RMAPanel(wxWindow *parent, int id, const wxPoint &pos, 
		   const wxSize &size, long style):
    wxPanel(parent, id, pos, size, wxTAB_TRAVERSAL)
{
  m_rmaOutputFileSizer_staticbox = new wxStaticBox(this, -1, 
						   wxT("Output File"));
  m_rmaOptionsSizer_staticbox = new wxStaticBox(this, -1, wxT("Options"));
  m_rmaBackgroundCheckBox = new wxCheckBox(this, -1, 
					   wxT("Use Background Correction"));

  const wxString m_rmaNormalizationOptionRadioBox_choices[] = 
  {
    wxT("Quantile"),
    wxT("None"),
    wxT("Mean")
  };
  
  m_rmaNormalizationOptionRadioBox = new wxRadioBox(this, 
						    -1, 
						    wxT("Normalization Method"),
						    wxDefaultPosition, 
						    wxDefaultSize, 
						    3, 
						    m_rmaNormalizationOptionRadioBox_choices, 0
   , wxRA_SPECIFY_ROWS);
  m_rmaNormalizeAFFXProbes = new wxCheckBox(this, -1, 
					    wxT("Normalize AFFX Probesets"));
  m_rmaOutputFileTextCtrl = new wxTextCtrl(this, -1, wxT(""));
  m_selectRMAOutputFileButton = new wxButton(this, 
					     ID_ChooseRMAOutputFile, 
					     wxT("Select File..."));
  m_executeRMAButton = new wxButton(this, ID_ExecuteRMA, wxT("Execute RMA"));
  
  set_properties();
  do_layout();
 
 // Set all fields to the appropriate defaults
  ri_init();
  m_rmaBackgroundCheckBox->SetValue(ri_getBackground());
  m_rmaNormalizationOptionRadioBox->SetStringSelection(wxString(ri_getNormalization(), wxConvUTF8));
  m_rmaOutputFileTextCtrl->SetValue(wxString(ri_getOutputFile(), wxConvUTF8));
  m_rmaNormalizeAFFXProbes->SetValue(ri_isAFFXProbeNormalization());
}

void RMAPanel::set_properties()
{
  m_rmaBackgroundCheckBox->SetToolTip(wxT("Use RMA background correction (default=yes)"));
  m_rmaBackgroundCheckBox->SetValue(1);
  m_rmaNormalizationOptionRadioBox->SetToolTip(wxT("Select the normalization method (quantile is default)"));
  m_rmaNormalizationOptionRadioBox->SetSelection(0);
  m_rmaNormalizeAFFXProbes->SetToolTip(wxT("Include AFFX Probesets in quantile normalization."));
  m_rmaNormalizeAFFXProbes->SetValue(1);
  m_executeRMAButton->SetDefault();
}

void RMAPanel::do_layout()
{
  wxBoxSizer       *m_rmaPanelSizer = new wxBoxSizer(wxVERTICAL);
  wxBoxSizer       *m_rmaInputSizer = new wxBoxSizer(wxVERTICAL);
  wxStaticBoxSizer *m_rmaOutputFileSizer = new wxStaticBoxSizer(m_rmaOutputFileSizer_staticbox, wxHORIZONTAL);
  wxStaticBoxSizer* m_rmaOptionsSizer = new wxStaticBoxSizer(m_rmaOptionsSizer_staticbox, wxVERTICAL);
  m_rmaOptionsSizer->Add(m_rmaBackgroundCheckBox, 0, wxALL|wxALIGN_CENTER_HORIZONTAL, 15);
  m_rmaOptionsSizer->Add(m_rmaNormalizationOptionRadioBox, 0, wxALL|wxALIGN_BOTTOM|wxALIGN_CENTER_HORIZONTAL, 15);
  m_rmaOptionsSizer->Add(m_rmaNormalizeAFFXProbes, 0, wxALIGN_CENTER_HORIZONTAL|wxADJUST_MINSIZE, 0);
  m_rmaInputSizer->Add(m_rmaOptionsSizer, 1, wxEXPAND, 0);
  m_rmaOutputFileSizer->Add(m_rmaOutputFileTextCtrl, 1, 0, 0);
  m_rmaOutputFileSizer->Add(m_selectRMAOutputFileButton, 0, wxALIGN_RIGHT, 0);
  m_rmaInputSizer->Add(m_rmaOutputFileSizer, 0, wxEXPAND, 0);
  m_rmaPanelSizer->Add(m_rmaInputSizer, 1, wxEXPAND, 3);
  m_rmaPanelSizer->Add(m_executeRMAButton, 0, wxALIGN_BOTTOM|wxALIGN_CENTER_HORIZONTAL, 0);
  SetAutoLayout(true);
  SetSizer(m_rmaPanelSizer);
  m_rmaPanelSizer->Fit(this);
  m_rmaPanelSizer->SetSizeHints(this);
}

void RMAPanel::ChooseOutputFile(wxCommandEvent& WXUNUSED(event))
{
  // Ask for file
  wxFileDialog dialog
  (
   this,
   wxT("Save Expressions..."),
   wxT(""),
   wxT("rma-exprs.txt"),
   wxT(""),
   wxSAVE
   );
  
  if (dialog.ShowModal() == wxID_OK)
    m_rmaOutputFileTextCtrl->SetValue(dialog.GetPath());
}

// Bean-type options for the field
int RMAPanel::GetBackground()
{
  return (m_rmaBackgroundCheckBox->GetValue());
}
int RMAPanel::GetAFFXProbeNormalization()
{
  return (m_rmaNormalizeAFFXProbes->GetValue());
}

wxString RMAPanel::GetNormalization()
{
  return (m_rmaNormalizationOptionRadioBox->GetStringSelection());
}
wxString RMAPanel::GetOutputFile()
{
  return (m_rmaOutputFileTextCtrl->GetValue());
}
