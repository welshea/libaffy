// -*- C++ -*- 

#include "gene.h"

BEGIN_EVENT_TABLE(PreferencesDialog, wxDialog)
EVT_BUTTON(ID_PreferencesOK, PreferencesDialog::OnOK)
EVT_BUTTON(ID_PreferencesCancel, PreferencesDialog::OnCancel)
EVT_BUTTON(ID_ChooseDefaultCDFDirectory, 
	   PreferencesDialog::OnChooseCDFDirectory)
END_EVENT_TABLE()

PreferencesDialog::PreferencesDialog(wxWindow* parent, 
				     int id, 
				     const wxString& title, 
				     const wxPoint& pos, 
				     const wxSize& size, 
				     long style):
    wxDialog(parent, id, title, pos, size, 
	     wxDIALOG_MODAL|wxCAPTION|wxSYSTEM_MENU)
{
  label_1 = new wxStaticText(this, 
			     -1, 
			     wxT("Default CDF Directory"), 
			     wxDefaultPosition, 
			     wxDefaultSize, 
			     wxALIGN_CENTRE);

  m_defaultCDFDirectoryTextCtrl = new wxTextCtrl(this, -1, wxT(""));
  m_defaultCDFDirectoryButton = new wxButton(this, 
					     ID_ChooseDefaultCDFDirectory, 
					     wxT("Choose Directory..."));
  m_OKButton = new wxButton(this, ID_PreferencesOK, wxT(" OK "));
  m_CancelButton = new wxButton(this, ID_PreferencesCancel, wxT("Cancel"));
  
  set_properties();
  do_layout();
    
  // Get the current default
  m_config = (wxConfig *)wxConfig::Get(false);
  m_config->Read(DEFAULT_CDF_DIR_KEY, &m_defaultCDFDirectory, 
		 DEFAULT_CDF_DIR_VALUE);
  
  m_defaultCDFDirectoryTextCtrl->SetValue(m_defaultCDFDirectory);
}

void PreferencesDialog::set_properties()
{
  SetTitle(wxT("Preferences"));
  m_defaultCDFDirectoryButton->SetDefault();
}

void PreferencesDialog::do_layout()
{
  wxBoxSizer *m_PreferencesSizer         = new wxBoxSizer(wxVERTICAL);
  wxBoxSizer *m_buttonPanelSizer         = new wxBoxSizer(wxHORIZONTAL);
  wxBoxSizer *m_defaultCDFDirectorySizer = new wxBoxSizer(wxHORIZONTAL);
  m_defaultCDFDirectorySizer->Add(label_1, 0, wxALL, 5);
  m_defaultCDFDirectorySizer->Add(m_defaultCDFDirectoryTextCtrl, 0, wxALL, 5);
  m_defaultCDFDirectorySizer->Add(m_defaultCDFDirectoryButton, 0, wxALL, 5);
  m_PreferencesSizer->Add(m_defaultCDFDirectorySizer, 1, 
			  wxALIGN_CENTER_HORIZONTAL, 0);
  m_buttonPanelSizer->Add(m_OKButton, 0, 0, 0);
  m_buttonPanelSizer->Add(m_CancelButton, 0, 0, 0);
  m_PreferencesSizer->Add(m_buttonPanelSizer, 
			  1, 
			  wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 
			  0);
  SetAutoLayout(true);
  SetSizer(m_PreferencesSizer);
  m_PreferencesSizer->Fit(this);
  m_PreferencesSizer->SetSizeHints(this);
  Layout();
  Centre();
}

void PreferencesDialog::OnCancel(wxCommandEvent& WXUNUSED(event))
{
  this->EndModal(0);
}

void PreferencesDialog::OnOK(wxCommandEvent& WXUNUSED(event))
{
  m_config->Write(DEFAULT_CDF_DIR_KEY, m_defaultCDFDirectory);
  
  this->EndModal(1);
}

void PreferencesDialog::OnChooseCDFDirectory(wxCommandEvent& WXUNUSED(event))
{
  // Ask for files
  wxDirDialog dialog(this, wxT("Choose Default CDF Directory ..."),
		     m_defaultCDFDirectory);
  
  if (dialog.ShowModal() == wxID_OK) 
  {
    m_defaultCDFDirectory = dialog.GetPath();
    m_defaultCDFDirectoryTextCtrl->SetValue(m_defaultCDFDirectory);
  }
}
