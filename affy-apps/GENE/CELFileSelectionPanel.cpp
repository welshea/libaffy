// -*- C++ -*-

#include "gene.h"

void wx_handle_error(AFFY_ERROR *err);

BEGIN_EVENT_TABLE(CELFileSelectionPanel, wxPanel)
EVT_BUTTON(ID_CELFileSelectionPanel_AddFiles, CELFileSelectionPanel::AddFiles)
EVT_BUTTON(ID_CELFileSelectionPanel_Reset, CELFileSelectionPanel::ResetList)
EVT_BUTTON(ID_CELFileSelectionPanel_ChooseCDFFile, 
	   CELFileSelectionPanel::ChooseCDFFile)
END_EVENT_TABLE()


CELFileSelectionPanel::CELFileSelectionPanel(wxWindow* parent, 
					     int id, 
					     const wxPoint &pos, 
					     const wxSize &size, 
					     long style):
  wxPanel(parent, id, pos, size, wxTAB_TRAVERSAL)
{
  m_celFileSelectionList = new CELFileListCtrl(this, 
					       -1, 
					       wxDefaultPosition,
					       wxDefaultSize, 
					       wxLC_REPORT|wxSUNKEN_BORDER);
  m_celFileAddFilesButton = new wxButton(this, 
					 ID_CELFileSelectionPanel_AddFiles, 
					 wxT("Add Files..."));
  m_celFileResetButton = new wxButton(this, 
				      ID_CELFileSelectionPanel_Reset, 
				      wxT("Reset List"));
  m_cdfFileButton = new wxButton(this, 
				 ID_CELFileSelectionPanel_ChooseCDFFile, 
				 wxT("Choose CDF File..."));
  m_cdfFileTextCtrl = new wxTextCtrl(this, -1, wxT(""));

  set_properties();
  do_layout();
}

void CELFileSelectionPanel::set_properties()
{
  m_celFileAddFilesButton->SetDefault();
  m_cdfFileButton->SetDefault();
  m_cdfFileTextCtrl->SetMinSize(wxSize(200, 21));
}

void CELFileSelectionPanel::do_layout()
{
  wxFlexGridSizer *m_celFileSelectionSizer = new wxFlexGridSizer(4, 1, 14, 0);
  wxFlexGridSizer *m_cdfLocationlSizer = new wxFlexGridSizer(1, 2, 0, 20);
  wxBoxSizer      *sizer_2 = new wxBoxSizer(wxHORIZONTAL);

  m_celFileSelectionSizer->Add(m_celFileSelectionList, 1, wxEXPAND, 0);
  sizer_2->Add(m_celFileAddFilesButton, 0, 
	       wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 0);
  sizer_2->Add(m_celFileResetButton, 0, 
	       wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 0);
  m_celFileSelectionSizer->Add(sizer_2, 1, wxALIGN_CENTER_HORIZONTAL, 0);
  m_cdfLocationlSizer->Add(m_cdfFileButton, 0, 
			   wxALIGN_RIGHT|wxALIGN_CENTER_VERTICAL, 0);
  m_cdfLocationlSizer->Add(m_cdfFileTextCtrl, 1, wxALIGN_CENTER_VERTICAL, 0);
  m_cdfLocationlSizer->AddGrowableRow(0);
  m_cdfLocationlSizer->AddGrowableCol(0);
  m_cdfLocationlSizer->AddGrowableCol(1);
  m_celFileSelectionSizer->Add(m_cdfLocationlSizer, 1, wxEXPAND, 0);
  SetAutoLayout(true);
  SetSizer(m_celFileSelectionSizer);
  m_celFileSelectionSizer->Fit(this);
  m_celFileSelectionSizer->SetSizeHints(this);
  m_celFileSelectionSizer->AddGrowableRow(0);
  m_celFileSelectionSizer->AddGrowableCol(0);
}

void CELFileSelectionPanel::ResetList(wxCommandEvent& WXUNUSED(event))
{
  m_celFileSelectionList->DeleteAllItems();
  m_celfiles.Clear();
}

void CELFileSelectionPanel::AddFiles(wxCommandEvent& WXUNUSED(event))
{
  AFFY_ERROR *err;
  // Ask for files
  wxFileDialog dialog
    (
     this,
     wxT("Open Files"),
     wxT(""),
     wxT(""),
     wxT("CEL files (*.CEL)|*.CEL|GZip Files (*.gz)|*.GZ|All Files (*.*)|*.*"),
     wxOPEN|wxMULTIPLE|wxCHANGE_DIR
     );

  if (dialog.ShowModal() == wxID_OK)
  {
    // Get selected filenames/paths.
    wxArrayString filenames;
    wxArrayString fullpaths;

    dialog.GetFilenames(filenames);
    dialog.GetPaths(fullpaths);
			
    err=affy_get_default_error();
    err->handler=&wx_handle_error;

    // Add each to the selection list
    for (int i=0; i < fullpaths.GetCount(); i++) 
    {
      // Call to C library for cel type
      wxString x=wxString(affy_get_cdf_name_from_cel(fullpaths[i].mb_str(wxConvUTF8),err), wxConvUTF8);
      m_celFileSelectionList->AddEntry(filenames[i], fullpaths[i], x);
      m_celfiles.Add(fullpaths[i]);
    }
  }
}

void CELFileSelectionPanel::ChooseCDFFile(wxCommandEvent& WXUNUSED(event))
{
  wxString  defaultDir;
  wxConfig *p_cfg;

  // Get the default CDF file location.
  p_cfg = (wxConfig *)wxConfig::Get(false);
  p_cfg->Read(DEFAULT_CDF_DIR_KEY, &defaultDir, DEFAULT_CDF_DIR_VALUE);

  // Ask for files
  wxFileDialog dialog
    (
     this,
     wxT("Choose CDF File..."),
     defaultDir,
     wxT(""),
     wxT("CDF files (*.CDF)|*.CDF|GZip Files (*.gz)|*.GZ|All Files (*.*)|*.*"),
     wxOPEN
     );
  
  if (dialog.ShowModal() == wxID_OK)
  {
    cdfFile = dialog.GetPath();
    m_cdfFileTextCtrl->SetValue(cdfFile);
  }
}

wxArrayString CELFileSelectionPanel::GetCELFileList()
{
  return (m_celfiles);
}

wxString CELFileSelectionPanel::GetCDFFile()
{
  return (m_cdfFileTextCtrl->GetValue());
}
