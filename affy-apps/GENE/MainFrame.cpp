// -*- C++ -*- 

#include "gene.h"
extern wxApp *wxapp;
/**
 * The C interface (RMA/MAS) must be defined here as 
 * C-code.
 * This section handles setter functions to handle
 * the state of static variables to pass to RMA.
 */
extern "C" {
  void mi_setBackground(int);
  void mi_setMeanNormalization(int);
  void mi_setMeanNormalizationValue(int);
  void mi_setScaleProbesets(int);
  void mi_setScaleProbesetsValue(int);
  void mi_setQuantileNormalization(int);
  void mi_setBioconductorCompatability(int);
  
  void ri_setBackground(int);
  void ri_setNormalization(const char *);
  void ri_setAFFXProbeNormalization(int);
  
  void mi_setOutputFile(const char *);
  void ri_setOutputFile(const char *);
  void mi_setCDFDirectory(const char *);
  void ri_setCDFDirectory(const char *);
  void mi_callMAS(char **files);
  void ri_callRMA(char **files);
};

BEGIN_EVENT_TABLE(MainFrame, wxFrame)
EVT_MENU(ID_FileAddCELFilesMenu, MainFrame::OnFileAddCELFiles)
EVT_MENU(wxID_EXIT, MainFrame::OnExit)
EVT_MENU(ID_EditPreferencesMenu, MainFrame::OnEditPreferences)
EVT_MENU(wxID_ABOUT, MainFrame::OnAbout)
EVT_MENU(ID_HelpHelpMenu, MainFrame::OnHelp)
EVT_BUTTON(ID_ExecuteRMA, MainFrame::ExecuteRMA)
EVT_BUTTON(ID_ExecuteMAS, MainFrame::ExecuteMAS)
END_EVENT_TABLE()

MainFrame::MainFrame(wxWindow* parent, int id, const wxString& title, 
		     const wxPoint& pos, const wxSize& size, long style):
  wxFrame(parent, id, title, pos, size, wxDEFAULT_FRAME_STYLE),
  m_helpController(NULL), m_preferencesDialog(NULL)
{
  // begin wxGlade: MainFrame::MainFrame
  panel_1    = new wxPanel(this, -1);
  sizer_3_staticbox = new wxStaticBox(panel_1, -1, wxT("Log Window"));
  m_notebook = new wxNotebook(this, -1, wxDefaultPosition, wxDefaultSize, 0);
  m_menuBar  = new wxMenuBar();

  // File menu.
  wxMenu *fileMenu = new wxMenu();
  fileMenu->Append(ID_FileAddCELFilesMenu, 
		   wxT("Add CEL Files..."), 
		   wxT("Add CEL files to list"), 
		   wxITEM_NORMAL);
  fileMenu->Append(wxID_EXIT, 
		   wxT("Exit"), 
		   wxT("Exit Program"), 
		   wxITEM_NORMAL);
  m_menuBar->Append(fileMenu, wxT("File"));

  // Edit menu.
  wxMenu *editMenu = new wxMenu();
  editMenu->Append(ID_EditPreferencesMenu, 
		   wxT("Preferences"), 
		   wxT(""), 
		   wxITEM_NORMAL);
  m_menuBar->Append(editMenu, wxT("Edit"));
#ifdef __WXMAC__
  wxapp->s_macPreferencesMenuItemId = ID_EditPreferencesMenu;
#endif

  wxMenu *helpMenu = new wxMenu();

  helpMenu->Append(wxID_ABOUT, 
		   wxT("About GENE"), 
		   wxT(""), 
		   wxITEM_NORMAL);
  helpMenu->Append(ID_HelpHelpMenu, 
		   wxT("View Documentation"), 
		   wxT(""), 
		   wxITEM_NORMAL);

  m_menuBar->Append(helpMenu, wxT("Help"));
#ifdef __WXMAC__
  wxapp->s_macHelpMenuTitleName     = wxT("Help");
  //    wxapp->s_macAboutMenuItemId       = ID_HelpAboutMenu;
#endif    

  SetMenuBar(m_menuBar);

  m_statusBar = CreateStatusBar(1, 0);

  m_celFileSelectionPanel = new CELFileSelectionPanel(m_notebook, 
						      -1, 
						      wxDefaultPosition, 
						      wxDefaultSize, 
						      0);
  m_rmaPanel =  new RMAPanel(m_notebook, 
			     -1, 
			     wxDefaultPosition, 
			     wxDefaultSize, 
			     0);

  m_masPanel  = new MASPanel(m_notebook, -1, wxDefaultPosition);

  m_logWindow = new wxTextCtrl(panel_1, 
			       -1, 
			       wxT(""), 
			       wxDefaultPosition, 
			       wxDefaultSize, 
			       wxTE_MULTILINE|wxTE_READONLY);

  set_properties();
  do_layout();

  // Initialize the config object; for Win32, this should be
  // HKEY_CURRENT_USER\HLMoffitt\GENE
  m_config = new wxConfig(wxT("GENE"), wxT("HLMoffitt"));
  if (m_config == NULL)
  {
    wxLogMessage(wxT("WARNING: couldn't access configuration data.\n"));
  }

  wxConfig::Set(m_config);

  // Initialize the logging window.
  m_logger = new LogTextCtrl(m_logWindow);
  m_logOld = wxLog::SetActiveTarget(m_logger);
  wxLogMessage(wxT("GENE initialized.\n"));
}

void MainFrame::set_properties()
{
  // begin wxGlade: MainFrame::set_properties
  SetTitle(wxT("GENE"));
  SetSize(wxSize(566, 687));
  int m_statusBar_widths[] = { -1 };
  m_statusBar->SetStatusWidths(1, m_statusBar_widths);
  const wxString m_statusBar_fields[] = {
    wxT("")
  };
  for(int i = 0; i < m_statusBar->GetFieldsCount(); ++i) {
    m_statusBar->SetStatusText(m_statusBar_fields[i], i);
  }
  m_logWindow->SetMinSize(wxSize(543, 73));
  // end wxGlade
}

void MainFrame::do_layout()
{
  wxFlexGridSizer  *m_frameSizer = new wxFlexGridSizer(2, 1, 0, 0);
  wxStaticBoxSizer *sizer_3      = new wxStaticBoxSizer(sizer_3_staticbox, 
							wxHORIZONTAL);

  m_notebook->AddPage(m_celFileSelectionPanel, wxT("CEL Files"));
  m_notebook->AddPage(m_rmaPanel, wxT("RMA"));
  m_notebook->AddPage(m_masPanel, wxT("MAS 5.0"));
  m_frameSizer->Add(m_notebook, 1, wxEXPAND, 0);
  sizer_3->Add(m_logWindow, 1, wxEXPAND, 0);
  panel_1->SetAutoLayout(true);
  panel_1->SetSizer(sizer_3);
  sizer_3->Fit(panel_1);
  sizer_3->SetSizeHints(panel_1);
  m_frameSizer->Add(panel_1, 1, wxEXPAND, 0);
  SetAutoLayout(true);
  SetSizer(m_frameSizer);
  m_frameSizer->AddGrowableRow(0);
  m_frameSizer->AddGrowableCol(0);
  Layout();
}

void MainFrame::OnFileAddCELFiles(wxCommandEvent& event)
{
  m_notebook->SetSelection(0);
  m_celFileSelectionPanel->AddFiles(event);
}

void MainFrame::OnExit(wxCommandEvent& WXUNUSED(event))
{
  if (m_helpController != NULL)
  {
    m_helpController->Quit();
    delete m_helpController;
  }

  if (m_config != NULL)
  {
    delete m_config;
  }

  this->Destroy();
}


void MainFrame::ExecuteMAS(wxCommandEvent& WXUNUSED(event))
{
  int           i;
  char        **celfiles;
  wxArrayString files;
  int           num_files;
  
  files = m_celFileSelectionPanel->GetCELFileList();
  num_files = files.GetCount();
  if (num_files == 0) 
  {
    wxMessageBox(wxT("No files selected"),
		 wxT("Warning"),
		 wxOK|wxICON_EXCLAMATION,
		 this);
    return;
  }
  
  celfiles = (char **)calloc((num_files+1),sizeof(char *));
  for (i=0; i < num_files; i++) 
  {
    celfiles[i] = strdup(files[i].mb_str(wxConvUTF8));
  }
  
  // Then set all options through appropriate beans
  mi_setBackground(m_masPanel->GetBackground());
  mi_setMeanNormalization(m_masPanel->GetMeanNormalization());
  mi_setMeanNormalizationValue(m_masPanel->GetMeanNormalizationValue());
  mi_setQuantileNormalization(m_masPanel->GetQuantileNormalization());
  mi_setScaleProbesets(m_masPanel->GetScaleProbesets());
  mi_setScaleProbesetsValue(m_masPanel->GetScaleProbesetsValue());
  mi_setBioconductorCompatability(m_masPanel->GetBioconductorCompatability());
  
  mi_setOutputFile(m_masPanel->GetOutputFile().mb_str());

  // CDF file can either be the file (explicitly) or
  // the directory (if other is empty).
  const char *cdf = m_celFileSelectionPanel->GetCDFFile().mb_str();
  if (strlen(cdf) != 0) 
  {
    mi_setCDFDirectory(cdf);
  } 
  else 
  {
    wxString defaultCDFDirectory;

    m_config->Read(DEFAULT_CDF_DIR_KEY, &defaultCDFDirectory, 
		   DEFAULT_CDF_DIR_VALUE);

    if (!defaultCDFDirectory.IsEmpty()) 
    {
      mi_setCDFDirectory(defaultCDFDirectory.mb_str(wxConvUTF8));
    }
  }

  // Disable timestamps during processing (cleaner)
  wxLogMessage(wxT("Starting MAS5.0\n"));
  const wxChar *timestamp = m_logger->GetTimestamp();
  m_logger->SetTimestamp(NULL);

  // Call the MAS routine
  mi_callMAS(celfiles);
 
 // Restore the timestamp
  m_logger->SetTimestamp(timestamp);
  wxLogMessage(wxT("MAS5.0 Finished\n"));
}

void MainFrame::ExecuteRMA(wxCommandEvent& WXUNUSED(event))
{
  int           i;
  char        **celfiles;
  wxArrayString files;
  int           num_files;
  
  files     = m_celFileSelectionPanel->GetCELFileList();
  num_files = files.GetCount();

  if (num_files == 0) 
  {
    wxMessageBox(wxT("No files selected"),
		 wxT("Warning"),
		 wxOK|wxICON_EXCLAMATION,
		 this);
    return;
  }
  
  celfiles = (char **)calloc((num_files+1),sizeof(char *));
  for (i=0; i < num_files; i++) 
  {
    celfiles[i] = strdup(files[i].mb_str(wxConvUTF8));
  }
  
  // Then set all options through appropriate beans
  ri_setBackground(m_rmaPanel->GetBackground());
  ri_setAFFXProbeNormalization(m_rmaPanel->GetAFFXProbeNormalization());
  ri_setNormalization(m_rmaPanel->GetNormalization().mb_str(wxConvUTF8));
  ri_setOutputFile(m_rmaPanel->GetOutputFile().mb_str(wxConvUTF8));
  
  // CDF file can either be the file (explicitly) or
  // the directory (if other is empty).
  const char *cdf = m_celFileSelectionPanel->GetCDFFile().mb_str(wxConvUTF8);
  if (strlen(cdf) != 0) 
  {
    ri_setCDFDirectory(cdf);
  } 
  else 
  {
    wxString defaultCDFDirectory;

    m_config->Read(DEFAULT_CDF_DIR_KEY, 
		   &defaultCDFDirectory, 
		   DEFAULT_CDF_DIR_VALUE);

    if (!defaultCDFDirectory.IsEmpty()) 
    {
      ri_setCDFDirectory(defaultCDFDirectory.mb_str(wxConvUTF8));
    }
  }
  
  // Disable timestamps during processing (cleaner)
  wxLogMessage(wxT("Starting RMA\n"));
  const wxChar *timestamp = m_logger->GetTimestamp();
  m_logger->SetTimestamp(NULL);
  // Call the process
  ri_callRMA(celfiles);	
  // Restore timestamps
  m_logger->SetTimestamp(timestamp);
  wxLogMessage(wxT("RMA Finished\n"));
}

// Have a load deafults and set defaults metohd
//bool Read(const wxString& key, wxString* str, const wxString& defaultVal) const
void MainFrame::OnEditPreferences(wxCommandEvent& event)
{
  if (m_preferencesDialog == NULL) 
  {
    m_preferencesDialog = new PreferencesDialog(this, -1, wxT("Preferences"));
  }

  int result = m_preferencesDialog->ShowModal();
}

void MainFrame::OnAbout(wxCommandEvent &event)
{
  wxMessageBox(wxT("GENE: A Gene Expression and Normalization Engine\n"
		   "Steven Eschrich and Andrew Hoerter\nv2.0, 2009"),
	       wxT("About GENE"),wxOK|wxICON_INFORMATION);
}

void MainFrame::OnHelp(wxCommandEvent &event)
{
  wxString  path;
 
  if (m_config->HasEntry(wxT("/InstallPath")))
    m_config->Read(wxT("/InstallPath"),path);
  else
    path = wxPathOnly(wxapp->argv[0]);

  path += wxT("/doc/usermanual.zip");
  if ( m_helpController == NULL ) {
    m_helpController = new wxHtmlHelpController();
    wxFileSystem::AddHandler(new wxZipFSHandler);
    m_helpController->AddBook(wxFileName(path),true);
  }

  m_helpController->Display(wxT("index.html"));
}
