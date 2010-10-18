
#ifndef GEXPRESS_H
#define GEXPRESS_H

// wxWidgets includes.
#include <wx/wx.h>
#include <wx/image.h>
#include <wx/html/helpctrl.h>
#include <wx/config.h>
#include <wx/listctrl.h>
#include <wx/notebook.h>
#include <wx/event.h>
#include <wx/fs_zip.h>

// libaffy
#include "affy.h"

// libutils
#include "utils.h"
#include "util_log.h"

extern wxArrayString celFiles;
extern wxArrayString celTypes;
extern wxArrayString celFullPaths;
extern wxString      defaultCDFDirectory;
extern wxString      cdfFile;

#define DEFAULT_CDF_DIR_KEY wxT("defaultCDFDirectory")
#define DEFAULT_CDF_DIR_VALUE wxT("")

enum {
	ID_Quit = 1,
	ID_About,
	ID_AddCELFiles,
	ID_CELFileSelectionPanel_AddFiles,
	ID_CELFileSelectionPanel_Reset,
	ID_CELFileSelectionPanel_ChooseCDFFile,
	ID_ChooseMASOutputFile,
	ID_ChooseRMAOutputFile,
	ID_ExecuteMAS,
	ID_ExecuteRMA,
	ID_FileAddCELFilesMenu,
	ID_EditPreferencesMenu,
	ID_MASPanel_MeanNormalizationCheckBox,
	ID_MASPanel_ScaleProbesetValuesCheckBox,
	ID_PreferencesOK,
	ID_PreferencesCancel,
	ID_ChooseDefaultCDFDirectory,
	ID_HelpHelpMenu,
	ID_MASPanel_QuantileNormalizationCheckBox
};

class CELFileListCtrl: public wxListCtrl
{
public:
  CELFileListCtrl(wxWindow *parent, wxWindowID id, const wxPoint& pos, 
		  const wxSize& size, long style);

  void AddEntry(wxString& filename, wxString &path, wxString& type);
};

class CELFileSelectionPanel: public wxPanel 
{
public:
  CELFileSelectionPanel(wxWindow* parent, int id, 
			const wxPoint& pos = wxDefaultPosition, 
			const wxSize& size=wxDefaultSize, long style=0);

  void AddFiles(wxCommandEvent& event);
  void ResetList(wxCommandEvent& event);
  void ChooseCDFFile(wxCommandEvent& WXUNUSED(event) );
  wxArrayString GetCELFileList();
  wxString GetCDFFile();

private:
  void set_properties();
  void do_layout();
  
  // Use events for this class
 DECLARE_EVENT_TABLE()
   
 protected:
  CELFileListCtrl *m_celFileSelectionList;
  wxButton        *m_celFileAddFilesButton;
  wxButton        *m_celFileResetButton;
  wxButton        *m_cdfFileButton;
  wxTextCtrl      *m_cdfFileTextCtrl;
  wxArrayString    m_celfiles;
};

// log everything to a text window (GUI only of course)
class LogTextCtrl : public wxLog
{
public:
  LogTextCtrl(wxTextCtrl *pTextCtrl);

private:
  // implement sink function
  virtual void DoLogString(const wxChar *szString, time_t t);

  // the control we use
  wxTextCtrl *m_pTextCtrl;
};

class MASPanel: public wxPanel 
{
public:

  MASPanel(wxWindow* parent, int id, const wxPoint& pos = wxDefaultPosition, 
	   const wxSize &size = wxDefaultSize, long style = 0);

  // Bean-style getters for parameters set in this panel.
  wxString GetOutputFile();
  int GetBackground();
  int GetMeanNormalization();
  int GetMeanNormalizationValue();
  int GetScaleProbesets();
  int GetScaleProbesetsValue();
  int GetQuantileNormalization();
  int GetBioconductorCompatability();

private:
  void set_properties();
  void do_layout();

  void MASPanel::ChooseOutputFile(wxCommandEvent& WXUNUSED(event));
  void MASPanel::OnScaleProbesetsUpdate(wxCommandEvent& WXUNUSED(event));
  void MASPanel::OnMeanNormalizationUpdate(wxCommandEvent& WXUNUSED(event));
		
  // Use events for this class
 DECLARE_EVENT_TABLE()

protected:
  wxStaticBox *m_masOutputFileSizer_staticbox;
  wxStaticBox *m_masOptionsSizer_staticbox;
  wxStaticBox *sizer_1_staticbox;
  wxCheckBox  *m_masBackgroundCheckBox;
  wxCheckBox  *m_masQuantileNormalizationCheckBox;
  wxCheckBox  *m_useMeanNormalizationCheckBox;
  wxTextCtrl  *m_useMeanNormalizationTextCtrl;
  wxCheckBox  *m_scaleProbesetsCheckBox;
  wxTextCtrl  *m_scaleProbesetsTextCtrl;
  wxCheckBox  *m_masBioconductorCompatability;
  wxTextCtrl  *m_masOutputFileTextCtrl;
  wxButton    *m_selectMASOutputFileButton;
  wxButton    *m_executeMASButton;
};

class RMAPanel: public wxPanel 
{
public:
  RMAPanel(wxWindow* parent, int id, const wxPoint& pos = wxDefaultPosition, 
	   const wxSize &size = wxDefaultSize, long style = 0);

  int      GetBackground();
  wxString GetNormalization();
  wxString GetOutputFile();
  int      GetAFFXProbeNormalization();

private:
  void set_properties();
  void do_layout();
  
  void ChooseOutputFile(wxCommandEvent& event);
  void ExecuteRMA(wxCommandEvent &event);
    
  // Use events for this class
 DECLARE_EVENT_TABLE()

protected:
  wxStaticBox *m_rmaOutputFileSizer_staticbox;
  wxStaticBox *m_rmaOptionsSizer_staticbox;
  wxCheckBox  *m_rmaBackgroundCheckBox;
  wxRadioBox  *m_rmaNormalizationOptionRadioBox;
  wxCheckBox  *m_rmaNormalizeAFFXProbes;
  wxTextCtrl  *m_rmaOutputFileTextCtrl;
  wxButton    *m_selectRMAOutputFileButton;
  wxButton    *m_executeRMAButton;
};

class PreferencesDialog: public wxDialog 
{
public:
  PreferencesDialog(wxWindow* parent, int id, const wxString& title, 
		    const wxPoint& pos = wxDefaultPosition, 
		    const wxSize& size = wxDefaultSize, 
		    long style = wxDEFAULT_DIALOG_STYLE);

  // Event Handlers
  void OnOK(wxCommandEvent& WXUNUSED(event));
  void OnCancel(wxCommandEvent& WXUNUSED(event));
  void OnChooseCDFDirectory(wxCommandEvent& WXUNUSED(event));
  
  wxString m_defaultCDFDirectory;

private:
  void set_properties();
  void do_layout();

  wxConfig     *m_config;
  
  DECLARE_EVENT_TABLE()
  
protected:
  wxStaticText *label_1;
  wxTextCtrl   *m_defaultCDFDirectoryTextCtrl;
  wxButton     *m_defaultCDFDirectoryButton;
  wxButton     *m_OKButton;
  wxButton     *m_CancelButton;
};

class MainFrame: public wxFrame 
{
public:
  MainFrame(wxWindow* parent, int id, const wxString& title, 
	    const wxPoint& pos = wxDefaultPosition, 
	    const wxSize& size = wxDefaultSize, 
	    long style = wxDEFAULT_FRAME_STYLE);
  
  // Event Handlers
  void OnEditPreferences(wxCommandEvent& event);
  void OnFileAddCELFiles(wxCommandEvent& event);
  void OnExit(wxCommandEvent& WXUNUSED(event));
  void OnAbout(wxCommandEvent& event);
  void OnHelp(wxCommandEvent &event);
  void ExecuteRMA(wxCommandEvent& WXUNUSED(event));
  void ExecuteMAS(wxCommandEvent& WXUNUSED(event));

private:
  void set_properties();
  void do_layout();
		
  // For remembering the old log output.
  wxLog       *m_logOld;
  LogTextCtrl *m_logger;

  // Use events for this class
  DECLARE_EVENT_TABLE()

protected:
  wxStaticBox           *sizer_3_staticbox;
  wxMenuBar             *m_menuBar;
  wxStatusBar           *m_statusBar;
  CELFileSelectionPanel *m_celFileSelectionPanel;
  RMAPanel              *m_rmaPanel;
  MASPanel              *m_masPanel;
  wxNotebook            *m_notebook;
  wxTextCtrl            *m_logWindow;
  wxPanel               *panel_1;
  PreferencesDialog     *m_preferencesDialog;
  wxHtmlHelpController  *m_helpController;
  wxConfig              *m_config;
};

#endif // GEXPRESS_H
