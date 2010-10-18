// -*- C++ -*-

#include "gene.h"

/**
 * The MAS C interface must be defined here as C-code.
 * This section handles getter functions, the setters
 * are in the parent (MainFrame).
 */
extern "C" {
  void mi_init();
  int mi_getBackground();
  int mi_getMeanNormalization();
  int mi_getMeanNormalizationValue();
  int mi_getScaleProbesets();
  int mi_getScaleProbesetsValue();
  int mi_getQuantileNormalization();
  int mi_getBioconductorCompatability();
		
  char *mi_getOutputFile();
  void mi_callRMA();
};

BEGIN_EVENT_TABLE(MASPanel, wxPanel)
EVT_CHECKBOX(ID_MASPanel_MeanNormalizationCheckBox, 
	     MASPanel::OnMeanNormalizationUpdate)
EVT_CHECKBOX(ID_MASPanel_ScaleProbesetValuesCheckBox, 
	     MASPanel::OnScaleProbesetsUpdate)
EVT_BUTTON(ID_ChooseMASOutputFile, MASPanel::ChooseOutputFile)
END_EVENT_TABLE()

MASPanel::MASPanel(wxWindow *parent, int id, const wxPoint &pos, 
		   const wxSize &size, long style):
wxPanel(parent, id, pos, size, wxTAB_TRAVERSAL)
{
  m_masOptionsSizer_staticbox = new wxStaticBox(this, -1, wxT("Options"));
  m_masOutputFileSizer_staticbox = new wxStaticBox(this, -1, 
						   wxT("Output File"));
  sizer_1_staticbox = new wxStaticBox(this, -1, 
				      wxT("Probe-level Normalization"));
  m_masBackgroundCheckBox = new wxCheckBox(this, -1, 
					   wxT("Use Background Correction"));
  m_masQuantileNormalizationCheckBox = new wxCheckBox(this, 
						      ID_MASPanel_QuantileNormalizationCheckBox, wxT("Use Quantile Normalization"));
  m_useMeanNormalizationCheckBox = new wxCheckBox(this, ID_MASPanel_MeanNormalizationCheckBox, wxT("Use Mean Normalization"));
  m_useMeanNormalizationTextCtrl = new wxTextCtrl(this, -1, wxT(""));
  m_scaleProbesetsCheckBox = new wxCheckBox(this, ID_MASPanel_ScaleProbesetValuesCheckBox, wxT("Scale Probeset Values"));
  m_scaleProbesetsTextCtrl = new wxTextCtrl(this, -1, wxT("500"));
  m_masBioconductorCompatability = new wxCheckBox(this, -1, wxT("Use Bioconductor Compatability Mode"));
  m_masOutputFileTextCtrl = new wxTextCtrl(this, -1, wxT(""));
  m_selectMASOutputFileButton = new wxButton(this, ID_ChooseMASOutputFile, wxT("Select File..."));
  m_executeMASButton = new wxButton(this, ID_ExecuteMAS, wxT("Execute MAS"));

  set_properties();
  do_layout();

  // Set all fields to the appropriate defaults
  mi_init();
  m_masBackgroundCheckBox->SetValue(mi_getBackground());
  m_useMeanNormalizationCheckBox->SetValue(mi_getMeanNormalization());
  wxString s;
  s.Printf(wxT("%d"), mi_getMeanNormalizationValue());
  m_useMeanNormalizationTextCtrl->SetValue(s);
  m_useMeanNormalizationTextCtrl->Enable(mi_getMeanNormalization());
  m_scaleProbesetsCheckBox->SetValue(mi_getScaleProbesets());
  s.Printf(wxT("%d"), mi_getScaleProbesetsValue());
  m_scaleProbesetsTextCtrl->SetValue(s);
  m_scaleProbesetsTextCtrl->Enable(mi_getScaleProbesets());      
  m_masOutputFileTextCtrl->SetValue(wxString(mi_getOutputFile(), wxConvUTF8));
  m_masQuantileNormalizationCheckBox->SetValue(mi_getQuantileNormalization());
  m_masBioconductorCompatability->SetValue(mi_getBioconductorCompatability());
}

void MASPanel::set_properties()
{
  m_masBackgroundCheckBox->SetToolTip(wxT("Background correct the data before processing."));
  m_masBackgroundCheckBox->SetValue(1);
  m_masQuantileNormalizationCheckBox->SetToolTip(wxT("Perform full PM/MM Quantile Normalization prior to MAS5 algorithm."));
  m_useMeanNormalizationCheckBox->SetToolTip(wxT("Normalize the data at the probe level before calculating signal."));
  m_useMeanNormalizationTextCtrl->SetToolTip(wxT("Enter the constant value to normalize raw probes to, prior to summarization."));
  m_scaleProbesetsCheckBox->SetToolTip(wxT("Scale probesets to a constant value so results are comparable across chips."));
  m_scaleProbesetsCheckBox->SetValue(1);
  m_scaleProbesetsTextCtrl->SetToolTip(wxT("Choose the target scaling value for each chip."));
  m_masBioconductorCompatability->SetToolTip(wxT("Include masked probesets and other bioconductor-specific issues."));
}

void MASPanel::do_layout()
{
  wxBoxSizer *m_masPanelSizer = new wxBoxSizer(wxVERTICAL);
  wxBoxSizer *m_masInputSizer = new wxBoxSizer(wxVERTICAL);
  wxStaticBoxSizer *m_masOutputFileSizer = new wxStaticBoxSizer(m_masOutputFileSizer_staticbox, wxHORIZONTAL);
  wxStaticBoxSizer *m_masOptionsSizer = new wxStaticBoxSizer(m_masOptionsSizer_staticbox, wxVERTICAL);
  wxBoxSizer *m_scaleProbesetsSizer = new wxBoxSizer(wxHORIZONTAL);
  wxStaticBoxSizer *sizer_1 = new wxStaticBoxSizer(sizer_1_staticbox, wxVERTICAL);
  wxBoxSizer *m_meanNormalizationSizer = new wxBoxSizer(wxHORIZONTAL);
  m_masOptionsSizer->Add(m_masBackgroundCheckBox, 0, wxALL|wxALIGN_CENTER_HORIZONTAL, 15);
  sizer_1->Add(m_masQuantileNormalizationCheckBox, 0, wxALL|wxALIGN_CENTER_HORIZONTAL, 0);
  m_meanNormalizationSizer->Add(m_useMeanNormalizationCheckBox, 0, wxRIGHT, 10);
  m_meanNormalizationSizer->Add(m_useMeanNormalizationTextCtrl, 0, 0, 0);
  sizer_1->Add(m_meanNormalizationSizer, 0, wxALL|wxALIGN_CENTER_HORIZONTAL, 15);
  m_masOptionsSizer->Add(sizer_1, 0, wxALIGN_CENTER_HORIZONTAL|wxSHAPED, 0);
  m_scaleProbesetsSizer->Add(m_scaleProbesetsCheckBox, 0, wxRIGHT, 15);
  m_scaleProbesetsSizer->Add(m_scaleProbesetsTextCtrl, 0, 0, 0);
  m_masOptionsSizer->Add(m_scaleProbesetsSizer, 0, wxALL|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 15);
  m_masOptionsSizer->Add(m_masBioconductorCompatability, 0, wxALIGN_CENTER_HORIZONTAL, 0);
  m_masInputSizer->Add(m_masOptionsSizer, 1, wxEXPAND, 0);
  m_masOutputFileSizer->Add(m_masOutputFileTextCtrl, 1, 0, 0);
  m_masOutputFileSizer->Add(m_selectMASOutputFileButton, 0, 0, 0);
  m_masInputSizer->Add(m_masOutputFileSizer, 0, wxEXPAND, 0);
  m_masPanelSizer->Add(m_masInputSizer, 1, wxEXPAND, 0);
  m_masPanelSizer->Add(m_executeMASButton, 0, wxTOP|wxBOTTOM|wxALIGN_BOTTOM|wxALIGN_CENTER_HORIZONTAL, 10);
  SetAutoLayout(true);
  SetSizer(m_masPanelSizer);
  m_masPanelSizer->Fit(this);
  m_masPanelSizer->SetSizeHints(this);
}

void MASPanel::OnMeanNormalizationUpdate(wxCommandEvent& WXUNUSED(event))
{
  bool editable = m_useMeanNormalizationCheckBox->IsChecked();
  m_useMeanNormalizationTextCtrl->Enable(editable);
}

void MASPanel::OnScaleProbesetsUpdate(wxCommandEvent& WXUNUSED(event))
{
  bool editable = m_scaleProbesetsCheckBox->IsChecked();
  m_scaleProbesetsTextCtrl->Enable(editable);
}

void MASPanel::ChooseOutputFile(wxCommandEvent& WXUNUSED(event))
{
  // Ask for file
  wxFileDialog dialog
    (
     this,
     wxT("Save Expressions..."),
     wxT(""),
     wxT("mas-exprs.txt"),
     wxT(""),
     wxSAVE
     );

  if (dialog.ShowModal() == wxID_OK)
    m_masOutputFileTextCtrl->SetValue(dialog.GetPath());
}

wxString MASPanel::GetOutputFile()
{
  return (m_masOutputFileTextCtrl->GetValue());
}

int MASPanel::GetBackground()
{
  return (m_masBackgroundCheckBox->IsChecked());
}

int MASPanel::GetQuantileNormalization() 
{
  return (m_masQuantileNormalizationCheckBox->IsChecked());
}

int MASPanel::GetMeanNormalization()
{
  return (m_useMeanNormalizationCheckBox->IsChecked());
}
int MASPanel::GetMeanNormalizationValue()
{
  long val;

  m_useMeanNormalizationTextCtrl->GetValue().ToLong(&val);

  return (int)(val);
}
int MASPanel::GetScaleProbesets()
{
  return (m_scaleProbesetsCheckBox->IsChecked());
}

int MASPanel::GetScaleProbesetsValue()
{
  long val;

  m_scaleProbesetsTextCtrl->GetValue().ToLong(&val);

  return (int)(val);
}

int MASPanel::GetBioconductorCompatability()
{
  return (m_masBioconductorCompatability->IsChecked());
}
