
#include "gene.h"

CELFileListCtrl::CELFileListCtrl(wxWindow      *parent, 
				 wxWindowID     id, 
				 const wxPoint &pos, 
				 const wxSize  &size, 
				 long           style)
  :wxListCtrl(parent,id,pos,size,style)
{
  int width, height;
  wxListItem li;

  // Figure out overall size for list ctrl, to size last column.
  GetSize(&width, &height);

  // Column Headings
  li.m_image =- 1;
  li.m_mask  = wxLIST_MASK_TEXT;
  li.m_text  = wxT("Filename");
  InsertColumn(0, li);
  SetColumnWidth(0, wxLIST_AUTOSIZE_USEHEADER);

  li.m_text = wxT("CEL Type");
  InsertColumn(1, li);
  SetColumnWidth(1, wxLIST_AUTOSIZE_USEHEADER);

  li.m_text = wxT("Path");	
  InsertColumn(2, li);
  SetColumnWidth(2, width);
}

void CELFileListCtrl::AddEntry(wxString &filename, wxString &path, 
			       wxString &type)
{
  int  nextIndex = GetItemCount();
  long tmp       = InsertItem(nextIndex, filename,-1);

  SetItemData(tmp, nextIndex);
  SetItem(nextIndex, 1, type);
  SetItem(nextIndex, 2, path);

  SetColumnWidth(0, wxLIST_AUTOSIZE);
  SetColumnWidth(1, wxLIST_AUTOSIZE);
  SetColumnWidth(2, wxLIST_AUTOSIZE);
}

