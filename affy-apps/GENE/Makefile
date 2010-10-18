
WXCONFIG=wx-config

LDFLAGS= -laffy -lutils -lwxlog -L../../libutils -L../../libaffy $(shell $(WXCONFIG) --libs)

CFLAGS=-I../../libaffy/include -I../../libutils
CXXFLAGS=$(shell $(WXCONFIG) --cxxflags) $(CFLAGS)
PROG= GENE
OBJECTS=gene.o MainFrame.o \
	CELFileListCtrl.o CELFileSelectionPanel.o \
	RMAPanel.o \
	MASPanel.o \
	PreferencesDialog.o \
	LogTextCtrl.o \
	rmaInterface.o \
	masInterface.o \
	errorInterface.o

all: $(PROG)
clean:
	rm $(OBJECTS) GENE
$(PROG): $(OBJECTS)
	g++ $(OBJECTS) -o GENE $(LDFLAGS)


$(PROG).app: mac/Info.plist mac/wxmac.icns $(PROG)
	-mkdir $@
	-mkdir $@/Contents
	-mkdir $@/Contents/MacOS
	-mkdir $@/Contents/Resources
	-mkdir $@/Contents/Resources/English.lproj
	cp mac/Info.plist $@/Contents/
	echo -n 'APPL????' > $@/Contents/PkgInfo
	cp $(PROG) $@/Contents/MacOS
	cp mac/wxmac.icns $@/Contents/Resources/
