SUBSYSNAME := Util
MODULENAME := buildUtil

include make.initial

INSTROOT := $(DISTROOT)/$(SUBSYSNAME)/$(MODULENAME)/$(NVERS)
MODROOT := $(DEVROOT)/$(SUBSYSNAME)/$(MODULENAME)/$(NVERS)

PUBFILES = make.common getarch $(wildcard README*) buildDist config make.initial

INSTFILES := $(addprefix $(MODULENAME)/,$(PUBFILES))

setup : incs

incs : 
	- cd $(MODROOT); ln -s $(INSTFILES) .

install : installsetup
	cp -r $(PUBFILES) $(INSTROOT)
	chmod a+rx $(INSTROOT)/getarch
	chmod a+rx $(INSTROOT)/buildDist
	cp README.scorec $(INSTROOT)


installsetup : 
	mkdir -p $(INSTROOT)
