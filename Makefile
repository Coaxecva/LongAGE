###############################################################################
##
##     LongAGE -- Ultra-long Alignment with Gap Excision
##     Copyright (C)  Alexej Abyzov
##                
## This program is free software; you can redistribute it and/or
## modify it under the terms of the Creative Commons license
## (Attribution-NonCommerical).
## See terms at http://creativecommons.org/licenses/by-nc/2.5/legalcode
##
## Author: Alexej Abyzov
###############################################################################

VERSION = v0.1
DEFAULT_FLAGS = -DLONG_AGE_VERSION=\"$(VERSION)\" -DLONG_AGE_TIME
CXX	= g++ -std=c++11 $(DEFAULT_FLAGS)
OBJDIR = obj
OBJS   = $(OBJDIR)/LongAGE.o \
	 	 $(OBJDIR)/LongAGEaligner.o \
	 	 $(OBJDIR)/long_age_align.o

DISTRIBUTION    = $(PWD)/AGE_$(VERSION).zip
TESTFILE        = TestDecription_$(VERSION).txt
TMPDIR  	= /tmp
AGEDIR   	= age_$(VERSION)
MAINDIR		= $(TMPDIR)/$(AGEDIR)
SRCDIR   	= $(MAINDIR)/src
EXAMPLEDIR 	= $(MAINDIR)/examples

all: long_age_align

long_age_align: $(OBJS)
	$(CXX) -o $@ $(OBJS)

$(OBJDIR)/%.o: %.cpp
	@mkdir -p $(OBJDIR)
	$(CXX) -c $< -o $@
clean:
	rm -fR $(OBJDIR)
	rm -f long_age_align

distribution: clean all
	@echo Creating directory ...
	@rm -rf $(MAINDIR)
	@mkdir  $(MAINDIR)
	@mkdir  $(SRCDIR)
	@mkdir  $(EXAMPLEDIR)
	@echo Making tests ...
	@./test.sh > $(TESTFILE)
	@echo Copying files ...
	@cp *.cpp *.h  $(SRCDIR)
	@cp Makefile    $(SRCDIR)
	@cp README      $(MAINDIR)
	@cp CITATION    $(MAINDIR)
	@cp license.rtf $(MAINDIR)
	@cp test*.fa $(TESTFILE) $(EXAMPLEDIR)
	@echo Zipping ...
	@ln -s $(MAINDIR)
	@zip -qr $(DISTRIBUTION) $(AGEDIR)
	@rm $(AGEDIR)
