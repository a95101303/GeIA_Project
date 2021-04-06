NAME=main_L114_code
EXEC=exe_L114_ntuple

CXX=g++
CINT=rootcint
CFLAGS=-I/usr/include/root
LDFLAGS=-O -Wall `root-config --cflags --libs` -lfftw3 -lm
SOURCE=$(NAME).C

HEADER=txn_ntupe_class_L114.h read_rawdata_L114_file.h ntuple_variables_def_L114.h  lakhs_sgolay_200m.h tanh_fit.h avmw_filter_lakhs.h 
$(EXEC): $(SOURCE) $(HEADER)
	$(CXX) $(CFLAGS) $(LDFLAGS) $(SOURCE) -o $(EXEC)

clean:
	@rm -fv $(EXEC)
