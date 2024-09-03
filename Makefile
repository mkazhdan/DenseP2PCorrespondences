LIBRARY_SRC = Include/
SRC = ./

SOFT_SOURCE = \
	SOFT1.0/cospmls.cpp \
	SOFT1.0/csecond.cpp \
	SOFT1.0/FFTcode.cpp \
	SOFT1.0/fft_grids.cpp \
	SOFT1.0/fft_grids_so3.cpp \
	SOFT1.0/FST_semi_memo.cpp \
	SOFT1.0/FST_semi_memo_fftw.cpp \
	SOFT1.0/indextables.cpp \
	SOFT1.0/makeWigner.cpp \
	SOFT1.0/naive_synthesis.cpp \
	SOFT1.0/newFCT.cpp \
	SOFT1.0/oddweights.cpp \
	SOFT1.0/OURmods.cpp \
	SOFT1.0/OURperms.cpp \
	SOFT1.0/permroots.cpp \
	SOFT1.0/primitive.cpp \
	SOFT1.0/primitive_FST.cpp \
	SOFT1.0/rotate_so3.cpp \
	SOFT1.0/rotate_so3_mem.cpp \
	SOFT1.0/seminaive.cpp \
	SOFT1.0/seminaive_fftw.cpp \
	SOFT1.0/so3_correlate_fftw.cpp \
	SOFT1.0/so3_correlate_sym.cpp \
	SOFT1.0/soft.cpp \
	SOFT1.0/soft_fftw.cpp \
	SOFT1.0/soft_fftw_pc.cpp \
	SOFT1.0/soft_fftw_wo.cpp \
	SOFT1.0/soft_sym.cpp \
	SOFT1.0/utils_so3.cpp \
	SOFT1.0/utils_vec_cx.cpp \
	SOFT1.0/weights.cpp \
	SOFT1.0/wignerTransforms.cpp \
	SOFT1.0/wignerTransforms_fftw.cpp \
	SOFT1.0/wignerTransforms_sym.cpp

AUTHALIC_EVOLUTION_VIEWER_SOURCE=AuthalicEvolutionViewer/AuthalicEvolutionViewer.cpp
CMCF_VIEWER_SOURCE=CMCFViewer/CMCFViewer.cpp
OPTICAL_FLOW_VIEWER_SOURCE=OpticalFlowViewer/OpticalFlowViewer.cpp
GET_CORRESPONDENCES_SOURCE=GetCorrespondences/GetCorrespondences.cpp
GET_HKS_SOURCE=GetHKS/GetHKS.cpp

SOFT_TARGET=SOFT1.0
AUTHALIC_EVOLUTION_VIEWER_TARGET=AuthalicEvolutionViewer
CMCF_VIEWER_TARGET=CMCFViewer
OPTICAL_FLOW_VIEWER_TARGET=OpticalFlowViewer
GET_CORRESPONDENCES_TARGET=GetCorrespondences
GET_HKS_TARGET=GetHKS


CFLAGS += -fpermissive -fopenmp -Wno-deprecated -std=c++11 -DNO_OPEN_GL
LFLAGS += -lgomp -lfftw3 -lfftw3f

CFLAGS_DEBUG = -DDEBUG -g3
LFLAGS_DEBUG =

CFLAGS_RELEASE = -O3 -DRELEASE -funroll-loops -ffast-math -g
LFLAGS_RELEASE = -O3 -g

BIN = Bin/Linux/
OBJECTS = Bin/Linux/Objects/
INCLUDE = /usr/local/include/ -IInclude

CC=gcc
CXX=g++
MD=mkdir

SOFT_OBJECTS                     =$(addprefix $(OBJECTS), $(addsuffix .o, $(basename $(SOFT_SOURCE))))
AUTHALIC_EVOLUTION_VIEWER_OBJECTS=$(addprefix $(OBJECTS), $(addsuffix .o, $(basename $(AUTHALIC_EVOLUTION_VIEWER_SOURCE))))
CMCF_VIEWER_OBJECTS              =$(addprefix $(OBJECTS), $(addsuffix .o, $(basename $(CMCF_VIEWER_SOURCE))))
OPTICAL_FLOW_VIEWER_OBJECTS      =$(addprefix $(OBJECTS), $(addsuffix .o, $(basename $(OPTICAL_FLOW_VIEWER_SOURCE))))
GET_CORRESPONDENCES_OBJECTS      =$(addprefix $(OBJECTS), $(addsuffix .o, $(basename $(GET_CORRESPONDENCES_SOURCE))))
GET_HKS_OBJECTS                  =$(addprefix $(OBJECTS), $(addsuffix .o, $(basename $(GET_HKS_SOURCE))))

all: CFLAGS += $(CFLAGS_RELEASE)
all: LFLAGS += $(LFLAGS_RELEASE)
all: $(BIN)
all: $(BIN)$(SOFT_TARGET)
all: $(BIN)$(AUTHALIC_EVOLUTION_VIEWER_TARGET)
all: $(BIN)$(CMCF_VIEWER_TARGET)
all: $(BIN)$(OPTICAL_FLOW_VIEWER_TARGET)
all: $(BIN)$(GET_CORRESPONDENCES_TARGET)
all: $(BIN)$(GET_HKS_TARGET)

debug: CFLAGS += $(CFLAGS_DEBUG)
debug: LFLAGS += $(LFLAGS_DEBUG)
debug: $(BIN)
debug: $(BIN)$(SOFT_TARGET)
debug: $(BIN)$(AUTHALIC_EVOLUTION_VIEWER_TARGET)
debug: $(BIN)$(CMCF_VIEWER_TARGET)
debug: $(BIN)$(OPTICAL_FLOW_VIEWER_TARGET)
debug: $(BIN)$(GET_CORRESPONDENCES_TARGET)
debug: $(BIN)$(GET_HKS_TARGET)

clean:
	rm -r $(BIN)

$(BIN):
	$(MD) -p $(BIN)
	$(MD) -p $(OBJECTS)$(SOFT_TARGET)
	$(MD) -p $(OBJECTS)$(AUTHALIC_EVOLUTION_VIEWER_TARGET)
	$(MD) -p $(OBJECTS)$(CMCF_VIEWER_TARGET)
	$(MD) -p $(OBJECTS)$(OPTICAL_FLOW_VIEWER_TARGET)
	$(MD) -p $(OBJECTS)$(GET_CORRESPONDENCES_TARGET)
	$(MD) -p $(OBJECTS)$(GET_HKS_TARGET)

$(BIN)$(SOFT_TARGET): $(SOFT_OBJECTS)

$(BIN)$(AUTHALIC_EVOLUTION_VIEWER_TARGET): $(AUTHALIC_EVOLUTION_VIEWER_OBJECTS) $(SOFT_OBJECTS)
	$(CXX) -o $@ $(SOFT_OBJECTS) $(AUTHALIC_EVOLUTION_VIEWER_OBJECTS) $(LFLAGS)

$(BIN)$(CMCF_VIEWER_TARGET): $(CMCF_VIEWER_OBJECTS) $(SOFT_OBJECTS)
	$(CXX) -o $@ $(SOFT_OBJECTS) $(CMCF_VIEWER_OBJECTS) $(LFLAGS)

$(BIN)$(OPTICAL_FLOW_VIEWER_TARGET): $(OPTICAL_FLOW_VIEWER_OBJECTS) $(SOFT_OBJECTS)
	$(CXX) -o $@ $(SOFT_OBJECTS) $(OPTICAL_FLOW_VIEWER_OBJECTS) $(LFLAGS)

$(BIN)$(GET_CORRESPONDENCES_TARGET): $(GET_CORRESPONDENCES_OBJECTS) $(SOFT_OBJECTS)
	$(CXX) -o $@ $(SOFT_OBJECTS) $(GET_CORRESPONDENCES_OBJECTS) $(LFLAGS)

$(BIN)$(GET_HKS_TARGET): $(GET_HKS_OBJECTS) $(SOFT_OBJECTS)
	$(CXX) -o $@ $(SOFT_OBJECTS) $(GET_HKS_OBJECTS) $(LFLAGS)

$(BIN)%.o: $(LIBRARY_SRC)%.c
	$(CC) -c -o $@ $(CFLAGS) $(INCLUDE) $<

$(BIN)%.o: $(LIBRARY_SRC)%.cpp
	$(CC) -c -o $@ $(CFLAGS) $(INCLUDE) $<

$(BIN)%.o: $(SRC)%.cpp
	$(CXX) -c -o $@ $(CFLAGS) -I$(INCLUDE) $<

$(OBJECTS)%.o: $(LIBRARY_SRC)%.c
	$(CC) -c -o $@ $(CFLAGS) $(INCLUDE) $<

$(OBJECTS)%.o: $(LIBRARY_SRC)%.cpp
	$(CC) -c -o $@ $(CFLAGS) $(INCLUDE) $<

$(OBJECTS)%.o: $(SRC)%.cpp
	$(CXX) -c -o $@ $(CFLAGS) -I$(INCLUDE) $<
