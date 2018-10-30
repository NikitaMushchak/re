#define macros
EXECUTABLE_NAME = planar3D.exe
CPP_FLAGS = /Ox /favor:%PROCESSOR_ARCHITECTURE% \
	/D_TARGET_="%PROCESSOR_ARCHITECTURE%" \
	/D_TIMESTAMP_="%date% %time%" \
	/D_PLATFORM_="%OS%"
OPENMP_FLAGS = /DOPENMP -openmp
DIR_SRC = .\sources
DIR_DATA = $(DIR_SRC)\data
DIR_INCLUDE = .\include
DIR_BUILD = .\build
SRC_FILES= $(DIR_SRC)\*.cc
DATA_FILES= $(DIR_DATA)\*.txt

#SHELL = powershell.exe -NoLogo -NoProfile -NonInteractive -ExecutionPolicy Bypass -Command

all: clean create_dirs planar3D

openmp: clean create_dirs planar3D.openmp

clean:
	@if exist $(DIR_BUILD) rmdir /S /Q $(DIR_BUILD)

create_dirs:
	@if not exist $(DIR_BUILD) mkdir $(DIR_BUILD)
	@if not exist $(DIR_BUILD)\InitialConditions mkdir $(DIR_BUILD)\InitialConditions
	@if not exist $(DIR_BUILD)\Results mkdir $(DIR_BUILD)\Results
	@if not exist $(DIR_BUILD)\Results\Concentration mkdir $(DIR_BUILD)\Results\Concentration
	@if not exist $(DIR_BUILD)\Results\Opening mkdir $(DIR_BUILD)\Results\Opening
	@if exist $(DIR_DATA) if exist $(DATA_FILES) copy $(DATA_FILES) $(DIR_BUILD)\InitialConditions

$(EXECUTABLE_NAME) : $(SRC_FILES)
	cl /EHsc $(CPP_FLAGS) /Fe$(DIR_BUILD)\$(EXECUTABLE_NAME) /I$(DIR_INCLUDE) $(SRC_FILES)
	del *.obj
	
$(EXECUTABLE_NAME).openmp : $(SRC_FILES)
	cl /EHsc $(CPP_FLAGS) $(OPENMP_FLAGS) /Fe$(DIR_BUILD)\$(EXECUTABLE_NAME) /I$(DIR_INCLUDE) $(SRC_FILES)
	del *.obj
	
planar3D: $(EXECUTABLE_NAME)

planar3D.openmp: $(EXECUTABLE_NAME).openmp


#http://www.bojankomazec.com/2011/10/how-to-use-nmake-and-makefile.html
#https://msdn.microsoft.com/en-us/library/fwkeyyhe.aspx
