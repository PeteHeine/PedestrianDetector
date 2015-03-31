################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/PedestrianDetection.cpp \
../src/PedestrianDetector.cpp \
../src/mexFunctionsPiotrDollar.cpp 

OBJS += \
./src/PedestrianDetection.o \
./src/PedestrianDetector.o \
./src/mexFunctionsPiotrDollar.o 

CPP_DEPS += \
./src/PedestrianDetection.d \
./src/PedestrianDetector.d \
./src/mexFunctionsPiotrDollar.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/usr/local/include -I/home/pistol/Documents/matio-1.5.2/src -I/usr/include/opencv -I/usr/local/MATLAB/R2014b/extern/include -I../__GXX_EXPERIMENTAL_CXX0X__ -O0 -g3 -Wall -c -fmessage-length=0 -std=c++0x -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


