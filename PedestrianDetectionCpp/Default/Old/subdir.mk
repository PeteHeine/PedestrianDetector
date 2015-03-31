################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Old/acfDetect1.cpp \
../Old/convConst.cpp \
../Old/gradientMex.cpp \
../Old/imPadMex.cpp \
../Old/rgbConvertMex.cpp 

OBJS += \
./Old/acfDetect1.o \
./Old/convConst.o \
./Old/gradientMex.o \
./Old/imPadMex.o \
./Old/rgbConvertMex.o 

CPP_DEPS += \
./Old/acfDetect1.d \
./Old/convConst.d \
./Old/gradientMex.d \
./Old/imPadMex.d \
./Old/rgbConvertMex.d 


# Each subdirectory must supply rules for building sources it contributes
Old/%.o: ../Old/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


