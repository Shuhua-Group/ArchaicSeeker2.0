################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../EMExp.cpp \
../MultiWaver.cpp \
../ParamExp.cpp \
../Utils.cpp 

OBJS += \
./EMExp.o \
./MultiWaver.o \
./ParamExp.o \
./Utils.o 

CPP_DEPS += \
./EMExp.d \
./MultiWaver.d \
./ParamExp.d \
./Utils.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<" -fopenmp
	@echo 'Finished building: $<'
	@echo ' '


